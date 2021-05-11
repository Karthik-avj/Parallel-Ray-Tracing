#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 1920
#define M 1080
#define NUM_OBJ 4
#define OBJ_LEN 15
#define XL_MIN 4.5
#define XL_MAX 5.5
#define ZL_MIN 4.5
#define ZL_MAX 5.5
#define YL 5.0
#define L_RANDOM 7
#define MAX_DEPTH 2
#define AA 2

void ray_direction(float* origin, float* point, float* vector){
    float dr[3];
    dr[0] = point[0]-origin[0];
    dr[1] = point[1]-origin[1];
    dr[2] = point[2]-origin[2];

    float norm = sqrt(dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]);
    vector[0] = dr[0]/norm;
    vector[1] = dr[1]/norm;
    vector[2] = dr[2]/norm;
}

void reflected_direction(float* incoming, float* normal, float* reflected){
    float dr[3];
    float dp = incoming[0]*normal[0] + incoming[1]*normal[1] + incoming[2]*normal[2];
    dr[0] = incoming[0] - 2*dp*normal[0];
    dr[1] = incoming[1] - 2*dp*normal[1];
    dr[2] = incoming[2] - 2*dp*normal[2];

    float norm = sqrt(dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]);
    reflected[0] = dr[0]/norm;
    reflected[1] = dr[1]/norm;
    reflected[2] = dr[2]/norm;
}

void sphere_intersection(float* origin, float* ray_direction, float* center, float* radius, float* dist){
    float b, c, disc;
    float dr[3];
    dr[0] = origin[0]-center[0];
    dr[1] = origin[1]-center[1];
    dr[2] = origin[2]-center[2];

    b = 2 * (ray_direction[0]*dr[0]+
    ray_direction[1]*dr[1]+
    ray_direction[2]*dr[2]);

    c = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2] - *radius * *radius;

    disc = b*b - (4*c);
    if(disc<=0){
        *dist = -1.0;
        return;
    }
    else{
        float val1, val2;
        float sqd = sqrt(disc);
        val1 = (-b + sqd)/2;
        val2 = (-b - sqd)/2;
        if(val1>0 && val2>0){
            *dist = fmin(val1, val2);
            return;
        }
        else{
            *dist = -1.0;
            return;
        }
    }
}

void nearest_intersection_object(float *objects, float *origin, float *ray_direction, float *min_dist, int *object_idx){
	float distances[NUM_OBJ];
	for(int i=0; i<(NUM_OBJ); i++){
		float center[] = {objects[i*OBJ_LEN+0], objects[i*OBJ_LEN+1], objects[i*OBJ_LEN+2]};
		sphere_intersection(origin, ray_direction, center, &objects[i*OBJ_LEN+3], &distances[i]);
	}

	for(int i=0; i<(NUM_OBJ); i++){
		if((distances[i]<(*min_dist)) && (distances[i]>0)){
			*min_dist = distances[i];
			*object_idx = i;
		}
	}
}

void shadowed(int *is_shad, float *normal, float *light_dir, float *shifted_point, float *min_dist, float *origin, float *ray_dir,
              float *light_source, float *objects, int *object_idx){
    // line 43
    float intersection_point[3];
    for (int i=0; i<3; i++){
        intersection_point[i] = *min_dist * ray_dir[i] + origin[i];
    }

    // line 44
    float object_center[] = {objects[*object_idx * OBJ_LEN], objects[*object_idx * OBJ_LEN + 1], objects[*object_idx * OBJ_LEN + 2]};
    ray_direction(object_center, intersection_point, normal);


    // line 45
    for (int i=0; i<3; i++){
        shifted_point[i] = 0.0001 * normal[i] + intersection_point[i];
    }

    // line 46
    ray_direction(shifted_point, light_source, light_dir);

    // line 48
    float min_distance = __INT_MAX__;
    int useless = -1;
    nearest_intersection_object(objects, shifted_point, light_dir, &min_distance, &useless);

    // line 49
    float intersection_to_light_dist = sqrt((light_source[0]-intersection_point[0])*(light_source[0]-intersection_point[0]) + (light_source[1]-intersection_point[1])*(light_source[1]-intersection_point[1]) + (light_source[2]-intersection_point[2])*(light_source[2]-intersection_point[2]));

    // line 50
    if (min_distance < intersection_to_light_dist) {
        *is_shad = 1;
    }
    else{
        *is_shad = 0;
    }
}

void color(float* normal_surface, float* light_intersection, float* ray_dir, float* object, float* light, float* reflection, float* illumination){
    // float illumination[3] = {0, 0, 0};
    float ambient[3] = {object[4]*light[5], object[5]*light[6], object[6]*light[7]};

    float nl_dp = normal_surface[0] * light_intersection[0] +
                  normal_surface[1] * light_intersection[1] +
                  normal_surface[2] * light_intersection[2];

    float diffuse[3] = {object[7]*light[5]*nl_dp, object[8]*light[6]*nl_dp, object[9]*light[7]*nl_dp};

    float light_ray[3] = {light_intersection[0]-ray_dir[0], light_intersection[1]-ray_dir[1], light_intersection[2]-ray_dir[2]};
    float norm = sqrt(light_ray[0]*light_ray[0] + light_ray[1]*light_ray[1] + light_ray[2]*light_ray[2]);

    float nlr_dp = normal_surface[0] * light_ray[0] +
                   normal_surface[1] * light_ray[1] +
                   normal_surface[2] * light_ray[2];

    nlr_dp = nlr_dp / norm;
    nlr_dp = pow(nlr_dp, 0.25*object[13]);

    float specular[3] = {object[10]*light[8]*nlr_dp, object[11]*light[9]*nlr_dp, object[12]*light[10]*nlr_dp};

    illumination[0] += *reflection *(ambient[0] + diffuse[0] + specular[0]);
    illumination[1] += *reflection *(ambient[1] + diffuse[1] + specular[1]);
    illumination[2] += *reflection *(ambient[2] + diffuse[2] + specular[2]);
}

void single_pixel(float* objects ,float* lights, float* camera, float* illumination, float* single_object, float* point){
    float ray_dir[3];
    float origin[3];

    origin[0] = camera[0];
    origin[1] = camera[1];
    origin[2] = camera[2];

    ray_direction(origin, point, ray_dir);
    float reflection = 1.0;

    for (int k=0; k < MAX_DEPTH; k++){
        float min_dist = __INT_MAX__;
        int n_object_idx = -1;

        nearest_intersection_object(objects, origin, ray_dir, &min_dist, &n_object_idx);
    
        if (n_object_idx == -1){
            return;
        }
        
        int is_shad;
        float normal[3], light_dir[3], shifted_point[3];

        for (int i=0; i<OBJ_LEN; i++){
            single_object[i] = objects[n_object_idx*OBJ_LEN + i];
        }

        for (int l=0; l<L_RANDOM; l++){
            float l_reflection = reflection;
            float x_rand = (float)rand()/RAND_MAX;
            float z_rand = (float)rand()/RAND_MAX;

            float light_pos[] = {lights[0] + x_rand*(lights[1] - lights[0]), lights[4], lights[2] + z_rand*(lights[3] - lights[2])};
            
            shadowed(&is_shad, normal, light_dir, shifted_point, &min_dist, origin, ray_dir, light_pos, objects, &n_object_idx);

            if (is_shad == 1)
                continue;

            
            color(normal, light_dir, ray_dir, single_object, lights, &reflection, illumination);
        }

        reflection *= single_object[14];
        origin[0] = shifted_point[0];
        origin[1] = shifted_point[1];
        origin[2] = shifted_point[2];

        reflected_direction(ray_dir, normal, ray_dir);

    }
}


int main(){
    float objects[] = {-0.2, 0, -1, 0.7, 0.1, 0, 0, 0.7, 0, 0, 1, 1, 1, 100, 0.5,
                       0.1, -0.3, 0, 0.1, 0.1, 0, 0.1, 0.7, 0, 0.7, 1, 1, 1, 100, 0.5,
                       -0.3, 0, 0, 0.15, 0, 0.1, 0, 0, 0.6, 0, 1, 1, 1, 100, 0.5,
                       -0.2, -9000, -1, 9000-0.7, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 100, 0
                      };
    float light[] = {XL_MIN, XL_MAX, ZL_MIN, ZL_MAX, YL, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    float camera[] = {0, 0, 1};
    float single_object[OBJ_LEN];
    float screen[] = {-1.0, 1.0, -(float)M/N, (float)M/N};
    float dx = (screen[1] - screen[0]) / N;
    float dy = (screen[3] - screen[2]) / M;
    int *image;
    image = (int*)malloc(N*M*3*sizeof(int));

    for (int i=0; i<N;i++){
        for (int j=0; j<M;j++){
            float position[] = {screen[0] + dx * i, screen[2] + dy * j, 0};
            float position1[] = {screen[0] + dx * i - dx/2, screen[2] + dy * j - dy/2, 0};
            float position2[] = {screen[0] + dx * i - dx/2, screen[2] + dy * j + dy/2, 0};
            float position3[] = {screen[0] + dx * i + dx/2, screen[2] + dy * j - dy/2, 0};
            float position4[] = {screen[0] + dx * i + dx/2, screen[2] + dy * j + dy/2, 0};
            float position5[] = {screen[0] + dx * i, screen[2] + dy * j + dy/2, 0};
            float position6[] = {screen[0] + dx * i, screen[2] + dy * j - dy/2, 0};
            float position7[] = {screen[0] + dx * i + dx/2, screen[2] + dy * j, 0};
            float position8[] = {screen[0] + dx * i - dx/2, screen[2] + dy * j, 0};
            float illumination[] = {0, 0, 0};
            single_pixel(objects, light, camera, illumination, single_object, position);
            single_pixel(objects, light, camera, illumination, single_object, position1);
            single_pixel(objects, light, camera, illumination, single_object, position2);
            single_pixel(objects, light, camera, illumination, single_object, position3);
            single_pixel(objects, light, camera, illumination, single_object, position4);
            single_pixel(objects, light, camera, illumination, single_object, position5);
            single_pixel(objects, light, camera, illumination, single_object, position6);
            single_pixel(objects, light, camera, illumination, single_object, position7);
            single_pixel(objects, light, camera, illumination, single_object, position8);
            image[3*M*i+3*j+0] = fmin(fmax(0, sqrt(illumination[0]/(9*(L_RANDOM)))), 1)*255;
            image[3*M*i+3*j+1] = fmin(fmax(0, sqrt(illumination[1]/(9*(L_RANDOM)))), 1)*255;
            image[3*M*i+3*j+2] = fmin(fmax(0, sqrt(illumination[2]/(9*(L_RANDOM)))), 1)*255;
        }
    }

    for (int i=AA; i<N-AA; i++){
        for (int j=AA; j<M-AA; j++){

            int sum_red = 0;
            int sum_green = 0;
            int sum_blue = 0;

            for (int k=-AA; k<=AA; k++){
                for (int l=-AA; l<=AA; l++){
                    sum_red += image[3*M*(i+k)+3*(j+l)+0];
                    sum_green += image[3*M*(i+k)+3*(j+l)+1];
                    sum_blue += image[3*M*(i+k)+3*(j+l)+2];
                }
            }

            image[3*M*i+3*j+0] = sum_red/((2*AA+1)*(2*AA+1));
            image[3*M*i+3*j+1] = sum_green/((2*AA+1)*(2*AA+1));
            image[3*M*i+3*j+2] = sum_blue/((2*AA+1)*(2*AA+1));
        }
    }

    printf("P3\n");
    printf("%d %d\n", N, M);
    printf("255 \n");
    for(int j=M-1; j>=0;j--){
        for(int i=0; i<N; i++){
            printf("%d %d %d\n", image[3*M*i+3*j+0], image[3*M*i+3*j+1], image[3*M*i+3*j+2]);
        }
    }


    free(image);
}
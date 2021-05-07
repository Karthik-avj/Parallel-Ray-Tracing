#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 500

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

void sphere_intersection(float* origin, float* ray_direction,
    float* center, float* radius, float* dist){


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

void nearest_intersection_object(float *objects, float *origin, int *object_n, float *ray_direction, int *objects_len, float *min_dist, int *object_idx){
	float distances[4];
	// float distances[*objects_len];
	for(int i=0; i<(*objects_len); i++){
		float center[] = {objects[i*(*object_n)+0], objects[i*(*object_n)+1], objects[i*(*object_n)+2]};
		sphere_intersection(origin, ray_direction, center, &objects[i*(*object_n)+3], &distances[i]);
	}

	for(int i=0; i<(*objects_len); i++){
		if((distances[i]<(*min_dist)) && (distances[i]>0)){
			*min_dist = distances[i];
			*object_idx = i;
		}
	}
}

void shadowed(int *is_shad, float *normal, float *light_dir, float *min_dist, float *origin, float *ray_dir,
              float *light_source, float *objects, int *object_idx, int *object_n, int *objects_len){
    // line 43
    float intersection_point[3];
    for (int i=0; i<3; i++){
        intersection_point[i] = *min_dist * ray_dir[i] + origin[i];
    }

    // line 44
    float object_center[] = {objects[*object_idx * *object_n], objects[*object_idx * *object_n + 1], objects[*object_idx * *object_n + 2]};
    ray_direction(object_center, intersection_point, normal);


    // line 45
    float shifted_point[3];
    for (int i=0; i<3; i++){
        shifted_point[i] = 0.000001 * normal[i] + intersection_point[i];
    }

    // line 46
    ray_direction(shifted_point, light_source, light_dir);

    // line 48
    float min_distance = __INT_MAX__;
    int useless = -1;
    nearest_intersection_object(objects, shifted_point, object_n, light_dir, objects_len, &min_distance, &useless);

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

void color(float* normal_surface, float* light_intersection, float* ray_dir, float* object, float* light, float* illumination){
    // float illumination[3] = {0, 0, 0};
    float ambient[3] = {object[4]*light[3], object[5]*light[4], object[6]*light[5]};

    float nl_dp = normal_surface[0] * light_intersection[0] +
                  normal_surface[1] * light_intersection[1] +
                  normal_surface[2] * light_intersection[2];

    float diffuse[3] = {object[7]*light[3]*nl_dp, object[8]*light[4]*nl_dp, object[9]*light[5]*nl_dp};

    float light_ray[3] = {light_intersection[0]-ray_dir[0], light_intersection[1]-ray_dir[1], light_intersection[2]-ray_dir[2]};
    float norm = sqrt(light_ray[0]*light_ray[0] + light_ray[1]*light_ray[1] + light_ray[2]*light_ray[2]);

    float nlr_dp = normal_surface[0] * light_ray[0] +
                   normal_surface[1] * light_ray[1] +
                   normal_surface[2] * light_ray[2];

    nlr_dp = nlr_dp / norm;
    nlr_dp = pow(nlr_dp, 0.25*object[13]);

    float specular[3] = {object[10]*light[6]*nlr_dp, object[11]*light[7]*nlr_dp, object[12]*light[8]*nlr_dp};

    illumination[0] = ambient[0] + diffuse[0] + specular[0];
    illumination[1] = ambient[1] + diffuse[1] + specular[1];
    illumination[2] = ambient[2] + diffuse[2] + specular[2];
}

void single_pixel(float* objects, int* objects_len ,float* lights, float* camera, float* illumination, float* single_object, float* point){
    float ray_dir[3];
    ray_direction(camera, point, ray_dir);

    float min_dist = __INT_MAX__;
    int n_object_idx = -1, object_n = 14;

    nearest_intersection_object(objects, camera, &object_n, ray_dir, objects_len, &min_dist, &n_object_idx);

    if (n_object_idx == -1)
        return;
    
    int is_shad;
    float normal[3], light_dir[3];

    float light_pos[] = {lights[0], lights[1], lights[2]};
    shadowed(&is_shad, normal, light_dir, &min_dist, camera, ray_dir, light_pos, objects, &n_object_idx, &object_n, objects_len);

    if (is_shad == 1)
        return;

    for (int i=0; i<object_n; i++){
        single_object[i] = objects[n_object_idx*object_n + i];
    }
    color(normal, light_dir, ray_dir, single_object, lights, illumination);
}


int main(){
    float objects[] = {-0.2, 0, -1, 0.7, 0.1, 0, 0, 0.7, 0, 0, 1, 1, 1, 100, 
                       0.1, -0.3, 0, 0.1, 0.1, 0, 0.1, 0.7, 0, 0.7, 1, 1, 1, 100,
                       -0.3, 0, 0, 0.15, 0, 0.1, 0, 0, 0.6, 0, 1, 1, 1, 100,
                       -0.2, -9000.7, -1, 9000, 0.1, 0.1, 0.1, 0.7, 0, 0, 1, 1, 1, 100};
    int objects_len = 4;
    float light[] = {5, 5, 5, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    float camera[] = {0, 0, 1};
    float single_object[14];
    float dx = 2.0 / N;
    int image[N][N][3];

    for (int i=0; i<N;i++){
        for (int j=0; j<N;j++){
            float position[] = {-1 + dx * i, -1 + dx * j, 0};
            float illumination[] = {0, 0, 0};
            single_pixel(objects, &objects_len, light, camera, illumination, single_object, position);
            image[i][j][0] = fmin(fmax(0, illumination[0]), 1)*255;
            image[i][j][1] = fmin(fmax(0, illumination[1]), 1)*255;
            image[i][j][2] = fmin(fmax(0, illumination[2]), 1)*255;

        }
    }

    printf("P3\n");
    printf("%d %d\n", N, N);
    printf("255 \n");
    for(int j=N-1; j>=0;j--){
        for(int i=0; i<N; i++){
            printf("%d %d %d\n", image[i][j][0], image[i][j][1], image[i][j][2]);
        }
    }
}
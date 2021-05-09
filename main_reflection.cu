#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 1920
#define M 1080
#define NUM_OBJ 4
#define OBJ_LEN 15
#define MAX(x,y) ( ((x) > (y)) ? x : y )
#define MIN(x,y) ( ((x) < (y)) ? x : y )

__device__ void ray_direction(float* origin, float* point, float* vector){
    float dr[3];
    dr[0] = point[0]-origin[0];
    dr[1] = point[1]-origin[1];
    dr[2] = point[2]-origin[2];

    float norm = sqrt(dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]);
    vector[0] = dr[0]/norm;
    vector[1] = dr[1]/norm;
    vector[2] = dr[2]/norm;
}

__device__ void reflected_direction(float* incoming, float* normal, float* reflected, int i, int j){
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

__device__ void sphere_intersection(float* origin, float* ray_direction, float* center, float* radius, float* dist){
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

__device__ void nearest_intersection_object(float *objects, float *origin, float *ray_direction, float *min_dist, int *object_idx){
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

__device__ void shadowed(int *is_shad, float *normal, float *light_dir, float *shifted_point, float *min_dist, float *origin, float *ray_dir,
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
        shifted_point[i] = 0.000001 * normal[i] + intersection_point[i];
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

__device__ void color(float* normal_surface, float* light_intersection, float* ray_dir, float* object, float* light, float* reflection, float* illumination){
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
    
    illumination[0] += *reflection *(ambient[0] + diffuse[0] + specular[0]);
    illumination[1] += *reflection *(ambient[1] + diffuse[1] + specular[1]);
    illumination[2] += *reflection *(ambient[2] + diffuse[2] + specular[2]);
}

__global__ void single_pixel(float* objects, float* lights, float* camera, float* screen, int* image, int* max_depth){
    int row = blockIdx.x*blockDim.x+threadIdx.x;
    int col = blockIdx.y*blockDim.y+threadIdx.y;

    if((row > N) || (col > M))
        return;
    
    float ray_dir[3];
    float origin[3];
    float dx = (screen[1] - screen[0]) / N;
    float dy = (screen[3] - screen[2]) / M;
    float point[] = {screen[0] + dx * row, screen[2] + dy * col, 0};
    float illumination[3] = {};
    float single_object[OBJ_LEN];

    origin[0] = camera[0];
    origin[1] = camera[1];
    origin[2] = camera[2];

    ray_direction(origin, point, ray_dir);
    float reflection = 1.0;

    for (int k=0; k < *max_depth; k++){
        float min_dist = __INT_MAX__;
        int n_object_idx = -1;

        nearest_intersection_object(objects, origin, ray_dir, &min_dist, &n_object_idx);

        if (n_object_idx == -1){
            if (k == 0){
                illumination[0] = 0.52734;
                illumination[1] = 0.80468;
                illumination[2] = 0.91796;
            }
            break;
        }
        
        int is_shad;
        float normal[3], light_dir[3], shifted_point[3];

        float light_pos[] = {lights[0], lights[1], lights[2]};
        shadowed(&is_shad, normal, light_dir, shifted_point, &min_dist, origin, ray_dir, light_pos, objects, &n_object_idx);
        
        if (is_shad == 1)
            break;

        for (int i=0; i<OBJ_LEN; i++){
            single_object[i] = objects[n_object_idx*OBJ_LEN + i];
        }

        color(normal, light_dir, ray_dir, single_object, lights, &reflection, illumination);

        reflection *= single_object[14];
        origin[0] = shifted_point[0];
        origin[1] = shifted_point[1];
        origin[2] = shifted_point[2];

        reflected_direction(ray_dir, normal, ray_dir, row, col);


    image[3 * M * row + 3*col + 0] = MIN(MAX(0, illumination[0]), 1)*255;
    image[3 * M * row + 3*col + 1] = MIN(MAX(0, illumination[1]), 1)*255;
    image[3 * M * row + 3*col + 2] = MIN(MAX(0, illumination[2]), 1)*255;
}


int main(){
    float objects[] = {-0.2, 0, -1, 0.7, 0.1, 0, 0, 0.7, 0, 0, 1, 1, 1, 100, 0.5,
                       0.1, -0.3, 0, 0.1, 0.1, 0, 0.1, 0.7, 0, 0.7, 1, 1, 1, 100, 0.5,
                       -0.3, 0, 0, 0.15, 0, 0.1, 0, 0, 0.6, 0, 1, 1, 1, 100, 0.5,
                       -0.2, -9000, -1, 9000-0.7, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1, 1, 1, 100, 0
                      };
    float light[] = {5, 5, 5, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    float camera[] = {0, 0, 1};
    float screen[] = {-1.0, 1.0, -(float)M/N, (float)M/N};
    int max_depth = 2;
    int *image;
    int size_image = N*M*3*sizeof(int);
    image = (int*)malloc(size_image);

    int size_objects = NUM_OBJ*OBJ_LEN*sizeof(float);

    float *dev_objects, *dev_light, *dev_camera, *dev_screen;
    int *dev_image, *dev_max_depth;

    cudaMalloc((void**) &dev_objects, size_objects);
    cudaMalloc((void**) &dev_light, 12*sizeof(float));
    cudaMalloc((void**) &dev_camera, 3*sizeof(float));
    cudaMalloc((void**) &dev_image, size_image);
    cudaMalloc((void**) &dev_max_depth, sizeof(int));
    cudaMalloc((void**) &dev_screen, 4*sizeof(float));

    cudaMemcpy(dev_objects, objects, size_objects, cudaMemcpyHostToDevice);
    cudaMemcpy(dev_light, light, 12*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_camera, camera, 3*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_max_depth, &max_depth, sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_screen, screen, 4*sizeof(float), cudaMemcpyHostToDevice);

    dim3 blockDim(16, 16);
    int k, l;
    if (N%blockDim.x == 0) {
        k = N/blockDim.x;
    } else {
        k = N/blockDim.x + 1;
    }
    if (M%blockDim.x == 0) {
        l = M/blockDim.x;
    } else {
        l = M/blockDim.x + 1;
    }
    dim3 gridDim(k, l);

    single_pixel<<<gridDim, blockDim>>>(dev_objects, dev_light, dev_camera, dev_screen, dev_image, dev_max_depth);

    cudaMemcpy(image, dev_image, size_image, cudaMemcpyDeviceToHost);

    printf("P3\n");
    printf("%d %d\n", N, M);
    printf("255 \n");
    for(int j=M-1; j>=0;j--){
        for(int i=0; i<N; i++){
            printf("%d %d %d\n", image[3*M*i+3*j+0], image[3*M*i+3*j+1], image[3*M*i+3*j+2]);
        }
    }

    cudaFree(dev_objects);
    cudaFree(dev_light);
    cudaFree(dev_camera);
    cudaFree(dev_image);
    cudaFree(dev_screen);
    free(image);
}
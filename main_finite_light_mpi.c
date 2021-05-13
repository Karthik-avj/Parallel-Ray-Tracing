#include "stdio.h"
#include "math.h"
#include "stdlib.h"
#include "mpi.h"

// Image Dimensions
#define N 1920
#define M 1080

// Finite Light Dimensions
#define XL_MIN 4.5
#define XL_MAX 5.5
#define ZL_MIN 4.5
#define ZL_MAX 5.5
#define YL 5.0

// Shadow rays to Light
#define L_RANDOM 16

// Number of reflections
#define MAX_DEPTH 3

// Code constants
#define NUM_OBJ 4
#define OBJ_LEN 15
#define ROOT 0


// distributing the workload evenly among different processors
void distribution_calc(int number_rows, int npes, int *send_rows,
                       int*displs, int each_row){
    int rows_each = number_rows/npes;
    int rows_left = number_rows%npes;
    int rows_sum = 0;

    for (int i=0; i<npes; i++){
        send_rows[i] = rows_each*each_row;
        if (rows_left > 0){
            send_rows[i] += each_row;
            rows_left -= 1;
        }

        displs[i] = rows_sum;
        rows_sum += send_rows[i];
    }

}


// finding the ray direction for a given origin and destination
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


// finding reflected ray for given normal and ray direction
void reflected_direction(float* incoming, float* normal, 
                         float* reflected){
    float dr[3];
    float dp = incoming[0]*normal[0] + incoming[1]*normal[1]
               + incoming[2]*normal[2];
    dr[0] = incoming[0] - 2*dp*normal[0];
    dr[1] = incoming[1] - 2*dp*normal[1];
    dr[2] = incoming[2] - 2*dp*normal[2];

    float norm = sqrt(dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]);
    reflected[0] = dr[0]/norm;
    reflected[1] = dr[1]/norm;
    reflected[2] = dr[2]/norm;
}


// distance between origin and intersection of sphere if any
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


// Nearest object which intersects with light ray
void nearest_intersection_object(float *objects, float *origin,
                                 float *ray_direction, float *min_dist, 
                                 int *object_idx){
	float distances[NUM_OBJ];
	for(int i=0; i<(NUM_OBJ); i++){
		float center[] = {objects[i*OBJ_LEN+0], 
                          objects[i*OBJ_LEN+1],
                          objects[i*OBJ_LEN+2]};
		sphere_intersection(origin, ray_direction, center, 
                            &objects[i*OBJ_LEN+3], &distances[i]);
	}

	for(int i=0; i<(NUM_OBJ); i++){
		if((distances[i]<(*min_dist)) && (distances[i]>0)){
			*min_dist = distances[i];
			*object_idx = i;
		}
	}
}


// Checks if the point is in shadow of other objects
void shadowed(int *is_shad, float *normal, float *light_dir, 
              float *shifted_point, float *min_dist, float *origin, 
              float *ray_dir, float *light_source, float *objects, 
              int *object_idx){

    float intersection_point[3];
    for (int i=0; i<3; i++){
        intersection_point[i] = *min_dist * ray_dir[i] + origin[i];
    }

    float object_center[] = {objects[*object_idx * OBJ_LEN],
                             objects[*object_idx * OBJ_LEN + 1],
                             objects[*object_idx * OBJ_LEN + 2]};
    ray_direction(object_center, intersection_point, normal);


    for (int i=0; i<3; i++){
        shifted_point[i] = 0.0001 * normal[i] + intersection_point[i];
    }

    ray_direction(shifted_point, light_source, light_dir);

    float min_distance = __INT_MAX__;
    int useless = -1;
    nearest_intersection_object(objects, shifted_point, light_dir,
                                &min_distance, &useless);

    float light_dist[3];
    light_dist[0] = light_source[0]-intersection_point[0];
    light_dist[1] = light_source[1]-intersection_point[1];
    light_dist[2] = light_source[2]-intersection_point[2];
    float intersection_to_light_dist = sqrt(
        (light_dist[0])*(light_dist[0]) + 
        (light_dist[1])*(light_dist[1]) + 
        (light_dist[2])*(light_dist[2]));

    if (min_distance < intersection_to_light_dist) {
        *is_shad = 1;
    }
    else{
        *is_shad = 0;
    }
}


// Determines the color using Blinn-Phong reflection model
void color(float* normal_surface, float* light_intersection,
           float* ray_dir, float* object, float* light, 
           float* reflection, float* illumination){

    float ambient[3] = {object[4]*light[5],
                        object[5]*light[6],
                        object[6]*light[7]};

    float nl_dp = normal_surface[0] * light_intersection[0] +
                  normal_surface[1] * light_intersection[1] +
                  normal_surface[2] * light_intersection[2];

    float diffuse[3] = {object[7]*light[5]*nl_dp,
                        object[8]*light[6]*nl_dp,
                        object[9]*light[7]*nl_dp};

    float light_ray[3] = {light_intersection[0]-ray_dir[0],
                          light_intersection[1]-ray_dir[1],
                          light_intersection[2]-ray_dir[2]};

    float norm = sqrt(
        light_ray[0]*light_ray[0] +
        light_ray[1]*light_ray[1] +
        light_ray[2]*light_ray[2]);

    float nlr_dp = normal_surface[0] * light_ray[0] +
                   normal_surface[1] * light_ray[1] +
                   normal_surface[2] * light_ray[2];

    nlr_dp = nlr_dp / norm;
    nlr_dp = pow(nlr_dp, 0.25*object[13]);

    float specular[3] = {object[10]*light[8]*nlr_dp,
                         object[11]*light[9]*nlr_dp,
                         object[12]*light[10]*nlr_dp};

    illumination[0] += *reflection *(ambient[0] + diffuse[0] + specular[0]);
    illumination[1] += *reflection *(ambient[1] + diffuse[1] + specular[1]);
    illumination[2] += *reflection *(ambient[2] + diffuse[2] + specular[2]);
}


// Calculation for a single pixel
void single_pixel(float* objects ,float* lights, float* camera, 
                  float* illumination, float* single_object, 
                  float* point){
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

        nearest_intersection_object(objects, origin, ray_dir,
                                    &min_dist, &n_object_idx);
    
        if (n_object_idx == -1){
            break;
        }
        
        int is_shad;
        float normal[3], light_dir[3], shifted_point[3];

        for (int i=0; i<OBJ_LEN; i++){
            single_object[i] = objects[n_object_idx*OBJ_LEN + i];
        }

        for (int l=0; l<L_RANDOM; l++){
            float x_rand = (float)rand()/RAND_MAX;
            float z_rand = (float)rand()/RAND_MAX;

            float light_pos[] = {lights[0] + x_rand*(lights[1] - lights[0]),
                                 lights[4],
                                 lights[2] + z_rand*(lights[3] - lights[2])};
            
            shadowed(&is_shad, normal, light_dir, shifted_point, 
                     &min_dist, origin, ray_dir, light_pos, 
                     objects, &n_object_idx);

            if (is_shad == 1)
                continue;

            
            color(normal, light_dir, ray_dir, single_object, 
                  lights, &reflection, illumination);
        }

        reflection *= single_object[14];
        origin[0] = shifted_point[0];
        origin[1] = shifted_point[1];
        origin[2] = shifted_point[2];

        reflected_direction(ray_dir, normal, ray_dir);

    }
}


// MAIN Function
int main(int argc, char** argv){
    int mype, npes;

    // objects being used in image
    float objects[] = {
        -0.2, 0, -1, 0.7, 0.1, 0, 0, 0.7, 0, 0, 1, 1, 1, 100, 0.5,
        0.1, -0.3, 0, 0.1, 0.1, 0, 0.1, 0.7, 0, 0.7, 1, 1, 1, 100, 0.5,
        -0.3, 0, 0, 0.15, 0, 0.1, 0, 0, 0.6, 0, 1, 1, 1, 100, 0.5,
        -0.2, -9000, -1, 8999.3, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 100, 0
        };
    
    // global variables
    float light[] = {XL_MIN, XL_MAX, ZL_MIN, ZL_MAX, YL, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    float camera[] = {0, 0, 1};

    // 9 points in each pixel
    int N_big = 2*N+1;
    int M_big = 2*M+1;
    int each_row = M_big*3;
    float screen[] = {-1.0, 1.0, -(float)M_big/N_big, (float)M_big/N_big};
    float dx = (screen[1] - screen[0]) / N_big;
    float dy = (screen[3] - screen[2]) / M_big;

    // image_local is a local copy used for calculation in each processor
    // image final is the final image
    int *image, *image_local, *image_final;
    image = (int*)malloc(N_big*M_big*3*sizeof(int));
    image_local = (int*)malloc(N_big*M_big*3*sizeof(int));
    image_final = (int*)malloc(N*M*3*sizeof(int));

    MPI_Status status;


    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &npes);
    MPI_Comm_rank(MPI_COMM_WORLD, &mype);

    int *send_rows = malloc(sizeof(int)*npes);
    int *displs = malloc(sizeof(int)*npes);

    // distributing the workload
    distribution_calc(N_big, npes, send_rows, displs, each_row);

    // calculation for each pixel locally
    for (int i=0; i<send_rows[mype]/each_row;i++){
        for (int j=0; j<M_big;j++){
            float single_object[OBJ_LEN];
            float position[] = {screen[0] + dx * (i + displs[mype]/each_row),
                                screen[2] + dy * j,
                                0};
            float illumination[] = {0, 0, 0};
            // calculation
            single_pixel(objects, light, camera, illumination, single_object, position);
            // storing the color in image (clipping between 0 and 255) (PPM format)
            image_local[3*M_big*i+3*j+0] = fmin(fmax(0, sqrt(illumination[0]/(L_RANDOM))), 1)*255;
            image_local[3*M_big*i+3*j+1] = fmin(fmax(0, sqrt(illumination[1]/(L_RANDOM))), 1)*255;
            image_local[3*M_big*i+3*j+2] = fmin(fmax(0, sqrt(illumination[2]/(L_RANDOM))), 1)*255;
        }
    }

    // gathering the data from all processors
    MPI_Gatherv(image_local, send_rows[mype], MPI_INT, image, 
                send_rows, displs, MPI_INT, ROOT, MPI_COMM_WORLD);

    
    // final image generation in the ROOT processor
    if (mype == ROOT){
        // averaging the results of 9 points in each pixel
        for (int i=1; i<N_big-1; i+=2){
            for (int j=1; j<M_big-1; j+=2){

                int sum_red = 0;
                int sum_green = 0;
                int sum_blue = 0;

                for (int k=-1; k<=1; k++){
                    for (int l=-1; l<=1; l++){
                        sum_red += image[3*M_big*(i+k)+3*(j+l)+0];
                        sum_green += image[3*M_big*(i+k)+3*(j+l)+1];
                        sum_blue += image[3*M_big*(i+k)+3*(j+l)+2];
                    }
                }

                image_final[3*M*(i-1)/2 + 3*(j-1)/2 + 0] = sum_red/9;
                image_final[3*M*(i-1)/2 + 3*(j-1)/2 + 1] = sum_green/9;
                image_final[3*M*(i-1)/2 + 3*(j-1)/2 + 2] = sum_blue/9;
                
                }
        }

        // printing in PPM format
        printf("P3\n");
        printf("%d %d\n", N, M);
        printf("255 \n");
        for(int j=M-1; j>=0; j--){
            for(int i=0; i<N; i++){
                printf("%d %d %d\n", image_final[3*M*i+3*j+0],
                                     image_final[3*M*i+3*j+1],
                                     image_final[3*M*i+3*j+2]);
            }
        }
    }

    free(send_rows);
    free(displs);
    free(image);
    free(image_local);
    free(image_final);

    MPI_Finalize();

    return 0;

}

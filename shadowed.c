#include<stdio.h>

void shadowed(int *is_shad, float *normal, float *light_dir, float *min_dist, float *origin, float *ray_dir, float *light_source, float *objects, int *object_idx, int *object_n){
    // line 43
    float *intersection_point[3];
    for (int i=0; i<3; i++){
        intersection_point[i] = *min_dist * ray_dir[i] + origin[i];
    }
    
    // line 44
    float object_center[] = {objects[*object_idx * *object_n], objects[*object_idx * *object_n + 1], objects[*object_idx * *object_n + 2]};
    ray_direction(normal, object_center, intersection_point);
    
    // line 45
    float *shifted_point[3];
    for (int i=0; i<3; i++){
        shifted_point[i] = 0.000001 * normal[i] + intersection_point[i];
    }

    // line 46
    ray_direction(light_dir, shifted_point, light_source);

    // line 48
    float min_distance;
    nearest_intersection_object(min_distance, objects, shifted_point, light_dir);

    // line 49
    float intersection_to_light_dist = (light_source[0]-intersection_point[0])*(light_source[0]-intersection_point[0]) + (light_source[1]-intersection_point[1])*(light_source[1]-intersection_point[1]) + (light_source[2]-intersection_point[2])*(light_source[2]-intersection_point[2]);

    // line 50
    if min_distance < intersection_to_light_dist{
        *is_shad = 1;
    }
    else{
        *is_shad = 0;
    }
}
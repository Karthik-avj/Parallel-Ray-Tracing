#include <stdio.h>

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
            *dist = min(val1, val2);
            return;
        }
        else{
            *dist = -1.0;
            return;
        }

    }



}


void ray_direction(float* origin, float* point, float* vector){
    float dr[3];
    dr[0] = point[0]-origin[0];
    dr[1] = point[1]-origin[1];
    dr[2] = point[2]-origin[2];

    float norm = dr[0]*dr[0] + dr[1]*dr[1]+dr[2]*dr[2];
    vector[0] = dr[0]/norm;
    vector[1] = dr[1]/norm;
    vector[2] = dr[2]/norm;
}
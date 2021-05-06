#include<math.h>

void color(float* illumination, float* normal_surface, float* light_intersection, float* ray_dir, float* object, float* light){
    // float illumination[3] = {0, 0, 0};
    float ambient[3] = {object[4]*light[3], object[5]*light[4], object[6]*light[5]};

    float nl_dp = normal_surface[0] * light_intersection[0] + \
                  normal_surface[1] * light_intersection[1] + \
                  normal_surface[2] * light_intersection[2];

    float diffuse[3] = {object[7]*light[3]*nl_dp, object[8]*light[4]*nl_dp, object[9]*light[5]*nl_dp};

    float light_ray[3] = {light_intersection[0]-ray_dir[0], light_intersection[1]-ray_dir[1], light_intersection[2]-ray_dir[2]};
    float norm = sqrt(light_ray[0]*light_ray[0] + light_ray[1]*light_ray[1] + light_ray[2]*light_ray[2]);

    float nlr_dp = normal_surface[0] * light_ray[0] + \
                   normal_surface[1] * light_ray[1] + \
                   normal_surface[2] * light_ray[2];

    nlr_dp = nlr_dp / norm;
    nlr_dp = pow(nlr_dp, 0.25*object[13]);

    float specular[3] = {object[10]*light[6]*nlr_dp, object[11]*light[7]*nlr_dp, object[12]*light[8]*nlr_dp};

    illumination[0] = ambient[0] + diffuse[0] + specular[0];
    illumination[1] = ambient[1] + diffuse[1] + specular[1];
    illumination[2] = ambient[2] + diffuse[2] + specular[2];
}
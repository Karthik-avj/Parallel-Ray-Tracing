void nearest_intersection_object(float *objects, float *origin, float *ray_direction, int *objects_len, float *min_dist, int *object_idx){
	float distances[objects_len];
	for(int i=0; i<objects_len; i++){
		float *center = {objects[i][0], objects[i][1], objects[i][2]};
		sphere_intersection(origin, ray_direction, center, objects[i][3], distances[i]);	
	}

	for(int i=0; i<objects_len; i++){
		if((distances[i]<min_dist) && (distances[i]>0)){
			min_dist = distances[i];
			object_idx = i;
		}
	}
}
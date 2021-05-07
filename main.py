import numpy as np
import matplotlib.pyplot as plt

def ray_direction(origin, point):
    vector = point - origin
    return vector / np.linalg.norm(vector)
  
def sphere_intersection(origin, ray_direction, center, radius):
    a = 1 # norm(ray_direction)^{2}
    b = 2 * np.dot(ray_direction, origin - center)
    c = (np.linalg.norm(origin - center) ** 2) - (radius * radius)

    disc = (b * b) - (4 * a * c)

    if disc <= 0:
        return -1
    
    else:
        val1 = (-b + np.sqrt(disc)) / (2 * a)
        val2 = (-b - np.sqrt(disc)) / (2 * a)

        if val1 > 0 and val2 > 0:
            return min(val1, val2)
        else:
            return -1


def nearest_intersection_object(objects, origin, ray_direction):
    distances = [sphere_intersection(origin, ray_direction, object.center, object.radius) for object in objects]

    object_idx = None
    min_dist = np.inf

    for idx, dist in enumerate(distances):
        if dist < min_dist and dist > 0:
            min_dist = dist
            object_idx = idx

    return min_dist, object_idx


def shadowed(min_dist, origin, ray_dir, light_source, objects, object_idx):
    intersection_point = min_dist * ray_dir + origin
    normal_surface = ray_direction(objects[object_idx].center, intersection_point)
    shifted_point = 1e-5*normal_surface + intersection_point
    light_intersection = ray_direction(shifted_point, light_source)
    
    min_distance, _ = nearest_intersection_object(objects, shifted_point, light_intersection)
    intersection_to_light_dist = np.linalg.norm(light_source - intersection_point)
    is_shadowed = min_distance < intersection_to_light_dist
    
    return is_shadowed, normal_surface, light_intersection
    

def color(normal_surface, light_intersection, ray_dir, single_object, light):
	
    illumination = np.zeros(3)
    illumination += single_object.ambient * light.ambient
    illumination += single_object.diffuse * light.diffuse * np.dot(normal_surface, light_intersection)
    illumination += single_object.specular * light.specular * (np.dot(normal_surface, (light_intersection - ray_dir)/np.linalg.norm(light_intersection - ray_dir)))**(0.25*single_object.shininess)
    return illumination
    # return np.clip(illumination, 0, 1)



class Sphere:
  def __init__(self, c, r, ambient, diffuse, specular, shininess):
    self.center = np.array(c)
    self.radius = r
    self.ambient = np.array(ambient) 
    self.diffuse = np.array(diffuse)
    self.specular = np.array(specular)
    self.shininess = shininess

class Light:
  def __init__(self, position, ambient, diffuse, specular):
    self.position = np.array(position)
    self.ambient = np.array(ambient) 
    self.diffuse = np.array(diffuse)
    self.specular = np.array(specular)

nr = 20
nc = 20

pixels = np.zeros((nr, nc))

camera = np.array([0, 0, 1])
ratio = float(nc) / nr
# on the xy plane
screen = (-1, 1 / ratio, 1, -1 / ratio) # left, top, right, bottom
image = np.zeros([nr, nc, 3])

objects = np.array([
  					Sphere([-0.2, 0, -1], 0.7, [0.1, 0, 0], [0.7, 0, 0], [1, 1, 1], 100), 
                    Sphere([0.1, -0.3, 0], 0.1, [0.1, 0, 0.1], [0.7, 0, 0.7], [1, 1, 1], 100), 
                    Sphere([-0.3, 0, 0], 0.15, [0, 0.1, 0], [0, 0.6, 0], [1, 1, 1], 100),
  					Sphere([-0.2, -9000.7, -1], 9000, [0.1, 0.1, 0.1], [0.7, 0, 0], [1, 1, 1], 100), 
					])

lights = np.array([
                   Light([5, 5, 5], [1, 1, 1], [1, 1, 1], [1, 1, 1]),
                #    Light([-5, 5, 5], [1, 1, 1], [1, 1, 1], [1, 1, 1]),
                ])

light = Light([5, 5, 5], [1, 1, 1], [1, 1, 1], [1, 1, 1])


for i, y in enumerate(np.linspace(screen[1], screen[3], nr)):
    for j, x in enumerate(np.linspace(screen[0], screen[2], nc)):
        pixel = np.array([x, y, 0])
        origin = camera
        ray_dir = ray_direction(origin, pixel)
        
        min_dist, nearest_object_idx = nearest_intersection_object(objects, origin, ray_dir)
        if nearest_object_idx is None:
        	continue
        
        for light in lights:
            is_shadowed, normal_surface, light_intersection = shadowed(min_dist, origin, ray_dir, light.position, objects, nearest_object_idx)
            print(is_shadowed)
            if is_shadowed:
                continue
            
            image[i, j] += color(normal_surface, light_intersection, ray_dir, objects[nearest_object_idx], light)
        
    # print("%d/%d" % (i + 1, nr))

image = np.clip(image, 0, 1)
plt.imsave('image.png', image)
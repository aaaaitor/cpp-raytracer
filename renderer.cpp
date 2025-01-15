#define _USE_MATH_DEFINES
#define STB_IMAGE_WRITE_IMPLEMENTATION
#define STB_IMAGE_IMPLEMENTATION
#include <limits>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include "geometry.h"
#include "stb_image_write.h"
#include "stb_image.h"

int envmap_width, envmap_height;
std::vector<Vec3f> envmap;

//Define a light source
struct Light {
    Light(const Vec3f& p, const float i) : position(p), intensity(i) {}
    Vec3f position;
    float intensity;
};

//Define sphere material
struct Material {
    Material(const float r, const Vec4f& a, const Vec3f& color, const float spec) : refractive_index(r), albedo(a), diffuse_color(color), specular_exponent(spec) {}
    Material() : refractive_index(1), albedo(1, 0, 0, 0), diffuse_color(), specular_exponent() {}
    float refractive_index;
    Vec4f albedo;
    Vec3f diffuse_color;
    float specular_exponent;
};

//Define a sphere based on a point on the world (center) and a radius
struct Sphere {
    Vec3f       center;
    float       radius;
    Material    material;

    Sphere(const Vec3f& c, const float r, const Material& m) : center(c), radius(r), material(m) {}

    //Calculate if a ray casted form the camera position intersects with the sphere
    bool ray_intersect(const Vec3f& origin, const Vec3f& direction, float& t0) const {
        Vec3f L = center - origin;
        float tca = L * direction;
        float d2 = L * L - tca * tca;

        if (d2 > radius * radius) return false;

        float thc = sqrtf(radius * radius - d2);

        t0 = tca - thc;
        float t1 = tca + thc;

        if (t0 < 0) t0 = t1;
        if (t0 < 0) return false;

        return true;
    }
};

//Calculate the light reflection
Vec3f reflect(const Vec3f& I, const Vec3f& N) {
    return I - N * 2.f * (I * N);
}

//Calculate the refraction of an object using Snell's law
Vec3f refract(const Vec3f& I, const Vec3f& N, const float eta_t, const float eta_i=1.f) {
    float cos   = -std::max(-1.f, std::min(1.f, I * N));
    
    //If a ray is inside an object, swap air and media values
    if (cos < 0) return refract(I, -N, eta_i, eta_t);

    float  eta = eta_i / eta_t;
    float k = 1 - eta * eta * (1 - cos * cos);

    //If k<0 the object has total reflection, meaning it should not refract
    return k < 0 ? Vec3f(1, 0, 0) : I * eta + N * (eta * cos - sqrtf(k));
}

//Calculate the intersections in the scene
bool scene_intersect(const Vec3f& origin, const Vec3f& direction, const std::vector<Sphere>& spheres, Vec3f& hit, Vec3f& N, Material& material) {
    float spheres_dist = std::numeric_limits<float>::max();
    for (size_t i = 0; i < spheres.size(); i++) {
        float dist_i;
        if (spheres[i].ray_intersect(origin, direction, dist_i) && dist_i < spheres_dist) {
            spheres_dist = dist_i;
            hit = origin + direction * dist_i;
            N = (hit - spheres[i].center).normalize();
            material = spheres[i].material;
        }
    }

    float checkerboard_dist = std::numeric_limits<float>::max();
    if (fabs(direction.y) > 1e-3) {
        float d = -(origin.y + 4) / direction.y;
        Vec3f pt = origin + direction * d;
        if (d > 0 && fabs(pt.x) < 10 && pt.z<-10 && pt.z>-30 && d < spheres_dist) {
            checkerboard_dist = d;
            hit = pt;
            N = Vec3f(0, 1, 0);
            material.diffuse_color = (int(.5 * hit.x + 1000) + int(.5 * hit.z)) & 1 ? Vec3f(.3, .3, .3) : Vec3f(.3, .2, .1);
        }
    }
    return std::min(spheres_dist, checkerboard_dist) < 1000;
}

Vec3f get_environment_color(const Vec3f& ray_direction) {

    float theta = acos(ray_direction.y);
    float phi = atan2(ray_direction.z, ray_direction.x);

    float u = (phi + M_PI) / (2 * M_PI);
    float v = (theta) / M_PI;

    int x = int(u * envmap_width);
    int y = int(v * envmap_height);

    x = std::min(std::max(x, 0), envmap_width - 1);
    y = std::min(std::max(y, 0), envmap_height - 1);

    return envmap[x + y * envmap_width];
}

//-------------
// RAY CASTING
//-------------

//Cast a ray using a origin point and a direction and return a color for the pixel
// depending if it intersects or not

const Vec3f background_color    = Vec3f(0.2, 0.7, 0.8);         //Background color for when a ray doesnt intersect with an object
const size_t max_depth          = 4;                            //Max depth used to calculate reflections

Vec3f cast_ray_ambient(const Vec3f& orig, const Vec3f& dir, const std::vector<Sphere>& spheres) {
    Vec3f point, N;
    Material material;
    if (!scene_intersect(orig, dir, spheres, point, N, material)) {
        return background_color;
    }
    return material.diffuse_color;
}

Vec3f cast_ray_diffuse(const Vec3f& orig, const Vec3f& dir, const std::vector<Sphere>& spheres, const std::vector<Light>& lights) {
    Vec3f point, N;
    Material material;

    if (!scene_intersect(orig, dir, spheres, point, N, material)) {
        return background_color;
    }

    float diffuse_light_intensity = 0;
    for (size_t i = 0; i < lights.size(); i++) {
        Vec3f light_direction = (lights[i].position - point).normalize();
        diffuse_light_intensity += lights[i].intensity * std::max(0.f, light_direction * N);
    }
    return material.diffuse_color * diffuse_light_intensity;
}

Vec3f cast_ray_specular(const Vec3f& orig, const Vec3f& dir, const std::vector<Sphere>& spheres, const std::vector<Light>& lights) {
    Vec3f point, N;
    Material material;

    if (!scene_intersect(orig, dir, spheres, point, N, material)) {
        return Vec3f(0.0, 0.0, 0.0);
    }

    float diffuse_light_intensity = 0, specular_light_intensity = 0;
    for (size_t i = 0; i < lights.size(); i++) {
        Vec3f light_direction = (lights[i].position - point).normalize();
        diffuse_light_intensity += lights[i].intensity * std::max(0.f, light_direction * N);
        specular_light_intensity += powf(std::max(0.f, -reflect(-light_direction, N) * dir), material.specular_exponent) * lights[i].intensity;
    }
    return Vec3f(0,0,0) * diffuse_light_intensity * material.albedo[0] + Vec3f(1., 1., 1.) * specular_light_intensity * material.albedo[1];
}

Vec3f cast_ray_final(const Vec3f& orig, const Vec3f& dir, const std::vector<Sphere>& spheres, const std::vector<Light>& lights, size_t depth = 0) {
    Vec3f point, N;
    Material material;

    if (depth > max_depth || !scene_intersect(orig, dir, spheres, point, N, material)) {
        return get_environment_color(dir);
    }

    Vec3f reflect_dir = reflect(dir, N).normalize();

    //Avoid occluding the object with itself
    Vec3f reflect_orig = reflect_dir * N < 0 ? point - N * 1e-3 : point + N * 1e-3;
    Vec3f reflect_color = cast_ray_final(reflect_orig, reflect_dir, spheres, lights, depth + 1);

    Vec3f refract_dir = refract(dir, N, material.refractive_index).normalize();
    Vec3f refract_orig = refract_dir * N < 0 ? point - N * 1e-3 : point + N * 1e-3;
    Vec3f refract_color = cast_ray_final(refract_orig, refract_dir, spheres, lights, depth + 1);

    float diffuse_light_intensity = 0, specular_light_intensity = 0;
    for (size_t i = 0; i < lights.size(); i++) {

        //Calculate light direction and distance
        Vec3f light_direction = (lights[i].position - point).normalize();
        float light_distance = (lights[i].position - point).norm();

        //Check if the ray between the current point and the light source intersect any other object,
        // if so, the point should not take into account the checked light source
        Vec3f shadow_orig = light_direction * N < 0 ? point - N * 1e-3 : point + N * 1e-3;
        Vec3f shadow_pt, shadow_N;
        Material tmpmaterial;
        if (scene_intersect(shadow_orig, light_direction, spheres, shadow_pt, shadow_N, tmpmaterial) && (shadow_pt - shadow_orig).norm() < light_distance)
            continue;

        diffuse_light_intensity += lights[i].intensity * std::max(0.f, light_direction * N);
        specular_light_intensity += powf(std::max(0.f, reflect(light_direction, N) * dir), material.specular_exponent) * lights[i].intensity;
    }
    return material.diffuse_color * diffuse_light_intensity * material.albedo[0] + Vec3f(1., 1., 1.) * specular_light_intensity * material.albedo[1] + reflect_color*material.albedo[2] + refract_color*material.albedo[3];
}

//-----------------
// IMAGE RENDERING
//-----------------

const int width     = 1024;                                                     //Image width
const int height    = 768;                                                      //Image height
const int channels  = 3;                                                        //Image channels (R,G,B)
const int quality   = 100;                                                      //Image quality
const int fov       = M_PI / 2.;                                                //viewpoint field of view

////Render base image
//void render_base() {
//    std::vector<Vec3f> framebuffer(width * height);                             //Pixel buffer
//    unsigned char* data = (unsigned char*)malloc(width * height * channels);    //Image data allocation
//
//    //Image generation
//    #pragma omp parallel for
//    for (size_t j = 0; j < height; j++) {
//        for (size_t i = 0; i < width; i++) {
//            framebuffer[i + j * width] = Vec3f(j / float(height), i / float(width), 0);
//        }
//    }
//
//    for (size_t i = 0; i < height * width; ++i) {
//        Vec3f& c = framebuffer[i];
//        float max = std::max(c[0], std::max(c[1], c[2]));
//        if (max > 1) c = c * (1. / max);
//        for (size_t j = 0; j < 3; j++) {
//            data[i * 3 + j] = (unsigned char)(255 * std::max(0.f, std::min(1.f, framebuffer[i][j])));
//        }
//    }
//
//    //Save image
//    stbi_write_jpg("./output/0_base.jpg", width, height, channels, data, quality);
//    free(data);
//}

void render_ambient(const std::vector<Sphere>& spheres) {
    std::vector<Vec3f> framebuffer(width * height);                             //Pixel buffer
    unsigned char* data = (unsigned char*)malloc(width * height * channels);    //Image data allocation

#pragma omp parallel for
    for (size_t j = 0; j < height; j++) {
        for (size_t i = 0; i < width; i++) {
            float x = (2 * (i + 0.5) / (float)width - 1) * tan(fov / 2.) * width / (float)height;
            float y = -(2 * (j + 0.5) / (float)height - 1) * tan(fov / 2.);
            Vec3f dir = Vec3f(x, y, -1).normalize();
            framebuffer[i + j * width] = cast_ray_ambient(Vec3f(0, 0, 0), dir, spheres);
        }
    }

    for (size_t i = 0; i < height * width; ++i) {
        Vec3f& c = framebuffer[i];
        float max = std::max(c[0], std::max(c[1], c[2]));
        if (max > 1) c = c * (1. / max);
        for (size_t j = 0; j < 3; j++) {
            assert(data != nullptr);
            data[i * 3 + j] = (unsigned char)(255 * std::max(0.f, std::min(1.f, framebuffer[i][j])));
        }
    }

    //Save image
    stbi_write_jpg("./output/1_ambient.jpg", width, height, channels, data, quality);
    free(data);
}

void render_diffuse(const std::vector<Sphere>& spheres, const std::vector<Light>& lights) {
    std::vector<Vec3f> framebuffer(width * height);                             //Pixel buffer
    unsigned char* data = (unsigned char*)malloc(width * height * channels);    //Image data allocation

    //Image generation
    #pragma omp parallel for
    for (size_t j = 0; j < height; j++) {
        for (size_t i = 0; i < width; i++) {
            framebuffer[i + j * width] = Vec3f(j / float(height), i / float(width), 0);
            float x = (2 * (i + 0.5) / (float)width - 1) * tan(fov / 2.) * width / (float)height;
            float y = -(2 * (j + 0.5) / (float)height - 1) * tan(fov / 2.);
            Vec3f dir = Vec3f(x, y, -1).normalize();
            framebuffer[i + j * width] = cast_ray_diffuse(Vec3f(0, 0, 0), dir, spheres, lights);
        }
    }

    for (size_t i = 0; i < height * width; ++i) {
        Vec3f& c = framebuffer[i];
        float max = std::max(c[0], std::max(c[1], c[2]));
        if (max > 1) c = c * (1. / max);
        for (size_t j = 0; j < 3; j++) {
            assert(data != nullptr);
            data[i * 3 + j] = (unsigned char)(255 * std::max(0.f, std::min(1.f, framebuffer[i][j])));
        }
    }

    //Save image
    stbi_write_jpg("./output/2_diffuse.jpg", width, height, channels, data, quality);
    free(data);
}

void render_specular(const std::vector<Sphere>& spheres, const std::vector<Light>& lights) {
    std::vector<Vec3f> framebuffer(width * height);                             //Pixel buffer
    unsigned char* data = (unsigned char*)malloc(width * height * channels);    //Image data allocation

    //Image generation
    #pragma omp parallel for
    for (size_t j = 0; j < height; j++) {
        for (size_t i = 0; i < width; i++) {

            //Calculate the camera direction
            float x = (2 * (i + .5) / (float)width - 1) * tan(fov / 2.) * width / (float)height;
            float y = -(2 * (j + .5) / (float)height - 1) * tan(fov / 2.);
            Vec3f direction = Vec3f(x, y, -1).normalize();

            framebuffer[i + j * width] = cast_ray_specular(Vec3f(0, 0, 0), direction, spheres, lights);
        }
    }

    //Save image
    for (size_t i = 0; i < height * width; ++i) {
        Vec3f& c = framebuffer[i];
        float max = std::max(c[0], std::max(c[1], c[2]));
        if (max > 1) c = c * (1. / max);
        for (size_t j = 0; j < 3; j++) {
            assert(data != nullptr);
            data[i * 3 + j] = (unsigned char)(255 * std::max(0.f, std::min(1.f, framebuffer[i][j])));
        }
    }

    //Save image
    stbi_write_jpg("./output/3_specular.jpg", width, height, channels, data, quality);
    free(data);
}

void render_final(const std::vector<Sphere>& spheres, const std::vector<Light>& lights) {
    std::vector<Vec3f> framebuffer(width * height);     //Pixel buffer

    //Image generation
    #pragma omp parallel for
    for (size_t j = 0; j < height; j++) {
        for (size_t i = 0; i < width; i++) {

            //Calculate the camera direction and flip it
            float direction_x =         (i + 0.5) - width / 2.;
            float direction_y =       -(j + 0.5) + height / 2.;
            float direction_z = -height / (2. * tan(fov / 2.));

            framebuffer[i + j * width] = cast_ray_final(Vec3f(0, 0, 0), Vec3f(direction_x, direction_y, direction_z).normalize(), spheres, lights);
        }
    }

    //Enviroment pixmap
    std::vector<unsigned char> pixmap(width*height*3);
    //Save image
    for (size_t i = 0; i < height * width; ++i) {
        Vec3f& c = framebuffer[i];
        float max = std::max(c[0], std::max(c[1], c[2]));
        if (max > 1) c = c * (1. / max);
        for (size_t j = 0; j < 3; j++) {
            pixmap[i * 3 + j] = (unsigned char)(255 * std::max(0.f, std::min(1.f, framebuffer[i][j])));
        }
    }

    //Save image
    stbi_write_jpg("./output/4_final.jpg", width, height, channels, pixmap.data(), quality);
}

int main() {
    std::vector<Sphere> spheres;
    std::vector<Light>  lights;

    std::cout << "[INFO/Resources/Images] Loading enviroment image." << std::endl;

    //Load and read enviroment .jpg
    int n = -1;
    unsigned char* pixmap = stbi_load("./resources/envmap.jpg", &envmap_width, &envmap_height, &n, 0);
    if (!pixmap || 3 != n) {
        std::cerr << "[Error] Can not load the environment map." << std::endl;
        return -1;
    }
    envmap = std::vector<Vec3f>(envmap_width * envmap_height);
    for (int j = envmap_height - 1; j >= 0; j--) {
        for (int i = 0; i < envmap_width; i++) {
            envmap[i + j * envmap_width] = Vec3f(pixmap[(i + j * envmap_width) * 3 + 0], pixmap[(i + j * envmap_width) * 3 + 1], pixmap[(i + j * envmap_width) * 3 + 2]) * (1 / 255.);
        }
    }
    stbi_image_free(pixmap);

    std::cout << "[INFO/Resources/Images] Enviroment image loaded." << std::endl;

    std::cout << "[INFO/Resources/Materials] Defining materials." << std::endl;

    //Define materials
    Material emerald_green(1.0, Vec4f(0.6,  0.6, 0.6, 0.0),  Vec3f(0.1, 0.7, 0.1),    70.);
    Material    red_rubber(1.0, Vec4f(0.9,  0.1, 0.0, 0.0),  Vec3f(0.3, 0.1, 0.1),    10.);
    Material  marble_white(1.0, Vec4f(0.6,  0.6, 0.6, 0.0),  Vec3f(0.4, 0.4, 0.3),    70.);
    Material        mirror(1.0, Vec4f(0.0, 10.0, 0.8, 0.0),  Vec3f(1.0, 1.0, 1.0),  1425.);
    Material         glass(1.5, Vec4f(0.0,  0.5, 0.1, 0.8),  Vec3f(0.6, 0.7, 0.8),   125.);

    std::cout << "[INFO/Resources/Materials] Materials defined." << std::endl;

    std::cout << "[INFO/Scene/Meshes] Creating meshes." << std::endl;

    //Add spheres
    spheres.push_back(Sphere(Vec3f(  -3,    0, -16), 2,  emerald_green));
    spheres.push_back(Sphere(Vec3f(-1.0, -1.5, -12), 2,          glass));
    spheres.push_back(Sphere(Vec3f( 1.5, -0.5, -18), 3,     red_rubber));
    spheres.push_back(Sphere(Vec3f(   7,    5, -18), 4,         mirror));
    spheres.push_back(Sphere(Vec3f(  -7, -1.5, -12), 2,   marble_white));

    std::cout << "[INFO/Scene/Meshes] Meshes created." << std::endl;

    std::cout << "[INFO/Scene/Lighting] Creating lights." << std::endl;

    //Add light sources
    lights.push_back(Light(Vec3f(-20, 20, 20), 1.5));
    lights.push_back(Light(Vec3f(30, 50, -25), 1.8));
    lights.push_back(Light( Vec3f(30, 20, 30), 1.7));

    std::cout << "[INFO/Scene/Lighting] Lights created." << std::endl;

    ////Render base image
    //render_base();
    
    std::cout << "[INFO/Rendering] Starting rendering queue." << std::endl;

    std::cout << "[INFO/Rendering/Ambient] Starting ambient image render." << std::endl;
    //Render ambient image
    render_ambient(spheres);
    std::cout << "[INFO/Rendering/Ambient] Finished ambient image render." << std::endl;

    std::cout << "[INFO/Rendering/Diffuse] Starting diffuse image render." << std::endl;
    //Render diffuse image
    render_diffuse(spheres, lights);
    std::cout << "[INFO/Rendering/Diffuse] Finished diffuse image render." << std::endl;

    std::cout << "[INFO/Rendering/Specular] Starting specular image render." << std::endl;
    //Render specular
    render_specular(spheres, lights);
    std::cout << "[INFO/Rendering/Specular] Finished specular image render." << std::endl;

    std::cout << "[INFO/Rendering/Specular] Starting final image render." << std::endl;
    //Render final
    render_final(spheres, lights);
    std::cout << "[INFO/Rendering/Specular] Finished final image render." << std::endl;

    std::cout << "[INFO/Rendering] Rendering queue finished." << std::endl;

    return 0;
}
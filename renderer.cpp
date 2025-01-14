#define _USE_MATH_DEFINES
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <limits>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include "geometry.h"
#include "stb_image_write.h"

//Define a light source
struct Light {
    Light(const Vec3f& p, const float& i) : position(p), intensity(i) {}
    Vec3f position;
    float intensity;
};

//Define sphere material
struct Material {
    Material(const Vec2f& a, const Vec3f& color, const float& spec) : albedo(a), diffuse_color(color), specular_exponent(spec) {}
    Material() : albedo(1, 0), diffuse_color(), specular_exponent() {}
    Vec2f albedo;
    Vec3f diffuse_color;
    float specular_exponent;
};

//Define a sphere based on a point on the world (center) and a radius
struct Sphere {
    Vec3f       center;
    float       radius;
    Material    material;

    Sphere(const Vec3f& c, const float& r, const Material& m) : center(c), radius(r), material(m) {}

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
    return spheres_dist < 1000;
}

//-------------
// RAY CASTING
//-------------

//Cast a ray using a origin point and a direction and return a color for the pixel
// depending if it intersects or not

Vec3f cast_ray_ambient(const Vec3f& orig, const Vec3f& dir, const std::vector<Sphere>& spheres) {
    Vec3f point, N;
    Material material;
    if (!scene_intersect(orig, dir, spheres, point, N, material)) {
        return Vec3f(0.2, 0.7, 0.8);
    }
    return material.diffuse_color;
}

Vec3f cast_ray_diffuse(const Vec3f& orig, const Vec3f& dir, const std::vector<Sphere>& spheres, const std::vector<Light>& lights) {
    Vec3f point, N;
    Material material;

    if (!scene_intersect(orig, dir, spheres, point, N, material)) {
        return Vec3f(0.2, 0.7, 0.8);
    }

    float diffuse_light_intensity = 0;
    for (size_t i = 0; i < lights.size(); i++) {
        Vec3f light_dir = (lights[i].position - point).normalize();
        diffuse_light_intensity += lights[i].intensity * std::max(0.f, light_dir * N);
    }
    return material.diffuse_color * diffuse_light_intensity;
}

Vec3f cast_ray_specular(const Vec3f& orig, const Vec3f& dir, const std::vector<Sphere>& spheres, const std::vector<Light>& lights) {
    Vec3f point, N;
    Material material;

    if (!scene_intersect(orig, dir, spheres, point, N, material)) {
        return Vec3f(0.2, 0.7, 0.8);
    }

    float diffuse_light_intensity = 0, specular_light_intensity = 0;
    for (size_t i = 0; i < lights.size(); i++) {
        Vec3f light_dir = (lights[i].position - point).normalize();
        diffuse_light_intensity += lights[i].intensity * std::max(0.f, light_dir * N);
        specular_light_intensity += powf(std::max(0.f, -reflect(-light_dir, N) * dir), material.specular_exponent) * lights[i].intensity;
    }
    return Vec3f(0.2, 0.7, 0.8) * diffuse_light_intensity * material.albedo[0] + Vec3f(1., 1., 1.) * specular_light_intensity * material.albedo[1];
}

Vec3f cast_ray_final(const Vec3f& orig, const Vec3f& dir, const std::vector<Sphere>& spheres, const std::vector<Light>& lights) {
    Vec3f point, N;
    Material material;

    if (!scene_intersect(orig, dir, spheres, point, N, material)) {
        return Vec3f(0.2, 0.7, 0.8);
    }

    float diffuse_light_intensity = 0, specular_light_intensity = 0;
    for (size_t i = 0; i < lights.size(); i++) {
        Vec3f light_dir = (lights[i].position - point).normalize();
        diffuse_light_intensity += lights[i].intensity * std::max(0.f, light_dir * N);
        specular_light_intensity += powf(std::max(0.f, -reflect(-light_dir, N) * dir), material.specular_exponent) * lights[i].intensity;
    }
    return material.diffuse_color * diffuse_light_intensity * material.albedo[0] + Vec3f(1., 1., 1.) * specular_light_intensity * material.albedo[1];
}

//-----------------
// IMAGE RENDERING
//-----------------

const int width     = 1024;                                                 //Image width
const int height    = 768;                                                  //Image height
const int channels  = 3;                                                    //Image channels (R,G,B)
const int quiality  = 100;                                                  //Image quality
const int fov       = M_PI / 2.;                                            //viewpoint field of view

//Render base image
void render_base() {
    std::vector<Vec3f> framebuffer(width * height);                             //Pixel buffer
    unsigned char* data = (unsigned char*)malloc(width * height * channels);    //Image data allocation

    //Image generation
    #pragma omp parallel for
    for (size_t j = 0; j < height; j++) {
        for (size_t i = 0; i < width; i++) {
            framebuffer[i + j * width] = Vec3f(j / float(height), i / float(width), 0);
        }
    }

    for (size_t i = 0; i < height * width; ++i) {
        Vec3f& c = framebuffer[i];
        float max = std::max(c[0], std::max(c[1], c[2]));
        if (max > 1) c = c * (1. / max);
        for (size_t j = 0; j < 3; j++) {
            data[i * 3 + j] = (unsigned char)(255 * std::max(0.f, std::min(1.f, framebuffer[i][j])));
        }
    }

    stbi_write_jpg("./output/0_base.jpg", width, height, channels, data, quiality);
    free(data);
}

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
            data[i * 3 + j] = (unsigned char)(255 * std::max(0.f, std::min(1.f, framebuffer[i][j])));
        }
    }

    stbi_write_jpg("./output/1_ambient.jpg", width, height, channels, data, quiality);
    free(data);
}

void render_diffuse(const std::vector<Sphere>& spheres) {
    const int width = 1024;                             //Image width
    const int height = 768;                             //Image height
    const int fov = M_PI / 2.;                          //Viewpoint field of view
    std::vector<Vec3f> framebuffer(width * height);     //Pixel buffer

    //Image generation
    #pragma omp parallel for
    for (size_t j = 0; j < height; j++) {
        for (size_t i = 0; i < width; i++) {
            framebuffer[i + j * width] = Vec3f(j / float(height), i / float(width), 0);
            float x = (2 * (i + 0.5) / (float)width - 1) * tan(fov / 2.) * width / (float)height;
            float y = -(2 * (j + 0.5) / (float)height - 1) * tan(fov / 2.);
            Vec3f dir = Vec3f(x, y, -1).normalize();
            framebuffer[i + j * width] = cast_ray_diffuse(Vec3f(0, 0, 0), dir, spheres, std::vector<Light> {});
        }
    }

    //Save the image
    std::ofstream ofs;
    ofs.open("./output/2_diffuse.ppm", std::ofstream::out, std::ofstream::binary);
    ofs << "P6\n" << width << " " << height << "\n255\n";
    for (size_t i = 0; i < height * width; ++i) {
        for (size_t j = 0; j < 3; j++) {
            ofs << (char)(255 * std::max(0.f, std::min(1.f, framebuffer[i][j])));
        }
    }
    ofs.close();
}

void render_specular(const std::vector<Sphere>& spheres, const std::vector<Light>& lights) {
    const int width = 1024;                             //Image width
    const int height = 768;                             //Image height
    const int fov = M_PI / 2.;                          //Viewpoint field of view

    std::vector<Vec3f> framebuffer(width * height);     //Pixel buffer

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

    //Save the image
    std::ofstream ofs ("./output/3_specular.ppm", std::ofstream::out, std::ofstream::binary);
    ofs << "P6\n" << width << " " << height << "\n255\n";
    for (size_t i = 0; i < height * width; ++i) {
        Vec3f& c = framebuffer[i];
        float max = std::max(c[0], std::max(c[1], c[2]));
        if (max > 1) c = c * (1. / max);
        for (size_t j = 0; j < 3; j++) {
            ofs << (char)(255 * std::max(0.f, std::min(1.f, framebuffer[i][j])));
        }
    }
    ofs.close();
}

void render_final(const std::vector<Sphere>& spheres, const std::vector<Light>& lights) {
    const int width     = 1024;                         //Image width
    const int height    = 768;                          //Image height
    const int fov       = M_PI / 2.;                    //Viewpoint field of view

    std::vector<Vec3f> framebuffer(width * height);     //Pixel buffer

    //Image generation
    #pragma omp parallel for
    for (size_t j = 0; j < height; j++) {
        for (size_t i = 0; i < width; i++) {

            //Calculate the camera direction
            float x = (2 * (i + .5) / (float)width - 1) * tan(fov / 2.) * width / (float)height;
            float y = -(2 * (j + .5) / (float)height - 1) * tan(fov / 2.);
            Vec3f direction = Vec3f(x, y, -1).normalize();

            framebuffer[i + j * width] = cast_ray_final(Vec3f(0, 0, 0), direction, spheres, lights);
        }
    }

    //Save the image
    std::ofstream ofs;
    ofs.open("./output/4_final.ppm", std::ofstream::out, std::ofstream::binary);
    ofs << "P6\n" << width << " " << height << "\n255\n";
    for (size_t i = 0; i < height * width; ++i) {
        Vec3f& c = framebuffer[i];
        float max = std::max(c[0], std::max(c[1], c[2]));
        if (max > 1) c = c * (1. / max);
        for (size_t j = 0; j < 3; j++) {
            ofs << (char)(255 * std::max(0.f, std::min(1.f, framebuffer[i][j])));
        }
    }
    ofs.close();
}

int main() {
    std::vector<Sphere> spheres;
    std::vector<Light>  lights;

    //Define materials
    Material      ivory(Vec2f(0.6, 0.3), Vec3f(0.4, 0.4, 0.3), 50.);
    Material red_rubber(Vec2f(0.9, 0.1), Vec3f(0.3, 0.1, 0.1), 10.);

    //Add spheres
    spheres.push_back(Sphere(Vec3f(-3,      0,   -16), 2, ivory));
    spheres.push_back(Sphere(Vec3f(-1.0, -1.5,   -12), 2, red_rubber));
    spheres.push_back(Sphere(Vec3f(1.5,  -0.5,   -18), 3, red_rubber));
    spheres.push_back(Sphere(Vec3f(7,       5,   -18), 4, ivory));

    //Add light sources
    lights.push_back(Light(Vec3f(-20, 20, 20), 1.5));
    lights.push_back(Light(Vec3f(30, 50, -25), 1.8));
    lights.push_back(Light(Vec3f(30, 20, 30), 1.7));


    //Render base image
    render_base();
    //Render ambient image
    render_ambient(spheres);
    //Render diffuse image
    render_diffuse(spheres);
    //Render specular
    render_specular(spheres, lights);
    //Render final
    render_final(spheres, lights);

    return 0;
}
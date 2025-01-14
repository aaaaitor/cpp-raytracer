#define _USE_MATH_DEFINES
#include <limits>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include "geometry.h"

//Define sphere material
struct Material {
    Material(const Vec3f& color) : diffuse_color(color) {}
    Material() : diffuse_color() {}
    Vec3f diffuse_color;
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

//Cast a ray using a origin point and a direction and return a color for the pixel
// depending if it intersects or not
Vec3f cast_ray(const Vec3f& origin, const Vec3f& direction, const std::vector<Sphere>& spheres) {
    Vec3f point, N;
    Material material;
    if (!scene_intersect(origin, direction, spheres, point, N, material)) {
        return Vec3f(.2, .7, .8);       //Not intersecting color
    }

    return material.diffuse_color;      //Intersecting color
}

void render(const std::vector<Sphere>& spheres) {
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

            framebuffer[i + j * width] = cast_ray(Vec3f(0, 0, 0), direction, spheres);
        }
    }

    //Save the image
    std::ofstream ofs;
    ofs.open("./out.ppm", std::ofstream::out, std::ofstream::binary);
    ofs << "P6\n" << width << " " << height << "\n255\n";
    for (size_t i = 0; i < height * width; ++i) {
        for (size_t j = 0; j < 3; j++) {
            ofs << (char)(255 * std::max(0.f, std::min(1.f, framebuffer[i][j])));
        }
    }
    ofs.close();
}

int main() {
    Material      ivory(Vec3f(0.4, 0.4, 0.3));
    Material red_rubber(Vec3f(0.3, 0.1, 0.1));
    std::vector<Sphere> spheres;

    spheres.push_back(Sphere(Vec3f(-3,      0,   -16), 2, ivory));
    spheres.push_back(Sphere(Vec3f(-1.0, -1.5,   -12), 2, red_rubber));
    spheres.push_back(Sphere(Vec3f(1.5,  -0.5,   -18), 3, red_rubber));
    spheres.push_back(Sphere(Vec3f(7,       5,   -18), 4, ivory));

    render(spheres);

    return 0;
}
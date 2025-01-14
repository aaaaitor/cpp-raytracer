#define _USE_MATH_DEFINES
#include <limits>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include "geometry.h"

//Define a sphere based on a point on the world (center) and a radius
struct Sphere {
    Vec3f   center;
    float   radius;

    Sphere(const Vec3f& c, const float& r) : center(c), radius(r) {}

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

//Cast a ray using a origin point and a direction and return a color for the pixel
// depending if it intersects or not
Vec3f cast_ray(const Vec3f& origin, Vec3f& direction, const Sphere& sphere) {
    float sphere_distance = std::numeric_limits<float>::max();

    if (!sphere.ray_intersect(origin, direction, sphere_distance)) {
        return Vec3f(.2, .7, .8);       //Not intersecting color
    }

    return Vec3f(.4, .4, .3);           //Intersecting color
}

void render(const Sphere& sphere) {
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

            framebuffer[i + j * width] = cast_ray(Vec3f(0, 0, 0), direction, sphere);
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
    render();
    return 0;
}
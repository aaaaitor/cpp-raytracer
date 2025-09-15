# C++ Ray Tracer

A small but powerful ray tracing program implemented in C++ that demonstrates various lighting and reflection effects. This project implements a ray tracer from scratch, capable of rendering scenes with multiple lighting effects including ambient, diffuse, specular reflection, shadows, and environment mapping.

## Features

- **Multiple Lighting Effects:**
  - Ambient lighting
  - Diffuse reflection
  - Specular reflection
  - Phong reflection model
  - Shadow casting
  - Reflections and refractions
  - Environment mapping

- **Image Processing:**
  - JPEG input/output support using STB library
  - High-quality rendering output
  - Configurable scene parameters

## Project Structure

```
cpp-raytracer/
├── include/
│   ├── geometry.h          # Vector mathematics and geometry operations
│   ├── stb_image.h        # Image loading library
│   └── stb_image_write.h  # Image writing library
├── src/
│   └── renderer.cpp       # Main rendering engine
├── resources/
│   └── envmap.jpg        # Environment map texture
└── output/               # Generated output images
```

## Prerequisites

- C++ compiler with C++11 support
- Visual Studio 2022 or later
- Windows SDK 10.0 or later

## Building the Project

1. Open the solution file `cpp-raytracer.sln` in Visual Studio
2. Select your preferred configuration (Debug/Release) and platform (x64/x86)
3. Build the solution (F7 or Ctrl+Shift+B)

## Running the Program

When running the program, it will generate five different images in the `output/` directory, each demonstrating different lighting and reflection effects:

*1. Ambient image:*

![1_ambient](https://github.com/user-attachments/assets/44b27dea-6c19-4738-a7e9-487a0f9d1e4f)

*2. Diffuse image:*

![2_diffuse](https://github.com/user-attachments/assets/32b33884-e3c6-43fc-9626-9a36ee34eb8f)

*3. Specular image:*

![3_specular](https://github.com/user-attachments/assets/61deeb99-30a9-43cf-b23f-11bbfffbfc68)

*4. Phong reflection - (Ambient + Diffuse + Specular):*

![4_phong](https://github.com/user-attachments/assets/f4cab853-6b70-44b4-83d6-9fc27fb7750d)

*5. Final image with environment - (Phong reflection + Shadows + Reflections + Refractions):*

![5_final](https://github.com/user-attachments/assets/a690fdce-3467-4674-bd20-e43c359b6372)

## Technical Implementation

The ray tracer implements several key computer graphics concepts:

- Ray-sphere intersection calculations
- Vector mathematics for light reflection
- Phong illumination model
- Environment mapping using spherical coordinates
- Fresnel equations for realistic reflection/refraction

## Dependencies

- [STB](https://github.com/nothings/stb) - Image loading/writing library
  - `stb_image.h` for image loading
  - `stb_image_write.h` for image writing

## License

This project is released under the MIT License. See the LICENSE file for details.
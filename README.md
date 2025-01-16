# Small c++ program to render an image using ray-tracing.

The program includes the main file `(renderer.cpp)`, a small header file to create and modify vectors `(geometry.h)` and to header files from the stb library to read and write JPEG images.


When running the `'renderer.cpp'` file the program should output 5 images in the `'./output/'` directory. The images should look as follows (if nothing is changed):

*1. Ambient image:*

![1_ambient](https://github.com/user-attachments/assets/44b27dea-6c19-4738-a7e9-487a0f9d1e4f)

*2. Diffuse image:*

![2_diffuse](https://github.com/user-attachments/assets/32b33884-e3c6-43fc-9626-9a36ee34eb8f)

*3. Specular image:*

![3_specular](https://github.com/user-attachments/assets/61deeb99-30a9-43cf-b23f-11bbfffbfc68)

*4. Phong reflection - (Ambient + Diffuse + Specular):*

![4_phong](https://github.com/user-attachments/assets/f4cab853-6b70-44b4-83d6-9fc27fb7750d)

*5. Final image with enviroment - (Phong reflection + Shadows + Reflections + Refractions):*

![5_final](https://github.com/user-attachments/assets/a690fdce-3467-4674-bd20-e43c359b6372)

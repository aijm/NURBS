# NURBS
A really simple NURBS library based on libigl. [libigl](https://github.com/libigl/libigl/) is a simple C++ geometry processing library, and I use it mainly for displaying the NURBS curve and surface.

Existing functions:
1. compute points on curve or surface with parameters 
2. display NURBS curve and surface by libigl

to be continued:
1. readIO and writeIO 
2. knot insertion and knot removal
3. degree elevation and degree reduction
4. NURBS curve and surface fitting
5. T-Spline (it's hard ...)


This project's cmake structure is from [libigl/libigl-example-project](https://github.com/libigl/libigl-example-project), which is a blank project example showing how to use libigl and cmake. 

## See the tutorial first

Then build, run and understand the [libigl
tutorial](http://libigl.github.io/libigl/tutorial/).

## Compile
After downloading and installing libigl, you need to add the Path of libigl to this project by modifying `FindLIBIGL.cmake`
```
find_path(LIBIGL_INCLUDE_DIR igl/readOBJ.h
    HINTS
        ENV LIBIGL
        ENV LIBIGLROOT
        ENV LIBIGL_ROOT
        ENV LIBIGL_DIR
    PATHS
        ${CMAKE_SOURCE_DIR}/../..
        ${CMAKE_SOURCE_DIR}/..
        ${CMAKE_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/libigl
        ${CMAKE_SOURCE_DIR}/../libigl
        ${CMAKE_SOURCE_DIR}/../../libigl
        D:/Program\ Files/libigl # you can add the path like this
        /usr
        /usr/local
        /usr/local/igl/libigl
    PATH_SUFFIXES include
)

```
Compile this project using the standard cmake routine:

    mkdir build
    cd build
    cmake ..
    make

This should find and build the dependencies and create a `NURBS_bin` binary.

## Run

From within the `build` directory just issue:

    ./NURBS_bin

A glfw app should launch displaying a Torus surface.

## examples and results
1. Cylindrical surface
![cylinder.PNG](https://github.com/aijm/NURBS/blob/master/examples/cylinder.PNG)
2. Torus surface
![torus.PNG](https://github.com/aijm/NURBS/blob/master/examples/torus.PNG)



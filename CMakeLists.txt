cmake_minimum_required(VERSION 3.7)
project(DiplomSoldatovEigen)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES main.cpp square_triangle_grid/square_triangle_grid.cpp square_triangle_grid/square_triangle_grid.h constants/constants.h structures/structures.h fem_triangle_element/fem_triangle_element.cpp fem_triangle_element/fem_triangle_element.h Rocket/air_parameters.h variables.h debug.h)
add_executable(DiplomSoldatovEigen ${SOURCE_FILES})
include_directories(C:\\Users\\aksol\\CLionProjects\\DiplomSoldatovEigen\\Eigen)
cmake_minimum_required(VERSION 3.6)
project(3d_bem)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror -Wall -Wextra -pedantic")

set(SOURCE_FILES
        il/io/numpy.cpp
        main.cpp
        mesh_file_io.cpp
        matrix_asm.cpp
        tensor_oper.cpp
        ele_base.cpp
        elast_ker_int.cpp
        h_potential.cpp)

add_executable(3d_bem ${SOURCE_FILES})
target_include_directories(3d_bem PRIVATE ${CMAKE_SOURCE_DIR})
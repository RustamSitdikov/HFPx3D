cmake_minimum_required(VERSION 3.6)
project(HFPx3D_static)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror -Wall -Wextra -pedantic")
set(SOURCE_FILES
        il/io/numpy.cpp
        main.cpp
        mesh_file_io.cpp
        matrix_assembly.cpp
        tensor_utilities.cpp
        element_utilities.cpp
        elasticity_kernel_integration.cpp
        h_potential.cpp)
add_executable(HFPx3D_static ${SOURCE_FILES})
target_include_directories(HFPx3D_static PRIVATE ${CMAKE_SOURCE_DIR})
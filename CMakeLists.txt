cmake_minimum_required(VERSION 3.6)
project(3d_bem)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES
    il/io/numpy.cpp
    main.cpp)

add_executable(3d_bem ${SOURCE_FILES})
target_include_directories(3d_bem PRIVATE ${CMAKE_SOURCE_DIR})
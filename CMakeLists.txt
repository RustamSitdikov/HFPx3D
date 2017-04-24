cmake_minimum_required(VERSION 3.6)
project(HFPx3D_VC)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror -Wall -Wextra -pedantic -DIL_OPENBLAS -DIL_BLAS")
set(SOURCE_FILES_DIR ${PROJECT_SOURCE_DIR}/src)
#set(SOURCE_FILES_DIR ${PROJECT_SOURCE_DIR})
set(SOURCE_FILES
        il/io/numpy.cpp
        main.cpp
        ${SOURCE_FILES_DIR}/mesh_file_io.cpp
        #${SOURCE_FILES_DIR}/c_f_iteration.cpp
        ${SOURCE_FILES_DIR}/system_assembly.cpp
        ${SOURCE_FILES_DIR}/mesh_utilities.cpp
        ${SOURCE_FILES_DIR}/tensor_utilities.cpp
        ${SOURCE_FILES_DIR}/element_utilities.cpp
        ${SOURCE_FILES_DIR}/elasticity_kernel_integration.cpp
        ${SOURCE_FILES_DIR}/h_potential.cpp
        )
add_executable(HFPx3D_VC ${SOURCE_FILES})
target_include_directories(HFPx3D_VC PRIVATE ${CMAKE_SOURCE_DIR})
target_link_libraries(HFPx3D_VC openblas)
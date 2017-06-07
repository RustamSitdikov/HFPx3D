cmake_minimum_required(VERSION 3.6)
project(HFPx3D_VC)
set(CMAKE_CXX_STANDARD 11)

IF(CMAKE_BUILD_TYPE MATCHES DEBUG)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror -Wall -Wextra -pedantic")
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror -Wall -pedantic")
ENDIF(CMAKE_BUILD_TYPE MATCHES DEBUG)

set(SOURCE_FILES_DIR ${PROJECT_SOURCE_DIR}/src)

set(SOURCE_FILES
        il/io/numpy.cpp
        main.cpp
        ${SOURCE_FILES_DIR}/mesh_file_io.cpp
        ${SOURCE_FILES_DIR}/mesh_file_io.h
        #${SOURCE_FILES_DIR}/c_f_iteration.cpp
        ${SOURCE_FILES_DIR}/cohesion_friction.h
        ${SOURCE_FILES_DIR}/system_assembly.cpp
        ${SOURCE_FILES_DIR}/mesh_utilities.cpp
        ${SOURCE_FILES_DIR}/tensor_utilities.cpp
        ${SOURCE_FILES_DIR}/element_utilities.cpp
        ${SOURCE_FILES_DIR}/elasticity_kernel_integration.cpp
        ${SOURCE_FILES_DIR}/h_potential.cpp
        )

add_executable(HFPx3D_VC ${SOURCE_FILES})

####### IL_MKL case

link_directories($ENV{MKLROOT}/lib $ENV{MKLROOT}/../compiler/lib)
target_compile_definitions(HFPx3D_VC PRIVATE IL_MKL=1 IL_BLAS=1)
target_include_directories(HFPx3D_VC PRIVATE ${CMAKE_SOURCE_DIR} $ENV{MKLROOT}/include)
target_link_libraries(HFPx3D_VC mkl_intel_lp64 mkl_core mkl_intel_thread dl iomp5 pthread m)

add_custom_command(TARGET HFPx3D_VC
        POST_BUILD COMMAND
        /usr/bin/install_name_tool -add_rpath /opt/intel/compilers_and_libraries_2017.0.102/mac/mkl/lib
        $<TARGET_FILE:HFPx3D_VC>)

#[[
####### IL_OPENBLAS case

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DIL_OPENBLAS -DIL_BLAS")
target_include_directories(HFPx3D_VC PRIVATE ${CMAKE_SOURCE_DIR})
target_link_libraries(HFPx3D_VC openblas)
]]

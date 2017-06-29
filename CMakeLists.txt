cmake_minimum_required(VERSION 3.6)
project(HFPx3D)
set(CMAKE_CXX_STANDARD 11)

# some options related to IL
set(IL_OPENMP 0)
set(IL_OPENBLAS 1)
set(IL_MKL 0)

# treat warnings as errors
if(CMAKE_BUILD_TYPE MATCHES DEBUG)
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror -Wall -Wextra -pedantic")
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror -Wall -pedantic")
    if (UNIX)
        if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror -Wall -Wextra -pedantic")
        elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror -Wall -Wextra -pedantic")
        elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror -Wall -pedantic")
        endif()
    elseif (WIN32)
        if (MINGW)
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror -Wall -Wextra -pedantic")
        endif()
    endif()
endif()

# For OpenMP
if (IL_OPENMP)
    if (UNIX)
        if (UNIX)
            if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
                set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O3 -xHost -ansi-alias -DNDEBUG")
            endif()
        endif()
        if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
            #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
            #link_directories(/opt/clang-3.9.0/lib)
        elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
        elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qopenmp")
        endif()
    elseif (WIN32)
        if (MINGW)
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
        endif()
    endif()
endif()

set(SOURCE_FILES_DIR ${PROJECT_SOURCE_DIR}/src)

set(SOURCE_FILES
        main.cpp
        il/io/numpy/numpy.cpp
        il/io/toml/toml.cpp
        ${SOURCE_FILES_DIR}/mesh_file_io.cpp
        ${SOURCE_FILES_DIR}/mesh_file_io.h
        #${SOURCE_FILES_DIR}/c_f_iteration.cpp
        ${SOURCE_FILES_DIR}/cohesion_friction.h
        ${SOURCE_FILES_DIR}/system_assembly.cpp
        ${SOURCE_FILES_DIR}/surface_mesh_utilities.cpp
        ${SOURCE_FILES_DIR}/tensor_utilities.cpp
        ${SOURCE_FILES_DIR}/element_utilities.cpp
        ${SOURCE_FILES_DIR}/elasticity_kernel_integration.cpp
        ${SOURCE_FILES_DIR}/h_potential.cpp
        )

add_executable(HFPx3D ${SOURCE_FILES})

####### IL_MKL case
## ENV{MKLROOT} is smth. like /opt/intel/compilers_and_libraries_2017.0.102/mac/mkl
if (IL_MKL)
    #link_directories($ENV{MKLROOT}/lib)
    link_directories($ENV{MKLROOT}/lib $ENV{MKLROOT}/../compilers/lib)
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64 -DIL_MKL -DIL_BLAS")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DIL_MKL -DIL_BLAS")

    if (UNIX)
        if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64 -DIL_MKL -DIL_BLAS")
            set(CMAKE_MKL_LIBRARIES mkl_intel_lp64 mkl_core mkl_intel_thread iomp5 dl pthread m)
        elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64 -DIL_MKL -DIL_BLAS")
            set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--no-as-needed")
            set(CMAKE_MKL_LIBRARIES mkl_intel_lp64 mkl_core mkl_intel_thread iomp5 dl pthread m)
            link_directories($ENV{MKLROOT}/lib/intel64_lin)
        elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mkl=parallel -DIL_MKL -DIL_BLAS")
        endif()
    endif()

    target_compile_definitions(HFPx3D PRIVATE IL_MKL=1 IL_BLAS=1)
    target_include_directories(HFPx3D PRIVATE ${CMAKE_SOURCE_DIR} $ENV{MKLROOT}/include)

####### IL_OPENBLAS case
elseif (IL_OPENBLAS)
    if (UNIX)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DIL_OPENBLAS -DIL_BLAS")
        set(CMAKE_OPENBLAS_LIBRARIES "openblas")
    elseif (WIN32)
        #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DIL_OPENBLAS -DIL_BLAS")
        #set(CMAKE_OPENBLAS_LIBRARIES "openblas")
        if (MINGW)
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DIL_OPENBLAS -DIL_BLAS")
            set(CMAKE_OPENBLAS_LIBRARIES "openblas")
        endif()
    elseif (APPLE)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DIL_MKL -DIL_OPENBLAS -DIL_BLAS")
        set(CMAKE_OPENBLAS_LIBRARIES "openblas")
    endif()

    target_include_directories(HFPx3D PRIVATE ${CMAKE_SOURCE_DIR})
endif()

target_link_libraries(HFPx3D ${CMAKE_MKL_LIBRARIES} ${CMAKE_OPENBLAS_LIBRARIES})

if (IL_MKL)
    if (APPLE)
        add_custom_command(TARGET HFPx3D
                POST_BUILD COMMAND
                /usr/bin/install_name_tool -add_rpath $ENV{MKLROOT}/lib
                $<TARGET_FILE:HFPx3D>)
    endif()
endif()

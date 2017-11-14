cmake_minimum_required(VERSION 3.6)

project(HFPx3DUnitTest CXX)

project(HFPx3D)

set(CMAKE_CXX_STANDARD 11)

# some options related to IL
set(IL_OPENMP 1)
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
        if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
            #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
            #link_directories(/opt/clang-3.9.0/lib)
        elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
        elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qopenmp")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O3 -xHost -fp-model fast=2")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -L/opt/intel/vtune_amplifier_xe_2016.4.0.470476/lib64")
            #set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O3 -xHost -ansi-alias -DNDEBUG")
        endif()
    elseif (WIN32)
        if (MINGW)
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
        endif()
    endif()
endif()

set(SOURCE_FILES_DIR ${PROJECT_SOURCE_DIR}/src)

set(SOURCE_FILES
        il/io/numpy/numpy.cpp
        il/io/toml/toml.cpp
        ${SOURCE_FILES_DIR}/IO/config_file_io.cpp
        ${SOURCE_FILES_DIR}/IO/mesh_file_io.cpp
        ${SOURCE_FILES_DIR}/IO/data_file_io.h
        ${SOURCE_FILES_DIR}/IO/data_file_io.cpp
        ${SOURCE_FILES_DIR}/Core/surface_mesh_utilities.cpp
        ${SOURCE_FILES_DIR}/Core/tensor_utilities.cpp
        ${SOURCE_FILES_DIR}/Core/element_utilities.cpp
        ${SOURCE_FILES_DIR}/Elasticity/h_potential.cpp
        ${SOURCE_FILES_DIR}/Elasticity/elasticity_kernel_integration.cpp
        ${SOURCE_FILES_DIR}/Solvers/system_assembly.cpp
        ${SOURCE_FILES_DIR}/Solvers/test_radial_crack_static.cpp
        ${SOURCE_FILES_DIR}/Solvers/test_radial_crack_static.h
        ${SOURCE_FILES_DIR}/Development/cohesion_friction.cpp
        ${SOURCE_FILES_DIR}/Development/c_f_iteration.cpp
        )

add_executable(HFPx3D ${SOURCE_FILES} main.cpp)

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
            target_include_directories(HFPx3D PRIVATE ${CMAKE_SOURCE_DIR} $ENV{MKLROOT}/include)
        elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64 -DIL_MKL -DIL_BLAS")
            set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--no-as-needed")
            set(CMAKE_MKL_LIBRARIES mkl_intel_lp64 mkl_core mkl_intel_thread iomp5 dl pthread m)
            link_directories($ENV{MKLROOT}/lib/intel64_lin)
            target_include_directories(HFPx3D PRIVATE ${CMAKE_SOURCE_DIR} $ENV{MKLROOT}/include)
        elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mkl=parallel -DIL_MKL -DIL_BLAS")
            set(CMAKE_MKL_LIBRARIES mkl_intel_lp64 mkl_core mkl_intel_thread dl iomp5 pthread m ittnotify)
            target_include_directories(HFPx3D PRIVATE ${CMAKE_SOURCE_DIR} $ENV{MKLROOT}/include /opt/intel/vtune_amplifier_xe_2016.4.0.470476/include)
        endif()
    endif()

    target_compile_definitions(HFPx3D PRIVATE IL_MKL=1 IL_BLAS=1 IL_OPENMP=1)
    target_compile_definitions(HFPx3DUnitTest PRIVATE IL_MKL=1 IL_BLAS=1 IL_OPENMP=1)

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

################################################################################
# For Unit Tests  HFPx3DUnitTest
################################################################################

set(SR_FOLDER "${CMAKE_SOURCE_DIR}/src")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DIL_FOLDER='\"${SR_FOLDER}\"'")

set(UNIT_TEST_FILES
        ${SOURCE_FILES_DIR}/Solvers/_test/solver_test.cpp
        ${SOURCE_FILES_DIR}/Solvers/_test/matrix_assembly_test.cpp
        ${SOURCE_FILES_DIR}/Elasticity/_test/kernel_integration_test.cpp
        #${SOURCE_FILES_DIR}/Core/_test/dof_handle_test.cpp
        ${SOURCE_FILES_DIR}/Core/_test/element_shape_functions_test.cpp
        gtest/src/gtest-all.cc)

add_executable(HFPx3DUnitTest ${SOURCE_FILES} ${UNIT_TEST_FILES} test.cpp)

target_link_libraries(HFPx3DUnitTest ${CMAKE_MKL_LIBRARIES} ${CMAKE_OPENBLAS_LIBRARIES})
# target_link_libraries(HFPx3DUnitTest mkl_intel_lp64 mkl_core
# mkl_intel_thread dl iomp5 pthread m)

# should be the same -> make sure the MKLROOT is defined under linux
if(APPLE)
    target_include_directories(HFPx3DUnitTest PUBLIC ${CMAKE_SOURCE_DIR}
            ${CMAKE_SOURCE_DIR}/gtest $ENV{MKLROOT}/include)
elseif(UNIX)
    target_include_directories(HFPx3DUnitTest PUBLIC ${CMAKE_SOURCE_DIR}
            ${CMAKE_SOURCE_DIR}/gtest  /opt/intel/mkl/include)
else() #Windows
    target_include_directories(HFPx3DUnitTest PUBLIC ${CMAKE_SOURCE_DIR}
            ${CMAKE_SOURCE_DIR}/gtest)

endif()

# For unit tests: The precondition of our functions are checked with assert
# macros that terminate the program in debug mode. In order to test those macros
# in our unit tests, the flag IL_UNIT_TEST turns them into exceptions that can
# be caught and checked.
target_compile_options(HFPx3DUnitTest PRIVATE "-DIL_UNIT_TEST")

if(APPLE)
    # below for some under mac os sierra ;(
    add_custom_command(TARGET HFPx3DUnitTest
            POST_BUILD COMMAND
            /usr/bin/install_name_tool -add_rpath /opt/intel/compilers_and_libraries_2017.0.102/mac/mkl/lib
            $<TARGET_FILE:HFPx2DUnitTest>)
endif()

//
// This file is part of HFPx3D.
//
// Created by D. Nikolski on 11/11/2017.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//

#include <gtest/gtest/gtest.h>
#include <iostream>
#include <complex>
#include <il/norm.h>
#include <il/Array2D.h>
#include <il/StaticArray.h>
#include <il/StaticArray2D.h>
#include "src/Core/element_utilities.h"
#include "src/Elasticity/elasticity_kernel_integration.h"
#include "src/IO/data_file_io.h"

// UNIT TEST FOR ELASTIC KERNEL INTEGRATION
// Compares the results of elastic kernel integration
// with the reference matrices for a single triangular element
// and 3 different observation points: one on the plane, two off the plane
// The reference matrices are calculated with the original MATLAB code
// by D. Nikolski & S. Mogilevskaya, 2015-16

TEST(kernel_integration_1, t1) {
    // tolerance
    double tol = 1e-06;

    // elastic properties
    double G = 1.0, nu = 0.35;

    // connectivity
    il::StaticArray2D<il::int_t, 3, 1> mesh_conn;
    for (il::int_t k = 0; k < 3; ++k) {
        mesh_conn(k, 0) = k;
    }

    // vertex coordinates
    il::StaticArray2D<double, 3, 3> nodes_crd;
    nodes_crd(0, 0) = 0.0;
    nodes_crd(1, 0) = 0.1;
    nodes_crd(2, 0) = 0.0;

    nodes_crd(0, 1) = 1.8;
    nodes_crd(1, 1) = 0.0;
    nodes_crd(2, 1) = 0.0;

    nodes_crd(0, 2) = 1.2;
    nodes_crd(1, 2) = 1.8;
    nodes_crd(2, 2) = 0.0;

    // observation point #1 (on the plane)
    il::StaticArray<double, 3> observ_pt_1;
    observ_pt_1[0] = 1.0;
    observ_pt_1[1] = 1.0;
    observ_pt_1[2] = 0.0;

    // rotation tensor
    il::StaticArray2D<double, 3, 3> r_tensor;

    // shape function coefficient matrix
    il::StaticArray2D<std::complex<double>, 6, 6> sfm =
        hfp3d::make_el_sfm_uniform(nodes_crd, il::io, r_tensor);

    il::StaticArray<std::complex<double>, 3> tau =
        hfp3d::make_el_tau_crd(nodes_crd, r_tensor);

    hfp3d::HZ hz = hfp3d::make_el_pt_hz(nodes_crd, observ_pt_1, r_tensor);

    il::StaticArray2D<double, 6, 18> bem_matrix_1 =
        hfp3d::make_local_3dbem_submatrix(1, G, nu, hz.h, hz.z, tau, sfm);


    std::string work_directory
            {"C:/Users/nikolski/ClionProjects/HFPx3D/src/Elasticity/_test/"};
    std::string tf_name{"integration_test_triangle_pt1.csv"};

    // loading reference matrix (built with MATLAB code)
    il::Status status{};
    il::Array2D<double> bem_matrix_matlab =
            hfp3d::load_data_from_csv
                    (work_directory, tf_name, il::io, status);
    //status.abortOnError();

    // checking for equal size
    std::cout << "[          ] " << bem_matrix_matlab.size(0)
              << "x" << bem_matrix_matlab.size(1)
              << " matrix loaded" << std::endl;
    ASSERT_EQ(bem_matrix_matlab.size(0), bem_matrix_1.size(0));
    ASSERT_EQ(bem_matrix_matlab.size(1), bem_matrix_1.size(1));

    // calculating the difference
    double err = 0.0, diff_ij;
    il::Array2D <double> diff_matrix {6, 18, 0.0};
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 18; ++j) {
            diff_matrix(i, j) = bem_matrix_1(i, j) - bem_matrix_matlab(i, j);
            diff_ij = std::fabs(diff_matrix(i, j));
            if (diff_ij > tol) {
                std::cout << "[          ] err = " << diff_ij
                          <<" in (" << i << ", "<< j
                          << ") element" << std::endl;
            }
            if (diff_ij > err) {
                err = diff_ij;
            }
        }
    }

    std::cout << "[          ] Max |M1-M2|=" << err << std::endl;

    ASSERT_NEAR(0.0, err, std::sqrt(tol));
}

TEST(kernel_integration_2, t1) {
    // tolerance
    double tol = 1e-06;

    // elastic properties
    double G = 1.0, nu = 0.35;

    // connectivity
    il::StaticArray2D<il::int_t, 3, 1> mesh_conn;
    for (il::int_t k = 0; k < 3; ++k) {
        mesh_conn(k, 0) = k;
    }

    // vertex coordinates
    il::StaticArray2D<double, 3, 3> nodes_crd;
    nodes_crd(0, 0) = 0.0;
    nodes_crd(1, 0) = 0.1;
    nodes_crd(2, 0) = 0.0;

    nodes_crd(0, 1) = 1.8;
    nodes_crd(1, 1) = 0.0;
    nodes_crd(2, 1) = 0.0;

    nodes_crd(0, 2) = 1.2;
    nodes_crd(1, 2) = 1.8;
    nodes_crd(2, 2) = 0.0;

    // observation point #2 (on an edge extension)
    il::StaticArray<double, 3> observ_pt_2;
    observ_pt_2[0] = 1.2;
    observ_pt_2[1] = 1.8;
    observ_pt_2[2] = -1.21;

    // rotation tensor
    il::StaticArray2D<double, 3, 3> r_tensor;

    // shape function coefficient matrix
    il::StaticArray2D<std::complex<double>, 6, 6> sfm =
            hfp3d::make_el_sfm_uniform(nodes_crd, il::io, r_tensor);

    il::StaticArray<std::complex<double>, 3> tau =
            hfp3d::make_el_tau_crd(nodes_crd, r_tensor);

    hfp3d::HZ hz = hfp3d::make_el_pt_hz(nodes_crd, observ_pt_2, r_tensor);

    il::StaticArray2D<double, 6, 18> bem_matrix_2 =
            hfp3d::make_local_3dbem_submatrix(1, G, nu, hz.h, hz.z, tau, sfm);


    std::string work_directory
            {"C:/Users/nikolski/ClionProjects/HFPx3D/src/Elasticity/_test/"};
    std::string tf_name{"integration_test_triangle_pt2.csv"};

    // loading reference matrix (built with MATLAB code)
    il::Status status{};
    il::Array2D<double> bem_matrix_matlab =
            hfp3d::load_data_from_csv
                    (work_directory, tf_name, il::io, status);
    //status.abortOnError();

    // checking for equal size
    std::cout << "[          ] " << bem_matrix_matlab.size(0)
              << "x" << bem_matrix_matlab.size(1)
              << " matrix loaded" << std::endl;
    ASSERT_EQ(bem_matrix_matlab.size(0), bem_matrix_2.size(0));
    ASSERT_EQ(bem_matrix_matlab.size(1), bem_matrix_2.size(1));

    // calculating the difference
    double err = 0.0, diff_ij;
    il::Array2D <double> diff_matrix {6, 18, 0.0};
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 18; ++j) {
            diff_matrix(i, j) = bem_matrix_2(i, j) - bem_matrix_matlab(i, j);
            diff_ij = std::fabs(diff_matrix(i, j));
            if (diff_ij > tol) {
                std::cout << "[          ] err = " << diff_ij
                          <<" in (" << i << ", "<< j
                          << ") element" << std::endl;
            }
            if (diff_ij > err) {
                err = diff_ij;
            }
        }
    }

    std::cout << "[          ] Max |M1-M2|=" << err << std::endl;

    ASSERT_NEAR(0.0, err, std::sqrt(tol));
}

TEST(kernel_integration_3, t1) {
    // tolerance
    double tol = 1e-06;

    // elastic properties
    double G = 1.0, nu = 0.35;

    // connectivity
    il::StaticArray2D<il::int_t, 3, 1> mesh_conn;
    for (il::int_t k = 0; k < 3; ++k) {
        mesh_conn(k, 0) = k;
    }

    // vertex coordinates
    il::StaticArray2D<double, 3, 3> nodes_crd;
    nodes_crd(0, 0) = 0.0;
    nodes_crd(1, 0) = 0.1;
    nodes_crd(2, 0) = 0.0;

    nodes_crd(0, 1) = 1.8;
    nodes_crd(1, 1) = 0.0;
    nodes_crd(2, 1) = 0.0;

    nodes_crd(0, 2) = 1.2;
    nodes_crd(1, 2) = 1.8;
    nodes_crd(2, 2) = 0.0;

    // observation point 3 (off the plane)
    il::StaticArray<double, 3> observ_pt_3;
    observ_pt_3[0] = 1.0;
    observ_pt_3[1] = 1.0;
    observ_pt_3[2] = -1.21;

    // rotation tensor
    il::StaticArray2D<double, 3, 3> r_tensor;

    // shape function coefficient matrix
    il::StaticArray2D<std::complex<double>, 6, 6> sfm =
            hfp3d::make_el_sfm_uniform(nodes_crd, il::io, r_tensor);

    il::StaticArray<std::complex<double>, 3> tau =
            hfp3d::make_el_tau_crd(nodes_crd, r_tensor);

    // observation point #3 (off the plane)
    hfp3d::HZ hz = hfp3d::make_el_pt_hz(nodes_crd, observ_pt_3, r_tensor);

    il::StaticArray2D<double, 6, 18> bem_matrix_3 =
            hfp3d::make_local_3dbem_submatrix(1, G, nu, hz.h, hz.z, tau, sfm);

    std::string work_directory
            {"C:/Users/nikolski/ClionProjects/HFPx3D/src/Elasticity/_test/"};
    std::string tf_name{"integration_test_triangle_pt3.csv"};

    // loading reference matrix (built with MATLAB code)
    il::Status status{};
    il::Array2D<double> bem_matrix_matlab =
            hfp3d::load_data_from_csv
                    (work_directory, tf_name, il::io, status);
    //status.abortOnError();

    // checking for equal size
    std::cout << "[          ] " << bem_matrix_matlab.size(0)
              << "x" << bem_matrix_matlab.size(1)
              << " matrix loaded" << std::endl;
    ASSERT_EQ(bem_matrix_matlab.size(0), bem_matrix_3.size(0));
    ASSERT_EQ(bem_matrix_matlab.size(1), bem_matrix_3.size(1));

    // calculating the difference
    double err = 0.0, diff_ij;
    il::Array2D <double> diff_matrix {6, 18, 0.0};
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 18; ++j) {
            diff_matrix(i, j) = bem_matrix_3(i, j) - bem_matrix_matlab(i, j);
            diff_ij = std::fabs(diff_matrix(i, j));
            if (diff_ij > tol) {
                std::cout << "[          ] err = " << diff_ij
                          <<" in (" << i << ", "<< j
                          << ") element" << std::endl;
            }
            if (diff_ij > err) {
                err = diff_ij;
            }
        }
    }

    std::cout << "[          ] Max |M1-M2|=" << err << std::endl;

    ASSERT_NEAR(0.0, err, std::sqrt(tol));
}
//
// This file is part of HFPx3D.
//
// Created by nikolski on 11/14/2017.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//

#include <gtest/gtest/gtest.h>
#include <il/math.h>
#include <il/Status.h>
#include <il/Array2D.h>
#include "src/Solvers/system_assembly.h"
#include "src/Core/model_parameters.h"
#include "src/Core/surface_mesh_utilities.h"
#include "src/Core/element_utilities.h"
#include "src/Core/tensor_utilities.h"
#include "src/IO/mesh_file_io.h"
#include "src/IO/data_file_io.h"

//  UNIT TESTS FOR ELASTICITY INFLUENCE MATRIX / RHS ASSEMBLY

// Testing matrix assembly on one element
TEST(system_assembly_1, t1) {

// tolerance
double tol = 1e-06;

// initializing the mesh etc.
hfp3d::Mesh_Data_T mesh_data;
hfp3d::Properties_T mat_props;
hfp3d::Num_Param_T num_param;

il::Status status{};

    // elastic properties
    mat_props.shear_m[0] = 1.0;
    mat_props.poiss_r[0] = 0.35;

    // mesh (one element)
    mesh_data.mesh.conn.resize(3, 1);
    for (il::int_t k = 0; k < 3; ++k) {
        mesh_data.mesh.conn(k, 0) = k;
    }

    mesh_data.mesh.nods.resize(3, 3);
    mesh_data.mesh.nods(0, 0) = 0.0;
    mesh_data.mesh.nods(1, 0) = 0.1;
    mesh_data.mesh.nods(2, 0) = 0.0;

    mesh_data.mesh.nods(0, 1) = 1.8;
    mesh_data.mesh.nods(1, 1) = 0.0;
    mesh_data.mesh.nods(2, 1) = 0.0;

    mesh_data.mesh.nods(0, 2) = 1.2;
    mesh_data.mesh.nods(1, 2) = 1.8;
    mesh_data.mesh.nods(2, 2) = 0.0;

    // numerical model parameters
    num_param.tip_type = 0;
    num_param.is_dd_local = false;
    num_param.beta = 0.125;

    // initializing the DoF handle
    mesh_data.dof_h_dd = hfp3d::make_dof_h_crack
        (mesh_data.mesh, 2, num_param.tip_type);

    // assembly of the algebraic system (matrix + RHS)
    hfp3d::SAE_T sae;
    // matrix
    sae.matrix = hfp3d::make_3dbem_matrix_s
        (mat_props.shear_m[0], mat_props.poiss_r[0],
         mesh_data.mesh,
         num_param,
         il::io, mesh_data.dof_h_dd);

    // importing data (matrix made by the MATLAB code)
    std::string work_directory
        {"C:/Users/nikolski/ClionProjects/HFPx3D/src/Solvers/_test/"};
    std::string tf_name{"assembly_test_triangle.csv"};
    // std::string of_name{"test_assembly_1_ele.csv"};

    // loading reference matrix (built with MATLAB code)
    il::Array2D<double> bem_matrix_matlab =
        hfp3d::load_data_from_csv
                (work_directory, tf_name, il::io, status);
    //status.abortOnError();

    // checking for equal size
    std::cout << "[          ] " << bem_matrix_matlab.size(0)
    << "x" << bem_matrix_matlab.size(1)
    << " matrix loaded" << std::endl;
    ASSERT_EQ(bem_matrix_matlab.size(0), sae.matrix.size(0));
    ASSERT_EQ(bem_matrix_matlab.size(1), sae.matrix.size(1));

    // calculating the difference
    double err = 0.0, diff_ij;
    il::Array2D <double> diff_matrix {6, 18, 0.0};
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 18; ++j) {
            diff_matrix(i, j) = sae.matrix(i, j) - bem_matrix_matlab(i, j);
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

// Testing matrix assembly on two non-coplanar elements
TEST(system_assembly_2, t1) {

    // tolerance
    double tol = 1e-06;

    // initializing the mesh etc.
    hfp3d::Mesh_Data_T mesh_data;
    hfp3d::Properties_T mat_props;
    hfp3d::Num_Param_T num_param;
    // hfp3d::Load_T load_data;

    il::Status status{};

    // elastic properties
    mat_props.shear_m[0] = 1.0, mat_props.poiss_r[0] = 0.35;

    // mesh (one element)
    mesh_data.mesh.conn.resize(3, 2);
    for (il::int_t k = 0; k < 3; ++k) {
        mesh_data.mesh.conn(k, 0) = k;
        mesh_data.mesh.conn(k, 1) = k + 1;
    }

    mesh_data.mesh.nods.resize(3, 4);
    mesh_data.mesh.nods(0, 0) = 0.0;
    mesh_data.mesh.nods(1, 0) = 0.1;
    mesh_data.mesh.nods(2, 0) = 0.0;

    mesh_data.mesh.nods(0, 1) = 1.8;
    mesh_data.mesh.nods(1, 1) = 0.0;
    mesh_data.mesh.nods(2, 1) = 0.0;

    mesh_data.mesh.nods(0, 2) = 1.2;
    mesh_data.mesh.nods(1, 2) = 1.8;
    mesh_data.mesh.nods(2, 2) = 0.0;

    mesh_data.mesh.nods(0, 3) = 2.2;
    mesh_data.mesh.nods(1, 3) = 1.72;
    mesh_data.mesh.nods(2, 3) = -1.21;

    // numerical model parameters
    num_param.tip_type = 0;
    num_param.is_dd_local = false;
    num_param.beta = 0.125;

    // load (S_inf)

    // initializing the DoF handle
    mesh_data.dof_h_dd = hfp3d::make_dof_h_crack
        (mesh_data.mesh, 2, num_param.tip_type);

    // assembly of the algebraic system (matrix + RHS)
    hfp3d::SAE_T sae;
    // matrix
    sae.matrix = hfp3d::make_3dbem_matrix_s
        (mat_props.shear_m[0], mat_props.poiss_r[0],
         mesh_data.mesh,
         num_param,
         il::io, mesh_data.dof_h_dd);
//    // right-hand side
//    hfp3d::add_s_inf_to_3dbem_rhs
//            (mesh_data, load_data, il::io, sae);

//    hfp3d::save_data_to_csv(sae.matrix, work_directory, of_name);

    // importing data (matrix made by the MATLAB code)
    std::string work_directory
        {"C:/Users/nikolski/ClionProjects/HFPx3D/src/Solvers/_test/"};
    std::string tf_name{"assembly_test_bitriangle.csv"};
    // std::string of_name{"test_assembly_2_ele.csv"};

    // loading reference matrix (built with MATLAB code)
    il::Array2D<double> bem_matrix_matlab =
        hfp3d::load_data_from_csv
                (work_directory, tf_name, il::io, status);
    //status.abortOnError();

    // checking for equal size
    std::cout << "[          ] " << bem_matrix_matlab.size(0)
    << "x" << bem_matrix_matlab.size(1)
    << " matrix loaded" << std::endl;
    ASSERT_EQ(bem_matrix_matlab.size(0), sae.matrix.size(0));
    ASSERT_EQ(bem_matrix_matlab.size(1), sae.matrix.size(1));

    // calculating the difference
    double err = 0.0, diff_ij;
    il::Array2D <double> diff_matrix {6, 18, 0.0};
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 18; ++j) {
            diff_matrix(i, j) = sae.matrix(i, j) - bem_matrix_matlab(i, j);
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

// Testing matrix assembly on a radial crack with 24 elements
TEST(system_assembly_3, t1) {

    // tolerance
    double tol = 1e-06;

    // initializing the mesh etc.
    hfp3d::Mesh_Data_T mesh_data;
    hfp3d::Properties_T mat_props;
    hfp3d::Num_Param_T num_param;

    il::Status status{};

    // elastic properties
    mat_props.shear_m[0] = 1.0, mat_props.poiss_r[0] = 0.35;

    // mesh (importing from 2 numpy files)
    std::string work_directory
            {"C:/Users/nikolski/ClionProjects/HFPx3D/src/Solvers/_test/"};
    std::string m_c_f_name{"Elems_pennymesh24el_32.npy"};
    std::string m_n_f_name{"Nodes_pennymesh24el_32.npy"};
    hfp3d::load_mesh_from_numpy_32
            (work_directory, m_c_f_name, m_n_f_name, 1,
             il::io, mesh_data.mesh);

    // numerical model parameters
    num_param.tip_type = 0;
    num_param.is_dd_local = false;
    num_param.beta = 0.125;

    // initializing the DoF handle
    mesh_data.dof_h_dd = hfp3d::make_dof_h_crack
            (mesh_data.mesh, 2, num_param.tip_type);

    // assembly of the algebraic system (matrix + RHS)
    hfp3d::SAE_T sae;
    // matrix
    sae.matrix = hfp3d::make_3dbem_matrix_s
            (mat_props.shear_m[0], mat_props.poiss_r[0],
             mesh_data.mesh,
             num_param,
             il::io, mesh_data.dof_h_dd);

    // importing data (matrix made by the MATLAB code)
    std::string tf_name{"assembly_test_radial_crack_24el.csv"};

    // loading reference matrix (built with MATLAB code)
    il::Array2D<double> bem_matrix_matlab =
            hfp3d::load_data_from_csv
                    (work_directory, tf_name, il::io, status);
    //status.abortOnError();

    // checking for equal size
    std::cout << "[          ] " << bem_matrix_matlab.size(0)
              << "x" << bem_matrix_matlab.size(1)
              << " matrix loaded" << std::endl;
    ASSERT_EQ(bem_matrix_matlab.size(0), sae.matrix.size(0));
    ASSERT_EQ(bem_matrix_matlab.size(1), sae.matrix.size(1));

    // calculating the difference
    double err = 0.0, diff_ij;
    il::Array2D <double> diff_matrix {6, 18, 0.0};
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 18; ++j) {
            diff_matrix(i, j) = sae.matrix(i, j) - bem_matrix_matlab(i, j);
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

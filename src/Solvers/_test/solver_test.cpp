//
// This file is part of HFPx3D.
//
// Created by nikolski on 11/11/2017.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//

#include <gtest/gtest/gtest.h>
#include "src/Solvers/test_radial_crack_static.h"

// UNIT TEST FOR SOLVERS & POST-PROCESSOR

// Testing static solver on a radial crack
TEST(static_solver_1, t1) {

    // tolerance
    double tol = 0.1;

    // mesh (importing from 2 numpy files)
    std::string work_directory
            {"C:/Users/nikolski/ClionProjects/HFPx3D/src/Solvers/_test/"};
    // connectivity and node coordinates
    // NOTE: doesn't work for read-only files
    std::string m_c_f_name{"Elems_pennymesh24el_32.npy"};
    std::string m_n_f_name{"Nodes_pennymesh24el_32.npy"};

    // numerical solution
    double err = hfp3d::test_dd_radial_crack_static
            (work_directory, m_c_f_name, m_n_f_name, 1, 0.1 * tol, true);

    std::cout << "[          ] Max |DD-DD_ref|=" << err << std::endl;

    ASSERT_NEAR(0.0, err, tol);
}

// Testing post-processing (stress field reconstruction) on a radial crack
TEST(postprocessing_1, t1) {
    // tolerance
    double tol = 0.1;

    // mesh (importing from 2 numpy files)
    std::string work_directory
            {"C:/Users/nikolski/ClionProjects/HFPx3D/src/Solvers/_test/"};
    // connectivity and node coordinates
    // NOTE: doesn't work for read-only files
    std::string m_c_f_name{"Elems_pennymesh24el_32.npy"};
    std::string m_n_f_name{"Nodes_pennymesh24el_32.npy"};

    double err = hfp3d::test_stresses_radial_crack_static
            (work_directory, m_c_f_name, m_n_f_name, 5, 0.1 * tol, true);

    std::cout << "[          ] Max |S-S_ref|=" << err << std::endl;

    ASSERT_NEAR(0.0, err, tol);
}
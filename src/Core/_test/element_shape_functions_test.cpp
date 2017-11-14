//
// This file is part of HFPx3D.
//
// Created by nikolski on 11/11/2017.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//

#include <gtest/gtest/gtest.h>
#include <complex>
#include <il/norm.h>
//#include <il/math.h>
#include <il/StaticArray.h>
#include <il/StaticArray2D.h>
#include "src/Core/element_utilities.h"

//  UNIT TEST FOR SHAPE FUNCTIONS, COORDINATE ROTATION $ TRANSLATION

TEST(sf_test_1, t1) {
    // tolerance
    double tol = 1e-06;

    // elastic properties
    double G = 1.0, poiss_r = 0.35;

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

    // observation point
    il::StaticArray<double, 3> observ_pt_crd;
    observ_pt_crd[0] = 1.0;
    observ_pt_crd[1] = 1.4;
    observ_pt_crd[2] = -1.21;

    // rotation tensor
    il::StaticArray2D<double, 3, 3> r_tensor;

    // shape function coefficient matrix
    il::StaticArray2D<std::complex<double>, 6, 6> sfm =
        hfp3d::make_el_sfm_uniform(nodes_crd, il::io, r_tensor);

    il::StaticArray<std::complex<double>, 3> tau =
        hfp3d::make_el_tau_crd(nodes_crd, r_tensor);

    hfp3d::HZ hz = hfp3d::make_el_pt_hz(nodes_crd, observ_pt_crd, r_tensor);

    il::StaticArray2D<std::complex<double>, 6, 6> reference_sfm, diff_m;

    reference_sfm(0, 0) = std::complex<double>{1.0, 0.0};
    reference_sfm(0, 1) = std::complex<double>{-0.832050294337844, 0.329680305303674};
    reference_sfm(0, 2) = std::complex<double>{-0.832050294337844, -0.329680305303674};
    reference_sfm(0, 3) = std::complex<double>{0.129693019689460, -0.121915820029028};
    reference_sfm(0, 4) = std::complex<double>{0.129693019689460, 0.121915820029028};
    reference_sfm(0, 5) = std::complex<double>{0.355998576005696, 0.0};

    reference_sfm(1, 0) = std::complex<double>{0.0, 0.0};
    reference_sfm(1, 1) = std::complex<double>{-0.277350098112615, -0.173561853850347};
    reference_sfm(1, 2) = std::complex<double>{-0.277350098112615, 0.173561853850347};
    reference_sfm(1, 3) = std::complex<double>{0.0935987196222155, 0.192549588776004};
    reference_sfm(1, 4) = std::complex<double>{0.0935987196222155, -0.192549588776004};
    reference_sfm(1, 5) = std::complex<double>{0.428187176140184, 0.0};

    reference_sfm(2, 0) = std::complex<double>{0.0, 0.0};
    reference_sfm(2, 1) = std::complex<double>{0.0, 0.283455288951571};
    reference_sfm(2, 2) = std::complex<double>{0.0, -0.283455288951571};
    reference_sfm(2, 3) = std::complex<double>{-0.160693801669238, 0.0};
    reference_sfm(2, 4) = std::complex<double>{-0.160693801669238, -0.0};
    reference_sfm(2, 5) = std::complex<double>{0.321387603338476, 0.0};

    reference_sfm(3, 0) = std::complex<double>{0.0, 0.0};
    reference_sfm(3, 1) = std::complex<double>{0.0, 0.0};
    reference_sfm(3, 2) = std::complex<double>{0.0, 0.0};
    reference_sfm(3, 3) = std::complex<double>{0.196788101736482, -0.314465408805031};
    reference_sfm(3, 4) = std::complex<double>{0.196788101736482, 0.314465408805031};
    reference_sfm(3, 5) = std::complex<double>{-0.393576203472964, 0.0};

    reference_sfm(4, 0) = std::complex<double>{0.0, 0.0};
    reference_sfm(4, 1) = std::complex<double>{0.0, -1.13382115580629};
    reference_sfm(4, 2) = std::complex<double>{0.0, 1.13382115580629};
    reference_sfm(4, 3) = std::complex<double>{0.124599501601994, 0.314465408805031};
    reference_sfm(4, 4) = std::complex<double>{0.124599501601994, -0.314465408805031};
    reference_sfm(4, 5) = std::complex<double>{-0.249199003203987, 0.0};

    reference_sfm(5, 0) = std::complex<double>{0.0, 0.0};
    reference_sfm(5, 1) = std::complex<double>{1.10940039245046, 0.694247415401387};
    reference_sfm(5, 2) = std::complex<double>{1.10940039245046, -0.694247415401387};
    reference_sfm(5, 3) = std::complex<double>{-0.383985540980913, -0.0706337687469762};
    reference_sfm(5, 4) = std::complex<double>{-0.383985540980913, 0.0706337687469762};
    reference_sfm(5, 5) = std::complex<double>{-0.462798148807405, 0.0};

    double err = 0.0, diff_ij;
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            diff_m(i, j) = sfm(i, j) - reference_sfm(i, j);
            diff_ij = std::abs(diff_m(i, j));
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

    ASSERT_NEAR(0.0, err, tol);

}

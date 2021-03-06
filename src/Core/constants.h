//
// This file is part of HFPx3D.
//
// Created by nikolski on 6/29/2017.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//

#ifndef INC_HFPX3D_CONSTANTS_H
#define INC_HFPX3D_CONSTANTS_H

#include <complex>

namespace hfp3d {

    const double pi = 3.1415926535897932385;
    const std::complex<double> ii = std::complex<double>{0.0, 1.0};
    // tolerance parameters
    const double h_tol = 2.221e-016;
    const double a_tol = 1.825e-008;

}

#endif //INC_HFPX3D_CONSTANTS_H

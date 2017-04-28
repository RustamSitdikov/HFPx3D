//
// This file is part of HFPx3D_VC.
//
// Created by D. Nikolski on 1/12/2017.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//

#ifndef INC_HFPX3D_ELAST_KER_INT_H
#define INC_HFPX3D_ELAST_KER_INT_H

#include <complex>
#include <il/StaticArray.h>
#include <il/StaticArray3D.h>
#include <il/StaticArray4D.h>

// Integration of a kernel of the elasticity equation
// over a part of a polygonal element (a sector associated with one edge)
// with 2nd order polynomial approximating (shape) functions.

namespace hfp3d {

// Coefficient matrices (rank 3) to be contracted with the vector of
// constituing functions defined below (via right multiplication)
// and with the vector of shape function coefficients
// associated with each node of the element (via left multiplication)
//
// Stress components (vs local Cartesian coordinate system of the element)
// combined as S11+S22, S11-S22+2*I*S12, S13+I*S23, S33

    il::StaticArray4D<std::complex<double>, 6, 4, 3, 9> s_integral_gen
                (const int ker,
                 double nu, std::complex<double> eix,
                 double h, std::complex<double> d);

    il::StaticArray4D<std::complex<double>, 6, 4, 3, 5> s_integral_red
                (const int kernel_id,
                 double nu, std::complex<double> eix,
                 double h);

    il::StaticArray3D<std::complex<double>, 6, 4, 3> s_integral_lim
                (const int ker,
                 double nu, std::complex<double> eix,
                 std::complex<double> d);


// Constituing functions for the integrals
// of any kernel of the elasticity equation
// over a part of a polygonal element.
// Example of usage:
// dot(S11_22H(nu, eix, h, d), integral_cst_fun(h, d, a, x, eix))
// dot(S11_22H_red(nu, eip, h, d), integral_cst_fun_red(h, d, a))
// dot(S13_23T(nu, eix, h, d), integral_cst_fun(h, d, a, x, eix))
// dot(S33T_red(nu, eip, h, d), integral_cst_fun_red(h, d, a))
// where eip = std::exp(I*std::arg(t-z));
// eix = std::exp(I*x); x = std::arg(t-z)-std::arg(d);
// a = std::fabs(t-z-d)*sign(x);

    il::StaticArray<std::complex<double>, 9> integral_cst_fun
            (double h, std::complex<double> d, double a,
             double x, std::complex<double> eix);

    il::StaticArray<std::complex<double>, 5> integral_cst_fun_red
            (double h, std::complex<double> d, double a);

}

#endif //INC_HFPX3D_ELAST_KER_INT_H

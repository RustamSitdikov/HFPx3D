//
// This file is part of 3d_bem.
//
// Created by nikolski on 1/12/2017.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//

#ifndef INC_3D_BEM_SIJK_H
#define INC_3D_BEM_SIJK_H

#endif //INC_3D_BEM_SIJK_H

//#include <il/StaticArray.h>
//#include <il/StaticArray2D.h>
//#include <il/StaticArray3D.h>
//#include <il/linear_algebra/dense/blas/blas.h>
//#include <il/linear_algebra/dense/blas/dot.h>
//#include <cmath>
//#include <complex>

// Integration of a kernel of the elasticity equation
// over a part of a polygonal element (a sector associated with one edge)
// with 2nd order polynomial approximating (shape) functions.
//
// Coefficient matrices (rank 3) to be contracted (via right multiplication)
// with the vector of constituing functions defined in ICFNS.cpp (via right multiplication)
// and with the vector of shape function coefficients
// associated with each node of the element (via left multiplication)
//
// Stress components (vs local Cartesian coordinate system of the element)
// combined as S11+S22, S11-S22+2*I*S12, S13+I*S23, S33

class Kernel_Integration{
public:
    il::StaticArray3D<std::complex<double>, 6, 3, 9> S11_22(double nu, std::complex<double> eix, double h, std::complex<double> d);
    il::StaticArray3D<std::complex<double>, 6, 3, 9> S11_22_12(double nu, std::complex<double> eix, double h, std::complex<double> d);
    il::StaticArray3D<std::complex<double>, 6, 3, 9> S13_23(double nu, std::complex<double> eix, double h, std::complex<double> d);
    il::StaticArray3D<std::complex<double>, 6, 3, 9> S33(double nu, std::complex<double> eix, double h, std::complex<double> d);

    il::StaticArray3D<std::complex<double>, 6, 3, 5> S11_22_red(double nu, std::complex<double> eix, double h, std::complex<double> d);
    il::StaticArray3D<std::complex<double>, 6, 3, 5> S11_22_12_red(double nu, std::complex<double> eix, double h, std::complex<double> d);
    il::StaticArray3D<std::complex<double>, 6, 3, 5> S13_23_red(double nu, std::complex<double> eix, double h, std::complex<double> d);
    il::StaticArray3D<std::complex<double>, 6, 3, 5> S33_red(double nu, std::complex<double> eix, double h, std::complex<double> d);

    il::StaticArray3D<std::complex<double>, 6, 4, 3> SijLim(double nu, std::complex<double> eix, std::complex<double> d);
};

//
// This file is part of 3d_bem.
//
// Created by D. Nikolski on 1/12/2017.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//

#ifndef INC_3D_BEM_SIJK_H
#define INC_3D_BEM_SIJK_H

#endif //INC_3D_BEM_SIJK_H

//#include <il/StaticArray2D.h>
#include <il/StaticArray.h>
#include <il/StaticArray3D.h>
#include <cmath>
#include <complex>

// Integration of a kernel of the elasticity equation
// over a part of a polygonal element (a sector associated with one edge)
// with 2nd order polynomial approximating (shape) functions.
//
// Coefficient matrices (rank 3) to be contracted (via right multiplication)
// with the vector of constituing functions defined below (via right multiplication)
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

// Constituing functions for the integrals
// of any kernel of the elasticity equation
// over a part of a polygonal element.
// Example of usage:
// dot(S11_22H(nu, eix, h, d), ICFns(h, d, a, x, eix))
// dot(S11_22H_red(nu, eip, h, d), ICFns_red(h, d, a))
// dot(S13_23T(nu, eix, h, d), ICFns(h, d, a, x, eix))
// dot(S33T_red(nu, eip, h, d), ICFns_red(h, d, a))
// where eip = std::exp(I*std::arg(t-z));
// eix = std::exp(I*x); x = std::arg(t-z)-std::arg(d);
// a = std::fabs(t-z-d)*sign(x);

il::StaticArray<std::complex<double>, 9> ICFns(double, std::complex<double>, double, double, std::complex<double>);
il::StaticArray<std::complex<double>, 5> ICFns_red(double, std::complex<double>, double);

// General case (h!=0, collocation point projected into or outside the element)
// powers of r, G0=arctan((ah)/(dr)), H0=arctanh(a/r) and its derivatives w.r. to h

il::StaticArray<std::complex<double>, 9> ICFns(double h, std::complex<double> d, double a, double x, std::complex<double> eix) {

    double D1 = std::abs(d), D2 = D1*D1, a2 = a*a,
            r = std::sqrt(h*h + a2 + D2),
            r2 = r*r, r3 = r2*r, r5 = r3*r2,
            ar = a/r, ar2 = ar*ar,
            hr = std::fabs(h/r),
            B = 1.0/(r2 - a2), B2=B*B, B3=B2*B;
    // evaluation of x and eix can be added here; then change arguments to pointer type
    double TanHi = std::imag(eix)/std::real(eix), tr = hr*TanHi,
            G0 = std::atan(tr), H0 = std::atanh(ar),
            H1 = -0.5*ar*B, H2 = 0.25*(3.0 - ar2)*ar*B2,
            H3 = -0.125*(15.0 - 10.0*ar2 + 3.0*ar2*ar2)*ar*B3;

    il::StaticArray<std::complex<double>, 9> BCE {0.0};
    BCE[0] = r; BCE[1] = 1.0/r; BCE[2] = 1.0/r3; BCE[3] = 1.0/r5;
    BCE[4] = G0-x; BCE[5] = H0; BCE[6] = H1; BCE[7] = H2; BCE[8] = H3;

    return BCE;
}

// Special case (reduced summation, collocation point projected onto the element contour) - additional terms

il::StaticArray<std::complex<double>, 5> ICFns_red(double h, std::complex<double> d, double a) {

    double h2 = h*h, h4 = h2*h2, h6 = h4*h2,
            D1 = std::abs(d), D2 = D1*D1, a2 = a*a,
            ro = std::sqrt(a2 + D2),
            r = std::sqrt(h2 + a2 + D2),
            rr = ro/r, rr2 = rr*rr, rr4 = rr2*rr2,
            L0 = std::atanh(rr), L1 = -0.5*rr/h2, L2 = 0.25*(3.0 - rr2)*rr/h4,
            L3 = -0.125*(15.0 - 10.0*rr2 + 3.0*rr4)*rr/h6;

    il::StaticArray<std::complex<double>, 5> BCE {0.0};
    BCE[0] = 1.0; BCE[1] = L0; BCE[2] = L1; BCE[3] = L2; BCE[4] = L3;

    return BCE;
}

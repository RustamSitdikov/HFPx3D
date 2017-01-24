//
// This file is part of 3d_bem.
//
// Created by D. Nikolski on 1/12/2017.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//

#ifndef INC_3D_BEM_ELAST_KER_INT_H
#define INC_3D_BEM_ELAST_KER_INT_H

#include <cmath>
#include <complex>
#include <il/StaticArray.h>
#include <il/StaticArray3D.h>

// Integration of a kernel of the elasticity equation
// over a part of a polygonal element (a sector associated with one edge)
// with 2nd order polynomial approximating (shape) functions.
//
// Coefficient matrices (rank 3) to be contracted with the vector of
// constituing functions defined below (via right multiplication)
// and with the vector of shape function coefficients
// associated with each node of the element (via left multiplication)
//
// Stress components (vs local Cartesian coordinate system of the element)
// combined as S11+S22, S11-S22+2*I*S12, S13+I*S23, S33

namespace hfp3d {

    class Kernel_Integration {
    public:
        il::StaticArray3D<std::complex<double>, 6, 3, 9> S11_22
                (double nu, std::complex<double> eix,
                 double h, std::complex<double> d);

        il::StaticArray3D<std::complex<double>, 6, 3, 9> S11_22_12
                (double nu, std::complex<double> eix,
                 double h, std::complex<double> d);

        il::StaticArray3D<std::complex<double>, 6, 3, 9> S13_23
                (double nu, std::complex<double> eix,
                 double h, std::complex<double> d);

        il::StaticArray3D<std::complex<double>, 6, 3, 9> S33
                (double nu, std::complex<double> eix,
                 double h, std::complex<double> d);

        il::StaticArray3D<std::complex<double>, 6, 3, 5> S11_22_red
                (double nu, std::complex<double> eix,
                 double h, std::complex<double> d);

        il::StaticArray3D<std::complex<double>, 6, 3, 5> S11_22_12_red
                (double nu, std::complex<double> eix,
                 double h, std::complex<double> d);

        il::StaticArray3D<std::complex<double>, 6, 3, 5> S13_23_red
                (double nu, std::complex<double> eix,
                 double h, std::complex<double> d);

        il::StaticArray3D<std::complex<double>, 6, 3, 5> S33_red
                (double nu, std::complex<double> eix,
                 double h, std::complex<double> d);

        il::StaticArray3D<std::complex<double>, 6, 4, 3> SijLim
                (double nu, std::complex<double> eix, std::complex<double> d);
    };

    // Hypersingular potential - for fracture boundaries
    class H_Potential : public Kernel_Integration {
    };
    // This if for non-fracture boundaries (e.g. borehole)
    // class T_Potential: public Kernel_Integration {};

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

    il::StaticArray<std::complex<double>, 9> ICFns
            (double, std::complex<double>, double, double,
             std::complex<double>);

    il::StaticArray<std::complex<double>, 5> ICFns_red
            (double, std::complex<double>, double);

}

#endif //INC_3D_BEM_ELAST_KER_INT_H

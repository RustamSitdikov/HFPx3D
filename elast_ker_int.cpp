//
// This file is part of 3d_bem.
//
// Created by D. Nikolski on 1/24/2017.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//

// Integration of a kernel of the elasticity equation
// over a part of a polygonal element (a sector associated with one edge)
// with 2nd order polynomial approximating (shape) functions.

#include <complex>
#include <il/StaticArray.h>
#include <il/StaticArray3D.h>
#include <il/StaticArray4D.h>
#include "elast_ker_int.h"
#include "h_potential.h"
//#include "t_potential.h"

namespace hfp3d {

// Coefficient matrices (rank 3) to be contracted with the vector of
// constituing functions defined below (via right multiplication)
// and with the vector of shape function coefficients
// associated with each node of the element (via left multiplication)
//
// Stress components (vs local Cartesian coordinate system of the element)
// combined as S11+S22, S11-S22+2*I*S12, S13+I*S23, S33

    il::StaticArray4D<std::complex<double>, 6, 4, 3, 9> s_integral_gen
            (const int kernel_id,
             double nu, std::complex<double> eix,
             double h, std::complex<double> d) {
        il::StaticArray4D<std::complex<double>, 6, 4, 3, 9> c;
        switch (kernel_id) {
            case 1:
                c = s_ij_gen_h(nu, eix, h, d);
                break;
            case 0:
                // c = s_ij_gen_t(nu, eix, h, d);
                break;
            default:break;
        }
        return c;
    }

    il::StaticArray4D<std::complex<double>, 6, 4, 3, 5> s_integral_red
            (const int kernel_id,
             double nu, std::complex<double> eix,
             double h) {
        il::StaticArray4D<std::complex<double>, 6, 4, 3, 5> c;
        switch (kernel_id) {
            case 1:
                c = s_ij_red_h(nu, eix, h);
                break;
            case 0:
                // c = s_ij_red_t(nu, eix, h, d);
                break;
            default:break;
        }
        return c;
    }

    il::StaticArray3D<std::complex<double>, 6, 4, 3> s_integral_lim
            (const int kernel_id,
             double nu, std::complex<double> eix,
             std::complex<double> d) {
        il::StaticArray3D<std::complex<double>, 6, 4, 3> c;
        switch (kernel_id) {
            case 1:
                c = s_ij_lim_h(nu, eix, d);
                break;
            case 0:
                // c = s_ij_lim_t(nu, eix, sgnh, d);
                break;
            default:break;
        }
        return c;
    }


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

// General case (h!=0, collocation point projected into or outside the element)
// powers of r, g0=arctan((ah)/(dr)),
// f0=arctanh(a/r) and its derivatives w.r. to h

    il::StaticArray<std::complex<double>, 9> integral_cst_fun
            (double h, std::complex<double> d, double a,
             double x, std::complex<double> eix) {

        double abs_d = std::abs(d), d2 = abs_d * abs_d, a2 = a * a,
                r = std::sqrt(h * h + a2 + d2),
                r2 = r * r, r3 = r2 * r, r5 = r3 * r2,
                ar = a / r, ar2 = ar * ar,
                hr = std::fabs(h / r),
                b = 1.0 / (r2 - a2), b2 = b * b, b3 = b2 * b;
        double tah_x = std::imag(eix) / std::real(eix), tr = hr * tah_x,
                g0 = std::atan(tr), f0 = std::atanh(ar),
                f1 = -0.5 * ar * b, f2 = 0.25 * (3.0 - ar2) * ar * b2,
                f3 = -0.125 * (15.0 - 10.0 * ar2 + 3.0 * ar2 * ar2) * ar * b3;

        il::StaticArray<std::complex<double>, 9> fun_list{0.0};
        fun_list[0] = r;
        fun_list[1] = 1.0 / r;
        fun_list[2] = 1.0 / r3;
        fun_list[3] = 1.0 / r5;
        fun_list[4] = g0 - x;
        fun_list[5] = f0;
        fun_list[6] = f1;
        fun_list[7] = f2;
        fun_list[8] = f3;

        return fun_list;
    }

// Special case (reduced summation,
// collocation point projected onto the element contour) - additional terms

    il::StaticArray<std::complex<double>, 5> integral_cst_fun_red
            (double h, std::complex<double> d, double a) {

        double h2 = h * h, h4 = h2 * h2, h6 = h4 * h2,
                abs_d = std::abs(d), d2 = abs_d * abs_d, a2 = a * a,
                ro = std::sqrt(a2 + d2),
                r = std::sqrt(h2 + a2 + d2),
                rr = ro / r, rr2 = rr * rr, rr4 = rr2 * rr2,
                f0 = std::atanh(rr), f1 = -0.5 * rr / h2, f2 =
                0.25 * (3.0 - rr2) * rr / h4,
                f3 = -0.125 * (15.0 - 10.0 * rr2 + 3.0 * rr4) * rr / h6;

        il::StaticArray<std::complex<double>, 5> fun_list{0.0};
        fun_list[0] = 1.0;
        fun_list[1] = f0;
        fun_list[2] = f1;
        fun_list[3] = f2;
        fun_list[4] = f3;

        return fun_list;
    }

}
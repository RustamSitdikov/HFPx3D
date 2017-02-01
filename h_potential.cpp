//
// This file is part of 3d_bem.
//
// Created by D. Nikolski on 1/24/2017.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//

// Integration of the hypersingular kernel of the elasticity equation
// over a part of a polygonal element (a sector associated with one edge)
// with 2nd order polynomial approximating (shape) functions.
//
// To be contracted (via right multiplication) with the vector of
// constituing functions defined in SijK.cpp
// and (via left multiplication) with the vector of
// shape function coefficients associated with each node of the element
//
// Stress components (vs local Cartesian coordinate system of the element)
// combined as S11+S22, S11-S22+2*I*S12, S13+I*S23, S33

#include <complex>
#include <il/StaticArray.h>
#include <il/StaticArray3D.h>
#include "h_potential.h"

namespace hfp3d {

// General case (h!=0, collocation point projected into or outside the element)

    il::StaticArray3D<std::complex<double>, 6, 3, 9> s_11_22_h
            (double nu, std::complex<double> eix,
             double h, std::complex<double> d) {
        const std::complex<double> I(0.0, 1.0);

        double c_1_nu = 1.0 + nu;
        double c_1_2nu = 1.0 + 2.0 * nu;
        double c_2_nu = 2.0 + nu;
        double c_3_2nu = 3.0 + 2.0 * nu;
        double c_4_nu = 4.0 + nu;
        double c_5_4nu = 5.0 + 4.0 * nu;
        double c_7_2nu = 7.0 + 2.0 * nu;
        double c_7_5nu = 7.0 + 5.0 * nu;
        double c_7_6nu = 7.0 + 6.0 * nu;
        double c_11_4nu = 11.0 + 4.0 * nu;
        double c_11_5nu = 11.0 + 5.0 * nu;
        double c_13_2nu = 13.0 + 2.0 * nu;
        double c_13_10nu = 13.0 + 10.0 * nu;

        double cos_x = std::real(eix);
        double tan_x = std::imag(eix) / cos_x;
        std::complex<double> tcos_x = cos_x * eix;

        double h2 = h * h;
        double h4 = h2 * h2;
        double sgh = ((h < 0) ? -1.0 : double((h > 0))); // sign(h)

        double abs_d = std::abs(d);
        double abs_d_2 = abs_d * abs_d;
        double abs_d_4 = abs_d_2 * abs_d_2;
        // std::complex<double> d_c = std::conj(d);
        std::complex<double> d_e = std::polar(1.0, std::arg(d)); //  = d/abs_d
        std::complex<double> d_e_2 = d_e * d_e; //  = d^2/abs_d^2

        std::complex<double> c_d_h = abs_d_2 + h2;
        std::complex<double> c_d_3h = abs_d_2 + 3.0 * h2;
        std::complex<double> c_d_m3h = abs_d_2 - 3.0 * h2;

        std::complex<double> p0, p1, p2;

        il::StaticArray3D<std::complex<double>, 6, 3, 9> c_array{0.0};

        c_array(0, 0, 2) = h * std::imag(d);
        c_array(0, 1, 2) = -h * std::real(d);
        p0 = 0.1875 * h;
        p1 = 3.0 * h2;
        p2 = abs_d_2 * tan_x;
        c_array(0, 0, 3) = -p0 * (p1 * std::imag(d) + p2 * std::real(d));
        c_array(0, 1, 3) = p0 * (p1 * std::real(d) - p2 * std::imag(d));
        p0 = c_7_2nu * h;
        c_array(0, 0, 6) = p0 * std::real(d_e);
        c_array(0, 1, 6) = p0 * std::imag(d_e);
        c_array(0, 2, 6) = -c_1_2nu * abs_d;
        p0 = 3.0 * h * c_d_3h;
        c_array(0, 0, 7) = p0 * std::real(d_e);
        c_array(0, 1, 7) = p0 * std::imag(d_e);
        c_array(0, 2, 7) = -2.0 * h2 * abs_d;
        p0 = -0.5 * h * c_d_m3h * c_d_h;
        c_array(0, 0, 8) = p0 * std::real(d_e);
        c_array(0, 1, 8) = p0 * std::imag(d_e);

        c_array(1, 1, 1) = 0.2 * c_11_5nu * h * d_e_2 * tcos_x;
        c_array(1, 0, 1) = I * c_array(1, 1, 1);
        p1 = 0.1 * (abs_d_2 * (7.0 + 2.0 * I * tan_x) + 
                16.0 * h2 * tcos_x) * d_e_2;
        p2 = c_7_5nu / 60.0 * abs_d_2 * tan_x;
        c_array(1, 0, 2) = -(I * p1 + p2) * h;
        c_array(1, 1, 2) = -(p1 + I * p2) * h;
        p1 = (abs_d_2 * h2 * (0.4 + 0.11875 * I * tan_x) + 
                0.09375 * I * abs_d_4 * tan_x +
              0.4 * h4 * tcos_x) * d_e_2;
        p2 = (-0.05625 * abs_d_2 + 0.11875 * h2) * abs_d_2 * tan_x;
        c_array(1, 0, 3) = (I * p1 + p2) * h;
        c_array(1, 1, 3) = (p1 + I * p2) * h;
        c_array(1, 0, 4) = -c_1_nu * sgh;
        c_array(1, 1, 4) = -I * c_1_nu * sgh;
        c_array(1, 2, 5) = 0.5 * c_1_2nu * d_e;
        p1 = 0.5 * c_7_2nu * d * d_e;
        p2 = abs_d * (0.3 + 4.0 / 3.0 * c_2_nu);
        c_array(1, 0, 6) = (p1 + p2) * h;
        c_array(1, 1, 6) = I * (-p1 + p2) * h;
        c_array(1, 2, 6) = 2.0 * c_2_nu * h2 * d_e;
        p1 = 1.5 * d * d_e * c_d_3h;
        p2 = abs_d * (1.0 / 6.0 * c_1_2nu * abs_d_2 + h2 * (43.0 / 30.0 + c_2_nu / 3.0));
        c_array(1, 0, 7) = (p1 + p2) * h;
        c_array(1, 1, 7) = I * (-p1 + p2) * h;
        c_array(1, 2, 7) = 2.0 * h4 * d_e;
        p0 = h * c_d_h;
        p1 = 0.15 * abs_d * abs_d_2 - 1.9 / 6.0 * abs_d * h2;
        p2 = 0.25 * d * d_e * c_d_m3h;
        c_array(1, 0, 8) = -p0 * (p1 + p2);
        c_array(1, 1, 8) = I * p0 * (-p1 + p2);

        for (int k = 0; k < c_array.size(2); ++k) {
            for (int j = 0; j < c_array.size(1); ++j) {
                c_array(2, j, k) = std::conj(c_array(1, j, k));
            }
        }

        c_array(3, 2, 0) = I * c_1_2nu * d_e_2 * tcos_x;
        p0 = d * h;
        p1 = 0.0625 * c_13_10nu;
        p2 = d_e_2 * (0.0625 * c_13_10nu + 0.5 * c_3_2nu * tcos_x);
        c_array(3, 0, 1) = -I * p0 * (p1 - p2);
        c_array(3, 1, 1) = p0 * (p1 + p2);
        c_array(3, 2, 1) = 2.0 * I * c_2_nu * d_e_2 * h2 * tcos_x;
        //p1 = ; p2 = ;
        c_array(3, 0, 2) = d * h *
                     (0.09375 * I * c_7_2nu * h2 - 
                             0.03125 * c_3_2nu * abs_d_2 * tan_x +
                      d_e_2 *
                      (abs_d_2 * (-1.0 / 12.0 * I * c_7_6nu + 
                              0.09375 * c_3_2nu * tan_x) -
                       I * h2 * (0.09375 * c_7_2nu + 0.25 * c_4_nu * tcos_x)));
        c_array(3, 1, 2) = d * h *
                     (-0.09375 * c_7_2nu * h2 - 
                             I * 0.03125 * c_3_2nu * abs_d_2 * tan_x -
                      d_e_2 *
                      (abs_d_2 * (1.0 / 12.0 * c_7_6nu + 
                              0.09375 * I * c_3_2nu * tan_x) +
                       h2 * (0.09375 * c_7_2nu + 0.25 * c_4_nu * tcos_x)));
        c_array(3, 2, 2) = -I * d_e_2 * h4 * tcos_x;
        //p0 = ; p1 = ; p2 = ;
        c_array(3, 0, 3) = d * h * 
                (abs_d_2 * tan_x * ((0.09375 - 0.28125 * d_e_2) * h2 - 
                        (0.046875 + 0.109375 * d_e_2) * abs_d_2) +
                              I * h2 * (0.625 * d_e_2 * abs_d_2 - 
                                      0.234375 * h2 + d_e_2 * h2 *
                                        (0.234375 + 0.3125 * tcos_x)));
        c_array(3, 1, 3) = d * h * 
                (I * abs_d_2 * tan_x * ((0.09375 + 0.28125 * d_e_2) * h2 + 
                        (-0.046875 + 0.109375 * d_e_2) * abs_d_2) + 
                        h2 * (0.625 * d_e_2 * abs_d_2 + 0.234375 * h2 +
                              d_e_2 * h2 * (0.234375 + 0.3125 * tcos_x)));
        c_array(3, 0, 5) = c_3_2nu * h * d_e * (0.25 * d_e_2 - 0.75);
        c_array(3, 1, 5) = -I * c_3_2nu * h * d_e * (0.25 * d_e_2 + 0.75);
        c_array(3, 2, 5) = 0.5 * c_1_2nu * d * d_e;
        p0 = h * d_e;
        p1 = d_e_2 * (0.75 * c_5_4nu * abs_d_2 + 0.25 * c_11_4nu * h2);
        p2 = 0.25 * c_5_4nu * abs_d_2 + 0.75 * c_11_4nu * h2;
        c_array(3, 0, 6) = h * d_e * (p1 - p2);
        c_array(3, 1, 6) = -I * h * d_e * (p1 + p2);
        c_array(3, 2, 6) = 2.0 * c_2_nu * h2 * d * d_e;
        //p1 = ; p2 = ;
        c_array(3, 0, 7) = h * d_e * 
                (0.125 * c_1_2nu * abs_d_4 - 0.25 * c_7_2nu * abs_d_2 * h2 - 
                        0.375 * c_13_2nu * h4 + 
                        d_e_2 * (0.625 * c_1_2nu * abs_d_4 + 
                                0.75 * c_7_2nu * abs_d_2 * h2 + 
                                0.125 * c_13_2nu * h4));
        c_array(3, 1, 7) = I * h * d_e * 
                (0.125 * c_1_2nu * abs_d_4 - 0.25 * c_7_2nu * abs_d_2 * h2 - 
                        0.375 * c_13_2nu * h4 - 
                        d_e_2 * (0.625 * c_1_2nu * abs_d_4 +
                                 0.75 * c_7_2nu * abs_d_2 * h2 +
                                 0.125 * c_13_2nu * h4));
        c_array(3, 2, 7) = 2.0 * h4 * d * d_e;
        p0 = h * d_e * c_d_h;
        p1 = d_e_2 * (7.0 / 24.0 * abs_d_4 - 
                11.0 / 12.0 * abs_d_2 * h2 - 5.0 / 24.0 * h4);
        p2 = 0.125 * abs_d_4 - 0.25 * abs_d_2 * h2 + 0.625 * h4;
        c_array(3, 0, 8) = -p0 * (p1 + p2);
        c_array(3, 1, 8) = I * p0 * (p1 - p2);

        for (int k = 0; k < c_array.size(2); ++k) {
            for (int j = 0; j < c_array.size(1); ++j) {
                c_array(4, j, k) = std::conj(c_array(3, j, k));
            }
        }

        p0 = 0.125 * c_13_10nu * h;
        c_array(5, 0, 1) = p0 * std::imag(d);
        c_array(5, 1, 1) = -p0 * std::real(d);
        c_array(5, 2, 1) = 1.0 / 12.0 * c_1_2nu * abs_d_2 * tan_x;
        p1 = 0.0625 * c_3_2nu * abs_d_2 * tan_x;
        p2 = 0.1875 * c_7_2nu * h2;
        c_array(5, 0, 2) = -h * (p1 * std::real(d) + p2 * std::imag(d));
        c_array(5, 1, 2) = -h * (p1 * std::imag(d) - p2 * std::real(d));
        c_array(5, 2, 2) = -1.0 / 6.0 * c_2_nu * abs_d_2 * h2 * tan_x;
        p1 = 0.09375 * abs_d_2 * tan_x * (2 * h2 - abs_d_2);
        p2 = 0.46875 * h4;
        c_array(5, 0, 3) = h * (p1 * std::real(d) + p2 * std::imag(d));
        c_array(5, 1, 3) = h * (p1 * std::imag(d) - p2 * std::real(d));
        c_array(5, 2, 3) = 0.25 * abs_d_2 * h4 * tan_x;
        c_array(5, 2, 4) = -4.0 * c_1_nu * std::fabs(h);
        p0 = -1.5 * c_3_2nu * h;
        c_array(5, 0, 5) = p0 * std::real(d_e);
        c_array(5, 1, 5) = p0 * std::imag(d_e);
        c_array(5, 2, 5) = -0.5 * c_1_2nu * abs_d;
        p0 = -h * (0.5 * c_5_4nu * abs_d_2 + 1.5 * c_11_4nu * h2);
        c_array(5, 0, 6) = p0 * std::real(d_e);
        c_array(5, 1, 6) = p0 * std::imag(d_e);
        c_array(5, 2, 6) = (1.0 / 6.0 * c_1_2nu * abs_d_2 + 
                1.5 * h2 * (5 + 2 * nu)) * abs_d;
        p0 = h * (0.25 * c_1_2nu * abs_d_4 - 
                0.5 * c_7_2nu * abs_d_2 * h2 - 0.75 * c_13_2nu * h4);
        c_array(5, 0, 7) = p0 * std::real(d_e);
        c_array(5, 1, 7) = p0 * std::imag(d_e);
        c_array(5, 2, 7) = 2.0 / 3.0 * h2 * (c_2_nu * abs_d_2 + 
                h2 * (7 + nu)) * abs_d;
        p0 = -h * c_d_h * (0.25 * abs_d_4 - 0.5 * abs_d_2 * h2 + 1.25 * h4);
        c_array(5, 0, 8) = p0 * std::real(d_e);
        c_array(5, 1, 8) = p0 * std::imag(d_e);
        c_array(5, 2, 8) = 2.0 / 3.0 * h4 * c_d_h * abs_d;

        return c_array;
    }

    il::StaticArray3D<std::complex<double>, 6, 3, 9> s_12_h
            (double nu, std::complex<double> eix,
             double h, std::complex<double> d) {
        const std::complex<double> I(0.0, 1.0);

        double c_1_mnu = 1.0 - nu;
        double c_1_m2nu = 1.0 - 2.0 * nu;
        double c_2_mnu = 2.0 - nu;
        double c_3_mnu = 3.0 - nu;
        double c_3_m4nu = 3.0 - 4.0 * nu;
        double c_5_m2nu = 5.0 - 2.0 * nu;
        double c_5_m4nu = 5.0 - 4.0 * nu;
        double c_6_m5nu = 6.0 - 5.0 * nu;
        double c_7_m2nu = 7.0 - 2.0 * nu;
        double c_8_m5nu = 8.0 - 5.0 * nu;
        double c_9_m2nu = 9.0 - 2.0 * nu;
        double c_9_m4nu = 9.0 - 4.0 * nu;
        double c_9_m8nu = 9.0 - 8.0 * nu;
        double c_13_m2nu = 13.0 - 2.0 * nu;
        double c_15_m4nu = 15.0 - 4.0 * nu;
        double c_15_m8nu = 15.0 - 8.0 * nu;
        double c_115_m38nu_80 = 1.4375 - 0.475 * nu;

        double cos_x = std::real(eix);
        double tan_x = std::imag(eix) / cos_x;
        std::complex<double> tcos_x = cos_x * eix;
        std::complex<double> tcos_x_m1 = tcos_x - 1.0;
        std::complex<double> c_3_4tcos = 3.0 + 4.0 * tcos_x;
        std::complex<double> e2x = eix * eix;
        std::complex<double> c_tcos_n1 = tcos_x_m1 * tcos_x;
        std::complex<double> w_c_tcos_n2 = 
                (13.0 + e2x - 10.0 * tcos_x) * tcos_x;

        double h2 = h * h;
        double h4 = h2 * h2;
        double sgh = ((h < 0) ? -1.0 : double((h > 0))); // sign(h)

        double abs_d = std::abs(d);
        double abs_d_2 = abs_d * abs_d;
        double abs_d_4 = abs_d_2 * abs_d_2;
        std::complex<double> d_c = std::conj(d);
        std::complex<double> d_e = std::polar(1.0, std::arg(d)); //  = d/abs_d
        std::complex<double> d_e_c = std::conj(d_e);
        std::complex<double> d_e_2 = d_e * d_e; //  = d^2/abs_d^2
        std::complex<double> d_e_3 = d_e * d_e_2; //  = d^3/abs_d^3
        std::complex<double> d_e_4 = d_e_2 * d_e_2; //  = d^4/abs_d^4

        std::complex<double> d2h2 = abs_d_2 * h2;

        std::complex<double> c_d_h = abs_d_2 + h2;
        std::complex<double> c_d_3h = abs_d_2 + 3.0 * h2;
        std::complex<double> c_3d_h = 3.0 * abs_d_2 + h2;
        std::complex<double> c_d_m3h = abs_d_2 - 3.0 * h2;

        std::complex<double> p0, p1, p2;

        il::StaticArray3D<std::complex<double>, 6, 3, 9> c_array{0.0};

        c_array(0, 2, 1) = -I * c_1_m2nu * d_e_2 * tcos_x;
        c_array(0, 0, 2) = I * d * h * (-0.5 + d_e_2 * (0.5 + 0.75 * tcos_x));
        c_array(0, 1, 2) = d * h * (0.5 + d_e_2 * (0.5 + 0.75 * tcos_x));
        c_array(0, 2, 2) = I * h2 * d_e_2 * tcos_x;
        //p0 = d*h; p1 = ; p2 = ;
        c_array(0, 0, 3) = d * h * 
                (0.28125 * I * h2 - 0.09375 * abs_d_2 * tan_x + 
                        d_e_2 * (abs_d_2 * (-0.75 * I + 0.28125 * tan_x) -
                                 I * h2 * (0.28125 + 0.375 * tcos_x)));
        c_array(0, 1, 3) = -d * h * 
                (0.28125 * h2 + 0.09375 * I * abs_d_2 * tan_x + 
                        d_e_2 * (abs_d_2 * (0.75 + 0.28125 * I * tan_x) +
                                 h2 * (0.28125 + 0.375 * tcos_x)));
        p0 = 0.5 * d_e * h;
        p1 = 3.0 * d_e_2;
        c_array(0, 0, 6) = p0 * (-p1 + c_9_m4nu);
        c_array(0, 1, 6) = I * p0 * (p1 + c_9_m4nu);
        c_array(0, 2, 6) = -c_1_m2nu * d_e * d;
        p0 = 1.5 * h * d_e;
        p2 = c_3d_h * d_e_2;
        c_array(0, 0, 7) = p0 * (c_d_3h - p2);
        c_array(0, 1, 7) = I * p0 * (c_d_3h + p2);
        c_array(0, 2, 7) = -2.0 * h2 * d_e * d;
        p0 = 0.25 * h * d_e * c_d_h;
        p2 = d_e_2 * (5.0 * abs_d_2 + h2);
        c_array(0, 0, 8) = -p0 * (c_d_m3h + p2);
        c_array(0, 1, 8) = I * p0 * (-c_d_m3h + p2);

        p0 = 0.4 * d_e_2 * h * tcos_x;
        p2 = 8.0 * d_e_2 * (-1.0 + tcos_x);
        c_array(1, 0, 1) = I * p0 * (c_8_m5nu + p2);
        c_array(1, 1, 1) = p0 * (-c_8_m5nu + p2);
        c_array(1, 2, 1) = -I * c_1_m2nu * d * d_e_2 * (0.625 + tcos_x);
        p0 = d_e_2 * h; //p1 = ; p2 = ;
        c_array(1, 0, 2) = p0 * (abs_d_2 * (-0.7 * I + 0.2 * tan_x) - 
                1.6 * I * h2 * tcos_x - 
                2.0 * d_e_2 * (0.8 * I * h2 * c_tcos_n1 + 
                        abs_d_2 * (0.1 * tan_x - I * (0.7 + 0.4 * tcos_x))));
        c_array(1, 1, 2) = p0 * (abs_d_2 * (0.7 + 0.2 * I * tan_x) + 
                1.6 * h2 * tcos_x + 
                2.0 * d_e_2 * (-0.8 * h2 * c_tcos_n1 + 
                        abs_d_2 * (0.7 + 0.1 * I * tan_x + 0.4 * tcos_x)));
        c_array(1, 2, 2) = d * d_e_2 * 
                (-abs_d_2 * c_1_m2nu * (-0.5 * I + 0.1875 * tan_x) + 
                        I * h2 * (1.1875 + 1.75 * tcos_x - 
                                0.125 * nu * c_3_4tcos));
        //p1 = ; p2 = ;
        c_array(1, 0, 3) = p0 * (d2h2 * (0.4 * I - 0.11875 * tan_x) - 
                0.09375 * abs_d_4 * tan_x + 0.4 * I * h4 * tcos_x + 
                d_e_2 * (3.0 * abs_d_4 * (-0.4 * I + 0.18125 * tan_x) + 
                        0.4 * I * h4 * c_tcos_n1 + 
                        d2h2 * (0.11875 * tan_x - 0.4 * I * (2.0 + tcos_x))));
        c_array(1, 1, 3) = -p0 * (d2h2 * (0.4 + 0.11875 * I * tan_x) + 
                0.09375 * I * abs_d_4 * tan_x + 0.4 * h4 * tcos_x + 
                d_e_2 * (3.0 * abs_d_4 * (0.4 + 0.18125 * I * tan_x) - 
                        0.4 * h4 * c_tcos_n1 + 
                        d2h2 * (0.11875 * I * tan_x + 0.4 * (2.0 + tcos_x))));
        c_array(1, 2, 3) = -I * d * d_e_2 * h2 * 
                (abs_d_2 * (1.5 + 0.5625 * I * tan_x) + 
                        0.1875 * h2 * c_3_4tcos);
        c_array(1, 2, 5) = -0.5 * c_1_m2nu * d_e_3;
        p0 = d * d_e * h;
        c_array(1, 0, 6) = -p0 * (4.5 * d_e_2 - 0.5 * c_9_m4nu);
        c_array(1, 1, 6) = I * p0 * (4.5 * d_e_2 + 0.5 * c_9_m4nu);
        c_array(1, 2, 6) = -d_e_3 * (2.0 * c_2_mnu * h2 + 
                3.0 * c_1_m2nu * abs_d_2);
        p1 = 1.5 * abs_d_2 + 4.5 * h2;
        p2 = d_e_2 * (7.5 * abs_d_2 + 4.5 * h2);
        c_array(1, 0, 7) = p0 * (p1 - p2);
        c_array(1, 1, 7) = I * p0 * (p1 + p2);
        c_array(1, 2, 7) = -0.25 * d_e_3 * 
                (5.0 * abs_d_4 * c_1_m2nu + 
                6.0 * d2h2 * c_7_m2nu + h4 * c_13_m2nu);
        p0 = 0.25 * p0 * c_d_h;
        p1 = abs_d_2 - 3.0 * h2;
        p2 = d_e_2 * (7.0 * abs_d_2 + 3.0 * h2);
        c_array(1, 0, 8) = -p0 * (p1 + p2);
        c_array(1, 1, 8) = I * p0 * (-p1 + p2);
        c_array(1, 2, 8) = -0.5 * d_e_3 * h2 * c_d_h * (5.0 * abs_d_2 + h2);

        c_array(2, 0, 1) = 3.2 * I * h * d_e_2 * tcos_x;
        c_array(2, 1, 1) = 3.2 * h * d_e_2 * tcos_x;
        c_array(2, 2, 1) = 0.625 * I * c_1_m2nu * d;
        p1 = 0.1 * d_e_2 * 
                (abs_d_2 * (7.0 + 2.0 * I * tan_x) + 16.0 * h2 * tcos_x);
        p2 = I * c_6_m5nu / 30.0 * abs_d_2 * tan_x;
        c_array(2, 0, 2) = -I * h * (p1 - p2);
        c_array(2, 1, 2) = -h * (p1 + p2);
        c_array(2, 2, 2) = d * 
                (0.0625 * c_1_m2nu * abs_d_2 * tan_x - 
                        I * h2 * (1.1875 - 0.375 * nu));
        p1 = abs_d_2 * tan_x * (-0.05625 * abs_d_2 + 0.11875 * h2);
        p2 = d_e_2 * (abs_d_2 * tan_x * (0.11875 * h2 + 0.09375 * abs_d_2) -
                    0.4 * I * h2 * (abs_d_2 + h2 * tcos_x));
        c_array(2, 0, 3) = h * (p1 - p2);
        c_array(2, 1, 3) = I * h * (p1 + p2);
        c_array(2, 2, 3) = d * h2 * 
                (-0.1875 * abs_d_2 * tan_x + 0.5625 * I * h2);
        c_array(2, 0, 4) = -2.0 * c_1_mnu * sgh;
        c_array(2, 1, 4) = I * c_array(2, 0, 4);
        c_array(2, 2, 5) = 1.5 * c_1_m2nu * d_e;
        p1 = 4.5 * d * d_e * h;
        p2 = abs_d * h * (4.3 - 8.0 / 3.0 * nu);
        c_array(2, 0, 6) = p1 + p2;
        c_array(2, 1, 6) = I * (-p1 + p2);
        c_array(2, 2, 6) = d_e * (6.0 * c_2_mnu * h2 + c_1_m2nu * abs_d_2);
        p1 = 1.5 * d * d_e * h * c_d_3h;
        p2 = 1.0 / 3.0 * abs_d * h * 
                ((2.0 * c_2_mnu + 3.3) * h2 + 0.5 * c_3_m4nu * abs_d_2);
        c_array(2, 0, 7) = p1 + p2;
        c_array(2, 1, 7) = I * (-p1 + p2);
        c_array(2, 2, 7) = d_e * 
                (3.0 * h2 * c_d_3h - 0.25 * c_1_m2nu * c_d_m3h * c_d_h);
        p0 = h * abs_d * c_d_h;
        p1 = 0.15 * abs_d_2 - 0.95 / 3.0 * h2;
        p2 = 0.25 * d_e_2 * c_d_m3h;
        c_array(2, 0, 8) = -p0 * (p1 + p2);
        c_array(2, 1, 8) = I * p0 * (-p1 + p2);
        c_array(2, 2, 8) = -0.5 * d_e * h2 * c_d_m3h * c_d_h;

        c_array(3, 2, 0) = 6.4 / 3.0 * I * d_e_4 * c_1_m2nu * c_tcos_n1;
        p0 = d_e_2 * d * h;
        p1 = d_e_2 * (1.4625 + 0.375 * w_c_tcos_n2);
        p2 = 1.25 * (0.15 + c_1_mnu) + 0.5 * c_5_m4nu * tcos_x;
        c_array(3, 0, 1) = -I * p0 * (p1 - p2);
        c_array(3, 1, 1) = -p0 * (p1 + p2);
        c_array(3, 2, 1) = I * d_e_4 * 
                (12.8 / 3.0 * c_2_mnu * c_tcos_n1 * h2 - 
                        c_1_m2nu / 3.0 * (5.2 + 0.725 * I * tan_x + 
                                3.2 * tcos_x) * abs_d_2);
        //p1 = ; p2 = ;
        c_array(3, 0, 2) = I * p0 * (d_e_2 * (abs_d_2 * (3.275 + 
                tan_x * (0.78125 * I + 0.025 * tan_x) + tcos_x) + 
                h2 * (0.86875 + 0.1875 * w_c_tcos_n2)) - 
                abs_d_2 * (c_1_mnu + 0.25 / 3.0 + 
                        0.09375 * I * c_5_m4nu * tan_x) - 
                h2 * (0.0625 * c_5_m2nu * c_3_4tcos - 0.09375));
        c_array(3, 1, 2) = p0 * (d_e_2 * (abs_d_2 * (3.275 + 
                tan_x * (0.78125 * I + 0.025 * tan_x) + tcos_x) + 
                h2 * (0.86875 + 0.1875 * w_c_tcos_n2)) + 
                abs_d_2 * (c_1_mnu + 0.25 / 3.0 + 
                        0.09375 * I * c_5_m4nu * tan_x) + 
                h2 * (0.0625 * c_5_m2nu * c_3_4tcos - 0.09375));
        c_array(3, 2, 2) = d_e_4 * 
                (-abs_d_4 * c_1_m2nu * (-0.8 * I + 0.3625 * tan_x) - 
                        0.8 / 3.0 * I * h4 * c_13_m2nu * c_tcos_n1 + 
                        d2h2 / 3.0 * (-c_115_m38nu_80 * tan_x + 
                                0.4 * I * (1.0 + 8.0 * c_3_mnu + 
                                        2.0 * c_7_m2nu * tcos_x)));
        //p1 = ; p2 = ;
        c_array(3, 0, 3) = p0 * ((d2h2 * (0.625 * I - 0.28125 * tan_x) - 
                0.109375 * abs_d_4 * tan_x + 
                0.078125 * I * h4 * c_3_4tcos) + 
                d_e_2 * (abs_d_4 * (-2.0 * I + 1.015625 * tan_x) - 
                         I * h4 * (0.234375 + 0.046875 * w_c_tcos_n2) + 
                         d2h2 * (0.46875 * tan_x - 
                                 0.125 * I * (15.0 + 4.0 * tcos_x))));
        c_array(3, 1, 3) = -p0 * ((d2h2 * (0.625 + 0.28125 * I * tan_x) + 
                0.109375 * I * abs_d_4 * tan_x + 0.078125 * h4 * c_3_4tcos) + 
                d_e_2 * (abs_d_4 * (2.0 + 1.015625 * I * tan_x) + 
                        h4 * (0.234375 + 0.046875 * w_c_tcos_n2) + 
                        d2h2 * (0.46875 * I * tan_x + 
                                0.125 * (15.0 + 4.0 * tcos_x))));
        c_array(3, 2, 3) = d_e_4 * h2 * 
                (3.0 * abs_d_4 * (-0.8 * I + 0.3625 * tan_x) + 
                        0.8 * I * h4 * c_tcos_n1 + 
                        d2h2 * (0.2375 * tan_x - 0.8 * I * (2.0 + tcos_x)));
        p0 = 0.25 * h * d_e_3;
        c_array(3, 0, 5) = p0 * (-3.0 * d_e_2 + c_5_m4nu);
        c_array(3, 1, 5) = I * p0 * (3.0 * d_e_2 + c_5_m4nu);
        c_array(3, 2, 5) = -1.5 * d * d_e_3 * c_1_m2nu;
        p1 = 3.0 * abs_d_2 * c_9_m8nu + h2 * c_15_m8nu;
        p2 = 9.0 * d_e_2 * (5.0 * abs_d_2 + h2);
        c_array(3, 0, 6) = p0 * (p1 - p2);
        c_array(3, 1, 6) = I * p0 * (p1 + p2);
        c_array(3, 2, 6) = -d * d_e_3 * 
                (6.0 * h2 * c_2_mnu + 5.0 * abs_d_2 * c_1_m2nu);
        p0 = 0.5 * p0;
        p1 = 5.0 * abs_d_4 * c_3_m4nu + 6.0 * d2h2 * c_9_m4nu + h4 * c_15_m4nu;
        p2 = 3.0 * d_e_2 * (35.0 * abs_d_4 + 30.0 * d2h2 + 3.0 * h4);
        c_array(3, 0, 7) = p0 * (p1 - p2);
        c_array(3, 1, 7) = I * p0 * (p1 + p2);
        c_array(3, 2, 7) = -d * d_e_3 *
                     (0.75 * h4 * c_13_m2nu + 2.5 * d2h2 * c_7_m2nu + 
                             1.75 * abs_d_4 * c_1_m2nu);
        p0 = p0 * c_d_h;
        p1 = (-7.0 * abs_d_4 + 22.0 * d2h2 + 5.0 * h4) / 3.0;
        p2 = d_e_2 * (21.0 * abs_d_4 + 14.0 * d2h2 + h4);
        c_array(3, 0, 8) = p0 * (p1 - p2);
        c_array(3, 1, 8) = I * p0 * (p1 + p2);
        c_array(3, 2, 8) = -d * d_e_3 * h2 * c_d_h * (3.5 * abs_d_2 + 1.5 * h2);

        p1 = 1.25 * nu * d_c;
        c_array(4, 0, 1) = h * (2.875 * std::imag(d) - I * p1);
        c_array(4, 1, 1) = h * (p1 - 2.875 * std::real(d));
        c_array(4, 2, 1) = 0.725 / 3.0 * c_1_m2nu * abs_d_2 * tan_x;
        p0 = 0.0625 * h;
        p1 = 6.0 * nu * I * h2 - c_5_m2nu * abs_d_2 * tan_x;
        p2 = I * (3.0 * c_9_m2nu * I * h2 - 2.0 * nu * abs_d_2 * tan_x);
        c_array(4, 0, 2) = p0 * (p1 * std::real(d) + p2 * std::imag(d));
        c_array(4, 1, 2) = p0 * (p1 * std::imag(d) - p2 * std::real(d));
        c_array(4, 2, 2) = abs_d_2 * tan_x * 
                (0.0375 * c_1_m2nu * abs_d_2 - c_115_m38nu_80 / 3.0 * h2);
        p0 = 0.09375 * h;
        p1 = 5.0 * h4;
        p2 = (-abs_d_4 + 2.0 * d2h2) * tan_x;
        c_array(4, 0, 3) = p0 * (p1 * std::imag(d) + p2 * std::real(d));
        c_array(4, 1, 3) = p0 * (-p1 * std::real(d) + p2 * std::imag(d));
        c_array(4, 2, 3) = abs_d_2 * h2 * tan_x * 
                (0.2375 * h2 - 0.1125 * abs_d_2);
        c_array(4, 2, 4) = -8.0 * c_1_mnu * std::fabs(h);
        p0 = 1.5 * h;
        p2 = 2.0 * I * nu;
        c_array(4, 0, 5) = p0 * 
                (-c_5_m2nu * std::real(d_e) - p2 * std::imag(d_e));
        c_array(4, 1, 5) = p0 * 
                (-c_5_m2nu * std::imag(d_e) + p2 * std::real(d_e));
        c_array(4, 2, 5) = -1.5 * c_1_m2nu * abs_d;
        p1 = 2.0 * h * c_d_3h * nu * d_e_c;
        p2 = 4.5 * h * (abs_d_2 + 5.0 * h2);
        c_array(4, 0, 6) = p1 - p2 * std::real(d_e);
        c_array(4, 1, 6) = I * p1 - p2 * std::imag(d_e);
        c_array(4, 2, 6) = abs_d * 
                (1.0 / 3.0 * c_1_m2nu * abs_d_2 + (3.6 * c_3_mnu - 0.4) * h2);
        p1 = 0.5 * nu * h * c_d_h * c_d_m3h * d_e_c;
        p2 = 0.75 * h * (abs_d_4 - 6.0 * d2h2 - 15.0 * h4);
        c_array(4, 0, 7) = -p1 + p2 * std::real(d_e);
        c_array(4, 1, 7) = -I * p1 + p2 * std::imag(d_e);
        c_array(4, 2, 7) = abs_d * (-0.15 * c_1_m2nu * abs_d_4 + 
                1.0 / 6.0 * c_7_m2nu * d2h2 + 
                (0.35 + 1.9 * (8 - nu)) / 3.0 * h4);
        p2 = -0.25 * h * c_d_h * (abs_d_4 - 2.0 * d2h2 + 5.0 * h4);
        c_array(4, 0, 8) = p2 * std::real(d_e);
        c_array(4, 1, 8) = p2 * std::imag(d_e);
        c_array(4, 2, 8) = abs_d * h2 * c_d_h * 
                (1.9 / 3.0 * h2 - 0.3 * abs_d_2);

        c_array(5, 2, 0) = 6.4 / 3.0 * I * c_1_m2nu * d_e_2 * tcos_x;
        p0 = d * h;
        p1 = 0.1875 + 1.25 * c_1_mnu;
        p2 = d_e_2 * (1.4375 + 2.5 * tcos_x);
        c_array(5, 0, 1) = I * p0 * (-p1 + p2);
        c_array(5, 1, 1) = p0 * (p1 + p2);
        c_array(5, 2, 1) = d_e_2 / 3.0 * 
                (c_1_m2nu * abs_d_2 * (2.6 * I - 0.725 * tan_x) + 
                        12.8 * c_2_mnu * I * h2 * tcos_x);
        //p1 = ; p2 = ;
        c_array(5, 0, 2) = p0 * (0.09375 * I * h2 * c_9_m4nu - 
                0.03125 * abs_d_2 * c_5_m4nu * tan_x + 
                d_e_2 * (abs_d_2 * (-3.25 / 3.0 * I + 0.46875 * tan_x) - 
                         0.3125 * I * h2 * (2.7 + 4.0 * tcos_x)));
        c_array(5, 1, 2) = -p0 * (0.09375 * h2 * c_9_m4nu + 
                0.03125 * I * abs_d_2 * c_5_m4nu * tan_x + 
                d_e_2 * (abs_d_2 * (3.25 / 3.0 + 0.46875 * I * tan_x) + 
                        0.3125 * h2 * (2.7 + 4.0 * tcos_x)));
        c_array(5, 2, 2) = d_e_2 * (0.0625 * c_1_m2nu * tan_x * abs_d_4 + 
                (I * (-5.0 + 1.6 * nu) + c_115_m38nu_80 * tan_x) / 3.0 * d2h2 - 
                0.8 / 3.0 * I * c_13_m2nu * tcos_x * h4);
        //p1 = ; p2 = ;
        c_array(5, 0, 3) = p0 * (-(0.234375 * I * h4 - 0.09375 * d2h2 * tan_x + 
                0.046875 * abs_d_4 * tan_x) + 
                d_e_2 * (0.078125 * I * h4 * c_3_4tcos + 
                        (0.625 * I - 0.28125 * tan_x) * d2h2 - 
                        0.109375 * abs_d_4 * tan_x));
        c_array(5, 1, 3) = p0 * ((0.234375 * h4 + 0.09375 * I * d2h2 * tan_x - 
                0.046875 * I * abs_d_4 * tan_x) + 
                d_e_2 * (0.078125 * h4 * c_3_4tcos + 
                        (0.625 + 0.28125 * I * tan_x) * d2h2 + 
                        0.109375 * I * abs_d_4 * tan_x));
        c_array(5, 2, 3) = d_e_2 * h2 * (-0.1875 * abs_d_4 * tan_x + 
                (0.8 * I - 0.2375 * tan_x) * d2h2 + 0.8 * I * h4 * tcos_x);
        p0 = d_e * h;
        p1 = 1.25 * d_e_2;
        p2 = 3.75 - 3.0 * nu;
        c_array(5, 0, 5) = p0 * (p1 - p2);
        c_array(5, 1, 5) = -I * p0 * (p1 + p2);
        c_array(5, 2, 5) = 1.5 * d * d_e * c_1_m2nu;
        p0 = 0.25 * p0;
        p1 = d_e_2 * (15.0 * h2 + 27.0 * abs_d_2);
        p2 = 3.0 * c_15_m8nu * h2 + c_9_m8nu * abs_d_2;
        c_array(5, 0, 6) = p0 * (p1 - p2);
        c_array(5, 1, 6) = -I * p0 * (p1 + p2);
        c_array(5, 2, 6) = d * d_e * (c_1_m2nu * abs_d_2 + 6.0 * c_2_mnu * h2);
        p0 = 0.5 * p0;
        p1 = c_3_m4nu * abs_d_4 - 2.0 * c_9_m4nu * d2h2 - 3.0 * c_15_m4nu * h4;
        p2 = 3.0 * d_e_2 * (5.0 * abs_d_4 + 18.0 * d2h2 + 5.0 * h4);
        c_array(5, 0, 7) = p0 * (p1 + p2);
        c_array(5, 1, 7) = -I * p0 * (-p1 + p2);
        c_array(5, 2, 7) = 0.25 * d * d_e *
                     (-c_1_m2nu * abs_d_4 + 2.0 * c_7_m2nu * d2h2 + 
                             3.0 * c_13_m2nu * h4);
        p0 = p0 * c_d_h;
        p1 = (-abs_d_4 + 2.0 * d2h2 - 5.0 * h4);
        p2 = d_e_2 / 3.0 * (-7.0 * abs_d_4 + 22.0 * d2h2 + 5.0 * h4);
        c_array(5, 0, 8) = p0 * (p1 + p2);
        c_array(5, 1, 8) = I * p0 * (p1 - p2);
        c_array(5, 2, 8) = -0.5 * d * d_e * h2 * c_d_h * c_d_m3h;

        return c_array;
    }

    il::StaticArray3D<std::complex<double>, 6, 3, 9> s_13_23_h
            (double nu, std::complex<double> eix,
             double h, std::complex<double> d) {
        const std::complex<double> I(0.0, 1.0);

        double c_1_nu = 1.0 + nu;
        double c_2_mnu = 2.0 - nu;
        double c_3_nu = 3.0 + nu;
        double c_3_mnu = 3.0 - nu;
        double c_5_mnu = 5.0 - nu;
        double c_6_nu = 6.0 + nu;
        double c_12_nu = 12.0 + nu;

        double cos_x = std::real(eix);
        double tan_x = std::imag(eix) / cos_x;
        std::complex<double> c_8_3i_tan_x = 8.0 + 3.0 * I * tan_x;
        std::complex<double> tcos_x = cos_x * eix;
        std::complex<double> tcos_x_c = std::conj(tcos_x);
        std::complex<double> tcos_x_m1 = tcos_x - 1.0;
        std::complex<double> c_3_4_tcos_x = 3.0 + 4.0 * tcos_x;
        std::complex<double> c_5_8_tcos_x = 5.0 + 8.0 * tcos_x;
        std::complex<double> c_tcos_n1 = tcos_x_m1 * tcos_x;
//  std::complex<double> e2x = eix*eix;
//  std::complex<double> w_c_tcos_n2 = (13.0+e2x-10.0*tcos_x)*tcos_x;

        double h2 = h * h;
        double h3 = h2 * h;
        double h4 = h2 * h2;
//  double sgh = (( h < 0 )? -1.0 : double(( h > 0 ))); // sign(h)

        double abs_d = std::abs(d);
        double abs_d_2 = abs_d * abs_d;
        double abs_d_4 = abs_d_2 * abs_d_2;
        std::complex<double> d_c = std::conj(d);
        std::complex<double> d_e = std::polar(1.0, std::arg(d)); // =d/abs_d
        std::complex<double> d_e_c = std::conj(d_e);
        std::complex<double> d_e_2 = d_e * d_e; // =d^2/abs_d^2
        std::complex<double> d_e_2_c = std::conj(d_e_2);
        std::complex<double> d_e_3 = d_e * d_e_2; // =d^3/abs_d^3
        std::complex<double> d_e_4 = d_e_2 * d_e_2; // =d^4/abs_d^4

        std::complex<double> d2h2 = abs_d_2 * h2;

        std::complex<double> c_d_h = abs_d_2 + h2;
//  std::complex<double> c_d_3h = abs_d_2+3.0*h2;
//  std::complex<double> c_3d_h = 3.0*abs_d_2+h2;
        std::complex<double> c_d_m3h = abs_d_2 - 3.0 * h2;
//  std::complex<double> D3mh1 = 3.0*abs_d_2-h2;

        std::complex<double> p0, p1, p2;

        il::StaticArray3D<std::complex<double>, 6, 3, 9> c_array{0.0};

        c_array(0, 1, 1) = -0.5 * nu * d_e_2 * tcos_x;
        c_array(0, 0, 1) = I * c_array(0, 1, 1);
        c_array(0, 1, 2) = 0.5 * d_e_2 * h2 * tcos_x;
        c_array(0, 0, 2) = I * c_array(0, 1, 2);
        c_array(0, 0, 6) = -0.5 * (c_2_mnu * abs_d + nu * d * d_e);
        c_array(0, 1, 6) = 0.5 * I * (-c_2_mnu * abs_d + nu * d * d_e);
        c_array(0, 2, 6) = -d_e * h;
        c_array(0, 0, 7) = -h2 * abs_d * (d_e_2 + 1.0);
        c_array(0, 1, 7) = I * h2 * abs_d * (d_e_2 - 1.0);
        c_array(0, 2, 7) = -2.0 * h3 * d_e;

        c_array(1, 1, 1) = -0.0625 * d * d_e_2 * nu * c_5_8_tcos_x;
        c_array(1, 0, 1) = I * c_array(1, 1, 1);
        c_array(1, 2, 1) = -I * d_e_2 * h * tcos_x;
        c_array(1, 1, 2) = d * d_e_2 * (0.03125 * nu * abs_d_2 * c_8_3i_tan_x +
                h2 * (0.5 + 0.09375 * nu + 0.125 * (6.0 + nu) * tcos_x));
        c_array(1, 0, 2) = I * c_array(1, 1, 2);
        c_array(1, 2, 2) = I * d_e_2 * h3 * tcos_x;
        c_array(1, 1, 3) = -0.09375 * d * d_e_2 * h2 *
                (abs_d_2 * c_8_3i_tan_x + h2 * c_3_4_tcos_x);
        c_array(1, 0, 3) = I * c_array(1, 1, 3);
        c_array(1, 0, 5) = 0.25 * d_e * (c_2_mnu - nu * d_e_2);
        c_array(1, 1, 5) = 0.25 * I * d_e * (c_2_mnu + nu * d_e_2);
        p2 = 3.0 * nu * abs_d_2 + c_3_nu * h2;
        c_array(1, 0, 6) = 0.5 * d_e * (c_5_mnu * h2 - d_e_2 * p2);
        c_array(1, 1, 6) = 0.5 * I * d_e * (c_5_mnu * h2 + d_e_2 * p2);
        c_array(1, 2, 6) = -d * d_e * h;
        p2 = d_e_2 * (0.625 * nu * abs_d_4 + 0.75 * c_6_nu * d2h2 +
                0.125 * c_12_nu * h4);
        c_array(1, 0, 7) = d_e * (h4 - p2);
        c_array(1, 1, 7) = I * d_e * (h4 + p2);
        c_array(1, 2, 7) = -2.0 * d * d_e * h3;
        c_array(1, 0, 8) = -0.25 * d_e_3 * h2 * c_d_h * (5.0 * abs_d_2 + h2);
        c_array(1, 1, 8) = -I * c_array(1, 0, 8);

        c_array(2, 1, 1) = 0.3125 * nu * d;
        c_array(2, 0, 1) = I * c_array(2, 1, 1);
        c_array(2, 1, 2) = -0.03125 * d * (h2 * (16 + 3.0 * nu) +
                I * nu * abs_d_2 * tan_x);
        c_array(2, 0, 2) = I * c_array(2, 1, 2);
        c_array(2, 2, 2) = abs_d_2 / 12.0 * h * tan_x;
        c_array(2, 1, 3) = 0.09375 * d * h2 * (3.0 * h2 + I * abs_d_2 * tan_x);
        c_array(2, 0, 3) = I * c_array(2, 1, 3);
        c_array(2, 2, 3) = -0.25 * abs_d_2 * h3 * tan_x;
        p1 = 3.0 * nu * d_e;
        p2 = c_2_mnu * d_e_c;
        c_array(2, 0, 5) = 0.25 * (p1 + p2);
        c_array(2, 1, 5) = -0.25 * I * (p1 - p2);
        p1 = nu * abs_d_2 + 3.0 * c_3_nu * h2;
        p2 = c_5_mnu * h2;
        c_array(2, 0, 6) = 0.5 * (p1 * d_e + p2 * d_e_c);
        c_array(2, 1, 6) = -0.5 * I * (p1 * d_e - p2 * d_e_c);
        c_array(2, 2, 6) = -10.0 / 3.0 * abs_d * h;
        p1 = 0.125 * (nu * abs_d_4 - 2.0 * c_6_nu * d2h2 - 3.0 * c_12_nu * h4);
        c_array(2, 0, 7) = -p1 * d_e + h4 * d_e_c;
        c_array(2, 1, 7) = I * (p1 * d_e + h4 * d_e_c);
        c_array(2, 2, 7) = -abs_d * h * (abs_d_2 + 11.0 * h2) / 3.0;
        c_array(2, 0, 8) = -0.25 * d_e * h2 * c_d_m3h * c_d_h;
        c_array(2, 1, 8) = -I * c_array(2, 0, 8);
        c_array(2, 2, 8) = -2.0 / 3.0 * abs_d * h3 * c_d_h;

        p1 = 0.5 * c_2_mnu * tcos_x;
        p2 = 3.2 / 3.0 * nu * d_e_2 * c_tcos_n1;
        c_array(3, 0, 0) = I * d_e_2 * (p1 + p2);
        c_array(3, 1, 0) = d_e_2 * (-p1 + p2);
        p1 = 0.5 * c_5_mnu * h2 * tcos_x;
        p2 = d_e_2 * (-3.2 / 3.0 * c_3_nu * h2 * c_tcos_n1 +
                    nu * abs_d_2 / 3.0 * (2.6 + 0.3625 * I * tan_x +
                            1.6 * tcos_x));
        c_array(3, 0, 1) = I * d_e_2 * (p1 - p2);
        c_array(3, 1, 1) = -d_e_2 * (p1 + p2);
        c_array(3, 2, 1) = -0.125 * I * d * d_e_2 * h * c_5_8_tcos_x;
        p1 = d_e_2 * (abs_d_4 * nu * (0.4 + 0.18125 * I * tan_x) -
                    h4 * 0.4 / 3.0 * c_12_nu * c_tcos_n1 +
                    d2h2 * (1.4 + 0.8 / 3.0 * nu +
                            0.4 / 3.0 * c_6_nu * tcos_x +
                            I * (0.2 + 0.11875 / 3.0 * nu) * tan_x));
        p2 = 0.5 * tcos_x * h4;
        c_array(3, 0, 2) = I * d_e_2 * (p1 - p2);
        c_array(3, 1, 2) = d_e_2 * (p1 + p2);
        c_array(3, 2, 2) = 0.0625 * I * d * d_e_2 * h *
                     (abs_d_2 * c_8_3i_tan_x + h2 * (19.0 + 28.0 * tcos_x));
        p0 = d_e_4 * h2; //p1 = ; p2 = ;
        c_array(3, 0, 3) = p0 * (abs_d_4 * (-1.2 * I + 0.54375 * tan_x) +
                0.4 * I * h4 * c_tcos_n1 +
                d2h2 * (0.11875 * tan_x - 0.4 * I * (2.0 + tcos_x)));
        c_array(3, 1, 3) = p0 * (abs_d_4 * (-1.2 - 0.54375 * I * tan_x) +
                0.4 * h4 * c_tcos_n1 +
                d2h2 * (-0.11875 * I * tan_x - 0.4 * (2.0 + tcos_x)));
        c_array(3, 2, 3) = -0.1875 * I * d * d_e_2 * h3 *
                (abs_d_2 * c_8_3i_tan_x + h2 * c_3_4_tcos_x);
        p1 = 0.25 * c_2_mnu * d * d_e;
        p2 = 0.75 * nu * d * d_e_3;
        c_array(3, 0, 5) = p1 - p2;
        c_array(3, 1, 5) = I * (p1 + p2);
        c_array(3, 2, 5) = -0.5 * d_e_3 * h;
        p1 = 0.5 * c_5_mnu * d * d_e * h2;
        p2 = 0.5 * d * d_e_3 * (5.0 * nu * abs_d_2 + 3.0 * c_3_nu * h2);
        c_array(3, 0, 6) = p1 - p2;
        c_array(3, 1, 6) = I * (p1 + p2);
        c_array(3, 2, 6) = -d_e_3 * h * (3.0 * abs_d_2 + 4.0 * h2);
        p1 = d * d_e * h4;
        p2 = 0.125 * d * d_e_3 *
             (7.0 * nu * abs_d_4 + 10.0 * c_6_nu * d2h2 + 3.0 * c_12_nu * h4);
        c_array(3, 0, 7) = p1 - p2;
        c_array(3, 1, 7) = I * (p1 + p2);
        c_array(3, 2, 7) = -0.25 * d_e_3 * h *
                (5.0 * abs_d_4 + 42.0 * d2h2 + 13.0 * h4);
        c_array(3, 0, 8) = -0.25 * d * d_e_3 * h2 * c_d_h *
                (7.0 * abs_d_2 + 3.0 * h2);
        c_array(3, 1, 8) = -I * c_array(3, 0, 8);
        c_array(3, 2, 8) = -0.5 * d_e_3 * h3 * c_d_h * (5.0 * abs_d_2 + h2);

        c_array(4, 1, 0) = 0.5 * c_2_mnu * d_e_2_c * tcos_x_c;
        c_array(4, 0, 0) = -I * c_array(4, 1, 0);
        p1 = 0.3625 / 3.0 * nu * abs_d_2 * tan_x;
        p2 = 0.5 * c_5_mnu * h2 * d_e_2_c * tcos_x_c;
        c_array(4, 0, 1) = p1 - I * p2;
        c_array(4, 1, 1) = -I * p1 + p2;
        p1 = (0.01875 * nu * abs_d_2 -
                h2 * (0.2 + 0.11875 / 3.0 * nu)) * abs_d_2 * tan_x;
        p2 = 0.5 * h4 * d_e_2_c * tcos_x_c;
        c_array(4, 0, 2) = p1 + I * p2;
        c_array(4, 1, 2) = -I * p1 - p2;
        c_array(4, 0, 3) = abs_d_2 * h2 *
                (-0.05625 * abs_d_2 + 0.11875 * h2) * tan_x;
        c_array(4, 1, 3) = -I * c_array(4, 0, 3);
        c_array(4, 0, 4) = -2.0 * c_1_nu * std::fabs(h);
        c_array(4, 1, 4) = -I * c_array(4, 0, 4);
        p1 = 0.75 * nu * abs_d;
        p2 = 0.25 * c_2_mnu * d_c * d_e_c;
        c_array(4, 0, 5) = -p1 + p2;
        c_array(4, 1, 5) = I * (p1 + p2);
        p1 = abs_d * (0.5 / 3.0 * nu * abs_d_2 + h2 * (4.3 + 0.9 * nu));
        p2 = 0.5 * c_5_mnu * h2 * d_c * d_e_c;
        c_array(4, 0, 6) = p1 + p2;
        c_array(4, 1, 6) = -I * (p1 - p2);
        p1 = abs_d * (-0.075 * nu * abs_d_4 + 0.25 / 3.0 * c_6_nu * d2h2 +
                   h4 / 3.0 * (7.3 + 0.475 * nu));
        p2 = h4 * d_c * d_e_c;
        c_array(4, 0, 7) = p1 + p2;
        c_array(4, 1, 7) = -I * (p1 - p2);
        c_array(4, 0, 8) = abs_d * h2 * c_d_h *
                (-0.15 * abs_d_2 + 0.95 / 3.0 * h2);
        c_array(4, 1, 8) = -I * c_array(4, 0, 8);

        c_array(5, 1, 0) = 3.2 / 3.0 * nu * d_e_2 * tcos_x;
        c_array(5, 0, 0) = I * c_array(5, 1, 0);
        p1 = 0.125 / 3.0 * c_2_mnu * abs_d_2 * tan_x;
        p2 = d_e_2 / 3.0 * (nu * abs_d_2 * (1.3 + 0.3625 * I * tan_x) +
                3.2 * c_3_nu * h2 * tcos_x);
        c_array(5, 0, 1) = p1 + I * p2;
        c_array(5, 1, 1) = I * p1 + p2;
        p1 = 0.125 / 3.0 * c_5_mnu * abs_d_2 * h2 * tan_x;
        p2 = d_e_2 * (0.03125 * nu * abs_d_4 * tan_x +
                d2h2 / 3.0 * (-I * (2.1 + 0.4 * nu) +
                        (0.6 + 0.11875 * nu) * tan_x) -
                0.4 / 3.0 * I * h4 * c_12_nu * tcos_x);
        c_array(5, 0, 2) = -p1 + p2;
        c_array(5, 1, 2) = -I * (p1 + p2);
        p1 = 0.125 * abs_d_2 * h4 * tan_x;
        p2 = d_e_2 * h2 * (0.4 * tcos_x * h4 +
                (0.4 + 0.11875 * I * tan_x) * d2h2 +
                0.09375 * I * tan_x * abs_d_4);
        c_array(5, 0, 3) = p1 + I * p2;
        c_array(5, 1, 3) = I * p1 + p2;
        c_array(5, 0, 4) = -c_3_mnu * std::fabs(h);
        c_array(5, 1, 4) = I * c_array(5, 0, 4);
        p1 = 0.25 * c_2_mnu * abs_d;
        p2 = 0.75 * nu * d * d_e;
        c_array(5, 0, 5) = -p1 + p2;
        c_array(5, 1, 5) = -I * (p1 + p2);
        p1 = abs_d * (0.75 * (6.0 - nu) * h2 + 0.25 / 3.0 * c_2_mnu * abs_d_2);
        p2 = 0.5 * d * d_e * (nu * abs_d_2 + 3.0 * c_3_nu * h2);
        c_array(5, 0, 6) = p1 + p2;
        c_array(5, 1, 6) = I * (p1 - p2);
        p1 = 1.0 / 6.0 * abs_d * h2 * ((15.0 - nu) * h2 + c_5_mnu * abs_d_2);
        p2 = 0.125 * d * d_e *
                (3.0 * c_12_nu * h4 + 2.0 * c_6_nu * d2h2 - nu * abs_d_4);
        c_array(5, 0, 7) = p1 + p2;
        c_array(5, 1, 7) = I * (p1 - p2);
        p0 = h2 * c_d_h;
        p1 = 1.0 / 3.0 * h2 * abs_d;
        p2 = 0.25 * d * d_e * c_d_m3h;
        c_array(5, 0, 8) = p0 * (p1 - p2);
        c_array(5, 1, 8) = I * p0 * (p1 + p2);

        c_array(5, 2, 1) = 0.625 * I * h * d;
        c_array(5, 2, 2) = 0.0625 * h * d * (-19.0 * I * h2 + abs_d_2 * tan_x);
        c_array(5, 2, 3) = 0.1875 * h3 * d * (3.0 * I * h2 - abs_d_2 * tan_x);
        c_array(5, 2, 5) = 1.5 * h * d_e;
        c_array(5, 2, 6) = (abs_d_2 + 12.0 * h2) * h * d_e;
        c_array(5, 2, 7) = 0.25 * (-abs_d_4 +
                14.0 * d2h2 + 39.0 * h4) * h * d_e;
        c_array(5, 2, 8) = -0.5 * d_e * h3 * c_d_h * c_d_m3h;

        for (int j = 0; j < c_array.size(2); ++j) {
            c_array(4, 2, j) = std::conj(c_array(5, 2, j));
        }

        return c_array;
    }

    il::StaticArray3D<std::complex<double>, 6, 3, 9> s_33_h
            (double nu, std::complex<double> eix,
             double h, std::complex<double> d) {
        const std::complex<double> I(0.0, 1.0);

        double cos_x = std::real(eix);
        double tan_x = std::imag(eix) / cos_x;
        std::complex<double> tcos_x = cos_x * eix;
        std::complex<double> c_3_4tcos = 3.0 + 4.0 * tcos_x;
        std::complex<double> c_8_3i_tan_x = 8.0 + 3.0 * I * tan_x;

        double h2 = h * h;
        double h3 = h * h2;
        double h4 = h2 * h2;
//  double sgh = (( h < 0 )? -1.0 : double(( h > 0 ))); // sign(h)

        double abs_d = std::abs(d);
        double abs_d_2 = abs_d * abs_d;
        double abs_d_4 = abs_d_2 * abs_d_2;
        // std::complex<double> d_c = std::conj(d);
        std::complex<double> d_e = std::polar(1.0, std::arg(d)); //  = d/abs_d
        std::complex<double> d_e_2 = d_e * d_e; //  = d^2/abs_d^2

        std::complex<double> c_d_h = abs_d_2 + h2;
//  std::complex<double> c_d_3h = abs_d_2+3.0*h2;
        std::complex<double> c_d_m3h = abs_d_2 - 3.0 * h2;

        std::complex<double> p0, p1, p2;

        il::StaticArray3D<std::complex<double>, 6, 3, 9> c_array{0.0};

        c_array(0, 0, 6) = -2.0 * h * std::real(d_e);
        c_array(0, 1, 6) = -2.0 * h * std::imag(d_e);
        c_array(0, 2, 6) = -2.0 * abs_d;
        c_array(0, 0, 7) = -4.0 * h3 * std::real(d_e);
        c_array(0, 1, 7) = -4.0 * h3 * std::imag(d_e);
        c_array(0, 2, 7) = 4.0 * h2 * abs_d;

        c_array(1, 1, 1) = -d_e_2 * h * tcos_x;
        c_array(1, 0, 1) = I * c_array(1, 1, 1);
        p1 = 1.0 / 12.0 * tan_x * abs_d_2;
        p2 = d_e_2 * h2 * tcos_x;
        c_array(1, 0, 2) = (p1 + I * p2) * h;
        c_array(1, 1, 2) = (I * p1 + p2) * h;
        c_array(1, 0, 3) = -0.25 * abs_d_2 * h3 * tan_x;
        c_array(1, 1, 3) = I * c_array(1, 0, 3);
        c_array(1, 2, 5) = d_e;
        p2 = 10.0 / 3.0 * abs_d;
        c_array(1, 0, 6) = (-d * d_e - p2) * h;
        c_array(1, 1, 6) = I * (d * d_e - p2) * h;
        c_array(1, 2, 6) = -4.0 * h2 * d_e;
        p0 = abs_d * h;
        p1 = (abs_d_2 + 11.0 * h2) / 3.0;
        p2 = 2.0 * d_e_2 * h2;
        c_array(1, 0, 7) = -(p1 + p2) * p0;
        c_array(1, 1, 7) = -I * (p1 - p2) * p0;
        c_array(1, 2, 7) = -4.0 * h4 * d_e;
        c_array(1, 0, 8) = -2.0 / 3.0 * h3 * abs_d * c_d_h;
        c_array(1, 1, 8) = I * c_array(1, 0, 8);

        for (int k = 0; k < c_array.size(2); ++k) {
            for (int j = 0; j < c_array.size(1); ++j) {
                c_array(2, j, k) = std::conj(c_array(1, j, k));
            }
        }

        c_array(3, 2, 0) = 2.0 * I * d_e_2 * tcos_x;
        p0 = d * h;
        p2 = d_e_2 * (0.625 + tcos_x);
        c_array(3, 0, 1) = I * p0 * (0.625 - p2);
        c_array(3, 1, 1) = -p0 * (0.625 + p2);
        c_array(3, 2, 1) = -4.0 * I * d_e_2 * h2 * tcos_x;
        p0 = 0.0625 * d * h;
        p1 = -19.0 * I * h2 + abs_d_2 * tan_x;
        p2 = d_e_2 * (abs_d_2 * c_8_3i_tan_x + h2 * (19.0 + 28.0 * tcos_x));
        c_array(3, 0, 2) = p0 * (p1 + I * p2);
        c_array(3, 1, 2) = p0 * (I * p1 + p2);
        c_array(3, 2, 2) = 2.0 * I * d_e_2 * h4 * tcos_x;
        p0 = 0.1875 * d * h3;
        p1 = 3.0 * I * h2 - abs_d_2 * tan_x;
        p2 = d_e_2 * (abs_d_2 * c_8_3i_tan_x + h2 * c_3_4tcos);
        c_array(3, 0, 3) = p0 * (p1 - I * p2);
        c_array(3, 1, 3) = p0 * (I * p1 - p2);
        p0 = h * d_e;
        c_array(3, 0, 5) = p0 * (-0.5 * d_e_2 + 1.5);
        c_array(3, 1, 5) = I * p0 * (0.5 * d_e_2 + 1.5);
        c_array(3, 2, 5) = d * d_e;
        p1 = abs_d_2 + 12.0 * h2;
        p2 = d_e_2 * (3.0 * abs_d_2 + 4.0 * h2);
        c_array(3, 0, 6) = p0 * (p1 - p2);
        c_array(3, 1, 6) = I * p0 * (p1 + p2);
        c_array(3, 2, 6) = -4.0 * h2 * d * d_e;
        p1 = -0.25 * abs_d_4 + 3.5 * abs_d_2 * h2 + 9.75 * h4;
        p2 = d_e_2 * (1.25 * abs_d_4 + 10.5 * abs_d_2 * h2 + 3.25 * h4);
        c_array(3, 0, 7) = p0 * (p1 - p2);
        c_array(3, 1, 7) = I * p0 * (p1 + p2);
        c_array(3, 2, 7) = -4.0 * h4 * d * d_e;
        p0 = h3 * d_e * c_d_h;
        p1 = 0.5 * c_d_m3h;
        p2 = d_e_2 * (2.5 * abs_d_2 + 0.5 * h2);
        c_array(3, 0, 8) = -p0 * (p1 + p2);
        c_array(3, 1, 8) = I * p0 * (-p1 + p2);

        for (int k = 0; k < c_array.size(2); ++k) {
            for (int j = 0; j < c_array.size(1); ++j) {
                c_array(4, j, k) = std::conj(c_array(3, j, k));
            }
        }

        c_array(5, 0, 1) = -1.25 * h * std::imag(d);
        c_array(5, 1, 1) = 1.25 * h * std::real(d);
        c_array(5, 2, 1) = 1.0 / 6.0 * abs_d_2 * tan_x;
        p1 = 0.125 * abs_d_2 * tan_x;
        p2 = 2.375 * h2;
        c_array(5, 0, 2) = h * (p1 * std::real(d) + p2 * std::imag(d));
        c_array(5, 1, 2) = h * (p1 * std::imag(d) - p2 * std::real(d));
        c_array(5, 2, 2) = 1.0 / 3.0 * abs_d_2 * h2 * tan_x;
        p1 = 3.0 * p1;
        p2 = 1.125 * h2;
        c_array(5, 0, 3) = -h3 * (p1 * std::real(d) + p2 * std::imag(d));
        c_array(5, 1, 3) = -h3 * (p1 * std::imag(d) - p2 * std::real(d));
        c_array(5, 2, 3) = -0.5 * abs_d_2 * h4 * tan_x;
        c_array(5, 0, 5) = 3.0 * h * std::real(d_e);
        c_array(5, 1, 5) = 3.0 * h * std::imag(d_e);
        c_array(5, 2, 5) = -abs_d;
        p0 = 2.0 * h * (abs_d_2 + 12.0 * h2);
        c_array(5, 0, 6) = p0 * std::real(d_e);
        c_array(5, 1, 6) = p0 * std::imag(d_e);
        c_array(5, 2, 6) = (1.0 / 3.0 * abs_d_2 - 9.0 * h2) * abs_d;
        p0 = h * (-0.5 * abs_d_4 + 7.0 * abs_d_2 * h2 + 19.5 * h4);
        c_array(5, 0, 7) = p0 * std::real(d_e);
        c_array(5, 1, 7) = p0 * std::imag(d_e);
        c_array(5, 2, 7) = -h2 * (4.0 / 3.0 * abs_d_2 + 8.0 * h2) * abs_d;
        p0 = h3 * c_d_h * c_d_m3h;
        c_array(5, 0, 8) = -p0 * std::real(d_e);
        c_array(5, 1, 8) = -p0 * std::imag(d_e);
        c_array(5, 2, 8) = -4.0 / 3.0 * h4 * c_d_h * abs_d;

        return c_array;
    }

// Additional terms for a special case:
// reduced summation; collocation point projected onto
// an edge line or a vertex of the element

    il::StaticArray3D<std::complex<double>, 6, 3, 5> s_11_22_red_h
            (double nu, std::complex<double> eix,
             double h, std::complex<double> d) {
        const std::complex<double> I(0.0, 1.0);

        double c_1_nu = 1.0 + nu;
        double c_1_2nu = 1.0 + 2.0 * nu;
        double c_2_nu = 2.0 + nu;
        double c_3_2nu = 3.0 + 2.0 * nu;
        double c_7_2nu = 7.0 + 2.0 * nu;
        double c_11_4nu = 11.0 + 4.0 * nu;
        double c_13_2nu = 13.0 + 2.0 * nu;

        double cos_x = std::real(eix);
        double sin_x = std::imag(eix);
        std::complex<double> e2x = eix * eix;
        std::complex<double> c_eix_3_1 = eix * (3.0 + e2x);
        std::complex<double> c_eix_3_m1 = eix * (3.0 - e2x);

        double h2 = h * h;
        double h3 = h2 * h;
        double h4 = h2 * h2;
        double h5 = h4 * h;
        double h7 = h5 * h2;
        double sgh = ((h < 0) ? -1.0 : static_cast<double>((h > 0))); // sign(h)

        std::complex<double> p1, p2, p3, p4;

        il::StaticArray3D<std::complex<double>, 6, 3, 5> c_array{0.0};

        c_array(0, 0, 2) = -c_7_2nu * h * sin_x;
        c_array(0, 1, 2) = c_7_2nu * h * cos_x;
        c_array(0, 0, 3) = -9.0 * h3 * sin_x;
        c_array(0, 1, 3) = 9.0 * h3 * cos_x;
        c_array(0, 0, 4) = -1.5 * h5 * sin_x;
        c_array(0, 1, 4) = 1.5 * h5 * cos_x;

        c_array(1, 0, 0) = -0.5 * I * c_1_nu * e2x * sgh;
        c_array(1, 1, 0) = -0.5 * c_1_nu * e2x * sgh;
        c_array(1, 2, 1) = 0.5 * I * c_1_2nu * eix;
        c_array(1, 2, 2) = 2.0 * I * c_2_nu * h2 * eix;
        c_array(1, 2, 3) = 2.0 * I * h4 * eix;

        for (int k = 0; k < c_array.size(2); ++k) {
            for (int j = 0; j < c_array.size(1); ++j) {
                c_array(2, j, k) = std::conj(c_array(1, j, k));
            }
        }

        p1 = 0.25 * c_3_2nu * h;
        p2 = 0.25 * c_11_4nu * h3;
        p3 = 0.125 * c_13_2nu * h5;
        p4 = 0.625 * h7;

        c_array(3, 2, 0) = -2.0 * I * c_1_nu * e2x * std::fabs(h);
        c_array(3, 0, 1) = -I * p1 * c_eix_3_1;
        c_array(3, 1, 1) = p1 * c_eix_3_m1;
        c_array(3, 0, 2) = -I * p2 * c_eix_3_1;
        c_array(3, 1, 2) = p2 * c_eix_3_m1;
        c_array(3, 0, 3) = -I * p3 * c_eix_3_1;
        c_array(3, 1, 3) = p3 * c_eix_3_m1;
        c_array(3, 0, 4) = -I * p4 / 3.0 * c_eix_3_1;
        c_array(3, 1, 4) = p4 / 3.0 * c_eix_3_m1;

        for (int k = 0; k < c_array.size(2); ++k) {
            for (int j = 0; j < c_array.size(1); ++j) {
                c_array(4, j, k) = std::conj(c_array(3, j, k));
            }
        }

        c_array(5, 0, 1) = 6.0 * p1 * sin_x;
        c_array(5, 1, 1) = -6.0 * p1 * cos_x;
        c_array(5, 0, 2) = 6.0 * p2 * sin_x;
        c_array(5, 1, 2) = -6.0 * p2 * cos_x;
        c_array(5, 0, 3) = 6.0 * p3 * sin_x;
        c_array(5, 1, 3) = -6.0 * p3 * cos_x;
        c_array(5, 0, 4) = 2.0 * p4 * sin_x;
        c_array(5, 1, 4) = -2.0 * p4 * cos_x;

        return c_array;
    }

    il::StaticArray3D<std::complex<double>, 6, 3, 5> s_12_red_h
            (double nu, std::complex<double> eix,
             double h, std::complex<double> d) {
        const std::complex<double> I(0.0, 1.0);

        double c_1_mnu = 1.0 - nu;
        double c_1_m2nu = 1.0 - 2.0 * nu;
        double c_2_mnu = 2.0 - nu;
        double c_5_m4nu = 5.0 - 4.0 * nu;
        double c_13_m2nu = 13.0 - 2.0 * nu;
        double c_15_m4nu = 15.0 - 4.0 * nu;
        double c_15_m8nu = 15.0 - 8.0 * nu;

        double cos_x = std::real(eix);
        double sin_x = std::imag(eix);
        std::complex<double> e2x = eix * eix;
        //std::complex<double> emx = std::conj(eix);
        //std::complex<double> Ce2x3_1 = 3.0+e2x;
        //std::complex<double> Ce2x3m1 = 3.0-e2x;
        std::complex<double> c_eix_3_1 = eix * (3.0 + e2x);
        std::complex<double> c_eix_3_m1 = eix * (3.0 - e2x);

        double h2 = h * h;
        double h3 = h2 * h;
        double h4 = h2 * h2;
        double h5 = h4 * h;
        double h7 = h5 * h2;
        double sgh = ((h < 0) ? -1.0 : static_cast<double>((h > 0))); // sign(h)
        double abh = std::fabs(h);

        std::complex<double> p1, p2, p3, p4;

        il::StaticArray3D<std::complex<double>, 6, 3, 5> c_array{0.0};

        c_array(0, 2, 0) = -I * nu * e2x / abh;
        c_array(0, 0, 2) = 0.5 * I * h * (3.0 * c_eix_3_1 - 4.0 * nu * eix);
        c_array(0, 1, 2) = -0.5 * h * (3.0 * c_eix_3_m1 - 4.0 * nu * eix);
        c_array(0, 0, 3) = 1.5 * I * h3 * c_eix_3_1;
        c_array(0, 1, 3) = -1.5 * h3 * c_eix_3_m1;
        c_array(0, 0, 4) = 0.25 * I * h5 * c_eix_3_1;
        c_array(0, 1, 4) = -0.25 * h5 * c_eix_3_m1;

        p1 = 0.5 * I * c_1_m2nu * eix;
        p2 = 2.0 * I * c_2_mnu * h2 * eix;
        p3 = 0.25 * I * c_13_m2nu * h4 * eix;
        p4 = 0.5 * I * h2 * h4 * eix;

        c_array(1, 0, 0) = -I * sgh * e2x * (c_1_mnu + 0.5 * e2x);
        c_array(1, 1, 0) = sgh * e2x * (c_1_mnu - 0.5 * e2x);
        c_array(1, 2, 1) = p1 * e2x;
        c_array(1, 2, 2) = p2 * e2x;
        c_array(1, 2, 3) = p3 * e2x;
        c_array(1, 2, 4) = p4 * e2x;

        c_array(2, 1, 0) = -sgh * e2x;
        c_array(2, 0, 0) = I * c_array(2, 1, 0);
        c_array(2, 2, 1) = 3.0 * p1;
        c_array(2, 2, 2) = 3.0 * p2;
        c_array(2, 2, 3) = 3.0 * p3;
        c_array(2, 2, 4) = 3.0 * p4;

        p1 = 0.25 * h * eix;
        p2 = 0.25 * h3 * eix;
        p3 = 0.125 * h5 * eix;
        p4 = 0.125 * h7 * eix;

        c_array(3, 2, 0) = -2.0 * I * c_1_mnu * abh * e2x * e2x;
        c_array(3, 0, 1) = -I * e2x * p1 * (c_5_m4nu + 3.0 * e2x);
        c_array(3, 0, 2) = -I * e2x * p2 * (c_15_m8nu + 9.0 * e2x);
        c_array(3, 0, 3) = -I * e2x * p3 * (c_15_m4nu + 9.0 * e2x);
        c_array(3, 0, 4) = -I * e2x * p4 * (5.0 / 3.0 + e2x);

        c_array(3, 1, 1) = e2x * p1 * (c_5_m4nu - 3.0 * e2x);
        c_array(3, 1, 2) = e2x * p2 * (c_15_m8nu - 9.0 * e2x);
        c_array(3, 1, 3) = e2x * p3 * (c_15_m4nu - 9.0 * e2x);
        c_array(3, 1, 4) = e2x * p4 * (5.0 / 3.0 - e2x);

        c_array(5, 2, 0) = -4.0 * I * c_1_mnu * abh * e2x;
        c_array(5, 0, 1) = -I * p1 * (3.0 * c_5_m4nu + 5.0 * e2x);
        c_array(5, 0, 2) = -3.0 * I * p2 * (c_15_m8nu + 5.0 * e2x);
        c_array(5, 0, 3) = -3.0 * I * p3 * (c_15_m4nu + 5.0 * e2x);
        c_array(5, 0, 4) = -5.0 * I * p4 * (1.0 + e2x / 3.0);

        c_array(5, 1, 1) = p1 * (3.0 * c_5_m4nu - 5.0 * e2x);
        c_array(5, 1, 2) = 3.0 * p2 * (c_15_m8nu - 5.0 * e2x);
        c_array(5, 1, 3) = 3.0 * p3 * (c_15_m4nu - 5.0 * e2x);
        c_array(5, 1, 4) = 5.0 * p4 * (1.0 - e2x / 3.0);

        c_array(4, 0, 1) = 3.0 * I * std::conj(p1) * (c_5_m4nu - 5.0 * e2x);
        //c_array(4, 0, 1) = -0.75*I*h*(5.0*eix-c_5_m4nu*emx);
        c_array(4, 1, 1) = -3.0 * std::conj(p1) * (c_5_m4nu + 5.0 * e2x);
        //c_array(4, 1, 1) = -0.75*h*(5.0*eix+c_5_m4nu*emx);
        c_array(4, 0, 2) = 3.0 * I * std::conj(p2) * (c_15_m8nu - 15.0 * e2x);
        //c_array(4, 0, 2) = -0.75*I*h3*(15.0*eix-c_15_m8nu*emx);
        c_array(4, 1, 2) = -3.0 * std::conj(p2) * (c_15_m8nu + 15.0 * e2x);
        //c_array(4, 1, 2) = -0.75*h3*(15.0*eix+c_15_m8nu*emx);
        c_array(4, 0, 3) = 3.0 * I * std::conj(p3) * (c_15_m4nu - 15.0 * e2x);
        //c_array(4, 0, 3) = -0.375*I*h5*(15.0*eix-c_15_m4nu*emx);
        c_array(4, 1, 3) = -3.0 * std::conj(p3) * (c_15_m4nu + 15.0 * e2x);
        //c_array(4, 1, 3) = -0.375*h5*(15.0*eix+c_15_m4nu*emx);
        c_array(4, 0, 4) = 1.25 * h7 * sin_x;
        c_array(4, 1, 4) = -1.25 * h7 * cos_x;

        return c_array;
    }

    il::StaticArray3D<std::complex<double>, 6, 3, 5> s_13_23_red_h
            (double nu, std::complex<double> eix,
             double h, std::complex<double> d) {
        const std::complex<double> I(0.0, 1.0);

        double c_1_nu = 1.0 + nu;
        double c_1_mnu = 1.0 - nu;
        double c_2_mnu = 2.0 - nu;
        double c_3_nu = 3.0 + nu;
        double c_3_mnu = 3.0 - nu;
        double c_5_mnu = 5.0 - nu;
        double c_12_nu = 12.0 + nu;

        // double cos_x = std::real(eix);
        // double sin_x = std::imag(eix);
        std::complex<double> e2x = eix * eix;
        std::complex<double> e3x = e2x * eix;
        std::complex<double> emx = std::conj(eix);
        std::complex<double> em2 = std::conj(e2x);
        std::complex<double> c_e3x_3emx = 3.0 * emx + e3x;
        //std::complex<double> c_eix_3_1 = eix * (3.0 + e2x);
        //std::complex<double> c_eix_3_m1 = eix * (3.0 - e2x);

        double h2 = h * h;
        double h3 = h2 * h;
        double h4 = h2 * h2;
        double h5 = h4 * h;
        double h6 = h4 * h2;
        double h7 = h5 * h2;
        // double sgh = (( h < 0 )? -1.0 : double(( h > 0 ))); // sign(h)
        double abh = std::fabs(h);

        std::complex<double> p1, p2, p3, p4;

        il::StaticArray3D<std::complex<double>, 6, 3, 5> c_array{0.0};

        c_array(0, 1, 0) = -0.25 * c_1_mnu * e2x / abh;
        c_array(0, 0, 0) = I * c_array(0, 1, 0);
        c_array(0, 2, 2) = -I * h * eix;
        c_array(0, 2, 3) = -2.0 * I * h3 * eix;

        p1 = 0.25 * nu * c_e3x_3emx;
        p2 = 0.5 * h2 * c_3_nu * c_e3x_3emx;
        p3 = 0.125 * h4 * c_12_nu * c_e3x_3emx;
        // p4 = 0.25*h6*c_e3x_3emx;

        c_array(1, 0, 1) = 0.25 * I * eix * (c_2_mnu + nu * e2x);
        c_array(1, 1, 1) = -0.25 * eix * (c_2_mnu - nu * e2x);
        c_array(1, 0, 2) = 0.5 * I * h2 * eix * (c_5_mnu + c_3_nu * e2x);
        c_array(1, 1, 2) = -0.5 * h2 * eix * (c_5_mnu - c_3_nu * e2x);
        c_array(1, 0, 3) = 0.125 * I * h4 * eix * (8.0 + c_12_nu * e2x);
        c_array(1, 1, 3) = -0.125 * h4 * eix * (8.0 - c_12_nu * e2x);
        c_array(1, 1, 4) = 0.25 * h6 * e3x;
        c_array(1, 0, 4) = I * c_array(1, 1, 4);

        c_array(2, 0, 1) = std::conj(c_array(1, 0, 1) - I * p1);
        c_array(2, 1, 1) = std::conj(p1 - c_array(1, 1, 1));
        c_array(2, 0, 2) = std::conj(c_array(1, 0, 2) - I * p2);
        c_array(2, 1, 2) = std::conj(p2 - c_array(1, 1, 2));
        c_array(2, 0, 3) = std::conj(c_array(1, 0, 3) - I * p3);
        c_array(2, 1, 3) = std::conj(p3 - c_array(1, 1, 3));
        c_array(2, 1, 4) = 0.75 * eix * h6;
        c_array(2, 0, 4) = I * c_array(2, 1, 4);

        p1 = 0.5 * I * h * eix;
        p2 = 4.0 * I * h3 * eix;
        p3 = 3.25 * I * h5 * eix;
        p4 = 0.5 * I * h7 * eix;

        c_array(5, 1, 0) = -c_1_nu * abh * e2x;
        c_array(5, 0, 0) = I * c_array(5, 1, 0);
        c_array(5, 2, 1) = 3.0 * p1;
        c_array(5, 2, 2) = 3.0 * p2;
        c_array(5, 2, 3) = 3.0 * p3;
        c_array(5, 2, 4) = 3.0 * p4;

        c_array(3, 1, 0) = 0.5 * (c_3_mnu - c_1_nu * e2x) * abh * e2x;
        c_array(3, 0, 0) = -0.5 * I * (c_3_mnu + c_1_nu * e2x) * abh * e2x;
        c_array(3, 2, 1) = e2x * p1;
        c_array(3, 2, 2) = e2x * p2;
        c_array(3, 2, 3) = e2x * p3;
        c_array(3, 2, 4) = e2x * p4;

        c_array(4, 1, 0) = -0.5 * c_3_mnu * abh * em2;
        c_array(4, 0, 0) = -I * c_array(4, 1, 0);
        c_array(4, 2, 1) = std::conj(c_array(5, 2, 1));
        c_array(4, 2, 2) = std::conj(c_array(5, 2, 2));
        c_array(4, 2, 3) = std::conj(c_array(5, 2, 3));
        c_array(4, 2, 4) = std::conj(c_array(5, 2, 4));

        return c_array;
    }

    il::StaticArray3D<std::complex<double>, 6, 3, 5> s_33_red_h
            (double nu, std::complex<double> eix,
             double h, std::complex<double> d) {
        const std::complex<double> I(0.0, 1.0);

        double cos_x = std::real(eix);
        double sin_x = std::imag(eix);
        std::complex<double> e2x = eix * eix;
        std::complex<double> c_eix_3_1 = eix * (3.0 + e2x);
        std::complex<double> c_eix_3_m1 = eix * (3.0 - e2x);

        double h2 = h * h;
        double h3 = h2 * h;
        double h4 = h2 * h2;
        double h5 = h4 * h;
        double h7 = h5 * h2;
        // double sgh = (( h < 0 )? -1.0 : double(( h > 0 ))); // sign(h)

        il::StaticArray3D<std::complex<double>, 6, 3, 5> c_array{0.0};

        c_array(0, 0, 2) = 2.0 * h * sin_x;
        c_array(0, 1, 2) = -2.0 * h * cos_x;
        c_array(0, 0, 3) = 4.0 * h3 * sin_x;
        c_array(0, 1, 3) = -4.0 * h3 * cos_x;

        c_array(1, 2, 1) = I * eix;
        c_array(1, 2, 2) = -4.0 * I * h2 * eix;
        c_array(1, 2, 3) = -4.0 * I * h4 * eix;

        for (int j = 1; j < c_array.size(2) - 1; ++j) {
            c_array(2, 2, j) = std::conj(c_array(1, 2, j));
        }

        c_array(3, 0, 1) = 0.5 * I * h * c_eix_3_1;
        c_array(3, 1, 1) = -0.5 * h * c_eix_3_m1;
        c_array(3, 0, 2) = 4.0 * I * h3 * c_eix_3_1;
        c_array(3, 1, 2) = -4.0 * h3 * c_eix_3_m1;
        c_array(3, 0, 3) = 3.25 * I * h5 * c_eix_3_1;
        c_array(3, 1, 3) = -3.25 * h5 * c_eix_3_m1;
        c_array(3, 0, 4) = 0.5 * I * h7 * c_eix_3_1;
        c_array(3, 1, 4) = -0.5 * h7 * c_eix_3_m1;

        for (int k = 1; k < c_array.size(2); ++k) {
            for (int j = 0; j < c_array.size(1) - 1; ++j) {
                c_array(4, j, k) = std::conj(c_array(3, j, k));
            }
        }

        c_array(5, 0, 1) = -3.0 * h * sin_x;
        c_array(5, 1, 1) = 3.0 * h * cos_x;
        c_array(5, 0, 2) = -24.0 * h3 * sin_x;
        c_array(5, 1, 2) = 24.0 * h3 * cos_x;
        c_array(5, 0, 3) = -19.5 * h5 * sin_x;
        c_array(5, 1, 3) = 19.5 * h5 * cos_x;
        c_array(5, 0, 4) = -3.0 * h7 * sin_x;
        c_array(5, 1, 4) = 3.0 * h7 * cos_x;

        return c_array;
    }

// Limit case (h==0, plane) - all stress components

    il::StaticArray3D<std::complex<double>, 6, 4, 3> s_ij_lim_h
            (double nu, std::complex<double> eix,
             double sgnh, std::complex<double> d) {
        const std::complex<double> I(0.0, 1.0);

        double c_1_2nu = 1.0 + 2.0 * nu;
        double c_1_m2nu = 1.0 - 2.0 * nu;
        double c_2_mnu = 2.0 - nu;

        double cos_x = std::real(eix);
        double sin_x = std::imag(eix);
        // double tan_x = sin_x/cos_x;
        std::complex<double> e2x = eix * eix;
        double h0_lim = std::atanh(sin_x); // std::complex<double> ?
        // double h0_lim = 0.5*(std::log(1.0+sin_x)-std::log(1.0-sin_x));

        double abs_d = std::abs(d); 
        // double abs_d_2 = abs_d*abs_d; double abs_d_4 = abs_d_2*abs_d_2;
        std::complex<double> d_e = std::polar(1.0, std::arg(d)), // = d/abs_d
                d_e_2 = d_e * d_e, 
                d_e_3 = d_e * d_e_2, 
                d_e_4 = d_e_2 * d_e_2,
                p1, p2;

        il::StaticArray3D<std::complex<double>, 6, 4, 3> c_array{0.0};

        il::StaticArray<std::complex<double>, 6> v1{0.0}, v2{0.0};

        v1[0] = sin_x / abs_d;
        v1[1] = h0_lim * d_e;
        v1[2] = std::conj(v1[1]);
        v1[3] = d * d_e * (h0_lim + 2.0 * I * eix);
        v1[4] = std::conj(v1[3]);
        v1[5] = -h0_lim * abs_d;

        v2[0] = 0.5 * d_e_2 / abs_d * (sin_x - 2.0 * I * eix * cos_x * cos_x);
        v2[1] = -0.125 * d_e_3 * (4.0 * h0_lim + I * eix * (e2x + 8.0));
        v2[2] = 0.125 * d_e * (12.0 * h0_lim + 5.0 * I * eix);
        v2[3] = -0.5 * d_e_4 * abs_d * (3.0 * h0_lim - 2.0 * I * eix * (e2x - 3.0));
        v2[4] = -1.5 * abs_d * h0_lim;
        v2[5] = 1.5 * d_e_2 * abs_d * (h0_lim + 2.0 * I * eix);

        for (int j = 0; j < c_array.size(0); ++j) {
            c_array(j, 0, 2) = 0.5 * c_1_2nu * v1[j];
            c_array(j, 1, 2) = c_1_m2nu * v2[j];
            p1 = 0.25 * c_2_mnu * v1[j];
            p2 = 0.5 * nu * v2[j];
            c_array(j, 2, 0) = p1 + p2;
            c_array(j, 2, 1) = I * (p1 - p2);
            c_array(j, 3, 2) = v1[j];
        }

        return c_array;
    }
}
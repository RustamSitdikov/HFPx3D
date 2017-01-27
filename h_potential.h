//
// This file is part of 3d_bem.
//
// Created by D. Nikolski on 1/24/2017.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//

//

#ifndef INC_3D_BEM_H_POTENTIAL_H
#define INC_3D_BEM_H_POTENTIAL_H

#include <complex>
#include <il/StaticArray3D.h>

namespace hfp3d {

    il::StaticArray3D<std::complex<double>, 6, 3, 9> s_11_22_h
            (double nu, std::complex<double> eix,
             double h, std::complex<double> d);

    il::StaticArray3D<std::complex<double>, 6, 3, 9> s_12_h
            (double nu, std::complex<double> eix,
             double h, std::complex<double> d);

    il::StaticArray3D<std::complex<double>, 6, 3, 9> s_13_23_h
            (double nu, std::complex<double> eix,
             double h, std::complex<double> d);

    il::StaticArray3D<std::complex<double>, 6, 3, 9> s_33_h
            (double nu, std::complex<double> eix,
             double h, std::complex<double> d);

    il::StaticArray3D<std::complex<double>, 6, 3, 5> s_11_22_red_h
            (double nu, std::complex<double> eix,
             double h, std::complex<double> d);

    il::StaticArray3D<std::complex<double>, 6, 3, 5> s_12_red_h
            (double nu, std::complex<double> eix,
             double h, std::complex<double> d);

    il::StaticArray3D<std::complex<double>, 6, 3, 5> s_13_23_red_h
            (double nu, std::complex<double> eix,
             double h, std::complex<double> d);

    il::StaticArray3D<std::complex<double>, 6, 3, 5> s_33_red_h
            (double nu, std::complex<double> eix,
             double h, std::complex<double> d);

    il::StaticArray3D<std::complex<double>, 6, 4, 3> s_ij_lim_h
            (double nu, std::complex<double> eix,
             double signh, std::complex<double> d);

}
#endif //INC_3D_BEM_H_POTENTIAL_H

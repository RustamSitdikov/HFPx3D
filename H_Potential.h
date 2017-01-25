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

    il::StaticArray3D<std::complex<double>, 6, 3, 9> S_11_22_H
            (double nu, std::complex<double> eix,
             double h, std::complex<double> d);

    il::StaticArray3D<std::complex<double>, 6, 3, 9> S_12_H
            (double nu, std::complex<double> eix,
             double h, std::complex<double> d);

    il::StaticArray3D<std::complex<double>, 6, 3, 9> S_13_23_H
            (double nu, std::complex<double> eix,
             double h, std::complex<double> d);

    il::StaticArray3D<std::complex<double>, 6, 3, 9> S_33_H
            (double nu, std::complex<double> eix,
             double h, std::complex<double> d);

    il::StaticArray3D<std::complex<double>, 6, 3, 5> S_11_22_Red_H
            (double nu, std::complex<double> eix,
             double h, std::complex<double> d);

    il::StaticArray3D<std::complex<double>, 6, 3, 5> S_12_Red_H
            (double nu, std::complex<double> eix,
             double h, std::complex<double> d);

    il::StaticArray3D<std::complex<double>, 6, 3, 5> S_13_23_Red_H
            (double nu, std::complex<double> eix,
             double h, std::complex<double> d);

    il::StaticArray3D<std::complex<double>, 6, 3, 5> S_33_Red_H
            (double nu, std::complex<double> eix,
             double h, std::complex<double> d);

    il::StaticArray3D<std::complex<double>, 6, 4, 3> S_ij_Lim_H
            (double nu, std::complex<double> eix,
             double signh, std::complex<double> d);

}
#endif //INC_3D_BEM_H_POTENTIAL_H

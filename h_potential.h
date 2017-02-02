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
#include <il/StaticArray4D.h>

namespace hfp3d {

    il::StaticArray4D<std::complex<double>, 6, 4, 3, 9> s_ij_gen_h
            (double nu, std::complex<double> eix,
             double h, std::complex<double> d);

    il::StaticArray4D<std::complex<double>, 6, 4, 3, 5> s_ij_red_h
            (double nu, std::complex<double> eix,
             double h);

    il::StaticArray3D<std::complex<double>, 6, 4, 3> s_ij_lim_h
            (double nu, std::complex<double> eix,
             std::complex<double> d);

}
#endif //INC_3D_BEM_H_POTENTIAL_H

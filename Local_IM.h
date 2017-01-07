//
// This file is part of 3d_bem.
//
// Created by D. Nikolski on 1/6/2017.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#ifndef INC_3D_BEM_LOCAL_IM_H
#define INC_3D_BEM_LOCAL_IM_H

#endif //INC_3D_BEM_LOCAL_IM_H

il::StaticArray<std::complex<double>, 9> ICFns(double, std::complex<double>, double, double, std::complex<double>);
il::StaticArray<std::complex<double>, 5> ICFns_red(double, std::complex<double>, double);

il::StaticArray2D<double, 6, 18> Local_IM_H(double, double, std::complex<double>, il::StaticArray<std::complex<double>,3>, il::StaticArray2D<std::complex<double>,6,6>);
//il::StaticArray2D<double, 6, 18> Local_IM_T(double, double, std::complex<double>, il::StaticArray<std::complex<double>,3>, il::StaticArray2D<std::complex<double>,6,6>);

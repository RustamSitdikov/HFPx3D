//
// This file is part of 3d_bem.
//
// Created by nikolski on 1/10/2017.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//

#ifndef INC_3D_BEM_ELE_BASE_H
#define INC_3D_BEM_ELE_BASE_H

#endif //INC_3D_BEM_ELE_BASE_H

struct el_x_cr {
    double h;
    std::complex<double> z;
};

void El_LB_RT(il::StaticArray2D<double,3,3>&, il::StaticArray2D<double,3,3>);
//il::StaticArray2D<double,3,3> El_LB_RT(il::StaticArray2D<double,3,3>);
void El_RT_Tr(il::StaticArray<std::complex<double>,3>&, il::StaticArray2D<double,3,3>&, il::StaticArray2D<double,3,3>&, il::StaticArray2D<double,3,3>);
void El_X_CR(el_x_cr&, il::StaticArray2D<double, 3,3>&, il::StaticArray2D<double,3,3>, il::StaticArray<double,3>);
il::StaticArray2D<std::complex<double>,2,2> El_CT(il::StaticArray2D<double,3,3>&, il::StaticArray2D<double,3,3>);
//il::StaticArray2D<std::complex<double>,2,2> El_CT(il::StaticArray2D<double,3,3>);
il::StaticArray2D<std::complex<double>,6,6> El_SFM_S(il::StaticArray2D<double,3,3>&, il::StaticArray2D<double,3,3>);
//il::StaticArray2D<std::complex<double>,6,6> El_SFM_S(il::StaticArray2D<double,3,3>);
il::StaticArray2D<std::complex<double>,6,6> El_SFM_N(il::StaticArray2D<double,3,3>&, il::StaticArray2D<double,3,3>, il::StaticArray<double,3>);
//il::StaticArray2D<std::complex<double>,6,6> El_SFM_N(il::StaticArray2D<double,3,3>, il::StaticArray<double,3>);
//il::StaticArray2D<std::complex<double>,6,6> El_SFM_C(il::StaticArray2D<double,3,3>&, il::StaticArray2D<double,3,3>, il::StaticArray<double,3>, double);
//il::StaticArray2D<std::complex<double>,6,6> El_SFM_C(il::StaticArray2D<double,3,3>, il::StaticArray<double,3>, double);
il::StaticArray<il::StaticArray<double, 3>, 6> El_CP_S(il::StaticArray2D<double,3,3>, double);
//il::StaticArray<il::StaticArray<double, 3>, 6> El_CP_N(il::StaticArray2D<double,3,3>, il::StaticArray<double,3>, double);

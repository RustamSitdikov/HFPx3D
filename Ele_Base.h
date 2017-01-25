//
// This file is part of 3d_bem.
//
// Created by D. Nikolski on 1/10/2017.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//

#ifndef INC_3D_BEM_ELE_BASE_H
#define INC_3D_BEM_ELE_BASE_H

namespace hfp3d {

    struct el_x_cr {
        double h;
        std::complex<double> z;
    };

    // Element's local coordinate system manipulations

    il::StaticArray2D<double, 3, 3> El_LB_RT
            (const il::StaticArray2D<double, 3, 3> &EV);

    il::StaticArray<std::complex<double>, 3> El_RT_Tr
            (const il::StaticArray2D<double, 3, 3> &EV,
             const il::StaticArray2D<double, 3, 3> &RT);

    el_x_cr El_X_CR
            (const il::StaticArray2D<double, 3, 3> &EV,
             const il::StaticArray<double, 3> &X0,
             const il::StaticArray2D<double, 3, 3> &RT);

    il::StaticArray2D<std::complex<double>, 2, 2> El_CT
            (const il::StaticArray2D<double, 3, 3> &EV,
             const il::StaticArray2D<double, 3, 3> &RT);

    // Element's basis (shape) functions

    il::StaticArray2D<std::complex<double>, 6, 6> El_SFM_S
            (const il::StaticArray2D<double, 3, 3> &EV,
             il::io_t, il::StaticArray2D<double, 3, 3> &RT);

    il::StaticArray2D<std::complex<double>, 6, 6> El_SFM_N
            (const il::StaticArray2D<double, 3, 3> &EV,
             const il::StaticArray<double, 3> &VW,
             il::io_t, il::StaticArray2D<double, 3, 3> &RT);

//il::StaticArray2D<std::complex<double>, 6, 6> El_SFM_C
// (const il::StaticArray2D<double, 3, 3> &EV,
// const il::StaticArray<double, 3> &VW,
// double beta,
// il::io_t, il::StaticArray2D<double, 3, 3> &RT);

    il::StaticArray2D<std::complex<double>, 6, 6> El_Shift_SFM
            (std::complex<double> z);

    // Collocation points

    il::StaticArray<il::StaticArray<double, 3>, 6> El_CP_S
            (const il::StaticArray2D<double, 3, 3> &EV, double beta);

    il::StaticArray<il::StaticArray<double, 3>, 6> El_CP_N
            (const il::StaticArray2D<double, 3, 3> &EV,
             const il::StaticArray<double, 3> &VW,
             double beta);

    // auxiliary functions (norm, cross product)

    double VNorm(const il::StaticArray<double, 3> &a);
    il::StaticArray<double, 3> normalize(const il::StaticArray<double, 3> &a);
    il::StaticArray<double, 3> cross(const il::StaticArray<double, 3> &a,
                                     const il::StaticArray<double, 3> &b);

}

#endif //INC_3D_BEM_ELE_BASE_H

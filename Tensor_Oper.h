//
// This file is part of 3d_bem.
//
// Created by D. Nikolski on 1/12/2017.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//

#ifndef INC_3D_BEM_TENSOR_OPER_H
#define INC_3D_BEM_TENSOR_OPER_H

#include <il/Array2D.h>
#include <il/StaticArray2D.h>

namespace hfp3d {

// Vector and triple tensor multiplication
// for stress stored as 6-component vector (or 6*N matrix)

    il::StaticArray2D<double, 3, 18> N_dot_SIM
            (const il::StaticArray<double, 3>& NV,
             const il::StaticArray2D<double, 6, 18>& SIM);

    il::StaticArray2D<double, 6, 18> SIM_P_R
            (const il::StaticArray2D<double, 3, 3>& RT_L,
             const il::StaticArray2D<double, 3, 3>& RT_R,
             const il::StaticArray2D<double, 6, 18>& SIM);

// Matrix-submatrix operations

    template<typename T_sub, typename T_A>
    T_sub get_submatrix
            (const T_A& A,
             il::int_t i0, il::int_t i1,
             il::int_t j0, il::int_t j1);

    template<typename T_sub, typename T_A>
    void set_submatrix
            (const T_sub& B,
             il::int_t i0, il::int_t i1,
             il::io_t, T_A& A);

    template <typename T_sub, typename T_A>
    void add_submatrix
            (const T_sub& B, double alpha,
             il::int_t i0, il::int_t i1,
             il::io_t, T_A& A);

}
#endif //INC_3D_BEM_TENSOR_OPER_H

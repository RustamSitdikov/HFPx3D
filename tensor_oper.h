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

    il::StaticArray2D<double, 3, 18> nv_dot_sim
            (const il::StaticArray<double, 3>& nv,
             const il::StaticArray2D<double, 6, 18>& sim);

    il::StaticArray2D<double, 6, 18> rotate_sim
            (const il::StaticArray2D<double, 3, 3>& rt,
             const il::StaticArray2D<double, 6, 18>& sim);

    il::StaticArray2D<double, 6, 18> rotate_sim_c
            (const il::StaticArray2D<double, 3, 3>& rt_l,
             const il::StaticArray2D<double, 3, 3>& rt_r,
             const il::StaticArray2D<double, 6, 18>& sim);

// Matrix-submatrix operations

    template<typename T_sub, typename T_A>
    T_sub get_submatrix
            (const T_A& a,
             il::int_t i0, il::int_t i1,
             il::int_t j0, il::int_t j1);

    template<typename T_sub, typename T_A>
    void set_submatrix
            (const T_sub& b,
             il::int_t i0, il::int_t i1,
             il::io_t, T_A& a);

    template <typename T_sub, typename T_A>
    void add_submatrix
            (const T_sub& b, double alpha,
             il::int_t i0, il::int_t i1,
             il::io_t, T_A& a);

}
#endif //INC_3D_BEM_TENSOR_OPER_H
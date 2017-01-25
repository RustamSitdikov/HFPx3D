//
// This file is part of 3d_bem.
//
// Created by D. Nikolski on 1/24/2017.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//

//#include <il/Array2D.h>
#include <il/StaticArray.h>
#include <il/StaticArray2D.h>
#include <il/linear_algebra/dense/blas/dot.h>
#include "Tensor_Oper.h"

namespace hfp3d {
// Vector and triple tensor multiplication
// for stress stored as 6-component vector (or 6*N matrix)

    il::StaticArray2D<double, 3, 18> N_dot_SIM
            (const il::StaticArray<double, 3> &NV,
             const il::StaticArray2D<double, 6, 18> &SIM) {
        // Normal vector (NV) multiplied by stress influence matrix (SIM, 6*18)
        il::StaticArray2D<double, 3, 6> NM{0.0};
        NM(0, 0) = NV[0];
        NM(1, 1) = NV[1];
        NM(2, 2) = NV[2];
        NM(0, 3) = NV[1];
        NM(1, 3) = NV[0];
        NM(0, 4) = NV[2];
        NM(2, 4) = NV[0];
        NM(1, 5) = NV[2];
        NM(2, 5) = NV[1];
        il::StaticArray2D<double, 3, 18> TIM = il::dot(NM, SIM);
        return TIM;
    };

    il::StaticArray2D<double, 6, 18> SIM_P_R
            (const il::StaticArray2D<double, 3, 3> &RT_L,
             const il::StaticArray2D<double, 3, 3> &RT_R,
             const il::StaticArray2D<double, 6, 18> &SIM) {
        // Triple product (RT_L dot S dot RT_R)
        // for stress influence matrix (SIM, 6*18)
        il::StaticArray2D<double, 3, 3> STM, STM_I, STM_R;
        il::StaticArray2D<double, 6, 18> SIM_R{0.0};
        int j, k, l, m, n;
        for (k = 0; k < SIM.size(1); ++k) {
            for (j = 0; j < 3; ++j) {
                l = (j + 1) % 3;
                m = (l + 1) % 3;
                n = 3 + m;
                STM(j, j) = SIM(j, k);
                STM(l, m) = SIM(n, k);
                STM(m, l) = STM(l, m);
            }
            STM_I = il::dot(STM, RT_R);
            STM_R = il::dot(RT_L, STM_I);
            for (j = 0; j < 3; ++j) {
                l = (j + 1) % 3;
                m = (l + 1) % 3;
                n = 3 + m;
                SIM_R(j, k) = STM_R(j, j);
                SIM_R(n, k) = STM_R(l, m); // STM_R(m, l);
            }
        }
        return SIM_R;
    };

// Matrix-submatrix operations

    template<typename T_sub, typename T_A>
    T_sub get_submatrix
            (const T_A &A,
             il::int_t i0, il::int_t i1,
             il::int_t j0, il::int_t j1) {
        T_sub sub;
        IL_ASSERT((i1 - i0 + 1) == sub.size(0));
        IL_ASSERT((j1 - j0 + 1) == sub.size(1));
        IL_ASSERT(i0 <= A.size(0));
        IL_ASSERT(j0 <= A.size(1));
        IL_ASSERT(i1 <= A.size(0));
        IL_ASSERT(j1 <= A.size(1));

        for (il::int_t i = i0; i <= i1; ++i) {
            for (il::int_t j = j0; j <= j1; ++j) {
                sub(i - i0, j - j0) = A(i, j);
            }
        }
        return sub;
    };

    template<typename T_sub, typename T_A>
    void set_submatrix
            (const T_sub &B,
             il::int_t i0, il::int_t i1,
             il::io_t, T_A &A) {
        IL_ASSERT(i0 + B.size(0) <= A.size(0));
        IL_ASSERT(i1 + B.size(1) <= A.size(1));

        for (il::int_t j1 = 0; j1 < B.size(1); ++j1) {
            for (il::int_t j0 = 0; j0 < B.size(0); ++j0) {
                A(i0 + j0, i1 + j1) = B(j0, j1);
            }
        }
    };

    template <typename T_sub, typename T_A>
    void add_submatrix
            (const T_sub& B, double alpha,
             il::int_t i0, il::int_t i1,
             il::io_t, T_A& A) {
        IL_ASSERT(i0 + B.size(0) <= A.size(0));
        IL_ASSERT(i1 + B.size(1) <= A.size(1));

        for (il::int_t j1 = 0; j1 < B.size(1); ++j1) {
            for (il::int_t j0 = 0; j0 < B.size(0); ++j0) {
                A(i0 + j0, i1 + j1) += alpha*B(j0, j1);
            }
        }
    };
}
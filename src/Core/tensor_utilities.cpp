//
// This file is part of HFPx3D.
//
// Created by D. Nikolski on 1/24/2017.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//

//#include <il/Array2D.h>
#include <il/StaticArray.h>
#include <il/StaticArray2D.h>
#include <il/linear_algebra.h>
// #include <il/linear_algebra/dense/blas/dot.h>
#include "tensor_utilities.h"

namespace hfp3d {
// Vector and triple tensor multiplication
// for stress stored as 6-component vector (or 6*N matrix)

    il::StaticArray2D<double, 3, 18> nv_dot_sim
            (const il::StaticArray<double, 3> &nv,
             const il::StaticArray2D<double, 6, 18> &sim) {
        // Normal vector (nv) multiplied by stress influence matrix (sim, 6*18)
        il::StaticArray2D<double, 3, 6> n_m{0.0};
        n_m(0, 0) = nv[0];
        n_m(1, 1) = nv[1];
        n_m(2, 2) = nv[2];
        n_m(0, 3) = nv[1];
        n_m(1, 3) = nv[0];
        n_m(0, 4) = nv[2];
        n_m(2, 4) = nv[0];
        n_m(1, 5) = nv[2];
        n_m(2, 5) = nv[1];
        il::StaticArray2D<double, 3, 18> tim = il::dot(n_m, sim);
        return tim;
    }

    il::StaticArray<double, 3> nv_dot_sim
            (const il::StaticArray<double, 3> &nv,
             const il::StaticArray<double, 6> &svf) {
        // Normal vector (nv) multiplied by stress in vector (svf)
        il::StaticArray2D<double, 3, 6> n_m{0.0};
        n_m(0, 0) = nv[0];
        n_m(1, 1) = nv[1];
        n_m(2, 2) = nv[2];
        n_m(0, 3) = nv[1];
        n_m(1, 3) = nv[0];
        n_m(0, 4) = nv[2];
        n_m(2, 4) = nv[0];
        n_m(1, 5) = nv[2];
        n_m(2, 5) = nv[1];
        il::StaticArray<double, 3> trac = il::dot(n_m, svf);
        return trac;
    }

    il::StaticArray2D<double, 6, 18> rotate_sim
            (const il::StaticArray2D<double, 3, 3> &rt,
             const il::StaticArray2D<double, 6, 18> &sim) {
        // Triple product (rt_left dot S dot rt_right)
        // for stress influence matrix (sim, 6*18)
        il::StaticArray2D<double, 3, 3> sm_3x3, sm_3x3_interm, sm_3x3_rotated;
        il::StaticArray2D<double, 6, 18> sim_rotated{0.0};
        for (int k = 0; k < sim.size(1); ++k) {
            for (int j = 0; j < 3; ++j) {
                int l = (j + 1) % 3;
                int m = (l + 1) % 3;
                int n = 3 + m;
                sm_3x3(j, j) = sim(j, k);
                sm_3x3(l, m) = sim(n, k);
                sm_3x3(m, l) = sm_3x3(l, m);
            }
            sm_3x3_interm = il::dot(sm_3x3, rt);
            sm_3x3_rotated = il::dot(rt, il::Blas::Transpose, sm_3x3_interm);
            // il::Blas::kTranspose
            for (int j = 0; j < 3; ++j) {
                int l = (j + 1) % 3;
                int m = (l + 1) % 3;
                int n = 3 + m;
                sim_rotated(j, k) = sm_3x3_rotated(j, j);
                sim_rotated(n, k) = sm_3x3_rotated(l, m); // sm_3x3_R(m, l);
            }
        }
        return sim_rotated;
    }

    il::StaticArray2D<double, 6, 18> rotate_sim_c
            (const il::StaticArray2D<double, 3, 3> &rt_left,
             const il::StaticArray2D<double, 3, 3> &rt_right,
             const il::StaticArray2D<double, 6, 18> &sim) {
        // Triple product (rt_left dot S dot rt_right)
        // for stress influence matrix (sim, 6*18)
        il::StaticArray2D<double, 3, 3> sm_3x3, sm_3x3_interm, sm_3x3_rotated;
        il::StaticArray2D<double, 6, 18> sim_rotated{0.0};
        for (int k = 0; k < sim.size(1); ++k) {
            for (int j = 0; j < 3; ++j) {
                int l = (j + 1) % 3;
                int m = (l + 1) % 3;
                int n = 3 + m;
                sm_3x3(j, j) = sim(j, k);
                sm_3x3(l, m) = sim(n, k);
                sm_3x3(m, l) = sm_3x3(l, m);
            }
            sm_3x3_interm = il::dot(sm_3x3, rt_right);
            sm_3x3_rotated = il::dot(rt_left, sm_3x3_interm);
            for (int j = 0; j < 3; ++j) {
                int l = (j + 1) % 3;
                int m = (l + 1) % 3;
                int n = 3 + m;
                sim_rotated(j, k) = sm_3x3_rotated(j, j);
                sim_rotated(n, k) = sm_3x3_rotated(l, m); // sm_3x3_R(m, l);
            }
        }
        return sim_rotated;
    }

// Matrix-submatrix operations

    template <typename T_sub, typename T_A>
    T_sub get_submatrix
            (const T_A &a,
             il::int_t i0, il::int_t i1,
             il::int_t j0, il::int_t j1) {
        T_sub sub;
        IL_EXPECT_FAST((i1 - i0 + 1) == sub.size(0));
        IL_EXPECT_FAST((j1 - j0 + 1) == sub.size(1));
        IL_EXPECT_FAST(i0 <= a.size(0));
        IL_EXPECT_FAST(j0 <= a.size(1));
        IL_EXPECT_FAST(i1 <= a.size(0));
        IL_EXPECT_FAST(j1 <= a.size(1));

        for (il::int_t i = i0; i <= i1; ++i) {
            for (il::int_t j = j0; j <= j1; ++j) {
                sub(i - i0, j - j0) = a(i, j);
            }
        }
        return sub;
    }

    template <typename T_sub, typename T_A>
    void set_submatrix
            (const T_sub &sub,
             il::int_t i0, il::int_t i1,
             il::io_t, T_A &a) {
        IL_EXPECT_FAST(i0 + sub.size(0) <= a.size(0));
        IL_EXPECT_FAST(i1 + sub.size(1) <= a.size(1));

        for (il::int_t j1 = 0; j1 < sub.size(1); ++j1) {
            for (il::int_t j0 = 0; j0 < sub.size(0); ++j0) {
                a(i0 + j0, i1 + j1) = sub(j0, j1);
            }
        }
    }

    template <typename T_sub, typename T_A>
    void add_submatrix
            (const T_sub &sub, double alpha,
             il::int_t i0, il::int_t i1,
             il::io_t, T_A &a) {
        IL_EXPECT_FAST(i0 + sub.size(0) <= a.size(0));
        IL_EXPECT_FAST(i1 + sub.size(1) <= a.size(1));

        for (il::int_t j1 = 0; j1 < sub.size(1); ++j1) {
            for (il::int_t j0 = 0; j0 < sub.size(0); ++j0) {
                a(i0 + j0, i1 + j1) += alpha*sub(j0, j1);
            }
        }
    }
}
//
// This file is part of 3d_bem.
//
// Created by nikolski on 1/12/2017.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//

// Some tensor operations in 6*N format: rotation, multiplication by normal
// For 3*N matrix format simply use il::dot(R, T)

#include <il/StaticArray.h>
#include <il/StaticArray2D.h>
#include <il/linear_algebra/dense/blas/blas.h>
#include <il/linear_algebra/dense/blas/dot.h>
//#include <il/norm.h>
#include <cmath>

il::StaticArray2D<double, 3, 18> N_dot_SIM(il::StaticArray<double, 3> NV, il::StaticArray2D<double, 6, 18> SIM) {
    // Normal vector (NV) multiplied by stress influence matrix (SIM, 6*18)
    il::StaticArray2D<double, 3, 6> NM{0.0};
    NM(0, 0) = NV[0]; NM(1, 1) = NV[1]; NM(2, 2) = NV[2];
    NM(0, 3) = NV[1]; NM(1, 5) = NV[0];
    NM(0, 4) = NV[2]; NM(2, 5) = NV[0];
    NM(1, 5) = NV[2]; NM(2, 5) = NV[1];
    il::StaticArray2D<double, 3, 18> TIM = il::dot(NM, SIM);
    return TIM;
}

il::StaticArray2D<double, 6, 18> SIM_P_R(il::StaticArray2D<double, 3, 3> RTl, il::StaticArray2D<double, 3, 3> RTr, il::StaticArray2D<double, 6, 18> SIM) {
    // Triple product (RTl dot S dot RTr) for stress influence matrix (SIM, 6*18)
    il::StaticArray2D<double, 3, 3> STM, STM_I, STM_R;
    il::StaticArray2D<double, 6, 18> SIM_R{0.0};
    for (int k=0; k<SIM.size(1); ++k){
        for (int j=0; j<3; ++j) {
            int l = (j+1)%3;
            int m = (l+1)%3;
            int n = 3+m;
            STM(j, j) = SIM(j, k);
            STM(l, m) = SIM(n, k);
            STM(m, l) = STM(l, m);
        }
        STM_I = il::dot(STM, RTr);
        STM_R = il::dot(RTl, STM_I);
        for (int j=0; j<3; ++j) {
            int l = (j+1)%3;
            int m = (l+1)%3;
            int n = 3+m;
            SIM_R(j, k) = STM_R(j, j);
            SIM_R(n, k) = STM_R(l, m); // STM_R(m, l);
        }
    }
    return SIM_R;
}

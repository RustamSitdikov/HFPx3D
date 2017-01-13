//
// This file is part of 3d_bem.
//
// Created by nikolski on 1/12/2017.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//

#ifndef INC_3D_BEM_TENSOR_OPER_H
#define INC_3D_BEM_TENSOR_OPER_H

#endif //INC_3D_BEM_TENSOR_OPER_H

il::StaticArray2D<double, 3, 18> N_dot_SIM(il::StaticArray<double,3>, il::StaticArray2D<double,6,18>);
il::StaticArray2D<double, 6, 18> SIM_P_R(il::StaticArray2D<double,3,3>, il::StaticArray2D<double,3,3>, il::StaticArray2D<double,6,18>);

il::StaticArray2D<double, 3, 18> N_dot_SIM(il::StaticArray<double, 3> NV, il::StaticArray2D<double, 6, 18> SIM) {
    // Normal vector (NV) multiplied by stress influence matrix (SIM, 6*18)
    il::StaticArray2D<double, 3, 6> NM{0.0};
    NM(0, 0) = NV[0]; NM(1, 1) = NV[1]; NM(2, 2) = NV[2];
    NM(0, 3) = NV[1]; NM(1, 3) = NV[0];
    NM(0, 4) = NV[2]; NM(2, 4) = NV[0];
    NM(1, 5) = NV[2]; NM(2, 5) = NV[1];
    il::StaticArray2D<double, 3, 18> TIM = il::dot(NM, SIM);
    return TIM;
}

il::StaticArray2D<double, 6, 18> SIM_P_R(il::StaticArray2D<double, 3, 3> RT_L, il::StaticArray2D<double, 3, 3> RT_R, il::StaticArray2D<double, 6, 18> SIM) {
    // Triple product (RT_L dot S dot RT_R) for stress influence matrix (SIM, 6*18)
    il::StaticArray2D<double, 3, 3> STM, STM_I, STM_R;
    il::StaticArray2D<double, 6, 18> SIM_R{0.0};
    int j, k, l, m, n;
    for (k=0; k<SIM.size(1); ++k){
        for (j=0; j<3; ++j) {
            l = (j+1)%3;
            m = (l+1)%3;
            n = 3+m;
            STM(j, j) = SIM(j, k);
            STM(l, m) = SIM(n, k);
            STM(m, l) = STM(l, m);
        }
        STM_I = il::dot(STM, RT_R);
        STM_R = il::dot(RT_L, STM_I);
        for (j=0; j<3; ++j) {
            l = (j+1)%3;
            m = (l+1)%3;
            n = 3+m;
            SIM_R(j, k) = STM_R(j, j);
            SIM_R(n, k) = STM_R(l, m); // STM_R(m, l);
        }
    }
    return SIM_R;
}

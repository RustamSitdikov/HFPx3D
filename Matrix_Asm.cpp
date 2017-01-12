//
// This file is part of 3d_bem.
//
// Created by nikolski on 1/12/2017.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//

// Matrix assembly for the hypersingular BEM (DDM)
// on a triangular boundary mesh with 2nd order
// polynomial approximation of unknowns

#include <il/Array2D.h>
#include <il/StaticArray.h>
#include <il/StaticArray2D.h>
#include <Ele_Base.h>
#include <Local_IM.h>
#include <Submatrix.h>
#include <Tensor_Oper.h>
#include <il/linear_algebra/dense/blas/dot.h>
#include <complex>

void BEMatrix_H_S(il::StaticArray2D& IM_H, double Mu, double Nu, double beta, il::Array2D& Conn_Mtr, il::Array2D& Node_Crd) {
    // Naive BEM matrix assembly from boundary mesh geometry data:
    // mesh connectivity (Conn_Mtr) and nodes' coordinates (Node_Crd)

    IL_ASSERT(Conn_Mtr.size(0) >= 3);
    IL_ASSERT(Conn_Mtr.size(1) >= 1); // at least 1 element
    IL_ASSERT(Node_Crd.size(0) >= 3);
    IL_ASSERT(Node_Crd.size(1) >= 3); // at least 3 nodes

    long Num_El = Conn_Mtr.size(1);

    IL_ASSERT(IM_H.size(0) == 18*Num_El);
    IL_ASSERT(IM_H.size(1) == 18*Num_El);

    //il::StaticArray2D<double, 18*Num_El, 18*Num_El> IM_H;
    //il::StaticArray<double, 18*Num_El> RHS;

    il::StaticArray2D<std::complex<double>, 6, 6> SFM{0.0};
    //il::StaticArray<double, 6> VW_S, VW_T;
    il::StaticArray<std::complex<double>, 3> tau;
    il::StaticArray2D<double, 3, 3> EV_S, EV_T;
    il::StaticArray2D<double, 3, 3> RT_S, RT_S_t, RT_T;
    il::StaticArray2D<double, 3, 3> TI_NN, TI_NN_I, TI_NN_G;
    il::StaticArray<il::StaticArray<double, 3>, 6> CP_T;
    il::StaticArray<double, 3> N_CP, N_CP_L;
    il::Array2D<double> V;
    el_x_cr hz;
    il::StaticArray2D<double, 6, 18> S_H_CP_L, S_H_CP_G;
    il::StaticArray2D<double, 3, 18> T_H_CP_L, T_H_CP_G;
    il::StaticArray2D<double, 18, 18> IM_H_L;
    //il::StaticArray<double, 18> IM_T_L;

    for (long S_El=0; S_El<Num_El; ++S_El) {
        // "Source" element
        for (int j=0; j<3; ++j) {
            get_submatrix(V,0,2,Conn_Mtr(j, S_El),Conn_Mtr(j, S_El),Node_Crd);
            set_submatrix_2_static(EV_S,0,j,V);
            // set VW_S[j]
        }
        // Rotational tensor and basis (shape) functions for the source element
        SFM = El_SFM_S(RT_S, EV_S);
        //SFM = El_SFM_N(RT, EV, VW_S);
        El_RT_Tr(tau, RT_S_t, RT_S, EV_S);
        for (long T_El=0; T_El<Num_El; ++T_El) {
            // "Target" element
            for (int j=0; j<3; ++j) {
                get_submatrix(V,0,2,Conn_Mtr(j, T_El),Conn_Mtr(j, T_El),Node_Crd);
                set_submatrix_2_static(EV_T,0,j,V);
                // set VW_T[j]
            }
            // Rotational tensor for the target element
            El_LB_RT(RT_T, EV_T);
            // Normal vector
            for (int k=0; k<3; ++k) {
                N_CP[k] = -RT_T(k,2);
            }
            // Collocation points' coordinates
            CP_T = El_CP_S(EV_T,beta);
            //CP_T = El_CP_N(EV_T,VW_T,beta);
            for (int N_T=0; N_T<6; ++N_T) {
                // Shifting to the N_T-th collocation pt
                El_X_CR(hz, RT_S_t, EV_S, CP_T[N_T]);
                // Calculating DD-to stress influence w.r. to the source element's local coordinate system
                S_H_CP_L = Local_IM_B_H(Mu, Nu, hz.h, hz.z, tau, SFM);
                // Multiplication by N_CP
                // Alternative 1: rotating stress at CP to the reference ("global") coordinate system
                //S_H_CP_G = SIM_P_R(RT_S, RT_S_t, S_H_CP_L);
                //T_H_CP_G = N_dot_SIM(N_CP, S_H_CP_G);
                // Alternative 2: rotating N_CP to the source element's local coordinate system
                N_CP_L = il::dot(RT_S_t, N_CP);
                T_H_CP_L = N_dot_SIM(N_CP_L, S_H_CP_L);
                T_H_CP_G = il::dot(RT_S, T_H_CP_L);
                for (int N_S=0; N_S<6; ++N_S) {
                    get_static_sub_from_static(TI_NN,0,3*N_S,T_H_CP_G);
                    // Re- to traction vs DD w.r. to the reference coordinate system
                    TI_NN_G = il::dot(TI_NN, RT_S_t);
                    // Adding the block to element-to-element influence sub-matrix
                    set_static_sub_2_static(IM_H_L,3*N_T,3*N_S,TI_NN_G);
                }
            }
            // Adding the block to the global influence matrix
            set_static_sub_2_static(IM_H,18*T_El,18*S_El,IM_H_L);
        }
    }
    //return &IM_H;
}

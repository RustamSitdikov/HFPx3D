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

#include <complex>
#include <il/Array.h>
#include <il/Array2D.h>
#include <il/StaticArray.h>
#include <il/StaticArray2D.h>
#include <il/linear_algebra/dense/blas/dot.h>
#include <Ele_Base.h>
#include <Submatrix.h>
#include <Tensor_Oper.h>
//#include <ICFns.h>

template <class Kernel, class T, typename M_array>
T BEMatrix_S(double Mu, double Nu, double beta, M_array Conn_Mtr, M_array Node_Crd) {
    // BEM matrix assembly from boundary mesh geometry data:
    // mesh connectivity (Conn_Mtr) and nodes' coordinates (Node_Crd)
    // Naive way: no parallelization, no ACA

//    IL_ASSERT(Conn_Mtr.size(0) >= 3);
//    IL_ASSERT(Conn_Mtr.size(1) >= 1); // at least 1 element
//    IL_ASSERT(Node_Crd.size(0) >= 3);
//    IL_ASSERT(Node_Crd.size(1) >= 3); // at least 3 nodes

    const il::int_t Num_El = Conn_Mtr.size(1);
    il::int_t S_N, T_N, S_El, T_El;
    int j, n_S, n_T;

    //IL_ASSERT(IM_H.size(0) == 18*Num_El);
    //IL_ASSERT(IM_H.size(1) == 18*Num_El);

    il::StaticArray2D<double, 18*Num_El, 18*Num_El> IM_H;
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

    for (S_El=0; S_El<Num_El; ++S_El) {
        // "Source" element
        for (j=0; j<3; ++j) {
            S_N = Conn_Mtr(j, S_El);
            get_submatrix<il::Array2D, il::Array2D>(V, 0, 2, S_N, S_N, Node_Crd);
            set_submatrix<il::Array2D, il::StaticArray2D<double, 3, 3>>(EV_S, il::int_t(0), il::int_t(j), V);
            // set VW_S[j]
        }
        // Rotational tensor and basis (shape) functions for the source element
        SFM = El_SFM_S(RT_S, EV_S);
        //SFM = El_SFM_N(RT, EV, VW_S);
        El_RT_Tr(tau, RT_S_t, RT_S, EV_S);
        for (T_El=0; T_El<Num_El; ++T_El) {
            // "Target" element
            for (j=0; j<3; ++j) {
                T_N = Conn_Mtr(j, T_El);
                get_submatrix<il::Array2D, il::Array2D>(V, 0, 2, T_N, T_N, Node_Crd);
                set_submatrix<il::Array2D, il::StaticArray2D<double, 3, 3>>(EV_T, 0, j, V);
                // set VW_T[j]
            }
            // Rotational tensor for the target element
            El_LB_RT(RT_T, EV_T);
            // Normal vector
            for (j=0; j<3; ++j) {
                N_CP[j] = -RT_T(j, 2);
            }
            // Collocation points' coordinates
            CP_T = El_CP_S(EV_T, beta);
            //CP_T = El_CP_N(EV_T,VW_T,beta);
            for (n_T=0; n_T<6; ++n_T) {
                // Shifting to the n_T-th collocation pt
                El_X_CR(hz, RT_S_t, EV_S, CP_T[n_T]);
                // Calculating DD-to stress influence w.r. to the source element's local coordinate system
                S_H_CP_L = Local_IM<Kernel>(Mu, Nu, hz.h, hz.z, tau, SFM);
                // Multiplication by N_CP
                // Alternative 1: rotating stress at CP to the reference ("global") coordinate system
                //S_H_CP_G = SIM_P_R(RT_S, RT_S_t, S_H_CP_L);
                //T_H_CP_G = N_dot_SIM(N_CP, S_H_CP_G);
                // Alternative 2: rotating N_CP to the source element's local coordinate system
                N_CP_L = il::dot(RT_S_t, N_CP);
                T_H_CP_L = N_dot_SIM(N_CP_L, S_H_CP_L);
                T_H_CP_G = il::dot(RT_S, T_H_CP_L);
                for (n_S=0; n_S<6; ++n_S) {
                    get_submatrix<il::StaticArray2D>(TI_NN, 0, 2, 3*n_S, 3*n_S+2, T_H_CP_G);
                    // Re- to traction vs DD w.r. to the reference coordinate system
                    TI_NN_G = il::dot(TI_NN, RT_S_t);
                    // Adding the block to element-to-element influence sub-matrix
                    set_submatrix<il::StaticArray2D>(IM_H_L, 3*n_T, 3*n_S, TI_NN_G);
                }
            }
            // Adding the block to the global influence matrix
            set_submatrix<il::StaticArray2D>(IM_H, 18*T_El, 18*S_El, IM_H_L);
        }
    }
    return T(IM_H);
}

template <class Kernel>
il::StaticArray2D<double, 6, 18> Local_IM
        (double mu, double nu, double h, std::complex<double> z,
         il::StaticArray<std::complex<double>, 3> tau, il::StaticArray2D<std::complex<double>, 6, 6> SFM) {
    // This function assembles a local "stiffness" sub-matrix
    // (influence of DD at the element nodes to stresses at the point z)
    // in terms of a triangular element's local coordinates
    //
    // tau (3) are coordinates of element's vertices and
    // the rows of SFM (6*6) are coefficients of shape functions
    // in terms of the element's own local coordinate system (tau-coordinates);
    // h and z define the position of the (collocation) point x in the same coordinates

    const std::complex<double> I(0.0,1.0);

    // scaling ("-" sign comes from traction Somigliana ID, H-term)
    double scale = -mu/(4.0*M_PI*(1.0-nu));
    // tolerance parameters
    const double HTol = 1.0E-9, DTol = 1.0E-8;
    // geometrical
    double an, am;
    std::complex<double> eixm, eixn, eipm, eipn, zc=std::conj(z), ntau2;
    // tz[m] and d[m] can be calculated here
    il::StaticArray<std::complex<double>,3> tz, d, dtau;
    int j, k, l, m, n;
    for (j=0; j<=2; ++j) {
        tz[j] = tau[j] - z;
        l = (j+1)%3;
        dtau[j]=tau[l]-tau[j];
        ntau2=dtau[j]/std::conj(dtau[j]);
        d[j]=0.5*(tz[j]-std::conj(tz[j])*ntau2);
    }
    // also, "shifted" SFM from z, tau[m], and local SFM
    il::StaticArray2D<std::complex<double>, 6, 6> ShiftZ {0.0};
    ShiftZ(0, 0) = 1.0;
    ShiftZ(1, 0) = z; ShiftZ(1, 1) = 1.0;
    ShiftZ(2, 0) = zc; ShiftZ(2, 2) = 1.0;
    ShiftZ(3, 0) = z*z; ShiftZ(3, 1) = 2.0*z; ShiftZ(3, 3) = 1.0;
    ShiftZ(4, 0) = zc*zc; ShiftZ(4, 2) =  2.0*zc; ShiftZ(4, 4) = 1.0;
    ShiftZ(5, 0) = z*zc; ShiftZ(5, 1) = zc; ShiftZ(5, 2) = z; ShiftZ(5, 5) = 1.0;
    il::StaticArray2D<std::complex<double>, 6, 6> SFMz = il::dot(SFM, ShiftZ);

    // constituents of the integrals
    il::StaticArray<std::complex<double>, 9> Fn, Fm;
    il::StaticArray3D<std::complex<double>, 6, 3, 9> Cn, Cm;
    il::StaticArray<std::complex<double>, 5> FnD, FmD;
    il::StaticArray3D<std::complex<double>, 6, 3, 5> CnD, CmD;
    // DD-to-stress influence
    // [(S11+S22)/2; (S11-S22)/2+i*S12; (S13+i*S23)/2; S33]
    // vs SF monomials (Sij_M) and nodal values (Sij_N)
    il::StaticArray2D<std::complex<double>, 6, 3> Sij_M_1 {0.0}, Sij_N_1 {0.0},
            Sij_M_2 {0.0}, Sij_N_2 {0.0},
            Sij_M_3 {0.0}, Sij_N_3 {0.0},
            Sij_M_4 {0.0}, Sij_N_4 {0.0};
    il::StaticArray3D<std::complex<double>, 6, 4, 3> SincLn {0.0}, SincLm {0.0};

    il::StaticArray2D<double, 6, 18> LIM {0.0};

    // searching for "degenerate" edges: point x (collocation pt) projects onto an edge or a vertex
    bool IsDegen = std::abs(d[0])<DTol || std::abs(d[1])<DTol || std::abs(d[2])<DTol; // (d[0]*d[1]*d[2]==0);

    // calculating angles (phi, psi, chi)
    il::StaticArray<double,3> phi {0.0}, psi {0.0};
    il::StaticArray2D<double, 2, 3> chi {0.0};
    for (j=0; j<=2; ++j) {
        phi[j] = std::arg(tz[j]);
        // make sure it's between -pi and pi (add or subtract 2*pi)
        // phi[j] = (std::fmod(phi[j]/M_PI+1.0, 2.0)-1.0)*M_PI;
        psi[j] = std::arg(d[j]);
    }
    for (j=0; j<=2; ++j) {
        for (k=0; k<=1; ++k) {
            l = (j+k)%3;
            chi(k,j) = phi[l]-psi[j];
            // make sure it's between -pi and pi (add or subtract 2*pi)
            //chi(k,j) = (std::fmod(chi(k,j)/M_PI+1.0, 2.0)-1.0)*M_PI;
            chi(k,j) = (chi(k,j)<=-M_PI)? chi(k,j)+M_PI : ((chi(k,j)>M_PI)? chi(k,j)-M_PI : chi(k,j));
            // reprooving for "degenerate" edges
            if (fabs(M_PI_2-std::fabs(chi(k,j)))<DTol) IsDegen = true;
        }
    }

    // summation over edges
    for (m=0; m<=2; ++m) {
        n = (m+1)%3;
        if (std::abs(d[m])>=DTol &&
            std::fabs(M_PI_2-std::fabs(chi(0,m)))>=DTol &&
            std::fabs(M_PI_2-std::fabs(chi(1,m)))>=DTol) {
            eixm = std::exp(I*chi(0,m)); eixn = std::exp(I*chi(1,m));
            if(std::fabs(h)<HTol) { // limit case (point x on the element's plane)
                SincLn = Kernel::SijLim(nu, eixn, d[m]);
                SincLm = Kernel::SijLim(nu, eixm, d[m]);
                for (j=0; j<=5; ++j) {
                    for (k=0; k<=2; ++k) {
                        Sij_M_1(j, k) += SincLn(j, 0, k) - SincLm(j, 0, k);
                        Sij_M_2(j, k) += SincLn(j, 1, k) - SincLm(j, 1, k);
                        Sij_M_3(j, k) += SincLn(j, 2, k) - SincLm(j, 2, k);
                        Sij_M_4(j, k) += SincLn(j, 3, k) - SincLm(j, 3, k);
                    }
                }
            }
            else { // out-of-plane case
                an = std::abs(tz[n]-d[m])*((chi(1,m)<0)? -1.0 : double((chi(1,m)>0)));
                am = std::abs(tz[m]-d[m])*((chi(0,m)<0)? -1.0 : double((chi(0,m)>0)));
                Fn = ICFns(h, d[m], an, chi(1,m), eixn);
                Fm = ICFns(h, d[m], am, chi(0,m), eixm);
                // S11+S22
                Cn = Kernel::S11_22(nu, eixn, h, d[m]);
                Cm = Kernel::S11_22(nu, eixm, h, d[m]);
                il::blas(1.0, Cn, Fn, 1.0, il::io, Sij_M_1);
                il::blas(-1.0, Cm, Fm, 1.0, il::io, Sij_M_1);
                if (IsDegen) { // "degenerate" case
                    eipm = std::exp(I*phi[n]); eipn = std::exp(I*phi[m]);
                    FnD = ICFns_red(h, d[m], an);
                    FmD = ICFns_red(h, d[m], am);
                    CnD = Kernel::S11_22_red(nu, eipn, h, d[m]);
                    CmD = Kernel::S11_22_red(nu, eipm, h, d[m]);
                    il::blas(1.0, CnD, FnD, 1.0, il::io, Sij_M_1);
                    il::blas(-1.0, CmD, FmD, 1.0, il::io, Sij_M_1);
                }
                // S11-S22+2*I*S12
                Cn = Kernel::S11_22_12(nu, eixn, h, d[m]);
                Cm = Kernel::S11_22_12(nu, eixm, h, d[m]);
                il::blas(1.0, Cn, Fn, 1.0, il::io, Sij_M_2);
                il::blas(-1.0, Cm, Fm, 1.0, il::io, Sij_M_2);
                if (IsDegen) { // "degenerate" case
                    CnD = Kernel::S11_22_12_red(nu, eipn, h, d[m]);
                    CmD = Kernel::S11_22_12_red(nu, eipm, h, d[m]);
                    il::blas(1.0, CnD, FnD, 1.0, il::io, Sij_M_2);
                    il::blas(-1.0, CmD, FmD, 1.0, il::io, Sij_M_2);
                }
                // S13+I*S23
                Cn = Kernel::S13_23(nu, eixn, h, d[m]);
                Cm = Kernel::S13_23(nu, eixm, h, d[m]);
                il::blas(1.0, Cn, Fn, 1.0, il::io, Sij_M_3);
                il::blas(-1.0, Cm, Fm, 1.0, il::io, Sij_M_3);
                if (IsDegen) { // "degenerate" case
                    CnD = Kernel::S13_23_red(nu, eipn, h, d[m]);
                    CmD = Kernel::S13_23_red(nu, eipm, h, d[m]);
                    il::blas(1.0, CnD, FnD, 1.0, il::io, Sij_M_3);
                    il::blas(-1.0, CmD, FmD, 1.0, il::io, Sij_M_3);
                }
                // S33
                Cn = Kernel::S33(nu, eixn, h, d[m]);
                Cm = Kernel::S33(nu, eixm, h, d[m]);
                il::blas(1.0, Cn, Fn, 1.0, il::io, Sij_M_4);
                il::blas(-1.0, Cm, Fm, 1.0, il::io, Sij_M_4);
                if (IsDegen) { // "degenerate" case
                    CnD = Kernel::S33_red(nu, eipn, h, d[m]);
                    CmD = Kernel::S33_red(nu, eipm, h, d[m]);
                    il::blas(1.0, CnD, FnD, 1.0, il::io, Sij_M_4);
                    il::blas(-1.0, CmD, FmD, 1.0, il::io, Sij_M_4);
                }
            }
        }
    }

    // here comes contraction with "shifted" SFM (left)
    Sij_N_1 = il::dot(SFMz, Sij_M_1);
    Sij_N_2 = il::dot(SFMz, Sij_M_2);
    Sij_N_3 = il::dot(SFMz, Sij_M_3);
    Sij_N_4 = il::dot(SFMz, Sij_M_4);

    // re-shaping of the resulting matrix
    // and scaling (comment out if not necessary)
    for (j=0; j<=5; ++j) {
        int q = j*3;
        for (k=0; k<=2; ++k) {
            // [S11; S22; S33; S12; S13; S23] vs \delta{u}_k at j-th node
            LIM(0,q+k) = scale*(std::real(Sij_N_1(j,k))+std::real(Sij_N_2(j,k))); // S11
            LIM(1,q+k) = scale*(std::real(Sij_N_1(j,k))-std::real(Sij_N_2(j,k))); // S22
            LIM(2,q+k) = scale*std::real(Sij_N_4(j,k)); // S33
            LIM(3,q+k) = scale*std::imag(Sij_N_2(j,k)); // S12
            LIM(4,q+k) = scale*2.0*std::real(Sij_N_3(j,k)); // S13
            LIM(5,q+k) = scale*2.0*std::imag(Sij_N_3(j,k)); // S23
        }
    }

    return LIM;
}

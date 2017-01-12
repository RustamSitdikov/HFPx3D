#include <cstdio>

#include <il/Array.h>
#include <il/Array2D.h>
#include <il/StaticArray.h>
#include <il/StaticArray2D.h>
#include <il/StaticArray3D.h>
#include <il/linear_algebra/dense/blas/dot.h>
//#include <il/linear_algebra/dense/blas/blas.h>
#include <il/io/numpy.h>
#include <complex>
//#include <cmath>
//#include <SijH.h>
//#include <ICFns.h>
#include <Ele_Base.h>
#include <Local_IM.h>
#include <Matrix_asm.h>

// main.cpp will be used for testing the code parts under development

int main() {
    std::string MeshDirectory{"C:/Users/nikolski/.spyder-py3/3DBEM"};
    std::string WorkDirectory{"C:/Users/nikolski/ClionProjects/3D-bem"};

    il::Status status{};
    il::Array2D<int> Conn_Mtr = il::load<il::Array2D<int>>(MeshDirectory + std::string{"/Elems_pennymesh24el.npy"}, il::io, status);
    status.abort_on_error();

    il::Array2D<double> Node_Crd = il::load<il::Array2D<double>>(MeshDirectory + std::string{"/Nodes_pennymesh24el.npy"}, il::io, status);
    status.abort_on_error();

    //const std::complex<double> I(0.0,1.0);

    double Mu = 1.0, Nu = 0.35, h0 = -1.21;
    std::complex<double> z0(1.0,0.4), z1(0.0,0.1), z2(1.8,0.0), z3(1.2,1.8);
    il::StaticArray<std::complex<double>, 3> tau;
    il::StaticArray<double, 3> X0; //, V0, X0r, NV;
    X0[0] = std::real(z0); X0[1] = std::imag(z0); X0[2] = -h0;
    il::StaticArray2D<double, 3, 3> EV, RT, RTt; //, EVr{0.0};
    //il::StaticArray<double, 3> VW;
    tau[0] = z1; tau[1] = z2; tau[2] = z3;
    for (int k = 0; k<=2; ++k) {
        EV(0, k) = std::real(tau[k]);
        EV(1, k) = std::imag(tau[k]);
        EV(2, k) = 0.0;
        //VW[k] = 1.0;
    }
    il::StaticArray2D<std::complex<double>, 6, 6> SFM_1{0.0};
    el_x_cr hz;
    il::StaticArray2D<double, 6, 18> L_IM_H_1;
    // tau, RT, RTt, and SFM are calculated for an element
    SFM_1 = El_SFM_S(RT, EV);
    //SFM_1 = El_SFM_N(RT, EV, VW);
    El_RT_Tr(tau, RTt, RT, EV);
    // the rest is recalculated for every CP
    El_X_CR(hz, RTt, EV, X0);
    L_IM_H_1 = Local_IM_H(Mu, Nu, hz.h, hz.z, tau, SFM_1);

    //for (int j=0; j<=2; ++j) {
    //    for (int k=0; k<=2; ++k) {
    //        RTt(k,j) = RT(j,k);
    //    }
    //}
    //for (int k = 0; k<=2; ++k) {
    //    V0[k] = EV(k, 0);
    //    NV[k] = -RT(k, 2);
    //    X0r[k] = X0[k] - V0[k];
    //    EVr(k, 1) = EV(k, 1) - V0[k];
    //    EVr(k, 2) = EV(k, 2) - V0[k];
    //}
    //EVr = il::dot(RTt,EVr);
    //X0r = il::dot(RTt,X0r);
    //z0 = std::complex<double>(X0r[0], X0r[1]);
    //for (int j=0; j<=2; ++j) {
    //    tau[j] = std::complex<double>(EVr(0,j), EVr(1,j));
    //}
    //L_IM_H_1 = Local_IM_B_H(Mu, Nu, h0, z0, tau, SFM_1);

    std::string path = WorkDirectory+std::string{"/Int_test_1_ele.csv"};
    FILE* of=std::fopen(path.c_str(),"w");
    for (int j=0; j<L_IM_H_1.size(0); ++j){
        for (int k=0; k<L_IM_H_1.size(1); ++k){
            std::fprintf(of,"%.12g",L_IM_H_1(j,k));
            if (k<L_IM_H_1.size(1)-1) std::fprintf(of,",");
        }
        std::fprintf(of,"\n");
    }
    std::fclose(of);

    path = WorkDirectory+std::string{"/SFM_test_1_ele.csv"};
    of=std::fopen(path.c_str(),"w");
    for (int j=0; j<SFM_1.size(0); ++j){
        for (int k=0; k<SFM_1.size(1); ++k){
            std::fprintf(of,"%.12g%+.12g*I",SFM_1(j,k));
            if (k<SFM_1.size(1)-1) std::fprintf(of,",");
        }
        std::fprintf(of,"\n");
    }
    std::fclose(of);

    //il::save(L_IM_H_1, "C:/Users/nikolski/.spyder-py3/3DBEM/matrix_1_el.npy", il::io, status);
    //status.abort_on_error();

    return 0;
}
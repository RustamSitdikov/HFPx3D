#include <cstdio>

//#include <il/Array.h>
#include <il/Array2D.h>
#include <il/StaticArray.h>
#include <il/StaticArray2D.h>
//#include <il/StaticArray3D.h>
//#include <il/linear_algebra/dense/blas/dot.h>
//#include <il/linear_algebra/dense/blas/blas.h>
#include <il/io/numpy.h>
//#include <complex>
//#include <cmath>
//#include <SijH.h>
//#include <Ele_Base.h>
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

    il::Array2D<int> Ele(3, 1); //Ele.size[0] = 3; Ele.size[1] = 1;
    il::Array2D<double> Nod(3, 3);

    for (int n=0; n<3; ++n) {
        Ele(n, 0) = n;
        Nod(2, n) = 0.0;
    }
    Nod(0, 0) = std::real(z1); Nod(0, 1) = std::real(z2); Nod(0, 2) = std::real(z3);
    Nod(1, 0) = std::imag(z1); Nod(1, 1) = std::imag(z2); Nod(1, 2) = z1.imag();

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
    H_Int H_I;
    L_IM_H_1 = Local_IM<H_Int>(H_I, Mu, Nu, hz.h, hz.z, tau, SFM_1);

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

    std::string path = WorkDirectory+std::string{"/test_integration_1_ele.csv"};
    FILE* of=std::fopen(path.c_str(),"w");
    for (int j=0; j<L_IM_H_1.size(0); ++j){
        for (int k=0; k<L_IM_H_1.size(1); ++k){
            std::fprintf(of,"%18.12g",L_IM_H_1(j,k));
            if (k<L_IM_H_1.size(1)-1) std::fprintf(of,",");
        }
        std::fprintf(of,"\n");
    }
    std::fclose(of);

    path = WorkDirectory+std::string{"/test_SFM_1_ele.csv"};
    of=std::fopen(path.c_str(),"w");
    for (int j=0; j<SFM_1.size(0); ++j){
        for (int k=0; k<SFM_1.size(1); ++k){
            std::fprintf(of,"%18.12g%+.12g*I",SFM_1(j,k));
            if (k<SFM_1.size(1)-1) std::fprintf(of,",");
        }
        std::fprintf(of,"\n");
    }
    std::fclose(of);

    //il::save(L_IM_H_1, "C:/Users/nikolski/.spyder-py3/3DBEM/matrix_1_el.npy", il::io, status);
    //status.abort_on_error();

    il::Array2D<double> IM_1(18, 18);
    IM_1 = BEMatrix_S<H_Int, il::Array2D<int>, il::Array2D<double>>(Mu, Nu, 0.25, Ele, Nod);

    path = WorkDirectory+std::string{"/test_assembly_1_ele.csv"};
    of=std::fopen(path.c_str(),"w");
    for (int j=0; j<IM_1.size(0); ++j){
        for (int k=0; k<IM_1.size(1); ++k){
            std::fprintf(of,"%18.12g",IM_1(j,k));
            if (k<IM_1.size(1)-1) std::fprintf(of,",");
        }
        std::fprintf(of,"\n");
    }
    std::fclose(of);

    //il::save(IM_1, "C:/Users/nikolski/.spyder-py3/3DBEM/matrix_1_el.npy", il::io, status);
    //status.abort_on_error();

    il::int_t N_El = Conn_Mtr.size(1), N_DOF = 18*N_El;
    // conversion from Matlab to C++
    for (int n=0; n<N_El; ++n) {
        for (int j=0; j<3; ++j) {
            Conn_Mtr(j, n) -=1;
        }
    }
    il::Array2D<double> IM_2(N_DOF, N_DOF);
    IM_2 = BEMatrix_S<H_Int, il::Array2D<int>, il::Array2D<double>>(Mu, Nu, 0.25, Conn_Mtr, Node_Crd);

    path = WorkDirectory+std::string{"/test_assembly_24_ele.csv"};
    of=std::fopen(path.c_str(),"w");
    for (int j=0; j<IM_2.size(0); ++j){
        for (int k=0; k<IM_2.size(1); ++k){
            std::fprintf(of,"%18.12g",IM_2(j,k));
            if (k<IM_2.size(1)-1) std::fprintf(of,",");
        }
        std::fprintf(of,"\n");
    }
    std::fclose(of);

    //il::save(IM_1, "C:/Users/nikolski/.spyder-py3/3DBEM/matrix_1_el.npy", il::io, status);
    //status.abort_on_error();

    return 0;
}
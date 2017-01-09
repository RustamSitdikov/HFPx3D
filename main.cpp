#include <cstdio>

#include <il/Array.h>
#include <il/Array2D.h>
#include <il/StaticArray.h>
#include <il/StaticArray2D.h>
#include <il/StaticArray3D.h>
#include <il/linear_algebra/dense/blas/dot.h>
#include <il/linear_algebra/dense/blas/blas.h>
#include <il/io/numpy.h>
#include <complex>
#include <cmath>
#include <SijH.h>
#include <ICFns.h>
#include <Local_IM.h>

// main.cpp will be used for testing the code parts under development

int main() {
    std::string MeshDirectory{"C:/Users/nikolski/.spyder-py3/3DBEM"};
    std::string WorkDirectory{"C:/Users/nikolski/ClionProjects/3D-bem"};

    il::Status status{};
    il::Array2D<int> C=il::load<il::Array2D<int>>(MeshDirectory + std::string{"/Elems_pennymesh24el.npy"}, il::io, status);
    status.abort_on_error();

    il::Array2D<double> D=il::load<il::Array2D<double>>(MeshDirectory + std::string{"/Nodes_pennymesh24el.npy"}, il::io, status);
    status.abort_on_error();

    il::save(D, "C:/Users/nikolski/.spyder-py3/3DBEM/matrix.npy", il::io, status);
    status.abort_on_error();

    const std::complex<double> I(0.0,1.0);
    il::StaticArray2D<double, 2, 3> A{0.0};
    il::StaticArray2D<std::complex<double>, 3, 3> B{0.0};
    il::StaticArray<double, 3> V{0.0};
    il::StaticArray3D<double, 2, 4, 3> W{0.0};

    A(0, 0) = 1.0;
    A(1, 1) = 2.0;
    A(1, 2) = 3.0;
    V[0] = 14.5; V[1] = -6.2; V[2] = 1.8;
    il::StaticArray<double, 2> Y = il::dot(A, V);
    il::blas(1.0, A, V, 1.0, il::io, Y);

    double Nu=0.35, h0=-1.21;
    std::complex<double> z0(1.0,0.4), z1(0.0,0.1), z2(1.8,0.0), z3(1.2,1.8);
    std::complex<double> d1 = 0.5*(z1-z0-(std::conj(z1)-std::conj(z0))*(z2-z1)/(std::conj(z2)-std::conj(z1)));
    std::complex<double> eip2 = std::exp(I*std::arg(z2-z0)), eix12 = std::exp(I*(std::arg(z2-z0)-std::arg(d1))),
          eip1 = std::exp(I*std::arg(z1-z0)), eix11 = std::exp(I*(std::arg(z1-z0)-std::arg(d1)));
    il::StaticArray3D<std::complex<double>,6,3,9> SccH0;
    il::StaticArray3D<std::complex<double>,6,3,5> SccHr;
    SccH0 = S13_23H(Nu, eix12, h0, d1);
    SccHr = S13_23H_red(Nu, eip2, h0, d1);

    //std::cout << d1 << "\n" <<std::arg(z1-z0) <<"\n";

    std::string path = WorkDirectory+std::string{"/inttest1.csv"};
    FILE* of=std::fopen(path.c_str(),"w");
    for (int s=0; s<=5; ++s){
        std::fprintf(of,"0+0*I,0+0*I,0+0*I\n");
        for (int n=0; n<=8; ++n){
            for (int k=0; k<=2; ++k){
                std::fprintf(of,"%.12g%+.12g*I",SccH0(s,k,n));
                if (k<2) std::fprintf(of,",");
            }
            std::fprintf(of,"\n");
        }
    }
    std::fclose(of);
    path = WorkDirectory+std::string{"/inttest2.csv"};
    of=std::fopen(path.c_str(),"w");
    for (int s=0; s<=5; ++s){
        for (int n=0; n<=4; ++n){
            for (int k=0; k<=2; ++k){
                std::fprintf(of,"%.12g%+.12g*I",SccHr(s,k,n));
                if (k<2) std::fprintf(of,",");
            }
            std::fprintf(of,"\n");
        }
    }
    std::fclose(of);

    return 0;
}
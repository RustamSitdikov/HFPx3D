#include <cstdio>
#include <il/Array2D.h>
#include <il/StaticArray.h>
#include <il/StaticArray2D.h>
#include <il/io/numpy.h>
#include "Matrix_Asm.h"
#include "Ele_Base.h"

// main.cpp will be used for testing the code parts under development

int main() {
    std::string SrcDirectory{"C:/Users/nikolski/.spyder-py3/3DBEM"};
    std::string WorkDirectory{"C:/Users/nikolski/ClionProjects/3D-bem/Test_Output"};

    il::Status status{};
    il::Array2D<il::int_t> Conn_Mtr = il::load<il::Array2D<il::int_t>>
            (SrcDirectory + std::string{"/Elems_pennymesh24el.npy"},
             il::io, status);
    status.abort_on_error();

    il::Array2D<double> Node_Crd = il::load<il::Array2D<double>>
            (SrcDirectory + std::string{"/Nodes_pennymesh24el.npy"},
             il::io, status);
    status.abort_on_error();

    il::int_t N_El = Conn_Mtr.size(1), N_DOF = 18*N_El;
    // conversion from Matlab to C++
    for (il::int_t n=0; n<N_El; ++n) {
        for (il::int_t j=0; j<3; ++j) {
            Conn_Mtr(j, n) -=1;
        }
    }

    double Mu = 1.0, Nu = 0.35;

    il::Array2D<double> IM_2(N_DOF, N_DOF);
    // <il::Array2D<il::int_t>, il::Array2D<double>>
    IM_2 = hfp3d::BEMatrix_S(Mu, Nu, 0.25, Conn_Mtr, Node_Crd);

    std::string path = WorkDirectory+std::string{"/test_assembly_24_ele.csv"};
    FILE* of=std::fopen(path.c_str(),"w");
    for (int j=0; j<IM_2.size(0); ++j){
        for (int k=0; k<IM_2.size(1); ++k){
            std::fprintf(of,"%21.16g",IM_2(j,k));
            if (k<IM_2.size(1)-1) std::fprintf(of,",");
        }
        std::fprintf(of,"\n");
    }
    std::fclose(of);

    //il::save(IM_1, "C:/Users/nikolski/.spyder-py3/3DBEM/matrix_24_el.npy",
    // il::io, status);
    //status.abort_on_error();

    return 0;
}
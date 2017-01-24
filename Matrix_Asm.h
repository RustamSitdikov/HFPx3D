//
// This file is part of 3d_bem.
//
// Created by D. Nikolski on 1/12/2017.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//

// Matrix assembly for the hypersingular BEM (DDM)
// on a triangular boundary mesh with 2nd order
// polynomial approximation of unknowns

#ifndef INC_3D_BEM_MATRIX_ASM_H
#define INC_3D_BEM_MATRIX_ASM_H

#include <complex>
#include <il/Array2D.h>
#include <il/StaticArray.h>
#include <il/StaticArray2D.h>

namespace hfp3d {

    // Element-to-point influence matrix (submatrix of the global one)

    template<class Kernel>
    il::StaticArray2D<double, 6, 18> Local_IM
            (Kernel& K_I, double mu, double nu,
             double h, std::complex<double> z,
             const il::StaticArray<std::complex<double>, 3>& tau,
             const il::StaticArray2D<std::complex<double>, 6, 6>& SFM);

    // "Global" matrix assembly

    template<typename C_array, typename N_array>
    il::Array2D<double> BEMatrix_S
            (double Mu, double Nu, double beta,
             C_array& Conn_Mtr,
             N_array& Node_Crd);

}

#endif //INC_3D_BEM_MATRIX_ASM_H

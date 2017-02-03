//
// This file is part of HFPx3D.
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

// "Global" matrix assembly
    il::Array2D<double> make_3dbem_matrix_s
            (double mu, double nu, double beta,
             const il::Array2D<il::int_t> &mesh_conn,
             const il::Array2D<double> &nodes_crd);

// Stress at given points (m_pts_crd) vs DD at nodal points (nodes_crd)
    il::Array2D<double> make_3dbem_stress_f_s
            (double mu, double nu,
             const il::Array2D<il::int_t> &mesh_conn,
             const il::Array2D<double> &nodes_crd,
             const il::Array2D<double> &m_pts_crd,
             //const il::Array2D<double> &m_pts_dsp,
             const bool is_in_glob);

// Element-to-point influence matrix (submatrix of the global one)
    il::StaticArray2D<double, 6, 18> make_local_3dbem_submatrix
            (const int kernel_id,
             double mu, double nu, double h, std::complex<double> z,
             const il::StaticArray<std::complex<double>, 3> &tau,
             const il::StaticArray2D<std::complex<double>, 6, 6> &sfm);

}

#endif //INC_3D_BEM_MATRIX_ASM_H

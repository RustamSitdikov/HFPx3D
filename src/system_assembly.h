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
#include "mesh_utilities.h"

namespace hfp3d {

    struct Alg_Sys_T {
        il::int_t n_dof = 0;
        il::Array2D<double> matrix{};
        il::Array<double> rhside{};
    };

// Element-to-point influence matrix (submatrix of the global one)
    il::StaticArray2D<double, 6, 18> make_local_3dbem_submatrix
            (const int kernel_id,
             double mu, double nu, double h, std::complex<double> z,
             const il::StaticArray<std::complex<double>, 3> &tau,
             const il::StaticArray2D<std::complex<double>, 6, 6> &sfm);

// Static matrix assembly
    il::Array2D<double> make_3dbem_matrix_s
            (double mu, double nu,
             const Mesh_Geom &mesh,
             const Num_Param &n_par,
             il::io_t, DoF_Handle_T &dof_hndl);

// Volume Control matrix assembly (additional row $ column)
    il::Array2D<double> make_3dbem_matrix_vc
            (double mu, double nu,
             const Mesh_Geom &mesh,
             const Num_Param &n_par,
             il::io_t, DoF_Handle_T &dof_hndl);

// Volume Control system modification (for DD increments)
    Alg_Sys_T mod_3dbem_system_vc
            (const il::Array2D<double> &orig_matrix,
             const DoF_Handle_T &orig_dof_hndl,
             const DoF_Handle_T &dof_hndl,
             const il::Array<double> &delta_t,
             const double delta_v);

// Stress at given points (m_pts_crd) vs DD at nodal points (nodes_crd)
    il::Array2D<double> make_3dbem_stress_f_s
            (double mu, double nu,
             const Mesh_Geom &mesh,
             const Num_Param &n_par,
             // const Mesh_Data &m_data,
             const il::Array2D<double> &m_pts_crd);

}

#endif //INC_3D_BEM_MATRIX_ASM_H

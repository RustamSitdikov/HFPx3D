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

#ifndef INC_HFPX3D_MATRIX_ASM_H
#define INC_HFPX3D_MATRIX_ASM_H

#include <complex>
#include <il/Array2D.h>
#include <il/StaticArray.h>
#include <il/StaticArray2D.h>
#include "model_utilities.h"
#include "surface_mesh_utilities.h"

namespace hfp3d {

    // structure for system of algebraic equations
    struct SAE_T {
        il::int_t n_dof = 0;
        il::Array2D<double> matrix{};
        il::Array<double> rhs_v{};
    };

/////// Elastostatics utilities ///////

    // Element-to-point influence matrix (submatrix of the global one)
    il::StaticArray2D<double, 6, 18> make_local_3dbem_submatrix
            (const int kernel_id,
             double mu, double nu, double h, std::complex<double> z,
             const il::StaticArray<std::complex<double>, 3> &tau,
             const il::StaticArray2D<std::complex<double>, 6, 6> &sfm);

    // Static matrix assembly
    il::Array2D<double> make_3dbem_matrix_s
            (double mu, double nu,
             const Mesh_Geom_T &mesh,
             const Num_Param_T &n_par,
             il::io_t, DoF_Handle_T &dof_hndl);

    // Add S_inf to the RHS
    void add_s_inf_to_3dbem_rhs
            (const Mesh_Data_T &mesh_data,
             const Load_T &load,
             il::io_t, SAE_T &sae);

    // Stress at given points (m_pts_crd) vs DD at nodal points (nodes_crd)
    il::Array2D<double> make_3dbem_stress_f_s
            (double mu, double nu,
             const Mesh_Data_T &m_data,
             const Load_T &load,
             const Num_Param_T &n_par,
             const il::Array2D<double> &m_pts_crd);

/////// Volume Control scheme utilities ///////

    // Volume Control matrix assembly (additional row $ column)
    il::Array2D<double> make_3dbem_matrix_vc
            (double mu, double nu,
             const Mesh_Geom_T &mesh,
             const Num_Param_T &n_par,
             il::io_t, DoF_Handle_T &dof_hndl);

    // Volume Control system modification (for DD increments)
    SAE_T mod_3dbem_system_vc
            (const il::Array2D<double> &orig_matrix,
             const DoF_Handle_T &orig_dof_hndl,
             const DoF_Handle_T &dof_hndl,
             const il::Array<double> &delta_t,
             const double delta_v);
}

#endif //INC_HFPX3D_MATRIX_ASM_H

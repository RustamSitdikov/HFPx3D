//
// This file is part of HFPx3D.
//
// Created by D. Nikolski on 4/20/2017.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//

#ifndef INC_HFPX3D_MESH_UTILITIES_H
#define INC_HFPX3D_MESH_UTILITIES_H

#include <cstdio>
#include <complex>
#include <il/Status.h>
#include <il/Array.h>
#include <il/Array2D.h>
#include <il/StaticArray.h>
#include <il/StaticArray2D.h>

namespace hfp3d {

    // mesh geometry structure
    // to be imported from .npy or .csv files
    struct Mesh_Geom_T {
        // nodes' coordinates
        il::Array2D<double> nods;

        // mesh connectivity
        il::Array2D<il::int_t> conn;
    };

    // DoF handle structure
    struct DoF_Handle_T {
        // No of DoF in use
        il::int_t n_dof = 0;

        // element-wise list of DoF
        il::Array2D<il::int_t> dof_h;
        // dof_h.size(0) = number of elements
        // dof_h.size(1) = number of degrees of freedom per element
        // dof_h(j, k) = -1 means fixed degree of freedom

        // mixed boundary conditions parameters
        // il::Array2D<double> bc_c;
        // bc_c(k, 0)*t + bc_c(k, 1)*DD = bc_c(k, 2)
    };

    // fault (crack) state, node-wise
    //todo: think about making these 2D arrays, element/local node-wise
    struct Frac_State_T {
        // "damage state", node-wise
        il::Array<double> mr_open; // max. reached relative opening (w/cr_open)
        il::Array<double> mr_slip; // max. reached relative slip (s/cr_slip)

        // instant friction & cohesion forces, node-wise
        il::Array<double> friction_coef; // slip (shear) friction coefficient
        // (tan of friction angle)
        il::Array<double> slip_cohesion; // slip (shear) cohesion
        il::Array<double> open_cohesion; // opening cohesion

        // dilatancy etc.
    };

    // solution state
    struct Mesh_Data_T {
        // link to the Mesh object
        Mesh_Geom_T mesh;

        // material IDs, for each element / local node
        il::Array2D<int> mat_id;

        // current time
        double time = 0;

        // set of "active" (slid or opened) elements
        il::Array<il::int_t> ae_set;
        // set of "filled" (w. fluid) elements
        il::Array<il::int_t> fe_set;
        // in both cases, -1 means intact element/node

        // list of next-to-tip elements & edges (to propagate from)
        il::Array2D<il::int_t> tip_set;
        // elem No; node a (1..3); node b (1..3); prev. edge No; next edge No

        // element-wise DoF handles for DD and pressure
        DoF_Handle_T dof_h_dd;
        DoF_Handle_T dof_h_pp;

        // displacement discontinuities, for all nodes
        il::Array2D<double> dd;
        // fluid pressure, for all nodes
        il::Array<double> pp;

        // dependent variables: rates, damage %, dilatancy, etc
        Frac_State_T frac_state;

        // convergence (discrepancy)

        // status
        il::Status status{};
    };

/////// some utilities ///////

    // Material ID initialization (trivial case, zeros)
    il::Array2D<int> make_mat_id_triv
            (const Mesh_Geom_T &mesh,
             int ap_order);

    // DoF handle initialization for an isolated crack
    // (fixed DoF at crack tip nodes defined by tip_type)
    DoF_Handle_T make_dof_h_crack
            (const Mesh_Geom_T &mesh,
             int ap_order,
             int tip_type);

    // mesh (solution) data initialization for an undisturbed fault
    Mesh_Data_T init_mesh_data_p_fault
            (const Mesh_Geom_T &mesh,
             int ap_order,
             //int tip_type,
             il::Array2D<il::int_t> inj_loc);

    // 2D (Nx3) to 1D array conversion for DD
    // according to DoF handles
    il::Array<double> get_dd_vector_from_md
            (const Mesh_Data_T &m_data,
             bool include_p);

    // 1D to 2D array conversion for DD
    // according to DoF handles
    void write_dd_vector_to_md
            (const il::Array<double> &sol_v,
             const DoF_Handle_T &dof_h_dd,
             bool include_p,
             const DoF_Handle_T &dof_h_pp,
             il::io_t,
             Mesh_Data_T &m_data);
}

#endif //INC_HFPX3D_MESH_UTILITIES_H

//
// This file is part of HFPx3D_VC.
//
// Created by D. Nikolski on 4/20/2017.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//

#ifndef HFPX3D_VC_MESH_UTILITIES_H
#define HFPX3D_VC_MESH_UTILITIES_H

#include <cstdio>
#include <complex>
#include <il/Array.h>
#include <il/Array2D.h>
#include <il/StaticArray.h>
#include <il/StaticArray2D.h>

namespace hfp3d {

    struct Mesh_Geom {
        il::Array2D<double> nods; // nodes' coordinates
        il::Array2D<il::int_t> conn; // mesh connectivity
    };

    struct Num_Param {
        double beta = 0.125; // relative collocation points' position
        int tip_type = 1; // how to enforce zero DD at the tip:
        // 0 - no enforcement;
        // 1 - only at vertex points;
        // 2 - at vertex and edge nodes
        // bool is_part_uniform = true; // how to partition edges
        bool is_dd_in_glob = true; // defines in what coordinate system
        // DD components are sought: true - global; false - local
    };

    struct DoF_Handle_T {
        il::int_t n_dof = 0;
        // dof_h.size(0) = number of elements
        // dof_h.size(1) = number of degrees of freedom per element
        // dof_h(j, k) = -1 means fixed degree of freedom
        il::Array2D<il::int_t> dof_h{};
        // il::Array2D<double> bc_c;
        // bc_c(k, 0)*t + bc_c(k, 1)*DD = bc_c(k, 2)
    };

    struct Mesh_Data {
        il::Array2D<il::int_t> tip_el; // list of next-to-tip elements
                                       // & edges (to propagate from)

        il::Array<il::int_t> ae_set; // set of "active" (slid or opened) elems

        il::Array2D<double> dd; // displacement discontinuities
        il::Array<double> pp; // fluid pressure
    };

    // DoF handle initialization for a crack (fixed DoF at crack tip nodes)
    DoF_Handle_T make_dof_h_triangular
            (const Mesh_Geom &mesh,
             int ap_order,
             int tip_type);

    // 2D to 1D array conversion for DD
    il::Array<double> get_dd_vector_from_md
            (const Mesh_Data &m_data,
             const DoF_Handle_T &dof_h,
             bool include_p,
             const il::Array<il::int_t> inj_pts);

    // 1D to 2D array conversion for DD
    void write_dd_vector_to_md
            (const il::Array<double> &dd_v,
             const DoF_Handle_T &dof_h,
             bool include_p,
             const il::Array<il::int_t> inj_pts,
             il::io_t,
             Mesh_Data & m_data);

}

#endif //HFPX3D_VC_MESH_UTILITIES_H

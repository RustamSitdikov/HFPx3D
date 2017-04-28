//
// This file is part of HFPx3D_VC.
//
// Created by nikolski on 2/27/2017.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//

#ifndef HFPX3D_VC_CF_ITERATION_H
#define HFPX3D_VC_CF_ITERATION_H

#include <il/Array.h>
#include <il/Array2D.h>
#include <il/StaticArray.h>
#include "system_assembly.h"
#include "cohesion_friction.h"

namespace hfp3d {

    double vc_cf_iteration
            (const Mesh_Geom_T &mesh, // triangulation data
             const Num_Param_T &n_par, // mesh preprocessing params: beta etc
             double mu, double nu, // shear modulus, Poisson ratio
             const il::StaticArray<double, 6> &s_inf, // stress at infinity
             const F_C_Model &cf_m, // friction-cohesion model
             const Alg_Sys_T &orig_vc_sys,
             const DoF_Handle_T &orig_dof_h,
             const Frac_State_T &prev_cp_state, // "damage state" @ prev time step
             double t_vol, // injected volume at the current time step
             il::io_t,
             Mesh_Data_T &m_data, // DD, pressure at nodal points
             DoF_Handle_T &dof_h,
             Frac_State_T &iter_cp_state); // "damage state" @ current time step

}

#endif //HFPX3D_VC_CF_ITERATION_H

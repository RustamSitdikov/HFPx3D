//
// This file is part of HFPx3D.
//
// Created by nikolski on 2/27/2017.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//

#ifndef HFPX3D_CF_ITERATION_H
#define HFPX3D_CF_ITERATION_H

#include <il/Array.h>
#include <il/Array2D.h>
#include <il/StaticArray.h>
#include "system_assembly.h"
#include "cohesion_friction.h"

namespace hfp3d {

    double vc_cf_iteration
            (const Mesh_Data_T &m_data_p, // mesh & solution @ prev time step
             const Properties_T &prop,

             const SAE_T &orig_vc_sys,
             const DoF_Handle_T &orig_dof_h,

             const F_C_Model &f_c_model, // friction-cohesion model
             const Num_Param_T &n_par, // simulation parameters
             const Load_T &load, // stress at infinity
             double d_vol, // injected volume between the time steps

             il::io_t,
             Mesh_Data_T &m_data); // mesh & solution @ current time step

}

#endif //HFPX3D_CF_ITERATION_H

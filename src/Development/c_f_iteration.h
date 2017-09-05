//
// This file is part of HFPx3D.
//
// Created by nikolski on 2/27/2017.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//

#ifndef INC_HFPX3D_CF_ITERATION_H
#define INC_HFPX3D_CF_ITERATION_H

#include <il/Array.h>
#include <il/Array2D.h>
#include <il/StaticArray.h>
#include "src/Core/model_parameters.h"
#include "src/Core/surface_mesh_utilities.h"
#include "src/Solvers/system_assembly.h"
#include "cohesion_friction.h"

namespace hfp3d {

    // one iteration step of the volume control scheme on a pre-existing mesh
    double vc_cf_iteration
            (const Mesh_Data_T &m_data_p, // mesh & solution @ prev time step

             const SAE_T &orig_vc_sys,
             const DoF_Handle_T &orig_dof_h,

             const Properties_T &prop, // material properties incl. contact
             const Num_Param_T &n_par, // simulation parameters

             const Load_T &load, // stress at infinity
             double new_vol, // total injected volume at the current time step

             il::io_t,
             Mesh_Data_T &m_data); // mesh & solution @ current time step

}

#endif //INC_HFPX3D_CF_ITERATION_H

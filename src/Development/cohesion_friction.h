//
// This file is part of HFPx3D.
//
// Created by nikolski on 4/21/2017.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//

#ifndef INC_HFPX3D_COHESION_FRICTION_H
#define INC_HFPX3D_COHESION_FRICTION_H

#include <il/Array.h>
#include <il/Array2D.h>
#include "src/Core/model_parameters.h"
#include "src/Core/surface_mesh_utilities.h"

namespace hfp3d {

    // Calculation of friction & cohesion forces
    // (matching them to current DD and "damage state")
    void match_f_c
            (Properties_T prop,       //
             il::Array<int> mat_id,
             il::Array2D<double> &dd, // current displacements
             bool is_dd_local,        // Note: dd must be in local coordinates!
             il::io_t,
             Frac_State_T &f_state    // "damage state" & friction-cohesion
            );


    // Calculation of max. relative (to critical) displacements
    // (called after finished the new time step)
    void update_mrd
            (Properties_T prop,
             il::Array<int> mat_id,
             il::Array2D<double> &dd, // current displacements
             bool is_dd_local,
             il::io_t,
             Frac_State_T &f_state    // "damage state" & friction-cohesion
            );

    //todo: think about Mesh_Data_T input
}

#endif //INC_HFPX3D_COHESION_FRICTION_H

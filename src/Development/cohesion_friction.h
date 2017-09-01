//
// This file is part of HFPx3D.
//
// Created by nikolski on 4/21/2017.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//

#ifndef INC_HFPX3D_COHESION_FRICTION_H
#define INC_HFPX3D_COHESION_FRICTION_H

#include <il/Array.h>
#include <il/Array2D.h>

namespace hfp3d {

    // fault (crack) state, node-wise
    //todo: (discussed) make it 2D arrays, element/local node-wise
    struct Frac_State_T {
        // "damage state", node-wise
        il::Array<double> mr_open; // max. reached relative opening (w/cr_open)
        il::Array<double> mr_slip; // max. reached relative slip (s/cr_slip)

        // friction & cohesion forces, node-wise
        il::Array<double> friction_coef; // slip (shear) friction coefficient
        // (tan of friction angle)
        il::Array<double> slip_cohesion; // slip (shear) cohesion
        il::Array<double> open_cohesion; // opening cohesion

        // dilatancy etc.
    };

    // friction & cohesion parameters
    //todo: make it possible for f_c_param to be array of size = number of mat_id's
    struct F_C_Param_T {
        double cr_open; // critical opening
        double peak_ts; // peak tensile stress
        double cr_slip; // critical slip
        double peak_sc; // peak shear cohesion
        double peak_sf; // peak shear friction (tangent of friction angle)
        double res_sf;  // residual friction (after critical slip)
    };

    // Calculation of friction & cohesion forces
    // (matching them to current DD and "damage state")
    void match_f_c
            (Properties_T prop,            //
             il::Array2D<int> mat_id,
             il::Array2D<double> &dd, // current displacements
             bool is_dd_local,        // Note: dd must be in local coordinates!
             il::io_t,
             Frac_State_T &f_state    // "damage state" & friction-cohesion
            );


    // Calculation of max. relative (to critical) displacements
    // (called after finished the new time step)
    void update_mrd
            (il::Array2D<double> &dd, // current displacements
             bool is_dd_local,
             il::io_t,
             Frac_State_T &f_state    // "damage state" & friction-cohesion
            );

}

#endif //INC_HFPX3D_COHESION_FRICTION_H

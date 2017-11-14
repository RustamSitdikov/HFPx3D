//
// This file is part of HFPx3D.
//
// Created by nikolski on 7/4/2017.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//

#ifndef INC_HFPX3D_MODEL_PARAMETERS_H
#define INC_HFPX3D_MODEL_PARAMETERS_H

#include <il/Array.h>
#include <il/Array2D.h>
#include <il/StaticArray.h>

namespace hfp3d {

    // contact model (friction / cohesion / dilatancy) parameters
    // to be read from a .toml file
    //todo: think about components as arrays of size = number of mat_id's
    struct F_C_Param_T {
        double cr_open; // critical opening
        double peak_ts; // peak tensile stress
        double cr_slip; // critical slip
        double peak_sc; // peak shear cohesion
        double peak_sf; // peak shear friction (tangent of friction angle)
        double res_sf;  // residual friction (after critical slip)
        // dilatancy etc.
    };

    // physical model parameters
    // to be read from a .toml file
    struct Properties_T {
        // number of materials (solids, liquids, surfaces)
        il::int_t n_solid = 1;
        il::int_t n_liquid = 1;
        il::int_t n_surface = 1;

        // shear moduli, for each solid ID
        il::Array<double> mu {n_solid, 1.0};
        // Poisson ratios, for each solid ID
        il::Array<double> nu {n_solid, 0.0};

//        // solid IDs, for each surface ID (for "plus" and "minus" sides)
//        // note: different moduli and Poisson ratios will affect the assembly!
//        il::Array<int> s_id_p;
//        il::Array<int> s_id_n;

        // contact model IDs, for each surface ID
        il::Array2D<int> c_model_id; // size = n_surface x 2

        // contact model parameters, for each contact model ID
        il::Array<F_C_Param_T> f_c_param;
    };

    // "global" load parameters
    // to be read from a .toml file
    struct Load_T {
        // stress at infinity
        il::StaticArray<double, 6> s_inf;

        // injection locations (elements, nodes)
        il::Array2D<il::int_t> inj_loc;

        // injection rate(s), for each location
        il::Array<double> inj_rate;

//        // liquid ID, for each location
//        il::Array<int> liq_id;
    };

    // numerical scheme parameters
    // to be read from a .toml file
    struct Num_Param_T {
        // order (0, 1, or 2) of approximating (shape) functions for DD
        // int approx_order = 2;

        // relative collocation points' position
        double beta = 0.125;

        // how to enforce zero DD at the tip:
        int tip_type = 1;
        // 0 -> no enforcement;
        // 1 -> only at vertex points;
        // 2 -> at vertex and edge nodes

        // in what coordinate system DD are sought
        bool is_dd_local = false;
        // true -> local (for each element); false -> global (reference)

        // how to partition edges
        // bool is_part_uniform = true;
    };

    // simulation parameters
    // to be read from a .toml file
    struct Sim_Param_T {
        bool do_save_matrix = false;
        bool do_save_solution = false;
        bool do_postprocess = false;
    };

}
#endif //INC_HFPX3D_MODEL_PARAMETERS_H

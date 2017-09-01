//
// This file is part of HFPx3D.
//
// Created by nikolski on 7/4/2017.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//

#ifndef HFPX3D_MODEL_UTILS_H
#define HFPX3D_MODEL_UTILS_H

#include "src/Development/cohesion_friction.h"

namespace hfp3d {

    // physical model parameters
    // to be read from a .toml file
    struct Properties_T {
        // number of materials (solids, liquids, surfaces)
        il::int_t n_solid = 1;
        il::int_t n_liquid = 1;
        il::int_t n_surface = 1;

        // shear moduli, for each solid ID
        il::Array<double> mu;
        // Poisson ratios, for each solid ID
        il::Array<double> nu;

//        // shear moduli, for each surface ID (for "plus" and "minus" sides)
//        // different moduli and Poisson ratios affect the assembly
//        il::Array<double> mu_p;
//        il::Array<double> mu_m;
//        // Poisson ratios, for each surface ID (for "plus" and "minus" sides)
//        il::Array<double> nu_p;
//        il::Array<double> nu_m;

        // contact model ID, for each surface ID (see cohesion_friction.h)
        il::Array<int> c_model_id;

        // contact model parameters, for each contact model ID
        // (see cohesion_friction.h)
        il::Array<F_C_Param_T> f_c_param;

//        // max friction coeff-s
//        il::Array<double> fr_c;
//        // max cohesion (for shear stress vs shear DD)
//        il::Array<double> sc_f;
//        // max cohesive forces (for opening mode)
//        il::Array<double> oc_f;
//        // critical opening displacement (normal DD)
//        il::Array<double> w_cr;
//        // critical slip displacement (in-plane DD)
//        il::Array<double> s_cr;
    };

    // "global" load parameters
    // to be read from a .toml file
    struct Load_T {
        // stress at infinity
        il::StaticArray<double, 6> s_inf;

        // injection locations (elements, nodes)
        il::Array2D<il::int_t> inj_loc;

        // injection rate(s)
        il::Array<double> inj_rate;
    };

    // numerical simulation parameters
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

}
#endif //HFPX3D_MODEL_UTILS_H

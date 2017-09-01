//
// This file is part of HFPx3D.
//
// Created by nikolski on 9/1/2017.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//

#include <il/math.h>
#include <il/Array.h>
#include <il/Array2D.h>
#include "src/Development/cohesion_friction.h"
#include "src/Core/element_utilities.h"

// Particular friction-cohesion models *************************************

// "Box" or linear function for cohesion; linear slip-weakening for friction
void match_f_c
            (Properties_T prop,
             il::Array2D<int> mat_id,
             il::Array2D<double> &dd,
             bool is_dd_local,
             il::io_t,
             Frac_State_T &f_state) {
    IL_EXPECT_FAST(dd.size(1) == 3);
    il::int_t n_nod = dd.size(0);
    //todo: (discussed) make it 2D arrays, element/local node-wise
    IL_EXPECT_FAST(n_nod == f_state.mr_open.size());
    IL_EXPECT_FAST(n_nod == f_state.mr_slip.size());
    IL_EXPECT_FAST(n_nod == f_state.friction_coef.size());
    IL_EXPECT_FAST(n_nod == f_state.slip_cohesion.size());
    IL_EXPECT_FAST(n_nod == f_state.open_cohesion.size());
    il::StaticArray<double, 3> nodal_dd{};
    //todo: loop over elements / local nodes
    for (il::int_t n = 0; n < n_nod; ++n) {
        // contact model for local mat_ID (surface type)
        int model_id = prop.c_model_id[mat_id[n]];
        // local friction-cohesion parameters
        F_C_Param_T f_c_p = prop.f_c_param[mat_id[n]];
        // local shear / opening DD
        for (int i = 0; i < 3; ++i) {
            nodal_dd[i] = dd(n, i);
        }
        if (!is_dd_local) {
            // rotate the DD vector to local coordinates
            // use /Core/element_utilities
        }
        // relative opening DD
        double dl_open = nodal_dd[2] / f_c_p.cr_open;
        // make sure dd(2, n) is positive
        // prev. opening (damage %)
        double pr_open = f_state.mr_open[n];
        // calculating correction for partial damage (unloading case)
        double r_s = 0.0;
        if (dl_open > 0.0 && dl_open < 1.0) {
            switch (model_id) {
                case 1: // linear weakening
                    r_s = 1.0 - dl_open;
                    break;
                case 0: // "box" (Dugdale)
                    r_s = 1.0;
                    break;
                default:break;
            }
        }
        if (dl_open > 0.0 && dl_open < pr_open) {
            r_s *= dl_open / pr_open; // linear unloading
        }

        // calculating forces
        // opening cohesion
        f_state.open_cohesion[n] = f_c_p.peak_ts * r_s;
        // shear cohesion
        f_state.slip_cohesion[n] = f_c_p.peak_sc * r_s;

        // calculating friction
            double r_f = 0.0; // zero if fully opened (normal DD > cr_open)
            if (dl_open < 1.0) {
                r_f = f_c_p.res_sf;
                // current sliding DD
                double dl_slip = nodal_dd[0] * nodal_dd[0] + nodal_dd[1] * nodal_dd[1];
                // the same relative to cr_slip
                dl_slip = std::sqrt(dl_slip) / f_c_p.cr_slip;
                // prev. relative sliding DD (damage %)
                double pr_slip = f_state.mr_slip[n];
                // comparing and calculating friction
                if (dl_slip < 1.0) {
                    double m_dl_slip = dl_slip;
                    if (dl_slip < pr_slip) {
                        m_dl_slip = pr_slip;
                    }
                    r_f = res_sf() +
                          (f_c_p.peak_sf - f_c_p.res_sf)*(1.0 - m_dl_slip);
                }
            }
            f_state.friction_coef[n] = r_f; // friction coefficient(s)
        }
    }

// Calculation of max. relative (to critical) displacements
// (called after finished the new time step)
void update_mrd
        (il::Array2D<double> &dd, // current displacements
         bool is_dd_local,
         il::io_t,
         Frac_State_T &f_state // "damage state" & friction-cohesion
        ) {
    IL_EXPECT_FAST(dd.size(1) == 3);
    il::int_t n_nod = dd.size(0);
    IL_EXPECT_FAST(n_nod == f_state.mr_open.size());
    IL_EXPECT_FAST(n_nod == f_state.mr_slip.size());
    IL_EXPECT_FAST(n_nod == f_state.friction_coef.size());
    IL_EXPECT_FAST(n_nod == f_state.slip_cohesion.size());
    IL_EXPECT_FAST(n_nod == f_state.open_cohesion.size());
    il::StaticArray<double, 3> nodal_dd{};
    for (il::int_t n = 0; n < n_nod; ++n) {
        for (int i = 0; i < 3; ++i) {
            nodal_dd[i] = dd(n, i);
        }
        if (!is_dd_local) {}
        // current opening DD relative to cr_open
        double dl_open = nodal_dd[2] / cr_open(); // relative opening DD
        // make sure dd(2, n) is positive!
        // current sliding DD
        double dl_slip = nodal_dd[0] * nodal_dd[0] + nodal_dd[1] * nodal_dd[1];
        // current sliding DD relative to cr_slip
        dl_slip = std::sqrt(dl_slip) / cr_slip();
        double pr_open = f_state.mr_open[n]; // prev. opening (damage %)
        double pr_slip = f_state.mr_slip[n]; // prev. slip (damage %)
        // updating the "damage %" (state)
        if (dl_slip > pr_slip) {
            f_state.mr_slip[n] = dl_slip;
        }
        if (dl_slip >= 1.0) {
            // f_state.mr_slip[n] = 1.0;
            // f_state.mr_open[n] = 1.0;
        }
        if (dl_open >= 1.0) {
            f_state.mr_open[n] = 1.0;
            // f_state.mr_slip[n] = 1.0;
        } else {
            if (dl_open > pr_open) {
                f_state.mr_open[n] = dl_open;
            }
        }
    }
}

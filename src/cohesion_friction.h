//
// This file is part of HFPx3D.
//
// Created by nikolski on 4/21/2017.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//

#ifndef INC_HFPX3D_COHESION_FRICTION_H
#define INC_HFPX3D_COHESION_FRICTION_H

#include <il/math.h>
#include <il/Array.h>
#include <il/Array2D.h>

namespace hfp3d {

    // fault (crack) state, node-wise
    struct Frac_State_T {
        // "damage state", node-wise
        il::Array<double> mr_open; // max. reached relative opening (w/cr_open)
        il::Array<double> mr_slip; // max. reached relative slip (s/cr_slip)

        // friction & cohesion forces, node-wise
        il::Array<double> friction_coef; // slip (shear) friction coefficient
        // (tan of friction angle)
        il::Array<double> slip_cohesion; // slip (shear) cohesion
        il::Array<double> open_cohesion; // opening cohesion
    };

    // friction & cohesion parameters
    struct F_C_Param_T {
        double cr_open; // critical opening
        double peak_ts; // peak tensile stress
        double cr_slip; // critical slip
        double peak_sc; // peak shear cohesion
        double peak_sf; // peak shear friction (tangent of friction angle)
        double res_sf; // residual friction (after critical slip)
    };

    // General form of friction-cohesion model
    class F_C_Model {
    private:
        F_C_Param_T f_c_param_;

    public:
        F_C_Model(F_C_Param_T f_c_param) {
            f_c_param_.cr_open = f_c_param.cr_open;
            f_c_param_.peak_ts = f_c_param.peak_ts;
            f_c_param_.cr_slip = f_c_param.cr_slip;
            f_c_param_.peak_sc = f_c_param.peak_sc;
            f_c_param_.peak_sf = f_c_param.peak_sf;
            f_c_param_.res_sf = f_c_param.res_sf;
        };

        F_C_Param_T f_c_param() { return f_c_param_; };
        double cr_open() { return f_c_param_.cr_open; };
        double cr_slip() { return f_c_param_.cr_slip; };
        double peak_ts() { return f_c_param_.peak_ts; };
        double peak_sc() { return f_c_param_.peak_sc; };
        double peak_sf() { return f_c_param_.peak_sf; };
        double res_sf() { return f_c_param_.res_sf; };

        // Calculation of friction & cohesion forces
        // (matching them to current DD and "damage state")
        virtual void match_f_c
                (il::Array2D<double> &dd, // current displacements
                 bool is_dd_local,
                 il::io_t,
                 Frac_State_T &f_state) // "damage state" & friction-cohesion
        = 0; // purely virtual in general
        // Note: dd must be in local coordinates!

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

        };
    };

    // Particular friction-cohesion models *************************************

    // "Box" function for cohesion; linear slip-weakening for friction
    class F_C_BFW: F_C_Model {

        F_C_BFW(F_C_Param_T f_c_param) :
                F_C_Model(f_c_param){}

        void match_f_c
                (il::Array2D<double> &dd,
                 bool is_dd_local,
                 il::io_t,
                 Frac_State_T &f_state) {
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
                double dl_open = nodal_dd[2] / cr_open(); // relative 
                // opening DD
                // make sure dd(2, n) is positive
                double pr_open = f_state.mr_open[n]; // prev. opening
                // (damage %)
                // calculating correction for partial damage (unloading case)
                double r_s = 0.0;
                if (dl_open > 0.0 && dl_open < 1.0) {
                    r_s = 1.0;
                }
                if (dl_open > 0.0 && dl_open < pr_open) {
                    r_s *= dl_open / pr_open;
                }
                // calculating forces
                f_state.open_cohesion[n] = peak_ts() * r_s; // opening cohesion
                f_state.slip_cohesion[n] = peak_sc() * r_s; // shear cohesion
                // calculating friction
                double r_f = 0.0; // zero if fully opened (normal DD > cr_open)
                if (dl_open < 1.0) {
                    r_f = res_sf();
                    // current sliding DD
                    double dl_slip = nodal_dd[0] * nodal_dd[0] + nodal_dd[1] * nodal_dd[1];
                    // the same relative to cr_slip
                    dl_slip = std::sqrt(dl_slip) / cr_slip();
                    // prev. relative sliding DD (damage %)
                    double pr_slip = f_state.mr_slip[n];
                    // comparing and calculating friction
                    if (dl_slip < 1.0) {
                        double m_dl_slip = dl_slip;
                        if (dl_slip < pr_slip) {
                            m_dl_slip = pr_slip;
                        }
                        r_f = res_sf() +
                              (peak_sf() - res_sf())*(1.0 - m_dl_slip);
                    }
                }
                f_state.friction_coef[n] = r_f; // friction coefficient(s)
            }
        }

    };

}

#endif //INC_HFPX3D_COHESION_FRICTION_H

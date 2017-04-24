//
// This file is part of HFPx3D_VC.
//
// Created by nikolski on 4/21/2017.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//

#ifndef HFPX3D_VC_COHESION_FRICTION_H
#define HFPX3D_VC_COHESION_FRICTION_H

#include <il/math.h>
#include <il/Array.h>
#include <il/Array2D.h>

namespace hfp3d {

    // fault (crack) "damage state", node-wise
    struct Frac_State {
        il::Array<double> mr_open; // max. reached relative opening (w/cr_open)
        il::Array<double> mr_slip; // max. reached relative slip (s/cr_slip)
    };

    // fault (crack) friction & cohesion forces, node-wise
    struct Frac_F_C {
        il::Array<double> slip_friction; // slip (shear) friction
        // (tan of friction angle)
        il::Array<double> slip_cohesion; // slip (shear) cohesion
        il::Array<double> open_cohesion; // opening cohesion
    };

    // friction & cohesion parameters
    struct F_C_Param {
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
        F_C_Param f_c_param_;

    public:
        F_C_Model(F_C_Param f_c_param) {
            f_c_param_.cr_open = f_c_param.cr_open;
            f_c_param_.peak_ts = f_c_param.peak_ts;
            f_c_param_.cr_slip = f_c_param.cr_slip;
            f_c_param_.peak_sc = f_c_param.peak_sc;
            f_c_param_.peak_sf = f_c_param.peak_sf;
            f_c_param_.res_sf = f_c_param.res_sf;
        };

        F_C_Param f_c_param() { return f_c_param_; };
        double cr_open() { return f_c_param_.cr_open; };
        double cr_slip() { return f_c_param_.cr_slip; };
        double peak_ts() { return f_c_param_.peak_ts; };
        double peak_sc() { return f_c_param_.peak_sc; };
        double peak_sf() { return f_c_param_.peak_sf; };
        double res_sf() { return f_c_param_.res_sf; };

        // Calculation of friction & cohesion forces
        virtual Frac_F_C c_f_forces // (node-wise)
                (const il::Array2D<double> &dd, // current displacements
                        // (must be in local coords!)
                 const Frac_State &p_state) // previous state (damage)
        = 0; //
    };

    // Particular friction-cohesion models

    // "Box" function for cohesion; linear slip-weakening for friction
    class F_C_BFW: F_C_Model {

        F_C_BFW(F_C_Param f_c_param) :
                F_C_Model(f_c_param){}

        Frac_F_C c_f_forces
                (const il::Array2D<double> &dd,
                 const Frac_State &p_state) {
            IL_EXPECT_FAST(dd.size(0) == 3);
            il::int_t n_nod = dd.size(1);
            IL_EXPECT_FAST(n_nod == p_state.mr_open.size());
            IL_EXPECT_FAST(n_nod == p_state.mr_slip.size());
            Frac_F_C f_c;
            f_c.slip_friction.resize(n_nod);
            f_c.slip_cohesion.resize(n_nod);
            f_c.open_cohesion.resize(n_nod);
            for (il::int_t n = 0; n < n_nod; ++n) {
                double dl_open = dd(2, n) / cr_open();
                double pr_open = p_state.mr_open[n];
                double r_s = 0.0;
                if (dl_open > 0.0 && dl_open < 1.0) {
                    r_s = 1.0;
                }
                if (dl_open > 0.0 && dl_open < pr_open) {
                    r_s *= dl_open / pr_open;
                }
                f_c.open_cohesion[n] = peak_ts() * r_s; // opening cohesion
                f_c.slip_cohesion[n] = peak_sc() * r_s; // shear cohesion
                double r_f = 0.0;
                if (dl_open < 1.0) {
                    r_f = res_sf();
                    double dl_slip = dd(0, n) * dd(0, n) + dd(1, n) * dd(1, n);
                    dl_slip = std::sqrt(dl_slip) / cr_slip();
                    double pr_slip = p_state.mr_slip[n];
                    if (dl_slip < 1.0) {
                        double m_dl_slip = dl_slip;
                        if (dl_slip < pr_slip) {
                            m_dl_slip = pr_slip;
                        }
                        r_f = res_sf() +
                              (peak_sf() - res_sf())*(1.0 - m_dl_slip);
                    }
                }
                f_c.slip_friction[n] = r_f; // friction
            }
            return f_c;
        }

    };

}

#endif //HFPX3D_VC_COHESION_FRICTION_H

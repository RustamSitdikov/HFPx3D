//
// This file is part of HFPx3D.
//
// Created by nikolski on 9/1/2017.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//

//#include <il/math.h>
#include <il/Array.h>
#include <il/Array2D.h>
#include "src/Core/model_parameters.h"
#include "src/Core/element_utilities.h"
#include "src/Core/surface_mesh_utilities.h"

namespace hfp3d {

// Particular friction-cohesion models *************************************

    // "Box" or linear function for cohesion; linear slip-weakening for friction
    void match_f_c
            (Properties_T prop,
             il::Array<int> mat_id,
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
#pragma omp parallel for
        for (il::int_t n = 0; n < n_nod; ++n) {
            // contact model for local mat_ID (surface type)
            int c_m_id = prop.c_model_id(mat_id[n],0);
            int f_m_id = prop.c_model_id(mat_id[n],1);
            // local friction-cohesion parameters
            F_C_Param_T f_c_p = prop.f_c_param[mat_id[n]];
            // local shear / opening DD
            for (int i = 0; i < 3; ++i) {
                nodal_dd[i] = dd(n, i);
            }
            if (!is_dd_local) {
                // todo: add rotation of the DD vector to local coordinates
                // use /Core/element_utilities
//                // Vertices' coordinates
//                il::StaticArray2D<double, 3, 3> el_vert;
//                for (il::int_t j = 0; j < 3; ++j) {
//                    il::int_t n = m_data.mesh.conn(j, el);
//                    for (il::int_t k = 0; k < 3; ++k) {
//                        el_vert(k, j) = m_data.mesh.nods(k, n);
//                    }
//                }
//                Element_Struct_T ele_s = set_ele_struct(el_vert, n_par.beta);
//                //il::blas(1.0, ele_s.r_tensor, nodal_dd, 0.0, il::io, nodal_dd);
//                nodal_dd = il::dot(ele_s.r_tensor, nodal_dd);
            }
            // relative opening DD
            double dl_open = nodal_dd[2] / f_c_p.cr_open;
            // make sure dd(2, n) is positive
            // prev. opening (damage %)
            double pr_open = f_state.mr_open[n];
            // current sliding DD
            double dl_slip = nodal_dd[0] * nodal_dd[0] +
                             nodal_dd[1] * nodal_dd[1];
            // the same relative to cr_slip
            dl_slip = std::sqrt(dl_slip) / f_c_p.cr_slip;
            // prev. relative sliding DD (damage %)
            double pr_slip = f_state.mr_slip[n];
            // calculating correction for partial damage
            double r_s = 0.0;
            if (dl_open > 0.0 && dl_open < 1.0) {
                switch (c_m_id) {
                    case 2: // ???
                        break;
                    case 1: // linear weakening
                        r_s = 1.0 - dl_open;
                        break;
                    case 0: // "box" (Dugdale)
                        r_s = 1.0;
                        break;
                    default:break;
                }
            }
            // calculating correction for partial damage (unloading case)
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
                    // comparing and calculating friction
                    switch (f_m_id) {
                        case 2: // exponential weakening
                            break;
                        case 1: // linear weakening
                            if (dl_slip < 1.0) {
                                double m_dl_slip = dl_slip;
                                if (dl_slip < pr_slip) {
                                    m_dl_slip = pr_slip;
                                }
                                r_f = f_c_p.res_sf +
                                      (f_c_p.peak_sf - f_c_p.res_sf) *
                                              (1.0 - m_dl_slip);
                            }
                            break;
                        case 0: // no weakening
                            r_f = f_c_p.peak_sf;
                            break;
                        default:break;
                    }
                }
                f_state.friction_coef[n] = r_f; // friction coefficient(s)
            }
        }

    // Calculation of max. relative (to critical) displacements
    // (called after finished the new time step)
    void update_mrd
        (Properties_T prop,
         il::Array<int> mat_id,
         il::Array2D<double> &dd, // current displacements
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
#pragma omp parallel for
        for (il::int_t n = 0; n < n_nod; ++n) {
            for (int i = 0; i < 3; ++i) {
                nodal_dd[i] = dd(n, i);
            }
            if (!is_dd_local) {
                // todo: add rotation of the DD vector to local coordinates
                // use /Core/element_utilities
//                //il::blas(1.0, ele_s.r_tensor, nodal_dd, 0.0, il::io, nodal_dd);
//                nodal_dd = il::dot(ele_s.r_tensor, nodal_dd);
            }
            // local friction-cohesion parameters
            F_C_Param_T f_c_p = prop.f_c_param[mat_id[n]];
            // current opening DD relative to cr_open
            double dl_open = nodal_dd[2] / f_c_p.cr_open;
            // make sure dd(2, n) is positive!
            // current sliding DD
            double dl_slip = nodal_dd[0] * nodal_dd[0] +
                    nodal_dd[1] * nodal_dd[1];
            // current sliding DD relative to cr_slip
            dl_slip = std::sqrt(dl_slip) / f_c_p.cr_slip;
            // prev. opening (damage %)
            double pr_open = f_state.mr_open[n];
            // prev. slip (damage %)
            double pr_slip = f_state.mr_slip[n];
            // updating the "damage %" (state)
            //todo: think about consistency of "damage" in opening & slip
            if (dl_slip >= 1.0) {
                // f_state.mr_slip[n] = 1.0;
                // f_state.mr_open[n] = 1.0;
            }
            if (dl_slip > pr_slip) {
                f_state.mr_slip[n] = dl_slip;
            }
            if (dl_open >= 1.0) {
                f_state.mr_open[n] = 1.0;
                // f_state.mr_slip[n] = 1.0;
            }
            if (dl_open > pr_open) {
                f_state.mr_open[n] = dl_open;
            }
        }
    }

}

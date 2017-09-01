//
// This file is part of HFPx3D.
//
// Created by nikolski on 2/27/2017.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//

#include <complex>
#include <cmath>
#include <il/Array.h>
#include <il/Array2D.h>
#include <il/StaticArray.h>
#include <il/linear_algebra.h>
#include <il/linear_algebra/dense/factorization/LU.h>
#include "src/Solvers/system_assembly.h"
#include "src/Core/element_utilities.h"
#include "src/Core/tensor_utilities.h"
#include "cohesion_friction.h"
#include "src/Core/surface_mesh_utilities.h"
#include "c_f_iteration.h"

namespace hfp3d {

    double vc_cf_iteration
    // This function performs one iteration step
    // of the volume control scheme on a pre-existing mesh
    // with ability to add new elements (this part is under development)
            (const Mesh_Data_T &m_data_p, // mesh & solution @ prev time step
             const Properties_T &prop,

             const SAE_T &orig_vc_sys, // for pre-meshed surface case
             const DoF_Handle_T &orig_dof_h,

             F_C_Model &f_c_model,
             const Num_Param_T &n_par, // simulation parameters
             const Load_T &load, // stress at infinity
             double new_vol, // total injected volume @ current time step

             il::io_t,
             Mesh_Data_T &m_data // mesh & solution @ current time step
            ) {

        //const il::int_t num_of_ele = orig_dof_h.dof_h.size(0);
        const il::int_t num_of_ele = m_data.mesh.conn.size(1);
        const il::int_t num_of_nod = m_data.mesh.nods.size(1);
        //const il::int_t ndpe = orig_dof_h.dof_h.size(1);
        const il::int_t ndpe = m_data.dof_h_dd.dof_h.size(1);
        const il::int_t nnpe = ndpe / 3;
        IL_EXPECT_FAST(num_of_ele > 0);
        //IL_EXPECT_FAST(num_of_ele == m_data.mesh.conn.size(1));
        const il::int_t full_ndof = num_of_ele * ndpe;
        const il::int_t orig_ndof = orig_dof_h.n_dof;
        const il::int_t prev_ndof = m_data_p.dof_h_dd.n_dof;
        const il::int_t curr_ndof = m_data.dof_h_dd.n_dof;
        IL_EXPECT_FAST(orig_ndof > 0 && orig_ndof <= full_ndof);
        IL_EXPECT_FAST(orig_ndof == orig_vc_sys.matrix.size(0));
        IL_EXPECT_FAST(orig_vc_sys.matrix.size(0) == orig_vc_sys.matrix.size(1));
        IL_EXPECT_FAST(m_data.mat_id.size(0) == num_of_ele);
        IL_EXPECT_FAST(m_data.mat_id.size(1) == nnpe);

        il::Array<double> dd_v = get_dd_vector_from_md(m_data, true);

        // make sure elastic_traction_a = il::dot(orig_vc_sys.matrix, dd_a);
        il::Array<double> e_t_v (dd_v.size());
        e_t_v = il::dot(orig_vc_sys.matrix, dd_v);

        // current volume
        il::Array<double> v_int{orig_ndof, 0.0};
        for (il::int_t i = 0; i < orig_ndof; ++i) {
            v_int[i] = orig_vc_sys.matrix(orig_ndof, i);
        }
        double curr_vol = il::dot(v_int, dd_v);

        double delta_v = new_vol - curr_vol;

        // injection point(s)
        il::int_t inj_pt_0 = load.inj_loc(0, 0) * 6 + load.inj_loc(0, 1);
        // current pressure
        double pressure = m_data.pp[inj_pt_0];

        // "fracture state", by node (.dd should be in local coordinates)
        // includes friction, shear cohesion & opening cohesion
        f_c_model.match_f_c
                (m_data.dd, n_par.is_dd_local, il::io, m_data.frac_state);

        il::Array<double> delta_trac{orig_ndof, 0.0};

        il::int_t dd_dof_dec = 0, pp_dof_inc = 0;

        for (il::int_t el = 0; el < num_of_ele;  ++el) {
            // vertices' coordinates
            il::StaticArray2D<double, 3, 3> el_vert;
            for (il::int_t j = 0; j < 3; ++j) {
                il::int_t n = m_data.mesh.conn(j, el);
                for (il::int_t k = 0; k < 3; ++k) {
                    el_vert(k, j) = m_data.mesh.nods(k, n);
                }
            }

            Element_Struct_T ele_s = set_ele_struct(el_vert, n_par.beta);

            // normal vector
            il::StaticArray<double, 3> nv_el;
            for (int j = 0; j < 3; ++j) {
                nv_el[j] = -ele_s.r_tensor(2, j);
            }

            // traction (negative) induced by stress at infinity
            il::StaticArray<double, 3> ti_el = nv_dot_sim(nv_el, load.s_inf);
            
            // converting traction to local coordinate system
            //il::blas(1.0, ele_s.r_tensor, ti_el, 0.0, il::io, ti_el);
            ti_el = il::dot(ele_s.r_tensor, ti_el);

            // nodal DD & pressure
            il::StaticArray2D<double, 3, 6> dd_el{0.0};
            il::StaticArray<double, 6> pr_el{0.0};
            for (il::int_t lnn = 0; lnn < 6; ++lnn) {
                il::int_t gnn = el * 6 + lnn;
                pr_el[lnn] = m_data.pp[gnn];
                for (int i = 0; i < 3; ++i) {
                    il::int_t el_dof = lnn * 3 + i;
                    il::int_t dof = orig_dof_h.dof_h(el, el_dof);
                    if (dof != -1) {
                        dd_el(i, lnn) = m_data.dd(gnn, i);
                    }
                }
            }

            for (int lnn = 0; lnn < nnpe; ++lnn) {
                // "global" node (or CP) number
                il::int_t gnn = el * nnpe + lnn;

                // mat ID of the node (CP)
                int mat_id_cp = m_data.mat_id(el, lnn);

                // status of CP (partial/total breakup)
                //double prev_cp_st = m_data_p.frac_state.mr_open[gnn];
                double cp_fr_state = m_data.frac_state.mr_open[gnn];
                double cp_sl_state = m_data.frac_state.mr_slip[gnn];
                // double iter_cp_st = il::max (cp_fr_state, cp_sl_state);
                // double iter_cp_st = cp_fr_state * cp_fr_state +
                // cp_sl_state * cp_sl_state;
                // iter_cp_st = sqrt(iter_cp_st);

                // DD at CP
                il::StaticArray<double, 3> dd_cp =
                        il::dot(dd_el, ele_s.sf_cp[lnn]);
                
                // traction due to DD at CP
                il::StaticArray<double, 3> tr_cp;
                for (int i = 0; i < 3; ++i) {
                    il::int_t el_dof = lnn * 3 + i;
                    il::int_t dof = orig_dof_h.dof_h(el, el_dof);
                    if (dof != -1) {
                        tr_cp[i] = e_t_v[dof];
                    }
                }
                
                // converting DD & traction to local coordinate system
                if (!n_par.is_dd_local) {
                    //il::blas(1.0, ele_s.r_tensor, dd_cp, 0.0, il::io, dd_cp);
                    dd_cp = il::dot(ele_s.r_tensor, dd_cp);
                }
                //il::blas(1.0, ele_s.r_tensor, tr_cp, 0.0, il::io, tr_cp);
                tr_cp = il::dot(ele_s.r_tensor, tr_cp);

                // pressure at CP
                double pr_cp = il::dot(pr_el, ele_s.sf_cp[lnn]);
                // il::StaticArray<double, 3> pf_cp {0.0};
                // for (int i = 0; i < 3; ++i) {
                //     pf_cp[i] = nv_el[i] * pr_cp;
                // }

                // total (effective) normal traction at CP
                double nt_cp = pr_cp - ti_el[2] - tr_cp[2];

                // total shear traction at CP
                il::StaticArray<double, 2> st_cp {0.0};
                for (int i = 0; i < 2; ++i) {
                    st_cp[i] += - ti_el[i] - tr_cp[i];
                }
                double sh_cp = st_cp[0] * st_cp[0] + st_cp[1] * st_cp[1];
                sh_cp = std::sqrt(sh_cp);

                // admissible normal & shear traction at CP
                double adm_nt = m_data.frac_state.open_cohesion[gnn];
                double adm_st = m_data.frac_state.slip_cohesion[gnn] -
                        m_data.frac_state.friction_coef[gnn] * nt_cp;
                //nt_cp -= adm_nt;
                // "-" comes from "positive tension" convention

                // shear direction
                il::StaticArray<double, 3> sh_dir {0.0};
                if (sh_cp > 0) {
                    for (int i = 0; i < 2; ++i) {
                        sh_dir[i] = st_cp[i] / sh_cp;
                    }
                }

                // adjustment of traction at CP for the next iteration
                il::StaticArray<double, 3> dt_cp {0.0};

                // traction admissibility check
                // and calculation of traction adjustments
                // normal traction admissibility check
                if (nt_cp > adm_nt) {
                    // full separation
                    // note: if previously fully open, adm_nt==0
                    dt_cp[2] = - nt_cp; // release all normal traction
                    //if (cp_fr_state >= 1.0) { // && nt_cp > 0.0
                        for (int i = 0; i < 2; ++i) {
                            dt_cp[i] = - st_cp[i]; // release all shear
                        }
                    //} else {}
                    for (int i = 0; i < 3; ++i) {
                        il::int_t el_dof = lnn * 3 + i;
                        il::int_t dd_dof = orig_dof_h.dof_h(el, el_dof);
                        m_data.dof_h_dd.dof_h(el, el_dof) = dd_dof - dd_dof_dec;
                    }
                    // il::int_t pp_dof = ;
                    // m_data.dof_h_pp.dof_h(el, lnn) = pp_dof + pp_dof_inc;
                } else { // 0 < prev_cp_st < 1; partial separation or slip
                    if (cp_fr_state < 1.0 && nt_cp > 0.0) {
                        // && nt_cp <= adm_nt
                        dt_cp[2] = adm_nt - nt_cp; // applying cohesion
                        il::int_t el_dof = lnn * 3 + 2;
                        il::int_t dof = orig_dof_h.dof_h(el, el_dof);
                        m_data.dof_h_dd.dof_h(el, el_dof) = dof - dd_dof_dec;
                    }
                    // il::int_t pp_dof = ;
                    // m_data.dof_h_pp.dof_h(el, lnn) = pp_dof + pp_dof_inc;

                    // shear traction admissibility check
                    if (cp_sl_state > 0 || sh_cp > adm_st) {
                        //&& nt_cp <= adm_nt
                        for (int i = 0; i < 2; ++i) {
                            dt_cp[i] = adm_st * sh_dir[i] - st_cp[i];
                            il::int_t el_dof = lnn * 3 + i;
                            il::int_t dof = orig_dof_h.dof_h(el, el_dof);
                            m_data.dof_h_dd.dof_h(el, el_dof) =
                                    dof - dd_dof_dec;
                        }
                    } else if (nt_cp <= 0.0) { // intact CP; no slip or opening
                        for (int i = 0; i < 3; ++i) {
                            il::int_t el_dof = lnn * 3 + i;
                            il::int_t dof = orig_dof_h.dof_h(el, el_dof);
                            m_data.dof_h_dd.dof_h(el, el_dof) = -1;
                            ++dd_dof_dec;
                        }
                    }
                }

                // rotation of traction adj. to the reference coordinate system
                dt_cp = il::dot(ele_s.r_tensor, il::Blas::kTranspose, dt_cp);

                // adding traction adjustments to the right hand side
                for (int i = 0; i < 2; ++i) {
                    il::int_t el_dof = lnn * 3 + i;
                    il::int_t dof = orig_dof_h.dof_h(el, el_dof);
                    if (dof != -1) {
                        delta_trac[dof] = dt_cp[i];
                    }
                }
            }
        }

        m_data.dof_h_dd.n_dof = orig_dof_h.n_dof - dd_dof_dec;
        // m_data.dof_h_pp.n_dof = orig_dof_h.n_dof + pp_dof_inc;

        // system of linear algebraic equations
        SAE_T trc_vc_sys;
        // check if the used matrix is smaller that the original matrix
        if(m_data.dof_h_dd.n_dof > 0 && m_data.dof_h_dd.n_dof <= orig_ndof) {
            // truncation of the algebraic system to only "active" nodes
            // + additional row for VC (volume vs DD & pressure)
            trc_vc_sys = mod_3dbem_system_vc
                    (orig_vc_sys.matrix, orig_dof_h,
                     m_data.dof_h_dd, delta_trac, delta_v);
        } else {
            // growing mesh
        }

        // solving the system
        il::Status status{};
        il::LU<il::Array2D<double>> lu_dc
                (trc_vc_sys.matrix, il::io, status);
        status.abortOnError();
        // double cnd = lu_dc.condition_number(il::Norm::L2, );
        // std::cout << cnd << std::endl;
        il::Array<double> dd_pp_incr_v = lu_dc.solve(trc_vc_sys.rhs_v);
        // dd_pp_incr_v = il::linear_solve
        // (trc_vc_sys.matrix, trc_vc_sys.rhside, il::io, status);
        
        // adjusting pressure
        double delta_p = dd_pp_incr_v[m_data.dof_h_dd.n_dof];
        pressure += delta_p;

        // DD, pressure, nodes' "state" (damage %) adjustment
        for (il::int_t el = 0; el < num_of_ele;  ++el) {
            // Vertices' coordinates
            il::StaticArray2D<double, 3, 3> el_vert;
            for (il::int_t j = 0; j < 3; ++j) {
                il::int_t n = m_data.mesh.conn(j, el);
                for (il::int_t k = 0; k < 3; ++k) {
                    el_vert(k, j) = m_data.mesh.nods(k, n);
                }
            }

            Element_Struct_T ele_s = set_ele_struct(el_vert, n_par.beta);

            // nodal DD over the element in local coordinates (initialization)
            il::StaticArray2D<double, 3, 6> dd_el{0.0};

            for (il::int_t lnn = 0; lnn < nnpe;  ++lnn) {
                // "global" node number
                il::int_t gnn = el * nnpe + lnn;

                // nodal DD (initialization)
                il::StaticArray<double, 3> dd_n;

                // adding calculated increments to the nodal DD
                for (int i = 0; i < 3; ++i) {
                    dd_n[i] = m_data.dd(gnn, i);
                    il::int_t el_dof = lnn * 3 + i;
                    il::int_t orig_dof = orig_dof_h.dof_h(el, el_dof);
                    il::int_t dof = m_data.dof_h_dd.dof_h(el, el_dof);
                    if (orig_dof != -1 && dof != -1) {
                        dd_n[i] += dd_pp_incr_v[dof];
                    }
                }

                // converting the nodal DD to local coordinate system
                if (!n_par.is_dd_local) {
                    //il::blas(1.0, ele_s.r_tensor, dd_n, 0.0, il::io, dd_n);
                    dd_n = il::dot(ele_s.r_tensor, dd_n);
                }
                
                // overlap (negative opening) check
                if (dd_n[2] < 0.0) {
                    dd_n[2] = 0.0;
                }

                // updating the element DD (in local coordinates)
                for (int i = 0; i < 3; ++i) {
                    dd_el(i, lnn) = dd_n[i];
                }

                // converting the nodal DD back to reference coordinate system
                if (!n_par.is_dd_local) {
                    dd_n = il::dot(ele_s.r_tensor, il::Blas::kTranspose, dd_n);
                }
                
                // updating the "global" DD array
                for (int i = 0; i < 3; ++i) {
                    il::int_t el_dof = lnn * 3 + i;
                    il::int_t orig_dof = orig_dof_h.dof_h(el, el_dof);
                    if (orig_dof != -1) {
                        m_data.dd(gnn, i) = dd_n[i];
                     }
                }

                // node state check (comment if checked only at CP)
                //Frac_State prev_n_st = prev_cp_state(el, lnn);
                //double op_st = dd_n[2] / cf_par.cr_open;
                //if (op_st > 1.0 || prev_n_st.mr_open >= 1.0) {
                //    op_st = 1.0;
                //}
                //iter_cp_state(el, lnn).mr_open = op_st;
                //double sl_st = dd_n[0] * dd_n[0] + dd_n[1] * dd_n[1];
                //sl_st = std::sqrt(sl_st) / cf_par.cr_slip;
                //if (sl_st > 1.0 || prev_n_st.mr_slip >= 1.0) {
                //    sl_st = 1.0;
                //}
                //iter_cp_state(el, lnn).mr_slip = sl_st;
            }

            for (il::int_t cpe = 0; cpe < nnpe;  ++cpe) {
                // "global" CP number
                il::int_t gnn = el * nnpe + cpe;

                // DD at CP (in local coordinates)
                il::StaticArray<double, 3> dd_cp =
                        il::dot(dd_el, ele_s.sf_cp[cpe]);

                // converting the nodal DD to local coordinate system
                if (!n_par.is_dd_local) {
                    //il::blas(1.0, ele_s.r_tensor, dd_cp, 0.0, il::io, dd_cp);
                    dd_cp = il::dot(ele_s.r_tensor, dd_cp);
                }

                // status of CP (partial/total breakup)
                // double prev_cp_st = m_data_p.frac_state.mr_open[gnn];
                double cp_fr_state = m_data.frac_state.mr_open[gnn];
                // double cp_sl_state = m_data.frac_state.mr_slip[gnn];
                // double iter_cp_st = il::max (cp_fr_state, cp_sl_state);
                // double iter_cp_st = cp_fr_state * cp_fr_state +
                // cp_sl_state * cp_sl_state;
                // iter_cp_st = sqrt(iter_cp_st);

                double cr_op_cp = f_c_model.cr_open();
                double cp_op_st = dd_cp[2] / cr_op_cp;
                if (cp_op_st >= 1.0 || cp_fr_state >= 1.0) {
                    cp_op_st = 1.0;
                }

                // adding calculated increment of pressure at opened CP
                if (cp_op_st >= 1.0 ||
                        (cp_op_st > 0.0 && cp_fr_state >= 1.0)) {
                    m_data.pp[gnn] = pressure;
                }
            }
        }

        // output (norm of delta_dd + delta_p; norm of delta_trac)
        double res = il::norm(dd_pp_incr_v, il::Norm::L1) +
                il::norm(delta_trac, il::Norm::L1);
        return res;
    }

}
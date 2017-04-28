//
// This file is part of HFPx3D_VC.
//
// Created by nikolski on 2/27/2017.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//

#include <complex>
#include <cmath>
#include <il/Array.h>
#include <il/Array2D.h>
#include <il/StaticArray.h>
#include <il/linear_algebra.h>
#include <il/linear_algebra/dense/factorization/LU.h>
#include "system_assembly.h"
#include "element_utilities.h"
#include "tensor_utilities.h"
#include "cohesion_friction.h"
//#include "mesh_utilities.h"
//#include "c_f_iteration.h"

namespace hfp3d {

    double vc_cf_iteration
    // This function performs one iteration step
    // of the volume control scheme on a pre-existing mesh
    // with ability to add new elements (this part is under development)
            (const Mesh_Geom_T &mesh, // triangulation data
             const Num_Param_T &n_par, // mesh preprocessing params: beta etc
             double mu, double nu, // shear modulus, Poisson ratio
             const il::StaticArray<double, 6> &s_inf, // stress at infinity
             const F_C_Model &cf_m, // friction-cohesion model
             const Alg_Sys_T &orig_vc_sys,
             const DoF_Handle_T &orig_dof_h,
             const Frac_State &prev_cp_state, // "damage state" @ prev time step
             double t_vol, // injected volume at the current time step
             il::io_t,
             Mesh_Data_T &m_data, // DD, pressure at nodal points
             DoF_Handle_T &dof_h,
             Frac_State &iter_cp_state // "damage state" @ current time step
            ) {

        IL_EXPECT_FAST(orig_vc_sys.matrix.size(0) == orig_vc_sys.matrix.size(1));
        const il::int_t num_of_ele = orig_dof_h.dof_h.size(0);
        const il::int_t ndpe = orig_dof_h.dof_h.size(1);
        const il::int_t nnpe = ndpe / 3;
        IL_EXPECT_FAST(num_of_ele > 0);
        IL_EXPECT_FAST(num_of_ele == mesh.conn.size(1));
        const il::int_t full_ndof = num_of_ele * ndpe;
        const il::int_t orig_ndof = orig_dof_h.n_dof;
        IL_EXPECT_FAST(orig_ndof > 0 && orig_ndof <= full_ndof);
        IL_EXPECT_FAST(orig_ndof == orig_vc_sys.matrix.size(0));

        // make sure elastic_traction_a = il::dot(orig_vc_sys.matrix, dd_a);
        elastic_traction_a = il::dot(orig_vc_sys.matrix, dd_a);

        il::Array<double> dd_a = convert_DD_array_to_vector(m_data);
        
        // current volume
        il::Array<double> v_int{orig_ndof, 0.0};
        for (il::int_t i = 0; i < orig_ndof; ++i) {
            v_int[i] = orig_vc_sys.matrix(orig_ndof, i);
        }
        double c_vol = il::dot(v_int, dd_a);
        double delta_v = t_vol - c_vol;
        double pressure = il::max(m_data.pp);

        il::Array<double> delta_t{orig_ndof, 0.0};

        il::int_t dof_dec = 0;

        for (il::int_t el = 0; el < num_of_ele;  ++el) {
            // vertices' coordinates
            il::StaticArray2D<double, 3, 3> el_vert;
            for (il::int_t j = 0; j < 3; ++j) {
                il::int_t n = mesh_conn(j, el);
                for (il::int_t k = 0; k < 3; ++k) {
                    el_vert(k, j) = nodes_crd(k, n);
                }
            }

            Ele_Struct ele_s = set_ele_struct(el_vert, beta);

            // normal vector
            il::StaticArray<double, 3> nv_el;
            for (int j = 0; j < 3; ++j) {
                nv_el[j] = -ele_s.r_tensor(2, j);
            }

            // traction (negative) induced by stress at infinity
            il::StaticArray<double, 3> ti_el = nv_dot_sim(nv_el, s_inf);
            
            // converting traction to local coordinate system
            //il::blas(1.0, ele_s.r_tensor, ti_el, 0.0, il::io, ti_el);
            ti_el = il::dot(ele_s.r_tensor, ti_el);

            // nodal DD & pressure
            il::StaticArray2D<double, 3, 6> dd_el{0.0};
            il::StaticArray<double, 6> pr_el{0.0};
            for (il::int_t lnn = 0; lnn < 6; ++lnn) {
                pr_el[lnn] = m_data.pp(el, lnn);
                for (int i = 0; i < 3; ++i) {
                    il::int_t el_dof = lnn * 3 + i;
                    il::int_t dof = orig_dof_hndl.dof_h(el, el_dof);
                    if (dof != -1) {
                        dd_el(i, lnn) = dd_a[dof];
                    }
                }
            }

            for (il::int_t lnn = 0; lnn < nnpe; ++lnn) {
                il::int_t n = el * nnpe + lnn;
                
                // status of CP (partial/total breakup)
                Frac_State prev_cp_st = prev_cp_state(el, lnn);
                Frac_State iter_cp_st = iter_cp_state(el, lnn);
                
                // DD at CP
                il::StaticArray<double, 3> dd_cp =
                        il::dot(dd_el, ele_s.sf_cp[lnn]);
                
                // traction due to DD at CP
                il::StaticArray<double, 3> tr_cp;
                for (int i = 0; i < 3; ++i) {
                    il::int_t el_dof = lnn * 3 + i;
                    il::int_t dof = orig_dof_hndl.dof_h(el, el_dof);
                    if (dof != -1) {
                        tr_cp[i] = elastic_traction_a[dof];
                    }
                }
                
                // converting DD & traction to local coordinate system
                if (is_dd_in_glob) {
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

                // friction, shear cohesion & opening cohesion at CP
                Frac_F_C f_c_cp = cfm.c_f(dd_cp, prev_cp_st);

                // total normal traction at CP
                double nt_cp = pr_cp - ti_el[2] - tr_cp[2];

                // total shear traction at CP
                il::StaticArray<double, 2> st_cp {0.0};
                for (int i = 0; i < 2; ++i) {
                    st_cp[i] += - ti_el[i] - tr_cp[i];
                }
                double sh_cp = st_cp[0] * st_cp[0] + st_cp[1] * st_cp[1];
                sh_cp = std::sqrt(sh_cp);

                // admissible normal & shear traction at CP
                double adm_nt = f_c_cp.open_cohesion;
                double adm_st =  f_c_cp.slip_cohesion - 
                        f_c_cp.slip_friction * nt_cp;
                //nt_cp -= adm_nt;

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
                    dt_cp[2] = - nt_cp; // release all normal traction
                    //if (prev_cp_st.mr_open >= 1.0) { // && nt_cp > 0.0
                        for (int i = 0; i < 2; ++i) {
                            dt_cp[i] = - st_cp[i]; // release all shear
                        }
                    //} else {}
                    for (int i = 0; i < 3; ++i) {
                        il::int_t el_dof = lnn * 3 + i;
                        il::int_t dof = orig_dof_hndl.dof_h(el, el_dof);
                        dof_hndl.dof_h(el, el_dof) = dof - dof_dec;
                    }
                } else { // 0 < prev_cp_st < 1; partial separation or slip
                    if (prev_cp_st.mr_open < 1.0 && nt_cp > 0.0) {
                        // && nt_cp <= adm_nt
                        dt_cp[2] = adm_nt - nt_cp; // applying cohesion
                        il::int_t el_dof = lnn * 3 + 2;
                        il::int_t dof = orig_dof_hndl.dof_h(el, el_dof);
                        dof_hndl.dof_h(el, el_dof) = dof - dof_dec;
                    }
                    // shear traction admissibility check
                    if (iter_cp_st.mr_slip > 0 || sh_cp > adm_st) {
                        //&& nt_cp <= adm_nt
                        for (int i = 0; i < 2; ++i) {
                            dt_cp[i] = adm_st * sh_dir[i] - st_cp[i];
                            il::int_t el_dof = lnn * 3 + i;
                            il::int_t dof = orig_dof_hndl.dof_h(el, el_dof);
                            dof_hndl.dof_h(el, el_dof) = dof - dof_dec;
                        }
                    } else if (nt_cp <= 0.0) { // intact CP; no slip or opening
                        for (int i = 0; i < 3; ++i) {
                            il::int_t el_dof = lnn * 3 + i;
                            il::int_t dof = orig_dof_hndl.dof_h(el, el_dof);
                            dof_hndl.dof_h(el, el_dof) = -1;
                            ++dof_dec;
                        }
                    }
                }

                // rotation of traction adj. to the reference coordinate system
                dt_cp = il::dot(ele_s.r_tensor, il::Blas::transpose, dt_cp);

                // adding traction adjustments to the right hand side
                for (int i = 0; i < 2; ++i) {
                    il::int_t el_dof = lnn * 3 + i;
                    il::int_t dof = orig_dof_hndl.dof_h(el, el_dof);
                    if (dof != -1) {
                        delta_t[dof] = dt_cp[i];
                    }
                }
            }
        }

        dof_hndl.n_dof = orig_dof_hndl.n_dof - dof_dec;

        Alg_Sys_T trc_vc_sys;
        const il::int_t used_ndof = dof_hndl.n_dof;
        // check if the used matrix is smaller that the original matrix
        if(used_ndof > 0 && used_ndof <= orig_ndof) {
            // truncation of the algebraic system to only "active" nodes
            Alg_Sys_T trc_vc_sys = mod_3dbem_system_vc
                    (orig_vc_sys.matrix, orig_dof_hndl,
                     dof_hndl, delta_t, delta_v);
        } else {
            // growing mesh
        }

        il::Status status{};
        il::LU<il::Array2D<double>> lu_dc
                (trc_vc_sys.matrix, il::io, status);
        status.abort_on_error();
        // double cnd = lu_dc.condition_number(il::Norm::L2, );
        // std::cout << cnd << std::endl;
        il::Array<double> trc_dd_v = lu_dc.solve(trc_vc_sys.rhside);
        // trc_dd_v = il::linear_solve
        // (trc_vc_sys.matrix, trc_vc_sys.rhside, il::io, status);
        
        double delta_p = trc_dd_v[used_ndof];
        pressure += delta_p;

        for (il::int_t el = 0; el < num_of_ele;  ++el) {
            // Vertices' coordinates
            il::StaticArray2D<double, 3, 3> el_vert;
            for (il::int_t j = 0; j < 3; ++j) {
                il::int_t n = mesh_conn(j, el);
                for (il::int_t k = 0; k < 3; ++k) {
                    el_vert(k, j) = nodes_crd(k, n);
                }
            }

            Ele_Struct ele_s = set_ele_struct(el_vert, beta);

            // nodal DD over the element in local coordinates (initialization)
            il::StaticArray2D<double, 3, 6> dd_el{0.0};

            for (il::int_t lnn = 0; lnn < nnpe;  ++lnn) {
                // nodal DD (initialization)
                il::StaticArray<double, 3> dd_n;

                // adding calculated increments to the nodal DD
                for (int i = 0; i < 3; ++i) {
                    il::int_t el_dof = lnn * 3 + i;
                    il::int_t orig_dof = orig_dof_hndl.dof_h(el, el_dof);
                    il::int_t dof = dof_hndl.dof_h(el, el_dof);
                    if (orig_dof != -1 && dof != -1) {
                        dd_n[i] = dd_a[orig_dof] + trc_dd_v[dof];
                    }
                }

                // converting the nodal DD to local coordinate system
                if (is_dd_in_glob) {
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
                if (is_dd_in_glob) {
                    dd_n = il::dot(ele_s.r_tensor, il::Blas::transpose, dd_n);
                }
                
                // updating the "global" DD array
                for (int i = 0; i < 3; ++i) {
                    il::int_t el_dof = lnn * 3 + i;
                    il::int_t orig_dof = orig_dof_hndl.dof_h(el, el_dof);
                    if (orig_dof != -1) {
                        dd_a[orig_dof] = dd_n[i];
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
                // DD at CP (in local coordinates)
                il::StaticArray<double, 3> dd_cp =
                        il::dot(dd_el, ele_s.sf_cp[cpe]);

                // CP state check
                Frac_State prev_n_st = prev_cp_state(el, cpe);
                double cropen = cfm.cr_open();
                double cp_op_st = dd_cp[2] / cropen;
                if (cp_op_st >= 1.0 || prev_n_st.mr_open >= 1.0) {
                    cp_op_st = 1.0;
                }
                iter_cp_state(el, cpe).mr_open = cp_op_st;
                double crslip = cfm.cr_slip();
                double cp_sl_st = dd_cp[0] * dd_cp[0] + dd_cp[1] * dd_cp[1];
                cp_sl_st = std::sqrt(cp_sl_st) / crslip;
                if (cp_sl_st > 1.0 || prev_n_st.mr_slip >= 1.0) {
                    cp_sl_st = 1.0;
                }
                iter_cp_state(el, cpe).mr_slip = cp_sl_st;

                // adding calculated increment of pressure at opened CP
                if (cp_op_st >= 1.0 ||
                        (cp_op_st > 0.0 && prev_n_st.mr_open >= 1.0)) {
                    m_data.pp(el, cpe) = pressure;
                }
            }
        }

        // adding calculated increments to the elastic traction
        elastic_traction_a = il::dot(orig_vc_sys.matrix, dd_a);

        // output (norm of delta_dd + delta_p; norm of delta_t)
        double res = il::norm(trc_dd_v, il::Norm::L1) +
                il::norm(delta_t, il::Norm::L1);
        return res;
    }

}
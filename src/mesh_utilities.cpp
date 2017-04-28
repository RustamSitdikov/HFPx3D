//
// This file is_n_a part of HFPx3D_VC.
//
// Created by D. Nikolski on 4/20/2017.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//

#include <il/Array2D.h>
#include <il/StaticArray.h>
#include "mesh_utilities.h"

namespace hfp3d {

    // DoF handle initialization for a crack (fixed DoF at crack tip nodes)
    DoF_Handle_T make_dof_h_crack
            (const Mesh_Geom_T &mesh,
             int ap_order,
             int tip_type) {
        // This function fills up the DoF handle matrix
        // assuming the surface mesh given by mesh connectivity (mesh.conn) 
        // and nodes' coordinates (mesh.nods) is an isolated crack
        // i.e. the edge mated with only one element belongs to the tip

        il::int_t n_ele = mesh.conn.size(1);
        bool is_ext_msh = ap_order > 0 &&
                          mesh.nods.size(0) >= 5 && mesh.conn.size(0) >= 5;
        IL_EXPECT_FAST(ap_order >= 0);
        // number of nodes per element (triangular)
        int nnpe = (ap_order + 1) * (ap_order + 2) / 2;
        int edge_nn = 0;
        if (ap_order > 1) {
            edge_nn = ap_order - 1;
        }
        int ndpe = 3 * nnpe;
        DoF_Handle_T d_h;
        d_h.n_dof = n_ele * ndpe;
        d_h.dof_h = il::Array2D<il::int_t> {n_ele, ndpe, 0};
        for (il::int_t el = 0; el < n_ele; ++el) {
            for (int v = 0; v < nnpe; ++v) {
                for (int l = 0; l < 3; ++l) {
                    int ldof = v * 3 + l;
                    il::int_t gdof = el * ndpe + ldof;
                    d_h.dof_h(el, ldof) = gdof;
                }
            }
        }
        if (nnpe >= 3) {
            il::int_t dof_dec = 0;
            for (il::int_t el = 0; el < n_ele; ++el) {
                // status of vertices 1...3 (at the tip or not)
                il::StaticArray<bool, 3> n_st{false};
                if (is_ext_msh && mesh.conn(3, el) == mesh.conn(4, el)) {
                    for (int v = 0; v < 3; ++v) {
                        il::int_t n = mesh.conn(v, el);
                        if (mesh.nods(3, n) == 1 && mesh.nods(4, n) == 1) {
                            n_st[v] = true;
                        }
                    }
                }
                for (int v = 0; v < nnpe; ++v) {
                    // check if the node (vertex) v (v < 3) is_n_a at the tip
                    if (tip_type >= 1 && v < 3 && n_st[v]) {
                        for (int l = 0; l < 3; ++l) {
                            int ldof = v * 3 + l;
                            --d_h.n_dof;
                            ++dof_dec;
                            d_h.dof_h(el, ldof) = -1;
                        }
                    } else if (ap_order > 1 && tip_type == 2 && v >= 3) {
                        // the vertex across the v-th node
                        int w = (v - 3) / edge_nn;
                        if (w < 3) {
                            // the edge containing the v-th node
                            int a = (w + 1) % 3;
                            int b = (a + 1) % 3;
                            // check if the edge ab is_n_a at the tip
                            if (n_st[a] && n_st[b]) {
                                for (int l = 0; l < 3; ++l) {
                                    int ldof = v * 3 + l;
                                    --d_h.n_dof;
                                    ++dof_dec;
                                    d_h.dof_h(el, ldof) = -1;
                                }
                            }
                        }
                    } else {
                        for (int l = 0; l < 3; ++l) {
                            int ldof = v * 3 + l;
                            d_h.dof_h(el, ldof) -= dof_dec;
                        }
                    }
                }
            }
        }
        return d_h;
    }

    // mesh (solution) data initialization for an undisturbed fault
    Mesh_Data_T init_mesh_data_p_fault
            (const Mesh_Geom_T &i_mesh,
             int ap_order,
             //int tip_type,
             il::Array2D<il::int_t> inj_loc) {
        // This function makes the DoF handles for DD and pressure
        // assuming the surface mesh given by mesh connectivity (mesh.conn)
        // and nodes' coordinates (mesh.nods) is a fault w. injection points
        // defined by inj_pts listing elements and local nodes (1..6) where
        // the pressure is applied

        IL_EXPECT_FAST(ap_order >= 0);
        // number of elements
        il::int_t n_el = i_mesh.conn.size(1);
        // number of nodes per element (triangular)
        int nnpe = (ap_order + 1) * (ap_order + 2) / 2;
        // number of DD DoF per element
        int ndpe = 3 * nnpe;
        // number of injection-affected elements
        il::int_t n_ie = inj_loc.size(0);
        IL_EXPECT_FAST(inj_loc.size(1) >= nnpe + 1);

        Mesh_Data_T m_data;
        // pass input mesh (i_mesh) handle
        m_data.mesh = i_mesh;

        // DD DoF initialization
        m_data.dof_h_dd.dof_h = il::Array2D<il::int_t>{n_el, ndpe, -1};
        // pressure DoF initialization (for discontinuous Galerkin scheme)
        m_data.dof_h_pp.dof_h = il::Array2D<il::int_t>{n_el, nnpe, -1};
        // pressure DoF initialization (for continuous scheme, e.g. FV)
        //m_data.dof_h_dd = il::Array2D<il::int_t>{?, ?, -1};

        // DD and pressure arrays initialization
        m_data.dd = il::Array2D<double> {n_el * nnpe, 3, 0.0};
        // for discontinuous Galerkin scheme, element-wise
        m_data.pp = il::Array<double> {n_el * nnpe, 0.0};
        // for continuous scheme, e.g. FV
        //m_data.pp = il::Array<double> {?, 0.0};

        // fluid and failure activated elements initialization
        m_data.fe_set = il::Array<il::int_t> {n_ie};
        m_data.ae_set = il::Array<il::int_t> {n_ie};

        // setting fe, ae, and DoF handles
        // limited by the injection "spot"
        il::int_t g_dd_dof = 0, g_pp_dof = 0;
        // loop over elements w. injection points
        for (il::int_t i_el = 0; i_el < n_ie; ++i_el) {
            // element No
            il::int_t el = inj_loc(i_el, 0);
            // adding element to "filled" and "active" sets
            m_data.fe_set[i_el] = el;
            m_data.ae_set[i_el] = el;
            // loop over nodal points of the element
            for (int n = 0; n < nnpe; ++n) {
                il::int_t i_n = inj_loc(i_el, n + 1);
                // look if pressure is applied at the node
                if (i_n != -1) {
                    ++g_pp_dof;
                    m_data.dof_h_pp.dof_h(el, n) = g_pp_dof;
                    for (int l = 0; l < 3; ++l) {
                        int ldof = n * 3 + l;
                        ++g_dd_dof;
                        m_data.dof_h_dd.dof_h(el, ldof) = g_dd_dof;
                    }
                }
            }
        }
        m_data.dof_h_dd.n_dof = g_dd_dof;
        m_data.dof_h_pp.n_dof = g_pp_dof;

        return m_data;
    }

    // 2D (n x 3) to 1D (RHS) array conversion for DD
    // according to DoF handles
    il::Array<double> get_dd_vector_from_md
            (const Mesh_Data_T &m_data,
             const DoF_Handle_T &dof_h_dd,
             bool include_p,
             const DoF_Handle_T &dof_h_pp) {
        // number of elements
        il::int_t n_el = dof_h_dd.dof_h.size(0);

        // number of DoF for DD
        il::int_t n_dd_dof = dof_h_dd.n_dof;

        // number of DoF for pressure
        il::int_t n_pp_dof = (include_p ? dof_h_pp.n_dof: 0);

        // make sure that dimensions match
        IL_EXPECT_FAST(dof_h_dd.dof_h.size(0) == n_el);
        IL_EXPECT_FAST(m_data.dd.size(0) == n_el * 6);
        IL_EXPECT_FAST(m_data.dd.size(1) == 3);
        if (include_p) {
            IL_EXPECT_FAST(dof_h_pp.dof_h.size(0) == n_el);
            IL_EXPECT_FAST(m_data.pp.size() == n_el);
        }

        // initialization of RHS vector (output)
        il::Array<double> rhs_v {n_dd_dof + n_pp_dof, 0.};

        // copying DD values
        // loop over elements
        for (il::int_t el = 0; el < n_el; ++el) {
            // loop over local nodes (1 .. 6)
            for (int en = 0; en < 6; ++en) {
                // node No
                il::int_t n = 6 * el + en;
                // loop over nodal DoF (1 .. 3)
                for (int i = 0; i < 3; ++i) {
                    // local DoF (1 .. 18)
                    il::int_t j = 3 * en + i;
                    // copying the DD value if listed in dof_h_dd
                    if (dof_h_dd.dof_h(el, j) != -1) {
                        rhs_v[dof_h_dd.dof_h(el, j)] = m_data.dd(n, i);
                    }
                }
            }
        }

        // copying pressure values if requested
        if (include_p) {
            // loop over elements
            for (il::int_t el = 0; el < n_el; ++el) {
                // loop over local nodes (1 .. 6)
                for (int en = 0; en < 6; ++en) {
                    // node No
                    il::int_t n = 6 * el + en;
                    // copy pressure value if listed in dof_h_pp
                    if (dof_h_pp.dof_h(el, n) != -1) {
                        rhs_v[dof_h_pp.dof_h(el, n)] = m_data.pp[n];
                    }
                }
            }
        }
        return rhs_v;
    }

    // 1D (RHS) to 2D (n x 3) array conversion for DD
    // according to DoF handles
    void write_dd_vector_to_md
            (const il::Array<double> &rhs_v,
             const DoF_Handle_T &dof_h_dd,
             bool include_p,
             const DoF_Handle_T &dof_h_pp,
             il::io_t,
             Mesh_Data_T &m_data) {
        // number of elements
        il::int_t n_el = dof_h_dd.dof_h.size(0);

        // number of DoF for DD
        il::int_t n_dd_dof = dof_h_dd.n_dof;

        // number of DoF for pressure
        il::int_t n_pp_dof = (include_p ? dof_h_pp.n_dof: 0);

        // make sure that dimensions match
        IL_EXPECT_FAST(dof_h_dd.dof_h.size(0) == n_el);
        IL_EXPECT_FAST(m_data.dd.size(0) == n_el * 6);
        IL_EXPECT_FAST(m_data.dd.size(1) == 3);
        if (include_p) {
            IL_EXPECT_FAST(dof_h_pp.dof_h.size(0) == n_el);
            IL_EXPECT_FAST(m_data.pp.size() == n_el);
        }
        IL_EXPECT_FAST(rhs_v.size() == n_dd_dof + n_pp_dof);

        // copying values
        // loop over elements
        for (il::int_t el = 0; el < n_el; ++el) {
            // loop over local nodes (1 .. 6)
            for (int en = 0; en < 6; ++en) {
                // node No
                il::int_t n = 6 * el + en;
                // loop over nodal DoF (1 .. 3)
                for (int i = 0; i < 3; ++i) {
                    // local DoF (1 .. 18)
                    il::int_t j = 3 * en + i;
                    // copying the DD value if listed in dof_h_dd
                    if (dof_h_dd.dof_h(el, j) != -1) {
                        m_data.dd(n, i) = rhs_v[dof_h_dd.dof_h(el, j)];
                    }
                }
                // copying the pressure value if requested
                if (include_p) {
                    if (dof_h_pp.dof_h(el, en) != -1) {
                        m_data.pp[n] = rhs_v[n_dd_dof + dof_h_pp.dof_h(el, en)];
                    }
                }
            }
        }
    }

}

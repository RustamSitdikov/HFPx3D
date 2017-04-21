//
// This file is part of HFPx3D_VC.
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
    DoF_Handle_T make_dof_h_triangular
            (const Mesh_Geom &mesh,
             int ap_order,
             int tip_type) {
        // This function fills up the DoF handle matrix
        // assuming the surface mesh given by mesh connectivity (mesh.conn) 
        // and nodes' coordinates (mesh.nods) is a crack
        // i.e. the edge mated with only one element belongs to the tip

        il::int_t n_ele = mesh.conn.size(1);
        bool is_ext_msh = ap_order > 0 &&
                          mesh.nods.size(0) >= 5 && mesh.conn.size(0) >= 5;
        IL_EXPECT_FAST(ap_order >= 0);
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
                    // check if the node (vertex) v (v < 3) is at the tip
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
                            // check if the edge ab is at the tip
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
    
}

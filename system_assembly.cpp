//
// This file is part of HFPx3D.
//
// Created by D. Nikolski on 1/24/2017.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//

#include <iostream>
#include <complex>
#include <il/math.h>
#include <il/Array2D.h>
#include <il/StaticArray.h>
#include <il/StaticArray2D.h>
#include <il/StaticArray3D.h>
#include <il/StaticArray4D.h>
#include <il/linear_algebra.h>
// #include <il/linear_algebra/dense/blas/dot.h>
// #include <il/linear_algebra/dense/blas/blas.h>
#include "system_assembly.h"
#include "tensor_utilities.h"
#include "element_utilities.h"
#include "elasticity_kernel_integration.h"

namespace hfp3d {

    Dof_Handle make_dof_handle
            (const il::Array2D<il::int_t> &mesh_conn,
             const il::Array2D<double> &nodes_crd,
             int tip_type) {
        il::int_t n_ele = mesh_conn.size(1);
        bool is_ext_msh = nodes_crd.size(0) >=5 && mesh_conn.size(0) >=5;
        Dof_Handle d_h;
        d_h.n_nod = n_ele * 6;
        d_h.n_dof = n_ele * 18;
        d_h.dof_h = il::Array2D<il::int_t>{d_h.n_nod, 3, 0};
        //d_h.h_dof = il::Array2D<il::int_t>{d_h.n_dof, 2, 0};
        //d_h.bc_c = il::Array2D<double>{d_h.n_dof, 3, 0.0};
        il::int_t dof = 0;
        for (il::int_t i = 0; i < n_ele; ++i) {
            il::StaticArray<bool, 3> n_st{false};
            if (is_ext_msh && mesh_conn(3, i) == mesh_conn(4, i)) {
                for (int v = 0; v < 3; ++v) {
                    il::int_t n = mesh_conn(v, i);
                    if (nodes_crd(3, n) == 1 && nodes_crd(4, n) == 1) {
                        n_st[v] = true;
                    }
                }
            }
            for (int v = 0; v < 3; ++v) {
                il::int_t n = i * 6 + v;
                // check if the vertex v is at the tip
                if (n_st[v] && (tip_type == 1 || tip_type == 2)) {
                    for (int l = 0; l < 3; ++l) {
                        --d_h.n_dof;
                        d_h.dof_h(n, l) = -1;
                    }
                } else {
                    for (int l = 0; l < 3; ++l) {
                        d_h.dof_h(n, l) = dof;
                        //d_h.h_dof(dof, 0) = n;
                        //d_h.h_dof(dof, 1) = l;
                        //d_h.bc_c(dof, 0) = 0.0;
                        //d_h.bc_c(dof, 1) = 1.0;
                        //d_h.bc_c(dof, 2) = 0.0;
                        ++dof;
                    }
                }
                // the edge across the v-th vertex
                int a = (v + 1) % 3;
                int b = (a + 1) % 3;
                // check if the edge ab is at the tip
                if (n_st[a] && n_st[b] && tip_type == 2) {
                    for (int l = 0; l < 3; ++l) {
                        --d_h.n_dof;
                        d_h.dof_h(n + 3, l) = -1;
                    }
                } else {
                    for (int l = 0; l < 3; ++l) {
                        d_h.dof_h(n + 3, l) = dof;
                        //d_h.h_dof(dof, 0) = n + 3;
                        //d_h.h_dof(dof, 1) = l;
                        //d_h.bc_c(dof, 0) = 0.0;
                        //d_h.bc_c(dof, 1) = 1.0;
                        //d_h.bc_c(dof, 2) = 0.0;
                        ++dof;
                    }
                }
            }
        }
        //d_h.h_dof.resize(d_h.n_dof, 2);
        return d_h;
    }

// "Global" matrix assembly

    il::Array2D<double> make_3dbem_matrix_s
            (double mu, double nu, double beta,
             const il::Array2D<il::int_t> &mesh_conn,
             const il::Array2D<double> &nodes_crd,
             int tip_type,
             il::io_t, Dof_Handle &dof_hndl) {
// BEM matrix assembly from boundary mesh geometry data:
// mesh connectivity (mesh_conn) and nodes' coordinates (nodes_crd)

// Naive way: no parallelization, no ACA

        IL_EXPECT_FAST(mesh_conn.size(0) >= 3);
        IL_EXPECT_FAST(mesh_conn.size(1) >= 1); // at least 1 element
        IL_EXPECT_FAST(nodes_crd.size(0) >= 3);
        IL_EXPECT_FAST(nodes_crd.size(1) >= 3); // at least 3 nodes

        dof_hndl = make_dof_handle(mesh_conn, nodes_crd, tip_type);

        const il::int_t num_of_elems = mesh_conn.size(1);
        const il::int_t num_of_dof = dof_hndl.n_dof;
        //const il::int_t num_of_dof = 18 * num_of_elems;

        //IL_EXPECT_FAST(global_matrix.size(0) == 18*num_of_elems);
        //IL_EXPECT_FAST(global_matrix.size(1) == 18*num_of_elems);

        il::Array2D<double> global_matrix(num_of_dof, num_of_dof);
        //il::StaticArray2D<double, num_of_dof, num_of_dof> global_matrix;
        //il::StaticArray<double, num_of_dof> right_hand_side;

        // Loop over "source" elements
#pragma omp parallel for
        for (il::int_t source_elem = 0; 
             source_elem < num_of_elems; ++source_elem) {
            // Vertices' coordinates
            il::StaticArray2D<double, 3, 3> el_vert_s;
            //il::StaticArray<double, 3> vert_wts_s;
            for (il::int_t j = 0; j < 3; ++j) {
                il::int_t n = mesh_conn(j, source_elem);
                for (il::int_t k = 0; k < 3; ++k) {
                    el_vert_s(k, j) = nodes_crd(k, n);
                }
                // set vert_wts_s[j]
            }

            // Basis (shape) functions and rotation tensor of the el-t
            il::StaticArray2D<double, 3, 3> r_tensor_s;
            il::StaticArray2D<std::complex<double>, 6, 6> sfm = 
                    make_el_sfm_uniform(el_vert_s, il::io, r_tensor_s);
            //il::StaticArray2D<std::complex<double>, 6, 6> sfm = 
            // make_el_sfm_nonuniform(r_tensor, el_vert, vert_wts_s);

            // Complex-valued positions of "source" element nodes
            il::StaticArray<std::complex<double>, 3> tau = 
                    make_el_tau_crd(el_vert_s, r_tensor_s);

            // Loop over "Target" elements
            for (il::int_t target_elem = 0; 
                 target_elem < num_of_elems; ++target_elem) {
                // Vertices' coordinates
                il::StaticArray2D<double, 3, 3> el_vert_t;
                //il::StaticArray<double, 3> vert_wts_t;
                for (il::int_t j = 0; j < 3; ++j) {
                    il::int_t n = mesh_conn(j, target_elem);
                    for (il::int_t k = 0; k < 3; ++k) {
                        el_vert_t(k, j) = nodes_crd(k, n);
                    }
                    // set vert_wts_t[j]
                }

                // Rotation tensor for the target element
                il::StaticArray2D<double, 3, 3> r_tensor_t = 
                        make_el_r_tensor(el_vert_t);

                // Normal vector at collocation point (x)
                il::StaticArray<double, 3> nrm_cp_glob;
                for (int j = 0; j < 3; ++j) {
                    nrm_cp_glob[j] = -r_tensor_t(2, j);
                }

                // Collocation points' coordinates
                il::StaticArray<il::StaticArray<double, 3>, 6> el_cp_crd = 
                        el_cp_uniform(el_vert_t, beta);
                //il::StaticArray<il::StaticArray<double, 3>, 6> el_cp_crd = 
                // el_cp_nonuniform(el_vert_t,vert_wts_t,beta);

                il::StaticArray2D<double, 18, 18> trac_infl_el2el;
                // Loop over nodes of the "target" element
                for (int n_t = 0; n_t < 6; ++n_t) {
                    // Shifting to the n_t-th collocation pt
                    HZ hz = make_el_pt_hz
                            (el_vert_s, el_cp_crd[n_t], r_tensor_s);

                    // Calculating DD-to stress influence
                    // w.r. to the source element's local coordinate system
                    il::StaticArray2D<double, 6, 18> stress_infl_el2p_loc_h = 
                            make_local_3dbem_submatrix
                                    (1, mu, nu, hz.h, hz.z, tau, sfm);
                    //stress_infl_el2p_loc_t = make_local_3dbem_submatrix
                    // (0, mu, nu, hz.h, hz.z, tau, sfm);

                    // Multiplication by nrm_cp_glob

                    // Alternative 1: rotating stress at coll. pt.
                    // to the reference ("global") coordinate system
                    //stress_infl_el2p_glob = hfp3d::rotate_sim
                    // (r_tensor_s, stress_infl_el2p_loc_h);
                    //trac_cp_glob = hfp3d::nv_dot_sim(nrm_cp_glob, SIM_CP_G);

                    // Alternative 2: rotating nrm_cp_glob to
                    // the source element's local coordinate system
                    il::StaticArray<double, 3> nrm_cp_loc = 
                            il::dot(r_tensor_s, nrm_cp_glob);
                    il::StaticArray2D<double, 3, 18> trac_el2p_loc = 
                            nv_dot_sim(nrm_cp_loc, stress_infl_el2p_loc_h);
                    il::StaticArray2D<double, 3, 18> trac_cp_glob = il::dot
                            (r_tensor_s, il::Blas::transpose, trac_el2p_loc);

                    // Alternative 3: keeping everything
                    // in terms of local coordinates
                    //trac_cp_x_loc = il::dot
                    // (r_tensor_t, il::Blas::transpose, trac_cp_glob);

                    // Re-relating DD-to traction influence to DD
                    // w.r. to the reference coordinate system
                    il::StaticArray2D<double, 3, 3> trac_infl_n2p, 
                            trac_infl_n2p_glob;
                    for (int n_s = 0; n_s < 6; ++n_s) {
                        // taking a block (one node of the "source" element)
                        for (int j = 0; j < 3; ++j) {
                            for (int k = 0; k < 3; ++k) {
                                trac_infl_n2p(k, j) = 
                                        trac_cp_glob(k, 3 * n_s + j);
                            }
                        }

                        // Coordinate rotation (for the unknown)
                        trac_infl_n2p_glob = il::dot(trac_infl_n2p, r_tensor_s);

                        // Adding the block to the element-to-element
                        // influence sub-matrix
                        for (int j = 0; j < 3; ++j) {
                            for (int k = 0; k < 3; ++k) {
                                trac_infl_el2el(3 * n_t + k, 3 * n_s + j) =
                                        trac_infl_n2p_glob(k, j);
                            }
                        }
                    }
                }

                // Adding the element-to-element influence sub-matrix
                // to the global influence matrix
                //IL_EXPECT_FAST
                // (18 * (target_elem + 1) <= global_matrix.size(0));
                //IL_EXPECT_FAST
                // (18 * (source_elem + 1) <= global_matrix.size(1));
                for (il::int_t j1 = 0; j1 < 6; ++j1) {
                    il::int_t i1 = 6 * source_elem + j1;
                    for (il::int_t j0 = 0; j0 < 6; ++j0) {
                        il::int_t i0 = 6 * target_elem + j0;
                        for (il::int_t l1 = 0; l1 < 3; ++l1) {
                            il::int_t k1 = dof_hndl.dof_h(i1, l1);
                            il::int_t m1 = 3 * j1 + l1;
                            for (il::int_t l0 = 0; l0 < 3; ++l0) {
                                il::int_t k0 = dof_hndl.dof_h(i0, l0);
                                il::int_t m0 = 3 * j0 + l0;
                                if (k0 != -1 && k1 != -1) {
                                    global_matrix(k0, k1) =
                                            trac_infl_el2el(m0, m1);
                                }
                            }
                        }
                    }
                }
            }
        }
        return global_matrix;
    }

// Stress at given points (m_pts_crd) vs DD at nodal points (nodes_crd)

    il::Array2D<double> make_3dbem_stress_f_s
            (double mu, double nu,
             const il::Array2D<il::int_t> &mesh_conn,
             const il::Array2D<double> &nodes_crd,
             const il::Array2D<double> &m_pts_crd,
             //const il::Array2D<double> &m_pts_dsp,
             const bool is_in_glob) {
// Stress at given points (m_pts_crd) vs DD at nodal points (nodes_crd)
// from boundary mesh geometry data:
// mesh connectivity (mesh_conn) and nodes' coordinates (nodes_crd)

// Naive way: no parallelization, no ACA

        IL_EXPECT_FAST(mesh_conn.size(0) >= 3);
        IL_EXPECT_FAST(mesh_conn.size(1) >= 1); // at least 1 element
        IL_EXPECT_FAST(nodes_crd.size(0) >= 3);
        IL_EXPECT_FAST(nodes_crd.size(1) >= 3); // at least 3 nodes

        const il::int_t num_of_elems = mesh_conn.size(1);
        const il::int_t num_of_dof = 18 * num_of_elems;
        const il::int_t num_of_m_pts = m_pts_crd.size(1);

        //IL_EXPECT_FAST(global_matrix.size(0) == 18*num_of_elems);
        //IL_EXPECT_FAST(global_matrix.size(1) == 18*num_of_elems);

        il::Array2D<double> stress_infl_matrix(6 * num_of_m_pts, num_of_dof);

        // Loop over elements
        for (il::int_t source_elem = 0; source_elem < num_of_elems; ++source_elem) {
            // Vertices' coordinates
            il::StaticArray2D<double, 3, 3> el_vert_s;
            //il::StaticArray<double, 3> vert_wts_t;
            for (il::int_t j = 0; j < 3; ++j) {
                il::int_t n = mesh_conn(j, source_elem);
                for (il::int_t k = 0; k < 3; ++k) {
                    el_vert_s(k, j) = nodes_crd(k, n);
                }
                // get vert_wts_s[j]
            }

            // Basis (shape) functions and rotation tensor (r_tensor_s) of the el-t
            il::StaticArray2D<double, 3, 3> r_tensor_s;
            il::StaticArray2D<std::complex<double>, 6, 6> sfm = 
                    make_el_sfm_uniform(el_vert_s, il::io, r_tensor_s);
            //il::StaticArray2D<std::complex<double>, 6, 6> sfm = 
            // make_el_sfm_nonuniform(r_tensor, el_vert, vert_wts_s);

            // Complex-valued positions of "source" element nodes
            il::StaticArray<std::complex<double>, 3> tau = 
                    make_el_tau_crd(el_vert_s, r_tensor_s);

            // Loop over monitoring points
            for (il::int_t m_pt = 0; m_pt < num_of_m_pts; ++m_pt) {
                // Monitoring points' coordinates
                il::StaticArray<double, 3> m_p_crd;
                for (il::int_t j = 0; j < 3; ++j) {
                    m_p_crd[j] = m_pts_crd(j, m_pt);
                }

                // Shifting to the monitoring point
                HZ hz = make_el_pt_hz(el_vert_s, m_p_crd, r_tensor_s);

                // Calculating DD-to stress influence
                // w.r. to the source element's local coordinate system
                il::StaticArray2D<double, 6, 18> stress_infl_el2p_loc_h = 
                        make_local_3dbem_submatrix
                                (1, mu, nu, hz.h, hz.z, tau, sfm);
                //il::StaticArray2D<double, 6, 18> stress_infl_el2p_loc_t = 
                // make_local_3dbem_submatrix
                // (0, mu, nu, hz.h, hz.z, tau, sfm);

                // Rotating stress at m_pt
                // to the reference ("global") coordinate system
                il::StaticArray2D<double, 6, 18> stress_infl_el2p_glob = 
                        rotate_sim(r_tensor_s, stress_infl_el2p_loc_h);

                if (is_in_glob) {
                    // Re-relating DD-to stress influence to DD
                    // w.r. to the reference coordinate system
                    il::StaticArray2D<double, 3, 3> stress_infl_n2p, 
                            stress_infl_n2p_glob;
                    for (int n_s = 0; n_s < 6; ++n_s) {
                        // taking a block (one node of the "source" element)
                        for (int j = 0; j < 3; ++j) {
                            for (int k = 0; k < 6; ++k) {
                                stress_infl_n2p(k, j) = 
                                        stress_infl_el2p_glob(k, 3 * n_s + j);
                            }
                        }

                        // Coordinate rotation (inverse)
                        stress_infl_n2p_glob = 
                                il::dot(stress_infl_n2p, r_tensor_s);

                        // Adding the block to the element-to-point
                        // influence sub-matrix
                        for (int j = 0; j < 3; ++j) {
                            for (int k = 0; k < 6; ++k) {
                                stress_infl_el2p_glob(k, 3 * n_s + j) = 
                                        stress_infl_n2p_glob(k, j);
                            }
                        }
                    }
                }

                // Adding the element-to-point influence sub-matrix
                // to the global stress matrix
                IL_EXPECT_FAST(6 * (m_pt + 1) <= stress_infl_matrix.size(0));
                IL_EXPECT_FAST(18 * (source_elem + 1) <= stress_infl_matrix.size(1));
                for (il::int_t j1 = 0; j1 < 18; ++j1) {
                    for (il::int_t j0 = 0; j0 < 6; ++j0) {
                        stress_infl_matrix
                                (6 * m_pt + j0, 18 * source_elem + j1) = 
                                stress_infl_el2p_glob(j0, j1);
                    }
                }
            }
        }
        return stress_infl_matrix;
        // il::Array<double> disp_vect = ...;
        // il::Array<double> stress_vect = il::dot(stress_infl_matrix, disp_vect);
        // il::Array2D<double> stress_array = ...;
        // return stress_array;
    }

// Element-to-point influence matrix (submatrix of the global one)

    il::StaticArray2D<double, 6, 18>
    make_local_3dbem_submatrix(const int kernel_id,
             double mu, double nu, double h, std::complex<double> z,
             const il::StaticArray<std::complex<double>, 3> &tau,
             const il::StaticArray2D<std::complex<double>, 6, 6> &sfm) {
        // This function assembles a local "stiffness" sub-matrix
        // (influence of DD at the element nodes to stresses at the point z)
        // in terms of a triangular element's local coordinates
        //
        // tau (3) are coordinates of element's vertices and
        // the rows of sfm (6*6) are coefficients of shape functions
        // in terms of the element's own local coordinate system (tau-coordinates);
        // h and z define the position of the (collocation) point x
        // in the same coordinates

        il::StaticArray2D<double, 6, 18> stress_el_2_el_infl{0.0};

        // const std::complex<double> I(0.0, 1.0);

        // scaling ("-" sign comes from traction Somigliana ID, H-term)
        double scale = -mu / (4.0 * M_PI * (1.0 - nu));
        // tolerance parameters
        const double h_tol = 1.0E-16, a_tol = 1.0E-8;

        // tz[m] and d[m] can be calculated here
        il::StaticArray<std::complex<double>, 3> tz, d, dtau;
        std::complex<double> ntau2;
        for (int j = 0; j < 3; ++j) {
            int q = (j + 1) % 3;
            tz[j] = tau[j] - z;
            dtau[j] = tau[q] - tau[j];
            ntau2 = dtau[j] / std::conj(dtau[j]);
            d[j] = 0.5 * (tz[j] - ntau2 * std::conj(tz[j]));
        }
        // also, "shifted" sfm from z, tau[m], and local sfm
        il::StaticArray2D<std::complex<double>, 6, 6> shft_z = shift_el_sfm(z);
        il::StaticArray2D<std::complex<double>, 6, 6> sfm_z =
                il::dot(sfm, shft_z);

        // searching for "degenerate" edges:
        // point x (collocation pt) projects onto an edge line or a vertex
        bool IsDegen = std::abs(d[0]) < h_tol || std::abs(d[1]) < h_tol ||
                       std::abs(d[2]) < h_tol; // (d[0]*d[1]*d[2]==0);
        il::StaticArray2D<bool, 2, 3> is_90_ang{false};

        // calculating angles (phi, psi, chi)
        il::StaticArray<double, 3> phi{0.0}, psi{0.0};
        il::StaticArray2D<double, 2, 3> chi{0.0};
        for (int j = 0; j < 3; ++j) {
            phi[j] = std::arg(tz[j]);
            psi[j] = std::arg(d[j]);
        }
        for (int j = 0; j < 3; ++j) {
            for (int k = 0; k < 2; ++k) {
                int q = (j + k) % 3;
                chi(k, j) = phi[q] - psi[j];
                // make sure it's between -pi and pi (add or subtract 2*pi)
                if (chi(k, j) <= -M_PI)
                    while (chi(k, j) <= -M_PI)
                        chi(k, j) += 2.0 * M_PI;
                else if (chi(k, j) > M_PI)
                    while (chi(k, j) > M_PI)
                        chi(k, j) -= 2.0 * M_PI;
                // reprooving for "degenerate" edges
                // (chi angles too close to 90 degrees)
                if (fabs(M_PI_2 - std::fabs(chi(k, j))) < a_tol) {
                    is_90_ang(k, j) = true;
                    IsDegen = true;
                }
            }
        }

        // DD-to-stress influence
        // [(S11+S22)/2; (S11-S22)/2+i*S12; (S13+i*S23)/2; S33]
        // vs SF monomials (s_ij_infl_mon) and nodal values (s_ij_infl_nod)
        il::StaticArray3D<std::complex<double>, 6, 4, 3> s_ij_infl_mon{0.0};

        // summation over edges
        for (int m = 0; m < 3; ++m) {
            int n = (m + 1) % 3;
            std::complex<double> dm = d[m];
            if (std::abs(dm) >= h_tol && ~is_90_ang(0, m) && ~is_90_ang(1, m)) {
                std::complex<double>
                // exp(I * chi(0, m))
                        eixm = std::exp(std::complex<double>(0.0, chi(0, m))),
                // exp(I * chi(1, m))
                        eixn = std::exp(std::complex<double>(0.0, chi(1, m)));
                // limit case (point x on the element's plane)
                if (std::fabs(h) < h_tol) {
                    il::StaticArray3D<std::complex<double>, 6, 4, 3>
                            s_incr_n =
                            s_integral_lim(kernel_id, nu, eixn, dm),
                            s_incr_m =
                            s_integral_lim(kernel_id, nu, eixm, dm);
                    for (int j = 0; j < 6; ++j) {
                        for (int k = 0; k < 4; ++k) {
                            for (int l = 0; l < 3; ++l) {
                                s_ij_infl_mon(j, k, l) += s_incr_n(j, k, l) -
                                        s_incr_m(j, k, l);
                            }
                        }
                    }
                    // il::blas(1.0, s_incr_n, 1.0, il::io, s_ij_infl_mon);
                    // il::blas(-1.0, s_incr_m, 1.0, il::io, s_ij_infl_mon);
                } else { // out-of-plane case
                    double an = std::abs(tz[n] - dm),
                            am = std::abs(tz[m] - dm);
                    an = (chi(1, m) < 0) ? -an : an;
                    am = (chi(0, m) < 0) ? -am : am;
                    // constituing functions of the integrals
                    il::StaticArray<std::complex<double>, 9>
                    f_n = integral_cst_fun(h, dm, an, chi(1, m), eixn),
                    f_m = integral_cst_fun(h, dm, am, chi(0, m), eixm);
                    // coefficients, by 2nd index:
                    // 0: S11+S22; 1: S11-S22+2*I*S12; 2: S13+S23; 3: S33
                    il::StaticArray4D<std::complex<double>, 6, 4, 3, 9>
                    c_n = s_integral_gen(kernel_id, nu, eixn, h, dm),
                    c_m = s_integral_gen(kernel_id, nu, eixm, h, dm);
                    // combining constituing functions & coefficients
                    il::blas(1.0, c_n, f_n, 1.0, il::io, s_ij_infl_mon);
                    il::blas(-1.0, c_m, f_m, 1.0, il::io, s_ij_infl_mon);
                    // additional terms for "degenerate" case
                    if (IsDegen) {
                        std::complex<double>
                        // exp(I * phi[n])
                                eipn = std::exp(std::complex<double>
                                                        (0.0, phi[n])),
                        // exp(I * phi[m])
                                eipm = std::exp(std::complex<double>
                                                        (0.0, phi[m]));
                        il::StaticArray<std::complex<double>, 5>
                        f_n_red = integral_cst_fun_red(h, dm, an),
                        f_m_red = integral_cst_fun_red(h, dm, am);
                        il::StaticArray4D<std::complex<double>, 6, 4, 3, 5>
                        c_n_red = s_integral_red(kernel_id, nu, eipn, h),
                        c_m_red = s_integral_red(kernel_id, nu, eipm, h);
                        il::blas(1.0, c_n_red, f_n_red, 1.0,
                                 il::io, s_ij_infl_mon);
                        il::blas(-1.0, c_m_red, f_m_red, 1.0,
                                 il::io, s_ij_infl_mon);
                    }
                }
            }
        }

        // contraction with "shifted" sfm (left)
        il::StaticArray3D<std::complex<double>, 6, 4, 3>
        s_ij_infl_nod = il::dot(sfm_z, s_ij_infl_mon);

        // re-shaping and scaling of the resulting matrix
        for (int j = 0; j < 6; ++j) {
            int q = j * 3;
            for (int k = 0; k < 3; ++k) {
                // [S11; S22; S33; S12; S13; S23] vs \delta{u}_k at j-th node
                stress_el_2_el_infl(0, q + k) = 
                        scale * (std::real(s_ij_infl_nod(j, 0, k)) +
                                std::real(s_ij_infl_nod(j, 1, k)));
                stress_el_2_el_infl(1, q + k) = 
                        scale * (std::real(s_ij_infl_nod(j, 0, k)) -
                                std::real(s_ij_infl_nod(j, 1, k)));
                stress_el_2_el_infl(2, q + k) = 
                        scale * std::real(s_ij_infl_nod(j, 3, k));
                stress_el_2_el_infl(3, q + k) = 
                        scale * std::imag(s_ij_infl_nod(j, 1, k));
                stress_el_2_el_infl(4, q + k) = 
                        scale * 2.0 * std::real(s_ij_infl_nod(j, 2, k));
                stress_el_2_el_infl(5, q + k) = 
                        scale * 2.0 * std::imag(s_ij_infl_nod(j, 2, k));
            }
        }
        return stress_el_2_el_infl;
    }

}
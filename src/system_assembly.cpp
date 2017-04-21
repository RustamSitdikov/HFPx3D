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

// Element-to-point influence matrix (submatrix of the global one)
    il::StaticArray2D<double, 6, 18>
    make_local_3dbem_submatrix
            (const int kernel_id,
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

// Static matrix assembly
    il::Array2D<double> make_3dbem_matrix_s
            (double mu, double nu,
             const Mesh_Geom &mesh,
             const Num_Param &n_par,
             il::io_t, DoF_Handle_T &dof_hndl) {
// This function performs BEM matrix assembly from boundary mesh geometry data:
// mesh connectivity (mesh.conn) and nodes' coordinates (mesh.nods)

// Naive way: no parallelization, no ACA

        IL_EXPECT_FAST(mesh.conn.size(0) >= 3);
        IL_EXPECT_FAST(mesh.conn.size(1) >= 1); // at least 1 element
        IL_EXPECT_FAST(mesh.nods.size(0) >= 3);
        IL_EXPECT_FAST(mesh.nods.size(1) >= 3); // at least 3 nodes

        dof_hndl = make_dof_h_triangular(mesh, 2, n_par.tip_type);

        const il::int_t num_ele = mesh.conn.size(1);
        const il::int_t num_dof = dof_hndl.n_dof;
        //const il::int_t num_dof = 18 * num_ele;
        const il::int_t ndpe = dof_hndl.dof_h.size(1);
        IL_EXPECT_FAST(ndpe == 18);

        //IL_EXPECT_FAST(global_matrix.size(0) == 18*num_ele);
        //IL_EXPECT_FAST(global_matrix.size(1) == 18*num_ele);

        il::Array2D<double> global_matrix {num_dof, num_dof, 0.0};
        //il::StaticArray2D<double, num_dof, num_dof> global_matrix;
        //il::StaticArray<double, num_dof> right_hand_side;

        // Loop over "source" elements
        for (il::int_t source_elem = 0;
             source_elem < num_ele; ++source_elem) {
            // Vertices' coordinates
            il::StaticArray2D<double, 3, 3> el_vert_s;
            //il::StaticArray<double, 3> vert_wts_s;
            for (il::int_t j = 0; j < 3; ++j) {
                il::int_t n = mesh.conn(j, source_elem);
                for (il::int_t k = 0; k < 3; ++k) {
                    el_vert_s(k, j) = mesh.nods(k, n);
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
                 target_elem < num_ele; ++target_elem) {
                // Vertices' coordinates
                il::StaticArray2D<double, 3, 3> el_vert_t;
                //il::StaticArray<double, 3> vert_wts_t;
                for (il::int_t j = 0; j < 3; ++j) {
                    il::int_t n = mesh.conn(j, target_elem);
                    for (il::int_t k = 0; k < 3; ++k) {
                        el_vert_t(k, j) = mesh.nods(k, n);
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
                        el_cp_uniform(el_vert_t, n_par.beta);
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
                    // (s_ele_s.r_tensor, stress_infl_el2p_loc_h);
                    //trac_cp_glob = hfp3d::nv_dot_sim(nrm_cp_glob, SIM_CP_G);

                    // Alternative 2: rotating nrm_cp_glob to
                    // the source element's local coordinate system
                    il::StaticArray<double, 3> nrm_cp_loc =
                            il::dot(r_tensor_s, nrm_cp_glob);
                    il::StaticArray2D<double, 3, 18> trac_el2p_loc =
                            nv_dot_sim(nrm_cp_loc, stress_infl_el2p_loc_h);
                    il::StaticArray2D<double, 3, 18> trac_cp_glob =
                            il::dot(r_tensor_s, il::Blas::transpose,
                                    trac_el2p_loc);

                    // Alternative 3: calculating traction
                    // in terms of local coordinates at CP
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
                        trac_infl_n2p_glob =
                                il::dot(trac_infl_n2p, r_tensor_s);

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
                // (ndpe * (target_elem + 1) <= global_matrix.size(0));
                //IL_EXPECT_FAST
                // (ndpe * (source_elem + 1) <= global_matrix.size(1));
                for (il::int_t i1 = 0; i1 < ndpe; ++i1) {
                    il::int_t j1 = dof_hndl.dof_h(source_elem, i1);
                    for (il::int_t i0 = 0; i0 < ndpe; ++i0) {
                        il::int_t j0 = dof_hndl.dof_h(target_elem, i0);
                        if (j0 >= 0 && j1 >= 0) {
                            global_matrix(j0, j1) +=
                                    trac_infl_el2el(i0, i1);
                        }
                    }
                }
            }
        }
        return global_matrix;
    }

// Volume Control matrix assembly (additional row $ column)
    il::Array2D<double> make_3dbem_matrix_vc
            (double mu, double nu,
             const Mesh_Geom &mesh,
             const Num_Param &n_par,
             il::io_t, DoF_Handle_T &dof_hndl) {
// This function performs Volume Control BEM matrix assembly
// from boundary mesh geometry data:
// mesh connectivity (mesh.conn) and nodes' coordinates (mesh.nods)

// Naive way: no parallelization, no ACA

        IL_EXPECT_FAST(mesh.conn.size(0) >= 3);
        IL_EXPECT_FAST(mesh.conn.size(1) >= 1); // at least 1 element
        IL_EXPECT_FAST(mesh.nods.size(0) >= 3);
        IL_EXPECT_FAST(mesh.nods.size(1) >= 3); // at least 3 nodes

        if (dof_hndl.n_dof == 0 || dof_hndl.dof_h.size(0) == 0) {
            dof_hndl = make_dof_h_triangular(mesh, 2, n_par.tip_type);
        }

        const il::int_t num_ele = mesh.conn.size(1);
        const il::int_t num_dof = dof_hndl.n_dof;
        const il::int_t ndpe = dof_hndl.dof_h.size(1);
        IL_EXPECT_FAST(ndpe == 18);

        il::Array2D<double> global_matrix {num_dof + 1, num_dof + 1, 0.0};
        //il::Array<double> right_hand_side {num_dof, 0.0};
        //Alg_Sys_T alg_system;
        //alg_sys.matrix = il::Array2D<double>{num_dof+1, num_dof+1, 0.0};
        //alg_sys.rhside = il::Array<double>{num_dof+1, 0.0};

        // Loop over "source" elements
        for (il::int_t source_elem = 0; 
             source_elem < num_ele; ++source_elem) {
            // Vertices' coordinates
            il::StaticArray2D<double, 3, 3> el_vert_s;
            //il::StaticArray<double, 3> vert_wts_s;
            for (il::int_t j = 0; j < 3; ++j) {
                il::int_t n = mesh.conn(j, source_elem);
                for (il::int_t k = 0; k < 3; ++k) {
                    el_vert_s(k, j) = mesh.nods(k, n);
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
                 target_elem < num_ele; ++target_elem) {
                // Vertices' coordinates
                il::StaticArray2D<double, 3, 3> el_vert_t;
                //il::StaticArray<double, 3> vert_wts_t;
                for (il::int_t j = 0; j < 3; ++j) {
                    il::int_t n = mesh.conn(j, target_elem);
                    for (il::int_t k = 0; k < 3; ++k) {
                        el_vert_t(k, j) = mesh.nods(k, n);
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
                        el_cp_uniform(el_vert_t, n_par.beta);
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

                    // Multiplication by normal at CP

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

                    // Alternative 3: traction vector
                    // in terms of local coordinates at CP
                    //trac_cp_x_loc = il::dot
                    // (r_tensor_t, il::Blas::transpose, trac_cp_glob);

                    if (n_par.is_dd_in_glob) {
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
                            trac_infl_n2p_glob = il::dot(trac_infl_n2p,
                                                         r_tensor_s);

                            // Adding the block to the element-to-element
                            // influence sub-matrix
                            for (int j = 0; j < 3; ++j) {
                                for (int k = 0; k < 3; ++k) {
                                    trac_infl_el2el(3 * n_t + k, 3 * n_s + j) =
                                            trac_infl_n2p_glob(k, j);
                                }
                            }
                        }
                    } else {
                        for (int dof_s = 0; dof_s < ndpe; ++dof_s) {
                            for (int n_t = 0; n_t < 6; ++n_t) {
                                for (int k = 0; k < 3; ++k)
                                trac_infl_el2el(3 * n_t + k, dof_s) =
                                        trac_cp_glob(k, dof_s);
                            }
                        }
                    }
                }

                // Adding the element-to-element influence sub-matrix
                // to the global influence matrix
                //IL_EXPECT_FAST
                // (ndpe * (target_elem + 1) <= global_matrix.size(0));
                //IL_EXPECT_FAST
                // (ndpe * (source_elem + 1) <= global_matrix.size(1));
                for (il::int_t i1 = 0; i1 < ndpe; ++i1) {
                    il::int_t j1 = dof_hndl.dof_h(source_elem, i1);
                    for (il::int_t i0 = 0; i0 < ndpe; ++i0) {
                        il::int_t j0 = dof_hndl.dof_h(target_elem, i0);
                        if (j0 >= 0 && j1 >= 0) {
                            global_matrix(j0, j1) +=
                                    trac_infl_el2el(i0, i1);
                        }
                    }
                }
            }

            // Influence of DD & pressure on tractions & volume
            il::StaticArray<double, 6> el_sf_integral =
                    el_p2_sf_integral(sfm, tau);
            for (int n_s = 0; n_s < 6; ++n_s) {
                // Integral of n_s-th shape function over the s-element
                double sf_integral = el_sf_integral[n_s];
                il::StaticArray<double, 3> sf_i_v {0.0};
                // Integral of normal DD (opening) over the element
                // for the n_s-th shape function
                if (n_par.is_dd_in_glob) {
                    // dot([0, 0, sf_integral], r_tensor_s)
                    for (int j = 0; j < 3; ++j) {
                        sf_i_v[j] = sf_integral * r_tensor_s(2, j);
                    }
                } else {
                    // [0, 0, sf_integral]
                    sf_i_v[2] = sf_integral;
                }
                for (int j = 0; j < 3; ++j) {
                    int l = n_s * 3 + j;
                    il::int_t s_dof = dof_hndl.dof_h(source_elem, l);
                    if (s_dof >= 0) {
                        // Volume vs DD
                        global_matrix(num_dof, s_dof) = sf_i_v[j];
                        // Tractions vs pressure
                        global_matrix(s_dof, num_dof) =
                                -r_tensor_s(2, j); // Normal at element
                    }
                }
            }
        }
        // global_matrix(num_dof, num_dof) = compressibility * volume
        return global_matrix;
    }

// Volume Control system modification (for DD increments)
    Alg_Sys_T mod_3dbem_system_vc
            (const il::Array2D<double> &orig_matrix,
             const DoF_Handle_T &orig_dof_hndl,
             const DoF_Handle_T &dof_hndl,
             const il::Array<double> &delta_t,
             const double delta_v) {
// Truncated matrix & RHS assembly from given original BEM matrix
// and original & "truncated" DoF handles (free & fixed degrees of freedom)
        IL_EXPECT_FAST(orig_matrix.size(0) == orig_matrix.size(1));
        const il::int_t num_of_ele = orig_dof_hndl.dof_h.size(0);
        const il::int_t ndpe = orig_dof_hndl.dof_h.size(1);
        IL_EXPECT_FAST(num_of_ele > 0);
        const il::int_t full_ndof = num_of_ele * ndpe;
        const il::int_t orig_ndof = orig_dof_hndl.n_dof;
        IL_EXPECT_FAST(orig_ndof > 0 && orig_ndof <= full_ndof);
        IL_EXPECT_FAST(orig_ndof == orig_matrix.size(0));
        const il::int_t used_ndof = dof_hndl.n_dof;
        // check if the used matrix is smaller that the original matrix
        IL_EXPECT_FAST(used_ndof > 0 && used_ndof <= orig_ndof);
        Alg_Sys_T alg_system;
        // RHS
        alg_system.rhside = il::Array<double>{used_ndof + 1};
        // (sought volume delta)
        alg_system.rhside[used_ndof] = delta_v;
        const il::int_t tsize = delta_t.size();
        IL_EXPECT_FAST( tsize == full_ndof ||
                        tsize == orig_ndof ||
                        tsize == used_ndof );
        alg_system.matrix = il::Array2D<double>{used_ndof + 1,
                                                used_ndof + 1, 0.0};
        // Loop over nodes
        for (il::int_t s_ele = 0; s_ele < num_of_ele; ++s_ele) {
            for (int j = 0; j < ndpe; ++j) {
                il::int_t s_dof = dof_hndl.dof_h(s_ele, j);
                if (s_dof >= 0) {
                    il::int_t o_s_dof = orig_dof_hndl.dof_h(s_ele, j);
                    for (il::int_t t_ele = 0; t_ele < num_of_ele; ++t_ele) {
                        for (int k = 0; k < ndpe; ++k) {
                            il::int_t t_dof = dof_hndl.dof_h(t_ele, k);
                            if (t_dof >= 0) {
                                il::int_t o_t_dof =
                                        orig_dof_hndl.dof_h(t_ele, k);
                                alg_system.matrix(t_dof, s_dof) +=
                                        orig_matrix(o_t_dof, o_s_dof);
                            }
                        }
                    }
                    // Volume vs DD
                    alg_system.matrix(used_ndof, s_dof) +=
                            orig_matrix(full_ndof, o_s_dof);
                    // Traction vs pressure
                    alg_system.matrix(s_dof, used_ndof) +=
                            orig_matrix(o_s_dof, full_ndof);
                    // RHS (sought traction delta)
                    if (tsize == full_ndof) {
                        il::int_t f_s_dof = s_ele * ndpe + j;
                        alg_system.rhside[s_dof] = delta_t[f_s_dof];
                    } else if (tsize == orig_ndof) {
                        alg_system.rhside[s_dof] = delta_t[o_s_dof];
                    } else {
                        alg_system.rhside[s_dof] = delta_t[s_dof];
                    }
                }
            }
        }
        return alg_system;
    }

// use std::function<double(double)> cohesion
// or std::function<il::StaticArray<double, 3>
// (il::StaticArray<double, 3>)> fric_cohesion
// [friction, shear cohesion, opening cohesion]
// vs [shear 1, shear 2, normal] traction

// Stress at given points (m_pts_crd) vs DD at nodal points (mesh.nods)
    il::Array2D<double> make_3dbem_stress_f_s
            (double mu, double nu,
             const Mesh_Geom &mesh,
             const Num_Param &n_par,
             // const Mesh_Data &m_data,
             const il::Array2D<double> &m_pts_crd) {
// This function calculates Stress at given points (m_pts_crd)
// vs DD (m_data.DD) at nodal points (mesh.nods)
// using boundary mesh geometry data:
// mesh connectivity (mesh.conn) and nodes' coordinates (mesh.nods)

// Naive way: no parallelization, no ACA

        IL_EXPECT_FAST(mesh.conn.size(0) >= 3);
        IL_EXPECT_FAST(mesh.conn.size(1) >= 1); // at least 1 element
        IL_EXPECT_FAST(mesh.nods.size(0) >= 3);
        IL_EXPECT_FAST(mesh.nods.size(1) >= 3); // at least 3 nodes

        const il::int_t num_ele = mesh.conn.size(1);
        const il::int_t num_dof = 18 * num_ele;
        const il::int_t num_of_m_pts = m_pts_crd.size(1);

        il::Array2D<double> stress_infl_matrix(6 * num_of_m_pts, num_dof);

        // Loop over elements
        for (il::int_t source_elem = 0; source_elem < num_ele; ++source_elem) {
            // Vertices' coordinates
            il::StaticArray2D<double, 3, 3> el_vert_s;
            //il::StaticArray<double, 3> vert_wts_t;
            for (il::int_t j = 0; j < 3; ++j) {
                il::int_t n = mesh.conn(j, source_elem);
                for (il::int_t k = 0; k < 3; ++k) {
                    el_vert_s(k, j) = mesh.nods(k, n);
                }
                // get vert_wts_s[j]
            }

            // Basis (shape) functions and
            // rotation tensor (r_tensor_s) of the element (source_elem)
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

                if (n_par.is_dd_in_glob) {
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
                IL_EXPECT_FAST(18 * (source_elem + 1) <=
                                       stress_infl_matrix.size(1));
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
        // il::Array<double> disp_vect{num_dof};
        // for (il::int_t k = 0; k < 6 * num_ele; ++k) {
        // for (il::int_t j = 0; j < 3; ++j) {
        // il::int_t l = 3 * k + j;
        // disp_vect[l] = m_data.DD(k, j);
        // }
        // }
        // il::Array<double> stress_vect{6 * num_of_m_pts} =
        // il::dot(stress_infl_matrix, disp_vect);
        // il::Array2D<double> stress_array{num_of_m_pts, 6};
        // for (il::int_t k = 0; k < num_of_m_pts; ++k) {
        // for (il::int_t j = 0; j < 6; ++j) {
        // il::int_t l = 6 * k + j;
        // stress_array(k, j) = stress_vec[l];
        // }
        // }
        // return stress_array;
    }

}
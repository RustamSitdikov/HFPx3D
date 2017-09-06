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
//#include <il/math.h>
#include <il/Array2D.h>
#include <il/StaticArray.h>
#include <il/StaticArray2D.h>
#include <il/StaticArray3D.h>
#include <il/StaticArray4D.h>
#include <il/linear_algebra.h>
// #include <il/linear_algebra/dense/blas/dot.h>
// #include <il/linear_algebra/dense/blas/blas.h>
#include "src/Core/constants.h"
#include "src/Core/tensor_utilities.h"
#include "src/Core/element_utilities.h"
#include "src/Elasticity/elasticity_kernel_integration.h"
#include "system_assembly.h"

namespace hfp3d {

    // Static matrix assembly
    il::Array2D<double> make_3dbem_matrix_s
            (double mu, double nu,
             const Mesh_Geom_T &mesh,
             const Num_Param_T &n_par,
             il::io_t, DoF_Handle_T &dof_hndl) {
    // This function performs BEM matrix assembly from boundary mesh geometry data:
    // mesh connectivity (mesh.conn) and nodes' coordinates (mesh.nods)

    // Naive way: no ACA. For parallel assembly uncomment line 236

        IL_EXPECT_FAST(mesh.conn.size(0) >= 3);
        IL_EXPECT_FAST(mesh.conn.size(1) >= 1); // at least 1 element
        IL_EXPECT_FAST(mesh.nods.size(0) >= 3);
        IL_EXPECT_FAST(mesh.nods.size(1) >= 3); // at least 3 nodes

        if (dof_hndl.n_dof == 0 || dof_hndl.dof_h.size(0) == 0) {
            dof_hndl = make_dof_h_crack(mesh, 2, n_par.tip_type);
        }

        const il::int_t num_ele = mesh.conn.size(1);
        const il::int_t num_dof = dof_hndl.n_dof;
        //const il::int_t num_dof = 18 * num_ele;
        const il::int_t ndpe = dof_hndl.dof_h.size(1);
        IL_EXPECT_FAST(ndpe == 18);

        il::Array2D<double> global_matrix {num_dof, num_dof, 0.0};
        //il::StaticArray2D<double, num_dof, num_dof> global_matrix;
        //il::StaticArray<double, num_dof> right_hand_side;

        // Loop over "source" elements

#pragma omp parallel for

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
                            il::dot(r_tensor_s, il::Blas::Transpose,
                                    trac_el2p_loc); // il::Blas::kTranspose

                    // Alternative 3: calculating traction
                    // in terms of local coordinates at CP
                    //if (!n_par.is_tr_local) {
                    //    trac_cp_x_loc = il::dot
                    //            (r_tensor_t, il::Blas::kTranspose,
                    //        trac_cp_glob);
                    //}

                    // in case if DDs are sought in local coordinates
                    if (!n_par.is_dd_local) {
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

                            // Coordinate rotation (for the unknown DD)
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
                        // adding DD-to traction influence "as is" (local)
                        // todo: may need revision
                        for (int dof_s = 0; dof_s < ndpe; ++dof_s) {
                            for (int k = 0; k < 3; ++k) {
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
        }
        //IL_EXPECT_FAST(global_matrix.size(0) == 18*num_ele);
        //IL_EXPECT_FAST(global_matrix.size(1) == 18*num_ele);
        return global_matrix;
    }

    // Add S_inf (induced tractions) to the RHS
    void add_s_inf_to_3dbem_rhs
            (const Mesh_Data_T &mesh_data,
             const Load_T &load,
             il::io_t, SAE_T &sae) {
        // This function adds the tractions induced by in-situ stress
        // defined in "load" to the RHS of the system "sae"
        il::int_t num_elems = mesh_data.mesh.conn.size(1);
        il::int_t num_dof = mesh_data.dof_h_dd.n_dof;
        if (sae.rhs_v.size() != num_dof) {
            sae.rhs_v = il::Array<double> {num_dof, 0.0};
        }

#pragma omp parallel for

        for (il::int_t el = 0; el < num_elems; ++el) {
            // element vertices' coordinates
            il::StaticArray2D<double, 3, 3> el_vert;
            for (il::int_t j = 0; j < 3; ++j) {
                il::int_t n = mesh_data.mesh.conn(j, el);
                for (il::int_t k = 0; k < 3; ++k) {
                    el_vert(k, j) = mesh_data.mesh.nods(k, n);
                }
            }

            // Rotation tensor for the element
            il::StaticArray2D<double, 3, 3> r_tensor =
                    hfp3d::make_el_r_tensor(el_vert);

            // Normal vector at collocation point (x)
            il::StaticArray<double, 3> norm_v;
            for (int j = 0; j < 3; ++j) {
                norm_v[j] = -r_tensor(2, j);
            }

            // induced traction (from in-situ stress)
            il::StaticArray<double, 3> trac_inf =
                    hfp3d::nv_dot_sim(norm_v, load.s_inf);

            // setting the RHS of the system (loads)
            for (int ln = 0; ln < 6; ++ln) {
                //il::int_t n = el * 6 + ln;
                for (int l = 0; l < 3; ++l) {
                    il::int_t ldof = ln * 3 + l;
                    //il::int_t dof = dof_handle.dof_h(el, ldof);
                    il::int_t dof = mesh_data.dof_h_dd.dof_h(el, ldof);
                    if (dof != -1) {
                        //rhs[dof] = (l != 1 ? 1.0 : 0.0);
                        sae.rhs_v[dof] -= trac_inf[l];
                    }
                }
            }
        }
    }


    // Stress at given points (m_pts_crd) vs DD (m_data.dd)
    // at nodal points (mesh.nods)
    il::Array2D<double> make_3dbem_stress_f_s
            (double mu, double nu,
             const Mesh_Data_T &m_data,
             const Load_T &load,
             const Num_Param_T &n_par,
             const il::Array2D<double> &m_pts_crd) {
    // This function calculates Stress at given points (m_pts_crd)
    // vs DD (m_data.DD) at nodal points (mesh.nods)
    // using boundary mesh geometry data:
    // mesh connectivity (mesh.conn) and nodes' coordinates (mesh.nods)

    // Naive way: no ACA. For parallel assembly, uncomment line 424

        IL_EXPECT_FAST(m_data.mesh.conn.size(0) >= 3);
        IL_EXPECT_FAST(m_data.mesh.conn.size(1) >= 1); // at least 1 element
        IL_EXPECT_FAST(m_data.mesh.nods.size(0) >= 3);
        IL_EXPECT_FAST(m_data.mesh.nods.size(1) >= 3); // at least 3 nodes

        const il::int_t num_ele = m_data.mesh.conn.size(1);
//        const il::int_t ndpe = 18;
        const il::int_t num_dof = 18 * num_ele;
        const il::int_t num_of_m_pts = m_pts_crd.size(1);

//        il::Array2D<double> stress_infl_matrix{6 * num_of_m_pts, num_dof, 0.0};
        il::Array2D<double> stress_array{num_of_m_pts, 6, 0.0};

        // Loop over elements

#pragma omp parallel for

        for (il::int_t source_elem = 0; source_elem < num_ele; ++source_elem) {
            // Vertices' coordinates
            il::StaticArray2D<double, 3, 3> el_vert_s;
            //il::StaticArray<double, 3> vert_wts_t;
            for (il::int_t j = 0; j < 3; ++j) {
                il::int_t n = m_data.mesh.conn(j, source_elem);
                for (il::int_t k = 0; k < 3; ++k) {
                    el_vert_s(k, j) = m_data.mesh.nods(k, n);
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

            // Loop over observation points
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

                // in cse if DDs are given in local coordinates of each element
                if (!n_par.is_dd_local) {
                    // Re-relating DD-to stress influence to DD
                    // w.r. to the reference coordinate system
                    il::StaticArray2D<double, 6, 3> stress_infl_n2p_loc,
                            stress_infl_n2p_glob;
                    for (int n_s = 0; n_s < 6; ++n_s) {
                        // taking a block (one node of the "source" element)
                        for (int j = 0; j < 3; ++j) {
                            for (int k = 0; k < 6; ++k) {
                                stress_infl_n2p_loc(k, j) =
                                        stress_infl_el2p_glob(k, 3 * n_s + j);
                            }
                        }

                        // Coordinate rotation (inverse) for DD
                        stress_infl_n2p_glob =
                                il::dot(stress_infl_n2p_loc, r_tensor_s);

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

//                // Adding the element-to-point influence sub-matrix
//                // to the global stress matrix
//                IL_EXPECT_FAST(6 * (m_pt + 1) <= stress_infl_matrix.size(0));
//                IL_EXPECT_FAST(18 * (source_elem + 1) <=
//                               stress_infl_matrix.size(1));
//                for (il::int_t j1 = 0; j1 < 18; ++j1) {
//                    for (il::int_t j0 = 0; j0 < 6; ++j0) {
//                        stress_infl_matrix
//                                (6 * m_pt + j0, ndpe * source_elem + j1) +=
//                                stress_infl_el2p_glob(j0, j1);
//                    }
//                }

                // re-arranging local DDs into a vector
                il::StaticArray<double, 18> ele_disp_vect;
                for (il::int_t k = 0; k < 6; ++k) {
                    for (il::int_t j = 0; j < 3; ++j) {
                        il::int_t l = 3 * k + j;
                        ele_disp_vect[l] = m_data.dd(k, j);
                    }
                }
                // Multiplying by element DD and adding the result to stresses
                il::StaticArray<double, 6> part_stress_vect{0.0};
                il::blas(1.0, stress_infl_el2p_glob, ele_disp_vect,
                         1.0, il::io, part_stress_vect);
                for (il::int_t j = 0; j < 6; ++j) {
                    stress_array(m_pt, j) += part_stress_vect[j];
                }
            }
        }

//        // re-arranging DDs into a vector
//        il::Array<double> disp_vect{num_dof};
//        for (il::int_t k = 0; k < 6 * num_ele; ++k) {
//            for (il::int_t j = 0; j < 3; ++j) {
//                il::int_t l = 3 * k + j;
//                disp_vect[l] = m_data.dd(k, j);
//            }
//        }
//        // multiplying stress_infl_matrix by disp_vect
//        il::Array<double> stress_vect{6 * num_of_m_pts};
//        IL_EXPECT_FAST(stress_vect.size() == stress_infl_matrix.size(0));
//        IL_EXPECT_FAST(disp_vect.size() == stress_infl_matrix.size(1));
//        stress_vect = il::dot(stress_infl_matrix, disp_vect);
//        // re-arranging stress components into num_of_m_pts*6 matrix
//        for (il::int_t k = 0; k < num_of_m_pts; ++k) {
//            for (il::int_t j = 0; j < 6; ++j) {
//                il::int_t l = 6 * k + j;
//                stress_array(k, j) = stress_vect[l];
//            }
//        }

        // Adding in-situ stress

#pragma omp parallel for

        // Loop over observation points
        for (il::int_t m_pt = 0; m_pt < num_of_m_pts; ++m_pt) {
            for (il::int_t j = 0; j < 6; ++j) {
                stress_array(m_pt, j) += load.s_inf[j];
            }
        }

        // return stress_infl_matrix;
        return stress_array;
    }

    // Volume Control matrix assembly (additional row $ column)
    il::Array2D<double> make_3dbem_matrix_vc
            (double mu, double nu,
             const Mesh_Geom_T &mesh,
             const Num_Param_T &n_par,
             il::io_t, DoF_Handle_T &dof_hndl) {
    // This function performs Volume Control BEM matrix assembly
    // from boundary mesh geometry data:
    // mesh connectivity (mesh.conn) and nodes' coordinates (mesh.nods)

    // Naive way: no ACA. For parallel assembly, uncomment line 580

        IL_EXPECT_FAST(mesh.conn.size(0) >= 3);
        IL_EXPECT_FAST(mesh.conn.size(1) >= 1); // at least 1 element
        IL_EXPECT_FAST(mesh.nods.size(0) >= 3);
        IL_EXPECT_FAST(mesh.nods.size(1) >= 3); // at least 3 nodes

        if (dof_hndl.n_dof == 0 || dof_hndl.dof_h.size(0) == 0) {
            dof_hndl = make_dof_h_crack(mesh, 2, n_par.tip_type);
        }

        const il::int_t num_ele = mesh.conn.size(1);
        const il::int_t num_dof = dof_hndl.n_dof;
        //const il::int_t ndpe = dof_hndl.dof_h.size(1);
        //IL_EXPECT_FAST(ndpe == 18);

        il::Array2D<double> global_matrix;
        global_matrix.reserve(num_dof + 1, num_dof + 1);
        global_matrix = make_3dbem_matrix_s
                        (mu, nu, mesh, n_par, il::io, dof_hndl);
        IL_EXPECT_FAST(global_matrix.size(0) == num_dof);
        IL_EXPECT_FAST(global_matrix.size(1) == num_dof);

        // adding a row and a column for Volume Control
        global_matrix.resize(num_dof + 1, num_dof + 1);

        // Loop over "source" elements

#pragma omp parallel for

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

            // Influence of DD & pressure on tractions & volume
            il::StaticArray<double, 6> el_sf_integral =
                    el_p2_sf_integral(sfm, tau);

            for (int n_s = 0; n_s < 6; ++n_s) {
                // Integral of n_s-th shape function over the s-element
                double sf_integral = el_sf_integral[n_s];
                il::StaticArray<double, 3> sf_i_v {0.0};
                // Integral of normal DD (opening) over the element
                // for the n_s-th shape function

                // in case if DDs are given in local coordinates of each element
                if (!n_par.is_dd_local) {
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
                        // Tractions vs pressure (normal at each CP)
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
    SAE_T mod_3dbem_system_vc
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
        SAE_T alg_system;
        // RHS
        alg_system.rhs_v = il::Array<double>{used_ndof + 1};
        // (sought volume delta)
        alg_system.rhs_v[used_ndof] = delta_v;
        const il::int_t tsize = delta_t.size();
        IL_EXPECT_FAST( tsize == full_ndof ||
                        tsize == orig_ndof ||
                        tsize == used_ndof );
        alg_system.matrix = il::Array2D<double>{used_ndof + 1,
                                                used_ndof + 1, 0.0};
        // Loop over elements & nodes

#pragma omp parallel for

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
                        alg_system.rhs_v[s_dof] = delta_t[f_s_dof];
                    } else if (tsize == orig_ndof) {
                        alg_system.rhs_v[s_dof] = delta_t[o_s_dof];
                    } else {
                        alg_system.rhs_v[s_dof] = delta_t[s_dof];
                    }
                }
            }
        }
        return alg_system;
    }

}
//
// This file is part of HFPx3D.
//
// Created by nikolski on 11/14/2017.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//

#include <cmath>
#include <random>
#include <il/math.h>
#include <il/Status.h>
#include <il/Array2D.h>
#include <il/linear_algebra.h>
#include <il/linear_algebra/dense/factorization/LU.h>
#include "src/Solvers/system_assembly.h"
#include "src/Core/model_parameters.h"
#include "src/Core/surface_mesh_utilities.h"
#include "src/Core/element_utilities.h"
#include "src/Core/constants.h"
#include "src/IO/mesh_file_io.h"
#include "test_radial_crack_static.h"

namespace hfp3d {

    double test_dd_radial_crack_static
            (std::string work_directory,
             std::string m_c_f_name,
             std::string m_n_f_name,
             int tt, double tol, bool echo) {

        // initializing the mesh etc.
        Mesh_Data_T mesh_data;
        hfp3d::Properties_T mat_props;
        hfp3d::Num_Param_T num_param;
        hfp3d::Load_T load_data;

        il::Status status{};

        // elastic properties
        mat_props.shear_m[0] = 1.0, mat_props.poiss_r[0] = 0.35;

        hfp3d::load_mesh_from_numpy_32
                (work_directory, m_c_f_name, m_n_f_name, 1,
                 il::io, mesh_data.mesh);

        // numerical model parameters
        num_param.tip_type = tt;
        num_param.is_dd_local = false;
        num_param.beta = 0.125;

        // load (S_inf)
        load_data.s_inf = il::StaticArray<double, 6> {0.0};
        load_data.s_inf[2] = 1.0;

        // initializing the DoF handle
        mesh_data.dof_h_dd = hfp3d::make_dof_h_crack
                (mesh_data.mesh, 2, num_param.tip_type);

        // assembly of the algebraic system (matrix + RHS)
        hfp3d::SAE_T sae;
        // matrix
        sae.matrix = hfp3d::make_3dbem_matrix_s
                (mat_props.shear_m[0], mat_props.poiss_r[0],
                 mesh_data.mesh,
                 num_param,
                 il::io, mesh_data.dof_h_dd);
        // right-hand side
        hfp3d::add_s_inf_to_3dbem_rhs
                (mesh_data, load_data, il::io, sae);

        // solving the system
        il::Array<double> dd_v{sae.n_dof};
        il::LU<il::Array2D<double>> lu_decomposition(sae.matrix, il::io, status);
        // if (!status.ok()) {
        //     // The matrix is singular to the machine precision.
        //     // You should deal with the error.
        // }
        // double cnd = lu_decomposition.conditionNumber(il::Norm::L2, );
        // std::cout << cnd << std::endl;
        status.abortOnError();
        dd_v = lu_decomposition.solve(sae.rhs_v);
//        dd_v = il::linearSolve(sae.matrix, sae.rhs_v, il::io, status);
//        status.abortOnError();

        // comparison with the reference solution
        il::int_t n_ex_n = 0;
        il::int_t num_el = mesh_data.mesh.conn.size(1);
        double err = 0.0;
        double r_tol = 0.2 * std::sqrt(hfp3d::pi / num_el / std::sqrt(3.0));
        for (il::int_t el = 0; el < num_el; ++el) {
            // Vertices' coordinates
            il::StaticArray2D<double, 3, 3> el_vert;
            for (il::int_t j = 0; j < 3; ++j) {
                il::int_t n = mesh_data.mesh.conn(j, el);
                for (il::int_t k = 0; k < 3; ++k) {
                    el_vert(k, j) = mesh_data.mesh.nods(k, n);
                }
            }
            // All 6 nodes' coordinates
            il::StaticArray<il::StaticArray<double, 3>, 6> el_np;
            el_np = hfp3d::el_cp_uniform(el_vert, 0.0);
            // DD vectors at 6 nodes
            for (int np = 0; np < 6; ++np) {
                double r;
                r = el_np[np][0] * el_np[np][0] + el_np[np][1] * el_np[np][1];
                r = std::sqrt(r);
                r = std::fabs(r - 1.0);
                // excluding nodes close to the tip
                if ((num_param.tip_type == 1 && np < 3 && r < r_tol)
                    || (num_param.tip_type == 0 && r < r_tol)
                    || (num_param.tip_type == 2 && r < r_tol)) {
                    ++n_ex_n;
                    if (echo) {
                        std::cout << "[          ] nod. " << np+1
                                  << " of el. " << el+1
                                  << " excluded "
                                  << std::endl;
                    }
                } else {
                    il::StaticArray<double, 3> np_dd;
                    for (int j = 0; j < 3; ++j) {
                        il::int_t local_dof, global_dof;
                        local_dof = np * 3 + j;
                        global_dof = mesh_data.dof_h_dd.dof_h(el, local_dof);
                        if (global_dof == -1) {
                            np_dd[j] = 0.0;
                        } else {
                            np_dd[j] = dd_v[global_dof];
                        }
                    }
                    // calculating the the reference solution
                    il::StaticArray<double, 3> ref_dd;
                    ref_dd = ref::get_dd_at_pt
                            (mat_props.shear_m[0],
                             mat_props.poiss_r[0],
                             1.0,
                             load_data.s_inf[2],
                             el_np[np]);
                    // calculating the difference with the reference solution
                    for (int j = 0; j < 3; ++j) {
                        double diff_ij = std::fabs(np_dd[j] - ref_dd[j]);
                        if (diff_ij > tol && echo) {
                            std::cout << "[          ] nod. " << np+1
                                      << " of el. "<< el+1
                                      << ", dof " << j+1
                                      << ", err = " << diff_ij
                                      << std::endl;
                        }
                        if (diff_ij > err) {
                            err = diff_ij;
                        }
                    }
                }
            }
        }
        //if (echo) {
            std::cout << "[          ] " << n_ex_n
                      << " of " << num_el * 6
                      << " nodes excluded "
                      << std::endl;
        //}
        return err;
    }

    double test_stresses_radial_crack_static
            (std::string work_directory,
             std::string m_c_f_name,
             std::string m_n_f_name,
             il::int_t n_m_p,
             double tol, bool echo) {

        // initializing the mesh etc.
        Mesh_Data_T mesh_data;
        hfp3d::Properties_T mat_props;
        hfp3d::Num_Param_T num_param;
        hfp3d::Load_T load_data;

        il::Status status{};

        // elastic properties
        mat_props.shear_m[0] = 1.0, mat_props.poiss_r[0] = 0.35;

        hfp3d::load_mesh_from_numpy_32
                (work_directory, m_c_f_name, m_n_f_name, 1,
                 il::io, mesh_data.mesh);

        // numerical model parameters
        num_param.tip_type = 1;
        num_param.is_dd_local = false;
        num_param.beta = 0.125;

        // load (S_inf)
        load_data.s_inf = il::StaticArray<double, 6> {0.0};
        load_data.s_inf[2] = 1.0;

        // substituting the reference solution
        il::int_t num_el = mesh_data.mesh.conn.size(1);
        mesh_data.dd = il::Array2D<double> {num_el * 6, 3, 0.0};
        for (il::int_t el = 0; el < num_el; ++el) {
            // Vertices' coordinates
            il::StaticArray2D<double, 3, 3> el_vert;
            for (int j = 0; j < 3; ++j) {
                il::int_t n = mesh_data.mesh.conn(j, el);
                for (il::int_t k = 0; k < 3; ++k) {
                    el_vert(k, j) = mesh_data.mesh.nods(k, n);
                }
            }
            // All 6 nodes' coordinates
            il::StaticArray<il::StaticArray<double, 3>, 6> el_n_c;
            el_n_c = hfp3d::el_cp_uniform(el_vert, 0.0);
            // DD vectors at 6 nodes
            for (int np = 0; np < 6; ++np) {
                // calculating the reference solution for dd
                il::StaticArray<double, 3> ref_dd;
                ref_dd = ref::get_dd_at_pt
                        (mat_props.shear_m[0],
                         mat_props.poiss_r[0],
                         1.0,
                         load_data.s_inf[2],
                         el_n_c[np]);
                // substituting the reference solution into dd of mesh_data
                il::int_t nn = el * 6 + np;
                for (int j = 0; j < 3; ++j) {
                    mesh_data.dd(nn, j) = ref_dd[j];
                }
            }
        }

        // defining points to monitor stresses
        std::random_device rd;  //Will be used to obtain a seed for the random number engine
        std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
        std::uniform_real_distribution<> distr(-2.0, 2.0);
        il::Array2D<double> m_pts_crd {3, n_m_p, 0.0};
        double rm;
        double r_tol = 0.02 * std::sqrt(hfp3d::pi / num_el / std::sqrt(3.0));
        il::int_t pt = 0;
        while (pt < n_m_p) {
            for (int j = 0; j < 3; ++j) {
                m_pts_crd(j, pt) = distr(gen);
            }
            rm = 1.0;
            for (il::int_t el = 0; el < num_el; ++el) {
                // Vertices' coordinates
                il::StaticArray2D<double, 3, 3> el_vert;
                for (int k = 0; k < 3; ++k) {
                    il::int_t n = mesh_data.mesh.conn(k, el);
                    double r2 = 0.0;
                    for (int j = 0; j < 3; ++j) {
                        el_vert(j, k) = mesh_data.mesh.nods(j, n);
                        r2 += std::pow(m_pts_crd(j, pt) - el_vert(j, k), 2);
                    }
                    r2 = std::sqrt(r2);
                    if (r2 < rm) {
                        rm = r2;
                    }
                }
            }
            if (rm >= r_tol) {
                if (echo) {
                    std::cout << "[          ] (" << m_pts_crd(0, pt)
                              << ", " << m_pts_crd(1, pt)
                              << ", " << m_pts_crd(2, pt)
                              << ")" << std::endl;
                }
                ++pt;
            }
        }

        // calculating stresses at m_pts_crd
        // (substituting the dd vector into stresses)
        il::Array2D<double> stress_m_pts(m_pts_crd.size(0), 6);
        stress_m_pts = hfp3d::make_3dbem_stress_f_s
                (mat_props.shear_m[0], mat_props.poiss_r[0],
                 mesh_data, load_data, num_param, m_pts_crd);

        double err = 0.0, diff_ij;
        for (il::int_t pt = 0; pt < n_m_p; ++pt) {
            il::StaticArray<double, 3> p_crd;
            for (int j = 0; j < 3; ++j) {
                p_crd[j] = m_pts_crd(j, pt);
            }
            il::StaticArray<double, 6> ref_stress;
            ref_stress = ref::get_stress_at_pt
                    (mat_props.shear_m[0],
                     mat_props.poiss_r[0],
                     1.0,
                     load_data.s_inf[2],
                     p_crd);
            for (int j = 0; j < 6; ++j) {
                diff_ij = std::fabs(stress_m_pts(pt, j) - ref_stress[j]);
                if (diff_ij > err) {
                    err = diff_ij;
                }
                if (diff_ij > tol && echo) {
                    std::cout << "[          ] pt. " << pt+1
                              << ", S_" << j+1
                              << ", err = " << diff_ij
                              << std::endl;
                }
            }
        }
        return err;
    }
}

// Reference solutions

namespace ref{

// analytical solution for normally loaded penny-shaped crack
// in infinite elastic solid medium

// DD on the surface
    il::StaticArray<double, 3> get_dd_at_pt
            (double shear_m,
             double poiss_r,
             double a,
             double p,
             il::StaticArray<double, 3> crd) {

        double scale = (4.0 * (1.0 - poiss_r) * a * p) / (shear_m * hfp3d::pi);

        double x2 = crd[0] * crd[0],
                y2 = crd[1] * crd[1];

        double r2 = (x2 + y2); // ignores the 3rd coordinate
        double a2 = a * a;

        il::StaticArray<double, 3> dd {0.0};
        if (r2 / a2 < 1.0) {
            dd[2] = scale * std::sqrt(1.0 - r2 / a2);
        }

        return dd;
    };

// Stress field
    il::StaticArray<double, 6> get_stress_at_pt
            (double shear_m,
             double poiss_r,
             double a,
             double p,
             il::StaticArray<double, 3> crd) {

        double scale = p / shear_m,
                fac = 2. / hfp3d::pi;

        double a2 = a * a;
        double x2 = crd[0] * crd[0],
                y2 = crd[1] * crd[1],
                xy = crd[0] * crd[1];
        double r = std::sqrt(x2 + y2),
                r2 = r * r,
                r3 = r * r2;
        double z = std::fabs(crd[2]),
                z2 = z * z;

        double L1 = 0.5 * (std::sqrt((a + r) * (a + r) + z2)
                           - std::sqrt((a - r) * (a - r) + z2));
        double L2 = 0.5 * (std::sqrt((a + r) * (a + r) + z2)
                           + std::sqrt((a - r) * (a - r) + z2));

        double g0 = std::asin(L1 / r);
        double g1 = L1 / L2,
                L14 = L1 * L1 * L1 * L1;
        double g2 = L2 * L2 - L1 * L1,
                g23 = g2 * g2 * g2;
        double g3 = std::sqrt(a * a - L1 * L1);
        double g4 = std::sqrt(L2 * L2 - a * a);
        double g5 = std::sqrt(r * r - L1 * L1);

        double Fr = 0.5 * r * g0 - 0.5 * g1 * g4;
        double Frr = 0.5 * g0 - 0.5 * a * (1. + g1 * g1) * g4 / g2;
        double Fxx = x2 / r2 * Frr + y2 / r3 * Fr;
        double Fyy = y2 / r2 * Frr + x2 / r3 * Fr;
        double Fxy = xy / r2 * Frr - xy / r3 * Fr;
        double Fzz = -g0 + a * g4 / g2;
        double Frz = -a * g1 * g3 / g2;
        double Frrz = a2 * z / (L2 * L2) /
                      g4 * (4. * r2 * g4 /(g23) - (a2 + z2) / (g2 * g2));
        double Fxxz = x2 / r2 * Frrz + y2 / r3 * Frz;
        double Fyyz = y2 / r2 * Frrz + x2 / r3 * Frz;
        double Fxyz = xy / r2 * Frrz - xy / r3 * Frz;
        double Frzz = a / L2 * g5 / (g23) * (a2 * (4. * L2 * L2 - 5. * r2) + L14);
        double Fxzz = crd[0] / r * Frzz;
        double Fyzz = crd[1] / r * Frzz;
        double Fzzz = g3 / (g23) * (L14 + a2 * (2.*a2 + 2.*z2 - 3.*r2));

        il::StaticArray<double, 6> sv;

        sv[0] = -scale * fac * (z * Fxxz + Fxx + 2. * poiss_r * Fyy);
        sv[1] = -scale * fac * (z * Fyyz + Fyy + 2. * poiss_r * Fxx);
        sv[2] = scale *(1. - fac * (z * Fzzz - Fzz));
        sv[3] = -scale * fac * (z * Fxyz + (1. - 2. * poiss_r) * Fxy);
        sv[4] = -scale * fac * (z * Fxzz);
        sv[5] = -scale * fac * (z * Fyzz);

        return sv;
    }

}

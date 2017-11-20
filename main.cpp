//
// This file is part of HFPx3D.
//
// Created by D. Nikolski.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

//todo: main.cpp will be used for testing the code parts under development

//#include <cstdio>
#include <iostream>
//#include <cassert>
//#include <sys/stat.h>
//#include <complex>
//#include <ittnotify.h>

#include <il/Timer.h>
#include <il/Toml.h>
//#include <il/io/toml/toml.h>
#include <il/Status.h>
#include <il/String.h>
#include <il/Array.h>
#include <il/Array2D.h>
//#include <il/StaticArray.h>
//#include <il/StaticArray2D.h>
#include <il/linear_algebra.h>
#include <il/linear_algebra/dense/factorization/LU.h>
#include <il/linear_algebra/dense/factorization/linearSolve.h>

#include "src/IO/config_file_io.h"
#include "src/IO/mesh_file_io.h"
#include "src/IO/data_file_io.h"
#include "src/Development/c_f_iteration.h"
#include "src/Development/cohesion_friction.h"
#include "src/Solvers/system_assembly.h"
#include "src/Core/model_parameters.h"
#include "src/Core/surface_mesh_utilities.h"
#include "src/Core/element_utilities.h"
#include "src/Core/tensor_utilities.h"


int main(int argc, char* argv[]) {

//todo: input config file (path) via argv
    // source files directory (containing main.cpp)
    std::string src_f = __FILE__;
    while (src_f.find("\\")!=std::string::npos) {
        src_f.replace(src_f.find("\\"),1,"/");
    }
    std::string src_dir = src_f.substr(0, src_f.rfind("/"));
    // std::string src_dir{"C:/Users/nikolski/ClionProjects/HFPx3D_VC"};
    // std::string src_dir{"/home/nikolski/Documents/HFPx3D"};
    // std::string src_dir{"/home/lecampio/Documents/HFPx3D"};
    // std::string c_f_name{src_dir + "/config.toml"};
    // std::string mesh_conn_fname{"Elems_pennymesh24el_32.npy"};
    // std::string nodes_crd_fname{"Nodes_pennymesh24el_32.npy"};

    std::string c_f_name = src_dir + "/config.toml";

    // initializing the mesh etc.
    hfp3d::Mesh_Data_T mesh_data;
    hfp3d::Properties_T mat_props;
    hfp3d::Load_T load_data;
    hfp3d::Num_Param_T num_param;
    hfp3d::Sim_Param_T sim_param;
    hfp3d::IO_param_T io_param;

// todo: move towards using il::Status type in output functions
    bool ok = true;

    il::Status status{};

    hfp3d::read_config(c_f_name,
                       il::io,
                       mesh_data, mat_props, load_data,
                       num_param, sim_param, io_param, status);

    il::String out_dir = il::toString(io_param.output_dir);
    il::String mf_name = il::toString(io_param.matr_f_name);
    il::String of_name = il::toString(io_param.out_f_name);
    il::String sf_name = il::toString(io_param.strs_f_name);

    // number of elements
    il::int_t num_elems = mesh_data.mesh.conn.size(1);

    // resetting the timer
    il::Timer timer{};
    std::cout << "Assembly";
    timer.start();

//    __itt_resume();

    // initializing the material ID array
    //mesh_data.mat_id = hfp3d::make_mat_id_triv(mesh_data.mesh, 2);

    // initializing the DoF handle
    mesh_data.dof_h_dd = hfp3d::make_dof_h_crack
            (mesh_data.mesh, 2, num_param.tip_type);

    // initializing the DD array
    //mesh_data.dd = il::Array2D<double> {num_elems * 6, 3, 0.0};

//    il::Array2D<il::int_t> ip(0, 7, 0);
//    hfp3d::Mesh_Data_T mdf = hfp3d::init_mesh_data_p_fault(mesh_data.mesh, 2, ip);

    // assembly of the algebraic system (matrix + RHS)
    hfp3d::SAE_T sae;
    // matrix
    sae.matrix = hfp3d::make_3dbem_matrix_s
            (mat_props.shear_m[0], mat_props.poiss_r[0],
             mesh_data.mesh,
             num_param,
             il::io, mesh_data.dof_h_dd);
    // number of DoF
    sae.n_dof = mesh_data.dof_h_dd.n_dof;
    // right-hand side
    hfp3d::add_s_inf_to_3dbem_rhs
            (mesh_data, load_data, il::io, sae);

    timer.stop();

    std::cout << ": " << timer.elapsed() << "s" << std::endl;
    std::cout << 18 * num_elems << " DoF Total" << std::endl;
    std::cout << 18 * num_elems - sae.n_dof << " Fixed DoF" << std::endl;

//    __itt_pause();

    // saving matrix to a .CSV file
// todo: add binary output
    if (sim_param.do_save_matrix) {
        if (io_param.out_f_format == "npy32") {
            //
        } else if (io_param.in_f_format == "npy64") {
            //
        } else { // treat as csv by default
            hfp3d::save_data_to_csv
                    (sae.matrix,
                     out_dir, mf_name,
                     il::io, ok);
        }
        if (ok) {
            std::cout << "Matrix saved to "
                      << io_param.output_dir << "/"
                      << io_param.matr_f_name << std::endl;
        }
        else {
            std::cout << "Cannot save the matrix" << std::endl;
            //status.abortOnError();
        }
    }

    std::cout << "Solution";
    timer.reset();
    timer.start();

    // solving the system
    il::Array<double> dd_v{sae.n_dof};

//    dd_v = il::linearSolve(sae.matrix, sae.rhs_v, il::io, status);
//    status.abortOnError();

    il::LU<il::Array2D<double>> lu_decomposition(sae.matrix, il::io, status);
    // if (!status.ok()) {
    //     // The matrix is singular to the machine precision.
    //     // You should deal with the error.
    // }
    status.abortOnError();
    // double cnd = lu_decomposition.conditionNumber(il::Norm::L2, );
    // std::cout << cnd << std::endl;
    dd_v = lu_decomposition.solve(sae.rhs_v);

    timer.stop();
    std::cout << ": " << timer.elapsed() << "s" << std::endl;

//    __itt_pause();

    hfp3d::write_dd_vector_to_md
            (dd_v, mesh_data.dof_h_dd,
             false, mesh_data.dof_h_pp,
             il::io, mesh_data);

    // saving the solution (nodes + DD) to a .CSV file
    if (sim_param.do_save_solution) {
        // the 2D array for nodal points' coordinates and DD - initialization
        il::Array2D<double> out_dd(6 * num_elems, 7);
        for (il::int_t j = 0; j < num_elems; ++j) {
            il::StaticArray2D<double, 3, 3> el_vert;
            for (il::int_t k = 0; k < 3; ++k) {
                il::int_t n = mesh_data.mesh.conn(k, j);
                for (il::int_t l = 0; l < 3; ++l) {
                    el_vert(l, k) = mesh_data.mesh.nods(l, n);
                }
            }
            // nodal point coordinates
            il::StaticArray<il::StaticArray<double, 3>, 6> el_np;
            el_np = hfp3d::el_cp_uniform(el_vert, 0.0);
            // el_np = hfp3d::el_cp_nonuniform(el_vert, v_wts, 0.0);
            for (il::int_t k = 0; k < 6; ++k) {
                il::int_t n = j * 6 + k;
                for (il::int_t l = 0; l < 3; ++l) {
                    //il::int_t dof = n * 3 + l;
//                il::int_t l_dof = k * 3 + l;
//                il::int_t g_dof = mesh_data.dof_h_dd.dof_h(j, l_dof);
                    out_dd(n, l) = el_np[k][l];
                    out_dd(n, l + 4) = mesh_data.dd(n, l);
//                if (g_dof != -1) {
//                    //out_dd(n, l) = dd_v[dof];
//                    out_dd(n, l + 4) = dd_v[g_dof];
//                } else {
//                    //out_dd(n, l) = 0.0;
//                    out_dd(n, l + 4) = 0.0;
//                }
                }
            }
        }
// todo: add binary output
        if (io_param.out_f_format == "npy32") {
//            il::save(out_dd, out_dir_name + of_name, il::io, status);
        } else if (io_param.in_f_format == "npy64") {
            //
        } else { // treat as csv by default
            hfp3d::save_data_to_csv
                    (out_dd,
                     out_dir, of_name,
                     il::io, ok);
        }
        if (ok) {
            std::cout << "Solution (nodes + DD) saved to "
                      << io_param.output_dir << "/" // .asCString()
                      << io_param.out_f_name << std::endl;
        }
        else {
            std::cout << "Cannot save the solution" << std::endl;
            //status.abortOnError();
        }
    }

// todo: add calculation of stresses (post-processing)
    if (sim_param.do_postprocess) {
        // define points to monitor stresses
        il::Array2D<double> m_pts_crd;

        // loading the mesh from files
        if (io_param.out_f_format == "csv") {
            m_pts_crd = hfp3d::load_crd_from_csv
                    (io_param.obs_p_dir, io_param.obs_p_f_name,
                     il::io, status);
        } else { // treat as numpy by default
            m_pts_crd = hfp3d::load_crd_from_numpy
                    (io_param.obs_p_dir, io_param.obs_p_f_name,
                     il::io, status);
        }
        //status.abortOnError();

        std::cout << "Postprocessing";
        timer.reset();
        timer.start();

        // calculate stresses at m_pts_crd
        il::Array2D<double> stress_m_pts(m_pts_crd.size(0), 6);
        stress_m_pts = hfp3d::make_3dbem_stress_f_s
            (mat_props.shear_m[0], mat_props.poiss_r[0],
             mesh_data, load_data, num_param, m_pts_crd);

        timer.stop();
        std::cout << ": " << timer.elapsed() << "s" << std::endl;

        // saving stresses to the file
// todo: add binary output
        if (io_param.out_f_format == "npy32") {
            //
        } else if (io_param.in_f_format == "npy64") {
            //
        } else { // treat as csv by default
            hfp3d::save_data_to_csv
                    (stress_m_pts,
                     out_dir, sf_name,
                     il::io, ok);
        }
        if (ok) {
            std::cout << "Stresses at given points saved to "
                      << io_param.output_dir << "/"
                      << io_param.strs_f_name << std::endl;
        }
        else {
            std::cout << "Cannot save the stresses" << std::endl;
            //status.abortOnError();
        }
    }

    return 0;
}
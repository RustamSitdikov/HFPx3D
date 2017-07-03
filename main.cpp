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
#include <il/toml.h>
#include <il/String.h>
#include <il/Array.h>
#include <il/Array2D.h>
//#include <il/StaticArray.h>
//#include <il/StaticArray2D.h>
#include <il/linear_algebra/dense/factorization/LU.h>
#include <il/linear_algebra.h>
//#include <il/linear_algebra/dense/factorization/linearSolve.h>

#include "src/mesh_file_io.h"
#include "src/system_assembly.h"
#include "src/surface_mesh_utilities.h"
#include "src/element_utilities.h"
#include "src/tensor_utilities.h"


int main() {

    // source files directory (containing main.cpp)
    std::string src_f = __FILE__;
    while (src_f.find("\\")!=std::string::npos) {
        src_f.replace(src_f.find("\\"),1,"/");
    }
    std::string src_dir = src_f.substr(0, src_f.rfind("/"));
    // std::string src_dir{"C:/Users/nikolski/ClionProjects/HFPx3D_VC"};
    // std::string src_dir{"/home/nikolski/Documents/HFPx3D"};
    // std::string src_dir{"/home/lecampio/Documents/HFPx3D"};

    std::string default_f_name = src_dir + "/config.toml";
    il::String config_f_name(default_f_name.c_str());
    // il::String config_f_name =
    // "C:/Users/nikolski/ClionProjects/HFPx3D_VC/config.toml";

    std::string default_input_dir{src_dir + "/Mesh_Files/"};
    // std::string mesh_conn_fname{"Elems_pennymesh24el_32.npy"};
    // std::string nodes_crd_fname{"Nodes_pennymesh24el_32.npy"};

    std::string default_output_dir{src_dir + "/Test_Output/"};

    //il::String src_f_name(src_dir.c_str());
    il::String in_dir_name((src_dir + "/").c_str());
    il::String out_dir_name((src_dir + "/").c_str());
    il::String m_c_f_name;
    il::String m_n_f_name;
    il::String i_f_format;
    il::int_t array_origin;

    il::Status status{};

    // reading configuration & parameters
    auto config =
            il::load<il::MapArray<il::String, il::Dynamic>>
                    (config_f_name, il::io, status);

    status.abortOnError();

    // reading input (triangulation) files
    il::int_t pos = config.search("input_directory");
    if (config.found(pos) && config.value(pos).isString()) {
        in_dir_name.append(config.value(pos).asString());
    } else {
        in_dir_name = il::String(default_input_dir.c_str());
    }
    in_dir_name.append("/");
    pos = config.search("mesh_conn_fname");
    if (config.found(pos) && config.value(pos).isString()) {
        m_c_f_name.append(config.value(pos).asString());
    } else {
        std::cout << "Can't find the input file" << std::endl;
        abort();
    }
    pos = config.search("nodes_crd_fname");
    if (config.found(pos) && config.value(pos).isString()) {
        m_n_f_name.append(config.value(pos).asString());
    } else {
        std::cout << "Can't find the input file" << std::endl;
        abort();
    }
    pos = config.search("input_format");
    if (config.found(pos) && config.value(pos).isString()) {
        i_f_format = config.value(pos).asString();
    } else {
        i_f_format = "npy32";
    }
    pos = config.search("array_origin");
    if (config.found(pos) && config.value(pos).isInteger()) {
        array_origin = config.value(pos).toInteger();
    } else {
        array_origin = 0;
    }

    // reading output target
    pos = config.search("output_directory");
    if (config.found(pos) && config.value(pos).isString()) {
        out_dir_name.append(config.value(pos).asString());
    } else {
        out_dir_name = il::String(default_output_dir.c_str());
    }
    il::String mf_name = "test_assembly";
    il::String of_name = "test_solution";
    // il::String sf_name = "test_stresses";
    pos = config.search("output_signature");
    if (config.found(pos) && config.value(pos).isString()) {
        mf_name.append(config.value(pos).asString());
        of_name.append(config.value(pos).asString());
    }
    mf_name.append(".csv");
    of_name.append(".csv");
    // std::string mf_name{"test_assembly_24_ele.csv"};
    // std::string of_name{"test_solution_24_ele.csv"};

    // material properties (default)
    //hfp3d::Properties_T properties_list;
    //properties_list.mu = il::Array<double>{1, 1.0};
    //properties_list.nu = il::Array<double>{1, 0.35};
    double mu = 1.0, nu = 0.35;

    // numerical simulation parameters
    hfp3d::Num_Param_T n_par; // default: beta = 0.125; tip_type = 1; DD in global
    n_par.beta = 0.125; n_par.tip_type = 1; n_par.is_dd_local = false;

    hfp3d::Load_T load;
    // stress at infinity (in-situ stress)
    load.s_inf = il::StaticArray<double, 6> {0.0};
    // load.s_inf[2] = 1.0; load.s_inf[4] = 1.0;
    pos = config.search("S_xx");
    if (config.found(pos) && config.value(pos).isDouble()) {
        load.s_inf[0] = (config.value(pos).toDouble());
    }
    pos = config.search("S_yy");
    if (config.found(pos) && config.value(pos).isDouble()) {
        load.s_inf[1] = (config.value(pos).toDouble());
    }
    pos = config.search("S_zz");
    if (config.found(pos) && config.value(pos).isDouble()) {
        load.s_inf[2] = (config.value(pos).toDouble());
    }
    pos = config.search("S_xy");
    if (config.found(pos) && config.value(pos).isDouble()) {
        load.s_inf[3] = (config.value(pos).toDouble());
    }
    pos = config.search("S_xz");
    if (config.found(pos) && config.value(pos).isDouble()) {
        load.s_inf[4] = (config.value(pos).toDouble());
    }
    pos = config.search("S_yz");
    if (config.found(pos) && config.value(pos).isDouble()) {
        load.s_inf[5] = (config.value(pos).toDouble());
    }

    // mesh
    hfp3d::Mesh_Data_T mesh_data;

    // loading the mesh from files
    if (i_f_format == "csv") {
        hfp3d::load_mesh_from_csv
                (in_dir_name, m_c_f_name, m_n_f_name, array_origin,
                 il::io, mesh_data.mesh);
    } else if (i_f_format == "npy64") {
        hfp3d::load_mesh_from_numpy_64
                (in_dir_name, m_c_f_name, m_n_f_name, array_origin,
                 il::io, mesh_data.mesh);
    } else { // treat as 32-bit numpy by default
        hfp3d::load_mesh_from_numpy_32
                (in_dir_name, m_c_f_name, m_n_f_name, array_origin,
                 il::io, mesh_data.mesh);
    }

    //hfp3d::load_mesh_from_numpy_32
    // (input_dir, mesh_conn_fname, nodes_crd_fname,
    // true, il::io, mesh_data.mesh);

    // number of elements
    il::int_t num_elems = mesh_data.mesh.conn.size(1);

    // resetting the timer
    il::Timer timer{};
    timer.start();

//    __itt_resume();

    // initializing the material ID array
    //mesh_data.mat_id = hfp3d::make_mat_id_triv(mesh_data.mesh, 2);

    // initializing the DoF handle
    mesh_data.dof_h_dd = hfp3d::make_dof_h_crack
            (mesh_data.mesh, 2, n_par.tip_type);

    // number of DoF
    il::int_t num_dof = mesh_data.dof_h_dd.n_dof;

    // initializing the DD array
    //mesh_data.dd = il::Array2D<double> {num_elems * 6, 3, 0.0};

//    il::Array2D<il::int_t> ip(0, 7, 0);
//    hfp3d::Mesh_Data_T mdf = hfp3d::init_mesh_data_p_fault(mesh_data.mesh, 2, ip);

    // assembly of the algebraic system (matyrix + RHS)
    hfp3d::SAE_T sae;
    sae.matrix = hfp3d::make_3dbem_matrix_s
            (mu, nu, mesh_data.mesh, n_par, il::io, mesh_data.dof_h_dd);
    hfp3d::add_s_inf_to_3dbem_rhs
            (mesh_data, load, il::io, sae);

//    __itt_pause();
    timer.stop();

    std::cout << "Assembly: " << timer.elapsed() << "s" << std::endl;
    std::cout << 18 * num_elems << " DoF Total" << std::endl;
    std::cout << 18 * num_elems - num_dof << " Fixed DoF" << std::endl;

    timer.reset();
    timer.start();

    // solving the system
    il::Array<double> dd_v;

    il::LU<il::Array2D<double>> lu_decomposition(sae.matrix, il::io, status);
    // if (!status.ok()) {
    //     // The matrix is singular to the machine precision. You should deal with
    //     // the error.
    // }
    status.abortOnError();
    // double cnd = lu_decomposition.conditionNumber(il::Norm::L2, );
    // std::cout << cnd << std::endl;
    dd_v = lu_decomposition.solve(sae.rhs_v);

//    dd_v = il::linearSolve(sae.matrix, sae.rhs_v, il::io, status);
//    status.abortOnError();

    timer.stop();
    std::cout << "Solution: " << timer.elapsed() << "s" << std::endl;

    hfp3d::write_dd_vector_to_md
            (dd_v, mesh_data.dof_h_dd,
             false, mesh_data.dof_h_pp,
             il::io, mesh_data);

// todo: calculation of stresses (post-processing)
//    // define points to monitor stresses
//    il::Array2D<double> m_pts_crd;
//
//    // calculate stresses at m_pts_crd
//    il::Array2D<double> stress_m_pts(m_pts_crd.size(0), 6);
//    stress_m_pts = hfp3d::make_3dbem_stress_f_s
//            (mu, nu, mesh_data, n_par, m_pts_crd);

// todo: move towards using il::Status type
    bool ok = true;
    // saving matrix to a .CSV file
//    hfp3d::save_data_to_csv(sae.matrix, out_dir_name, mf_name, il::io, ok);
//    if (ok) {
//        std::cout << "Matrix saved to "
//                  << out_dir_name.asCString() << "/"
//                  << mf_name.asCString() << std::endl;
//    }

    // the 2D array for nodal points' coordinates and DD - initialization
    il::Array2D<double> out_dd(6 * num_elems, 7);

    // saving the solution (nodes + DD) to a .CSV file
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
                il::int_t l_dof = k * 3 + l;
                il::int_t g_dof = mesh_data.dof_h_dd.dof_h(j, l_dof);
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

    hfp3d::save_data_to_csv(out_dd, out_dir_name, of_name, il::io, ok);
    if (ok) {
        std::cout << "Solution (nodes + DD) saved to "
                  << out_dir_name.asCString() << "/"
                  << of_name.asCString() << std::endl;
    }

/*
    std::string npy_of_name{"/test_solution_24_ele.npy"};
    // il::Status status{};
    il::save(out_dd, output_dir + npy_of_name, il::io, status);
    status.abortOnError();
*/

    return 0;
}
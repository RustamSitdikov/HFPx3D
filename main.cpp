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
#include <il/io/toml/toml.h>
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
#include "src/model_utilities.h"
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
    il::String d_name;
    il::String in_dir_name((src_dir + "/").c_str());
    il::String m_c_f_name;
    il::String m_n_f_name;
    il::String in_f_format;
    il::String obs_dir_name((src_dir + "/").c_str());
    il::String o_p_f_name;
    il::String obs_f_format;
    il::String out_dir_name((src_dir + "/").c_str());
    il::String out_f_format;
    il::int_t array_origin;
    bool do_save_matrix = false;
    bool do_save_solution = false;
    bool do_postprocess = false;
// todo: move towards using il::Status type in output functions
    bool ok = true;
    il::Status status{};

    // reading configuration & parameters
    auto config =
            il::load<il::MapArray<il::String, il::Dynamic>>
                    (config_f_name, il::io, status);

    status.abortOnError();

    // reading input (triangulation) files
    il::int_t pos = config.search("mesh_input_directory");
    if (config.found(pos) && config.value(pos).isString()) {
        d_name = config.value(pos).asString();
        in_dir_name.append(d_name);
        if (!d_name.hasSuffix("/"))
            in_dir_name.append("/");
    } else {
        in_dir_name = il::String(default_input_dir.c_str());
    }
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
    pos = config.search("mesh_input_format");
    if (config.found(pos) && config.value(pos).isString()) {
        in_f_format = config.value(pos).asString();
    } else {
        in_f_format = "npy32";
    }
    pos = config.search("array_origin");
    if (config.found(pos) && config.value(pos).isInteger()) {
        array_origin = config.value(pos).toInteger();
    } else {
        array_origin = 0;
    }

    // reading observation points
    pos = config.search("observ_crd_directory");
    if (config.found(pos) && config.value(pos).isString()) {
        d_name = config.value(pos).asString();
        obs_dir_name.append(d_name);
        if (!d_name.hasSuffix("/"))
            obs_dir_name.append("/");
    } else {
        obs_dir_name = il::String(default_input_dir.c_str());
    }
    pos = config.search("observ_crd_fname");
    if (config.found(pos) && config.value(pos).isString()) {
        o_p_f_name.append(config.value(pos).asString());
    } else {
        std::cout << "Can't find the input file" << std::endl;
        abort();
    }
    pos = config.search("observ_input_format");
    if (config.found(pos) && config.value(pos).isString()) {
        obs_f_format = config.value(pos).asString();
    } else {
        obs_f_format = "npy32";
    }

    // reading output permissions
    pos = config.search("do_save_matrix");
    if (config.found(pos) && config.value(pos).isBool()) {
        do_save_matrix = config.value(pos).toBool();
    }
    pos = config.search("do_save_solution");
    if (config.found(pos) && config.value(pos).isBool()) {
        do_save_solution = config.value(pos).toBool();
    }
    pos = config.search("do_postprocess");
    if (config.found(pos) && config.value(pos).isBool()) {
        do_postprocess = config.value(pos).toBool();
    }

    // reading output target
    pos = config.search("output_directory");
    if (config.found(pos) && config.value(pos).isString()) {
        d_name = config.value(pos).asString();
        out_dir_name.append(d_name);
        if (!d_name.hasSuffix("/"))
            out_dir_name.append("/");
    } else {
        out_dir_name = il::String(default_output_dir.c_str());
    }
    il::String mf_name = "test_assembly";
    il::String of_name = "test_solution";
    il::String sf_name = "test_stresses";
    pos = config.search("output_signature");
    if (config.found(pos) && config.value(pos).isString()) {
        mf_name.append(config.value(pos).asString());
        of_name.append(config.value(pos).asString());
        sf_name.append(config.value(pos).asString());
    }
    pos = config.search("output_format");
    if (config.found(pos) && config.value(pos).isString()) {
        out_f_format = config.value(pos).asString();
    } else {
        out_f_format = "csv";
    }
//todo: add binary output
    mf_name.append(".csv");
    of_name.append(".csv");
    sf_name.append(".csv");

    // reading material properties (default)
    hfp3d::Properties_T solid_properties;
    solid_properties.n_solid = 1; // default
    pos = config.search("number_solids");
    if (config.found(pos) && config.value(pos).isInteger()) {
        solid_properties.n_solid = config.value(pos).toInteger();
        if (solid_properties.n_solid == 0) solid_properties.n_solid = 1;
    }
    solid_properties.mu = il::Array<double>{solid_properties.n_solid};
    solid_properties.nu = il::Array<double>{solid_properties.n_solid};
    solid_properties.mu[0] = 1.0; // default
    solid_properties.nu[0] = 0.35; // default
    pos = config.search("solid");
    if (config.found(pos) && config.value(pos).isMapArray()) {
        const il::MapArray<il::String, il::Dynamic> &solid =
                config.value(pos).asMapArray();

        il::int_t j = solid.search("Poisson_ratio");
        if (solid.found(j) && solid.value(j).isDouble()) {
            solid_properties.nu[0] = solid.value(j).toDouble();
        }

        j = solid.search("Shear_modulus");
        if (solid.found(j) && solid.value(j).isDouble()) {
            solid_properties.mu[0] = solid.value(j).toDouble();
        } else {
            j = solid.search("Young_modulus");
            if (solid.found(j) && solid.value(j).isDouble()) {
                double ym = (solid.value(j).toDouble());
                solid_properties.mu[0] = ym / 2.0 /
                                         (1.0 + solid_properties.nu[0]);
            }
        }
    } else {
//todo: read moduli as Array<double>
        pos = config.search("Poisson_ratio");
        if (config.found(pos) && config.value(pos).isDouble()) {
            solid_properties.nu[0] = (config.value(pos).toDouble());
        }
        pos = config.search("Shear_modulus");
        if (config.found(pos)) {
            if(config.value(pos).isDouble()) {
                solid_properties.mu[0] = (config.value(pos).toDouble());
            }
        } else {
            pos = config.search("Young_modulus");
            if(config.found(pos) && config.value(pos).isDouble()) {
                double ym = (config.value(pos).toDouble());
                solid_properties.mu[0] = ym / 2.0 /
                                         (1.0 + solid_properties.nu[0]);
            }
        }
    }
    if (solid_properties.n_solid > 1) {
//    for (il::int_t k = 1; k < solid_properties.n_solid; ++k) {
//
//    }
    }

    // reading numerical simulation parameters
    hfp3d::Num_Param_T n_par; // default: beta = 0.125; tip_type = 1; DD in global
    n_par.beta = 0.125; n_par.tip_type = 1; n_par.is_dd_local = false;
    pos = config.search("cp_offset");
    if (config.found(pos) && config.value(pos).isDouble()) {
        n_par.beta = (config.value(pos).toDouble());
    }
    pos = config.search("fix_tip");
    if (config.found(pos) && config.value(pos).isInteger()) {
        n_par.tip_type = (int)(config.value(pos).toInteger());
    }
    pos = config.search("is_dd_local");
    if (config.found(pos) && config.value(pos).isBool()) {
        n_par.is_dd_local = (config.value(pos).toBool());
    }

    // reading load parameters
    hfp3d::Load_T load;
    // stress at infinity (in-situ stress)
    load.s_inf = il::StaticArray<double, 6> {0.0};
    // load.s_inf[2] = 1.0; load.s_inf[4] = 1.0;
    pos = config.search("S_xx");
    if (config.found(pos) && config.value(pos).isDouble()) {
        load.s_inf[0] = -(config.value(pos).toDouble());
    }
    pos = config.search("S_yy");
    if (config.found(pos) && config.value(pos).isDouble()) {
        load.s_inf[1] = -(config.value(pos).toDouble());
    }
    pos = config.search("S_zz");
    if (config.found(pos) && config.value(pos).isDouble()) {
        load.s_inf[2] = -(config.value(pos).toDouble());
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

    // initializing the mesh
    hfp3d::Mesh_Data_T mesh_data;

    // loading the mesh from files
    if (in_f_format == "csv") {
        hfp3d::load_mesh_from_csv
                (in_dir_name, m_c_f_name, m_n_f_name, array_origin,
                 il::io, mesh_data.mesh);
    } else if (in_f_format == "npy64") {
        hfp3d::load_mesh_from_numpy_64
                (in_dir_name, m_c_f_name, m_n_f_name, array_origin,
                 il::io, mesh_data.mesh);
    } else { // treat as 32-bit numpy by default
        hfp3d::load_mesh_from_numpy_32
                (in_dir_name, m_c_f_name, m_n_f_name, array_origin,
                 il::io, mesh_data.mesh);
    }

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
            (mesh_data.mesh, 2, n_par.tip_type);

    // initializing the DD array
    //mesh_data.dd = il::Array2D<double> {num_elems * 6, 3, 0.0};

//    il::Array2D<il::int_t> ip(0, 7, 0);
//    hfp3d::Mesh_Data_T mdf = hfp3d::init_mesh_data_p_fault(mesh_data.mesh, 2, ip);

    // assembly of the algebraic system (matrix + RHS)
    hfp3d::SAE_T sae;
    // matrix
    sae.matrix = hfp3d::make_3dbem_matrix_s
            (solid_properties.mu[0], solid_properties.nu[0],
             mesh_data.mesh,
             n_par, il::io, mesh_data.dof_h_dd);
    // number of DoF
    sae.n_dof = mesh_data.dof_h_dd.n_dof;
    // right-hand side
    hfp3d::add_s_inf_to_3dbem_rhs
            (mesh_data, load, il::io, sae);

    timer.stop();

    std::cout << ": " << timer.elapsed() << "s" << std::endl;
    std::cout << 18 * num_elems << " DoF Total" << std::endl;
    std::cout << 18 * num_elems - sae.n_dof << " Fixed DoF" << std::endl;

//    __itt_pause();

    // saving matrix to a .CSV file
    if (do_save_matrix) {
        // todo: add binary output
        if (out_f_format == "npy32") {
            //
        } else if (in_f_format == "npy64") {
            //
        } else { // treat as csv by default
            hfp3d::save_data_to_csv
                    (sae.matrix, out_dir_name, mf_name, il::io, ok);
        }
        if (ok) {
            std::cout << "Matrix saved to "
                      << out_dir_name.asCString()
                      << mf_name.asCString() << std::endl;
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
    std::cout << ": " << timer.elapsed() << "s" << std::endl;

//    __itt_pause();

    hfp3d::write_dd_vector_to_md
            (dd_v, mesh_data.dof_h_dd,
             false, mesh_data.dof_h_pp,
             il::io, mesh_data);

    // saving the solution (nodes + DD) to a .CSV file
    if (do_save_solution) {
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
        if (out_f_format == "npy32") {
//            il::save(out_dd, out_dir_name + of_name, il::io, status);
        } else if (in_f_format == "npy64") {
            //
        } else { // treat as csv by default
            hfp3d::save_data_to_csv
                    (out_dd, out_dir_name, of_name, il::io, ok);
        }
        if (ok) {
            std::cout << "Solution (nodes + DD) saved to "
                      << out_dir_name.asCString()
                      << of_name.asCString() << std::endl;
        }
        else {
            std::cout << "Cannot save the solution" << std::endl;
            //status.abortOnError();
        }
    }

// todo: add calculation of stresses (post-processing)
    if (do_postprocess) {
        // define points to monitor stresses
        il::Array2D<double> m_pts_crd;

        // loading the mesh from files
        if (out_f_format == "csv") {
            m_pts_crd = hfp3d::load_crd_from_csv
                    (obs_dir_name, o_p_f_name,
                     il::io, status);
        } else { // treat as numpy by default
            m_pts_crd = hfp3d::load_crd_from_numpy
                    (obs_dir_name, o_p_f_name,
                     il::io, status);
        }
        //status.abortOnError();

        std::cout << "Postprocessing";
        timer.reset();
        timer.start();

        // calculate stresses at m_pts_crd
        il::Array2D<double> stress_m_pts(m_pts_crd.size(0), 6);
        stress_m_pts = hfp3d::make_3dbem_stress_f_s
            (solid_properties.mu[0], solid_properties.nu[0],
             mesh_data, load, n_par, m_pts_crd);

        timer.stop();
        std::cout << ": " << timer.elapsed() << "s" << std::endl;

        // saving stresses to the file
        // todo: add binary output
        if (out_f_format == "npy32") {
            //
        } else if (in_f_format == "npy64") {
            //
        } else { // treat as csv by default
            hfp3d::save_data_to_csv
                    (stress_m_pts, out_dir_name, sf_name, il::io, ok);
        }
        if (ok) {
            std::cout << "Stresses at given points saved to "
                      << out_dir_name.asCString()
                      << sf_name.asCString() << std::endl;
        }
        else {
            std::cout << "Cannot save the stresses" << std::endl;
            //status.abortOnError();
        }
    }

    return 0;
}
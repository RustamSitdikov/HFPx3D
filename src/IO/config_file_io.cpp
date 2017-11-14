//
// This file is part of HFPx3D.
//
// Created by nikolski on 9/5/2017.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//

#include <il/Toml.h>
//#include <il/io/toml/toml.h>
#include <il/Status.h>
#include <il/String.h>
#include <il/Map.h>
#include <il/Array.h>
#include <il/Array2D.h>
#include "src/IO/config_file_io.h"
#include "src/IO/mesh_file_io.h"
#include "src/IO/data_file_io.h"
#include "src/Core/surface_mesh_utilities.h"
#include "src/Core/model_parameters.h"

namespace hfp3d{

    //
    IO_param_T set_default_paths (std::string src_path) {

        IO_param_T default_p;

        if (src_path.size() > 0) {
            default_p.src_dir = src_path;
        } else {
            // source files directory (containing main.cpp)
            std::string src_f = __FILE__;
            while (src_f.find("\\")!=std::string::npos) {
                src_f.replace(src_f.find("\\"),1,"/");
            }
            default_p.src_dir = src_f.substr(0, src_f.rfind("/"));
        }


        default_p.cf_name = default_p.src_dir + "/config.toml";

        default_p.input_dir = default_p.src_dir + "/Mesh_Files/";

        default_p.output_dir = default_p.src_dir + "/Test_Output/";

        default_p.matr_f_name = "test_assembly";
        default_p.out_f_name = "test_solution";
        default_p.strs_f_name = "test_stresses";

        return default_p;
    }

    //
    void read_config(
            std::string f_path,
            il::io_t,
            Mesh_Data_T &mesh_data,
            Properties_T &props,
            Load_T &load,
            Num_Param_T &num_param,
            Sim_Param_T &sim_param,
            IO_param_T &io_param,
            il::Status &status
    ) {
        std::string src_p = f_path;
        while (src_p.find("\\")!=std::string::npos) {
            src_p.replace(src_p.find("\\"),1,"/");
        }
        std::string src_dir = src_p.substr(0, src_p.rfind("/"));
        io_param = set_default_paths(src_dir);

        //il::String src_f_name(src_dir.c_str());
        il::String d_name;
        il::String in_dir_name = il::toString(src_dir); //.c_str()
        il::String m_c_f_name;
        il::String m_n_f_name;
        il::String in_f_format;
        il::String obs_dir_name = il::toString(src_dir);
        il::String o_p_f_name;
        il::String obs_f_format;
        il::String out_dir_name = il::toString(src_dir);
        il::String out_f_format;
        il::String mf_name = il::toString(io_param.matr_f_name);
        il::String of_name = il::toString(io_param.out_f_name);
        il::String sf_name = il::toString(io_param.strs_f_name);

        il::int_t array_origin;

// todo: move towards using il::Status type in output functions
        bool ok = true;

        // reading configuration & parameters
        il::String config_f_name;
        if (src_dir.size() == 0 && f_path.size() > 0) {
            config_f_name = il::toString(io_param.src_dir + src_p);
        } else if (src_dir.size() == f_path.size()) {
            config_f_name = il::toString(io_param.cf_name);
        } else {
            config_f_name = il::toString(src_p);
        }
        auto config =
                il::load<il::MapArray<il::String, il::Dynamic>>
                        (config_f_name, il::io, status);
        status.abortOnError();

        // reading input (triangulation) files
        il::int_t pos = config.search("mesh_input_directory");
        if (config.found(pos) && config.value(pos).isString()) {
            d_name = config.value(pos).asString();
            if (d_name.begin() != "/")
                in_dir_name.append('/');
            in_dir_name.append(d_name);
            if (d_name.end() != "/") // (!d_name.hasSuffix("/"))
                in_dir_name.append('/');
            io_param.input_dir = in_dir_name.asCString();
        } else {
            // use default directory
            in_dir_name = il::toString(io_param.input_dir);
        }
        pos = config.search("mesh_conn_fname");
        if (config.found(pos) && config.value(pos).isString()) {
            m_c_f_name.append(config.value(pos).asString());
            io_param.conn_f_name = m_c_f_name.asCString();
        } else {
            std::cout << "Can't find the input file" << std::endl;
            abort();
        }
        pos = config.search("nodes_crd_fname");
        if (config.found(pos) && config.value(pos).isString()) {
            m_n_f_name.append(config.value(pos).asString());
            io_param.node_f_name = m_n_f_name.asCString();
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
        io_param.in_f_format = in_f_format.asCString();

        // array origin - C++ (0) or Matlab/Mathematica (1) defaults
        pos = config.search("array_origin");
        if (config.found(pos) && config.value(pos).isInteger()) {
            array_origin = config.value(pos).toInteger();
        } else {
            array_origin = 0;
        }

        // reading output permissions
        sim_param.do_save_matrix = false;
        sim_param.do_save_solution = false;
        sim_param.do_postprocess = false;
        pos = config.search("do_save_matrix");
        if (config.found(pos) && config.value(pos).isBool()) {
            sim_param.do_save_matrix = config.value(pos).toBool();
        }
        pos = config.search("do_save_solution");
        if (config.found(pos) && config.value(pos).isBool()) {
            sim_param.do_save_solution = config.value(pos).toBool();
        }
        pos = config.search("do_postprocess");
        if (config.found(pos) && config.value(pos).isBool()) {
            sim_param.do_postprocess = config.value(pos).toBool();
        }

        // reading material properties (default)
        props.n_solid = 1; // default
        pos = config.search("number_solids");
        if (config.found(pos) && config.value(pos).isInteger()) {
            props.n_solid = config.value(pos).toInteger();
            if (props.n_solid == 0) props.n_solid = 1;
        }
        props.shear_m = il::Array<double>{props.n_solid};
        props.poiss_r = il::Array<double>{props.n_solid};
        props.shear_m[0] = 1.0; // default
        props.poiss_r[0] = 0.0; // default // 0.35
        pos = config.search("solid");
        if (config.found(pos) && config.value(pos).isMapArray()) {
            const il::MapArray<il::String, il::Dynamic> &solid =
                    config.value(pos).asMapArray();

            il::int_t j = solid.search("Poisson_ratio");
            if (solid.found(j) && solid.value(j).isDouble()) {
                props.poiss_r[0] = solid.value(j).toDouble();
            }

            j = solid.search("Shear_modulus");
            if (solid.found(j) && solid.value(j).isDouble()) {
                props.shear_m[0] = solid.value(j).toDouble();
            } else {
                j = solid.search("Young_modulus");
                if (solid.found(j) && solid.value(j).isDouble()) {
                    double ym = (solid.value(j).toDouble());
                    props.shear_m[0] = ym / 2.0 /
                                             (1.0 + props.poiss_r[0]);
                }
            }
        } else {
            pos = config.search("Poisson_ratio");
            if (config.found(pos) && config.value(pos).isDouble()) {
                props.poiss_r[0] = (config.value(pos).toDouble());
            }
            pos = config.search("Shear_modulus");
            if (config.found(pos)) {
                if(config.value(pos).isDouble()) {
                    props.shear_m[0] = (config.value(pos).toDouble());
                }
            } else {
                pos = config.search("Young_modulus");
                if(config.found(pos) && config.value(pos).isDouble()) {
                    double ym = (config.value(pos).toDouble());
                    props.shear_m[0] = ym / 2.0 /
                                             (1.0 + props.poiss_r[0]);
                }
            }
        }
//todo: read moduli as Array<double>
        if (props.n_solid > 1) {
//    for (il::int_t k = 1; k < solid_properties.n_solid; ++k) {
//
//    }
        }

        // reading simulation parameters

        // reading numerical scheme parameters
        // default: beta = 0.125; tip_type = 1; DD in global
        num_param.beta = 0.125; num_param.tip_type = 1;
        num_param.is_dd_local = false;
        pos = config.search("cp_offset");
        if (config.found(pos) && config.value(pos).isDouble()) {
            num_param.beta = (config.value(pos).toDouble());
        }
        pos = config.search("fix_tip");
        if (config.found(pos) && config.value(pos).isInteger()) {
            num_param.tip_type = (int)(config.value(pos).toInteger());
        }
        pos = config.search("is_dd_local");
        if (config.found(pos) && config.value(pos).isBool()) {
            num_param.is_dd_local = (config.value(pos).toBool());
        }

        // reading load parameters
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
//todo: add injection

        // loading the mesh from files
        if (in_f_format == "csv") {
            load_mesh_from_csv
                    (in_dir_name, m_c_f_name, m_n_f_name, array_origin,
                     il::io, mesh_data.mesh);
        } else if (in_f_format == "npy64") {
            load_mesh_from_numpy_64
                    (in_dir_name, m_c_f_name, m_n_f_name, array_origin,
                     il::io, mesh_data.mesh);
        } else { // treat as 32-bit numpy by default
            load_mesh_from_numpy_32
                    (in_dir_name, m_c_f_name, m_n_f_name, array_origin,
                     il::io, mesh_data.mesh);
        }

        // reading observation points
        pos = config.search("observ_crd_directory");
        if (config.found(pos) && config.value(pos).isString()) {
            d_name = config.value(pos).asString();
            if (d_name.begin() != "/")
                obs_dir_name.append('/');
            obs_dir_name.append(d_name);
            if (d_name.end() != "/") // (!d_name.hasSuffix("/"))
                obs_dir_name.append('/');
        } else {
            obs_dir_name = il::toString(io_param.input_dir);
        }
        io_param.obs_p_dir = obs_dir_name.asCString();
        pos = config.search("observ_crd_fname");
        if (config.found(pos) && config.value(pos).isString()) {
            o_p_f_name.append(config.value(pos).asString());
            io_param.obs_p_f_name = o_p_f_name.asCString();
        } else {
            std::cout << "Can't find the input file" << std::endl;
            abort();
        }
        pos = config.search("observ_input_format");
        if (config.found(pos) && config.value(pos).isString()) {
            obs_f_format = config.value(pos).asString();
        } else {
            obs_f_format = il::toString("npy32");
        }
        io_param.obs_p_f_format = obs_f_format.asCString();

        // reading output target
        pos = config.search("output_directory");
        if (config.found(pos) && config.value(pos).isString()) {
            d_name = config.value(pos).asString();
            if (d_name.begin() != "/")
                out_dir_name.append('/');
            out_dir_name.append(d_name);
            if (d_name.end() != "/") // (!d_name.hasSuffix("/"))
                out_dir_name.append('/');
        } else {
            out_dir_name = il::toString(io_param.output_dir);
        }
        mf_name = il::toString(io_param.matr_f_name);
        of_name = il::toString(io_param.out_f_name);
        sf_name = il::toString(io_param.strs_f_name);
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
            out_f_format = il::toString("csv");
        }
        io_param.out_f_format = out_f_format.asCString();
        if (true) {
            mf_name.append(il::toString(".csv"));
            of_name.append(il::toString(".csv"));
            sf_name.append(il::toString(".csv"));
        }
//todo: add binary output
        io_param.matr_f_name = mf_name.asCString();
        io_param.out_f_name = of_name.asCString();
        io_param.strs_f_name = sf_name.asCString();

    }
}

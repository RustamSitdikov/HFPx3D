//
// This file is part of HFPx3D.
//
// Created by nikolski on 9/5/2017.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//

#ifndef HFPX3D_CONFIG_FILE_IO_H
#define HFPX3D_CONFIG_FILE_IO_H

#include <il/String.h>
#include <il/Status.h>
#include "src/Core/model_parameters.h"
#include "src/Core/surface_mesh_utilities.h"

namespace hfp3d{

    struct IO_param_T {
        std::string src_dir;
        std::string cf_name; // .toml config file
        std::string input_dir; // mesh files directory
        std::string conn_f_name; // mesh connectivity file name
        std::string node_f_name; // mesh node coordinates file name
        std::string in_f_format; // choose from "npy32", "npy64", "csv"
        std::string obs_p_dir; // observation points directory
        std::string obs_p_f_name; // observation points file name
        std::string obs_p_f_format; // choose from "npy32", "npy64", "csv"
        std::string output_dir; // output (solution) directory
        std::string out_f_name; // output (solution, mesh DD) file name
        std::string out_f_format; // choose from "npy32", "npy64", "csv"
        std::string matr_f_name; // matrix output file name
        std::string strs_f_name; // stress (postprocessing) file name
    };

    //
    IO_param_T set_default_paths (std::string src_path);

    //
    void read_config (
            std::string f_path,
            il::io_t,
            Mesh_Data_T &mesh_data,
            Properties_T &props,
            Load_T &load,
            Num_Param_T &num_param,
            Sim_Param_T &sim_param,
            IO_param_T &io_param,
            il::Status &status
    );
}

#endif //HFPX3D_CONFIG_FILE_IO_H

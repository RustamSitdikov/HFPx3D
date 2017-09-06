//
// This file is part of HFPx3D.
//
// Created by D. Nikolski on 1/27/2017.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#ifndef INC_HFPX3D_MESH_FILE_IO_H
#define INC_HFPX3D_MESH_FILE_IO_H

#include <il/Array.h>
#include <il/Array2D.h>
#include <il/StaticArray.h>
#include <il/StaticArray2D.h>
#include <il/Status.h>
#include <il/String.h>
#include <complex>
#include <cstdio>
#include <fstream>
#include <iostream>
#include "src/Core/surface_mesh_utilities.h"

namespace hfp3d {

// loading data (point coordinates, il::Array2D<double>)

    il::Array2D<double> load_crd_from_numpy
            (const std::string &src_dir,
             const std::string &crd_f_name,
             il::io_t, il::Status &status);

    // overload for il::String
    il::Array2D<double> load_crd_from_numpy
            (const il::String &src_dir,
             const il::String &crd_f_name,
             il::io_t, il::Status &status);

    il::Array2D<double> load_crd_from_csv
            (const std::string &src_dir,
             const std::string &crd_f_name,
             il::io_t, il::Status &status);

    // overload for il::String
    il::Array2D<double> load_crd_from_csv
            (const il::String &src_dir,
             const il::String &crd_f_name,
             il::io_t, il::Status &status);

// loading data (mesh, two of il::Array2D<double>)

    void load_mesh_from_numpy_32
        (const std::string &src_dir,
         const std::string &conn_f_name,
         const std::string &node_f_name,
         // bool is_matlab,
         il::int_t origin,
         il::io_t, Mesh_Geom_T &mesh);

    // overload for il::String
    void load_mesh_from_numpy_32
            (const il::String &src_dir,
             const il::String &conn_f_name,
             const il::String &node_f_name,
             // bool is_matlab,
             il::int_t origin,
             il::io_t, Mesh_Geom_T &mesh);

    void load_mesh_from_numpy_64
        (const std::string &src_dir,
         const std::string &conn_f_name,
         const std::string &node_f_name,
         // bool is_matlab,
         il::int_t origin,
         il::io_t, Mesh_Geom_T &mesh);

    // overload for il::String
    void load_mesh_from_numpy_64
            (const il::String &src_dir,
             const il::String &conn_f_name,
             const il::String &node_f_name,
             // bool is_matlab,
             il::int_t origin,
             il::io_t, Mesh_Geom_T &mesh);

    void load_mesh_from_csv
        (const std::string &src_dir,
         const std::string &conn_f_name,
         const std::string &node_f_name,
         // bool is_matlab,
         il::int_t origin,
         il::io_t, Mesh_Geom_T &mesh);

    // overload for il::String
    void load_mesh_from_csv
            (const il::String &src_dir,
             const il::String &conn_f_name,
             const il::String &node_f_name,
             // bool is_matlab,
             il::int_t origin,
             il::io_t, Mesh_Geom_T &mesh);

// saving data (il::Array<double> or il::Array2D<double>)

}

#endif  // INC_HFPX3D_MESH_FILE_IO_H

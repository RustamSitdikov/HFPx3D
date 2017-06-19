//
// This file is part of HFPx3D.
//
// Created by D. Nikolski on 1/27/2017.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//

#include <cstdio>
#include <il/Array.h>
#include <il/Array2D.h>
#include <il/StaticArray.h>
#include <il/StaticArray2D.h>
#include <il/io/numpy/numpy.h>
#include "mesh_file_io.h"

namespace hfp3d {

    void load_mesh_from_numpy_32
            (const std::string &src_dir,
             const std::string &conn_f_name,
             const std::string &node_f_name,
             bool is_matlab,
             il::io_t, Mesh_Geom_T &mesh) {
// This function reads the mesh connectivity matrix (3*N_elements)
// (32-bit integer)
// and node coordinates matrix (3*N_nodes) from numpy binary files
        il::Status status{};
        std::string f_path = src_dir + conn_f_name;

        il::Array2D<int> mesh_conn_tmp = il::load<il::Array2D<int>>
                (f_path, il::io, status);
        status.abortOnError();
        mesh.conn.resize(mesh_conn_tmp.size(0), mesh_conn_tmp.size(1));
        for (il::int_t k = 0; k < mesh.conn.size(1); ++k) {
            for (il::int_t j = 0; j < mesh.conn.size(0); ++j) {
                mesh.conn(j, k) = mesh_conn_tmp(j, k);
            }
        }

        mesh.nods = il::load<il::Array2D<double>>
                (src_dir + node_f_name, il::io, status);
        status.abortOnError();

        if (is_matlab) {
            // conversion from Matlab (array numbering starting with 1)
            // to C++ standard (array numbering starting with 0)
            il::int_t n_r = (mesh.conn.size(0) >= 3) ? 3 : mesh.conn.size(0);
            for (il::int_t k = 0; k < mesh.conn.size(1); ++k) {
                for (il::int_t j = 0; j < n_r; ++j) {
                    mesh.conn(j, k) -=1;
                }
            }
        }
    }

    void load_mesh_from_numpy_64
            (const std::string &src_dir,
             const std::string &conn_f_name,
             const std::string &node_f_name,
             bool is_matlab,
             il::io_t, Mesh_Geom_T &mesh) {
// This function reads the mesh connectivity matrix (3*N_elements)
// (64-bit integer)
// and node coordinates matrix (3*N_nodes) from numpy binary files
        il::Status status{};
        std::string f_path = src_dir + conn_f_name;

        mesh.conn = il::load<il::Array2D<il::int_t>>
                (f_path, il::io, status);
        status.abortOnError();

        mesh.nods = il::load<il::Array2D<double>>
                (src_dir + node_f_name, il::io, status);
        status.abortOnError();

        if (is_matlab) {
            // conversion from Matlab (array numbering starting with 1)
            // to C++ standard (array numbering starting with 0)
            il::int_t n_r = (mesh.conn.size(0) >= 3) ? 3 : mesh.conn.size(0);
            for (il::int_t k = 0; k < mesh.conn.size(1); ++k) {
                for (il::int_t j = 0; j < n_r; ++j) {
                    mesh.conn(j, k) -=1;
                }
            }
        }
    }

}

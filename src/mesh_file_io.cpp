//
// This file is part of HFPx3D.
//
// Created by D. Nikolski on 1/27/2017.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//

#include <cstdio>
#include <iostream>
#include <cstring>
#include <fstream>
#include <sstream>
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
             // bool is_matlab,
             il::int_t origin,
             il::io_t, Mesh_Geom_T &mesh) {
// This function reads the mesh connectivity matrix (3*N_elements)
// (32-bit integer)
// and node coordinates matrix (3*N_nodes) from numpy binary files
        il::Status status{};
        std::string c_f_path = src_dir + conn_f_name;
        std::string n_f_path = src_dir + node_f_name;

        il::Array2D<int> mesh_conn_tmp = il::load<il::Array2D<int>>
                (c_f_path, il::io, status);
        status.abortOnError();
        mesh.conn.resize(mesh_conn_tmp.size(0), mesh_conn_tmp.size(1));
        for (il::int_t k = 0; k < mesh.conn.size(1); ++k) {
            for (il::int_t j = 0; j < mesh.conn.size(0); ++j) {
                mesh.conn(j, k) = mesh_conn_tmp(j, k);
            }
        }

        mesh.nods = il::load<il::Array2D<double>>
                (n_f_path, il::io, status);
        status.abortOnError();

        if (origin !=0) {
        // if (is_matlab) {
            // conversion from Matlab (array numbering starting with 1)
            // to C++ standard (array numbering starting with 0)
            il::int_t n_r = (mesh.conn.size(0) >= 3) ? 3 : mesh.conn.size(0);
            for (il::int_t k = 0; k < mesh.conn.size(1); ++k) {
                for (il::int_t j = 0; j < n_r; ++j) {
                    mesh.conn(j, k) -= origin;
                }
            }
        }
    }

    // overloaded for il::String
    void load_mesh_from_numpy_32
            (const il::String &src_dir,
             const il::String &conn_f_name,
             const il::String &node_f_name,
             il::int_t origin,
             // bool is_matlab,
             il::io_t, Mesh_Geom_T &mesh) {
// This function reads the mesh connectivity matrix (3*N_elements)
// (32-bit integer)
// and node coordinates matrix (3*N_nodes) from numpy binary files
        il::Status status{};
        il::String c_f_path(src_dir); c_f_path.append("/");
        c_f_path.append(conn_f_name);
        il::String n_f_path(src_dir); n_f_path.append("/");
        n_f_path.append(node_f_name);

        il::Array2D<int> mesh_conn_tmp = il::load<il::Array2D<int>>
                (c_f_path, il::io, status);
        status.abortOnError();
        mesh.conn.resize(mesh_conn_tmp.size(0), mesh_conn_tmp.size(1));
        for (il::int_t k = 0; k < mesh.conn.size(1); ++k) {
            for (il::int_t j = 0; j < mesh.conn.size(0); ++j) {
                mesh.conn(j, k) = mesh_conn_tmp(j, k);
            }
        }

        mesh.nods = il::load<il::Array2D<double>>
                (n_f_path, il::io, status);
        status.abortOnError();

        if (origin !=0) {
        // if (is_matlab) {
            // conversion from Matlab (array numbering starting with 1)
            // to C++ standard (array numbering starting with 0)
            il::int_t n_r = (mesh.conn.size(0) >= 3) ? 3 : mesh.conn.size(0);
            for (il::int_t k = 0; k < mesh.conn.size(1); ++k) {
                for (il::int_t j = 0; j < n_r; ++j) {
                    mesh.conn(j, k) -= origin;
                }
            }
        }
    }

    void load_mesh_from_numpy_64
            (const std::string &src_dir,
             const std::string &conn_f_name,
             const std::string &node_f_name,
             // bool is_matlab,
             il::int_t origin,
             il::io_t, Mesh_Geom_T &mesh) {
// This function reads the mesh connectivity matrix (3*N_elements)
// (64-bit integer)
// and node coordinates matrix (3*N_nodes) from numpy binary files
        il::Status status{};
        std::string c_f_path = src_dir + conn_f_name;
        std::string n_f_path = src_dir + node_f_name;

        mesh.conn = il::load<il::Array2D<il::int_t>>
                (c_f_path, il::io, status);
        status.abortOnError();

        mesh.nods = il::load<il::Array2D<double>>
                (n_f_path, il::io, status);
        status.abortOnError();

        if (origin !=0) {
        // if (is_matlab) {
            // conversion from Matlab (array numbering starting with 1)
            // to C++ standard (array numbering starting with 0)
            il::int_t n_r = (mesh.conn.size(0) >= 3) ? 3 : mesh.conn.size(0);
            for (il::int_t k = 0; k < mesh.conn.size(1); ++k) {
                for (il::int_t j = 0; j < n_r; ++j) {
                    mesh.conn(j, k) -= origin;
                }
            }
        }
    }

    // overloaded for il::String
    void load_mesh_from_numpy_64
            (const il::String &src_dir,
             const il::String &conn_f_name,
             const il::String &node_f_name,
             // bool is_matlab,
             il::int_t origin,
             il::io_t, Mesh_Geom_T &mesh) {
// This function reads the mesh connectivity matrix (3*N_elements)
// (64-bit integer)
// and node coordinates matrix (3*N_nodes) from numpy binary files
        il::Status status{};
        il::String c_f_path(src_dir); c_f_path.append(conn_f_name);
        il::String n_f_path(src_dir); n_f_path.append(node_f_name);

        mesh.conn = il::load<il::Array2D<il::int_t>>
                (c_f_path, il::io, status);
        status.abortOnError();

        mesh.nods = il::load<il::Array2D<double>>
                (n_f_path, il::io, status);
        status.abortOnError();

        if (origin !=0) {
        // if (is_matlab) {
            // conversion from Matlab (array numbering starting with 1)
            // to C++ standard (array numbering starting with 0)
            il::int_t n_r = (mesh.conn.size(0) >= 3) ? 3 : mesh.conn.size(0);
            for (il::int_t k = 0; k < mesh.conn.size(1); ++k) {
                for (il::int_t j = 0; j < n_r; ++j) {
                    mesh.conn(j, k) -= origin;
                }
            }
        }
    }

    void load_mesh_from_csv
            (const std::string &src_dir,
             const std::string &conn_f_name,
             const std::string &node_f_name,
             // bool is_matlab,
             il::int_t origin,
             il::io_t, Mesh_Geom_T &mesh) {
// This function reads the mesh connectivity matrix (3*N_elements)
// and node coordinates matrix (3*N_nodes) from CSV files.
//
        std::string c_f_path = src_dir + conn_f_name;
        std::string n_f_path = src_dir + node_f_name;

        std::string line = "", subline = "";
        char delim = ',';

        std::ifstream cf; // connectivity file stream
        cf.open(c_f_path.c_str());
        if(cf.is_open()) {
            // counting the array size
            il::int_t n_row = 0, n_col = 0;
            while(!cf.eof()) {
                std::getline(cf, line);
                if (line.length() == 0) {
                    std::cout << "Empty row" << std::endl;
                } else {
                    std::stringstream linestream(line);
                    ++n_row;
                    il::int_t n_col_t = 0;
                    do {
                        subline = "";
                        std::getline(linestream, subline, delim);
                        if (subline.length() > 0){
                            ++n_col_t;
                        }
                    } while(subline.length() > 0);
                    if (n_col_t != n_col && n_col > 0) {
                        std::cout << "Row size is not constant" << std::endl;
                    }
                    else n_col = n_col_t;
                }
            }
            // resizing the output array (mesh connectivity)
            mesh.conn.resize(n_col, n_row);
            // importing the array
            for (il::int_t k = 0; k < mesh.conn.size(1); ++k) {
                std::getline(cf, line);
                std::stringstream linestream(line);
                for (il::int_t j = 0; j < mesh.conn.size(0); ++j) {
                    subline = "";
                    std::getline(linestream, subline, delim);
                    if (subline.length() > 0){
                        il::int_t value = -1;
                        std::stringstream sublinestream(subline);
                        sublinestream >> value;
                        mesh.conn(j, k) = value;
                    }
                }
            }
            if (origin !=0) {
            // if (is_matlab) {
                // conversion from Matlab (array numbering starting with 1)
                // to C++ standard (array numbering starting with 0)
                il::int_t n_r = (mesh.conn.size(0) >= 3) ? 3 : mesh.conn.size(0);
                for (il::int_t k = 0; k < mesh.conn.size(1); ++k) {
                    for (il::int_t j = 0; j < n_r; ++j) {
                        mesh.conn(j, k) -= origin;
                    }
                }
            }
            cf.close();
        }
        else {
            std::cout << "Can't open the file" << std::endl;
        }

        std::ifstream nf; // node coordinates file stream
        nf.open(n_f_path.c_str());
        if(nf.is_open()) {
            // counting the array size
            il::int_t n_row = 0, n_col = 0;
            while(!nf.eof()) {
                std::getline(nf, line);
                if (line.length() == 0) {
                    std::cout << "Empty row" << std::endl;
                } else {
                    std::stringstream linestream(line);
                    ++n_row;
                    il::int_t n_col_t = 0;
                    do {
                        subline = "";
                        std::getline(linestream, subline, delim);
                        if (subline.length() > 0){
                            ++n_col_t;
                        }
                    } while(subline.length() > 0);
                    if (n_col_t != n_col && n_col > 0) {
                        std::cout << "Row size is not constant" << std::endl;
                    }
                    else n_col = n_col_t;
                }
            }
            // resizing the output array (node coordinates)
            mesh.nods.resize(n_col, n_row);
            // importing the array
            for (il::int_t k = 0; k < mesh.nods.size(1); ++k) {
                std::getline(nf, line);
                std::stringstream linestream(line);
                for (il::int_t j = 0; j < mesh.nods.size(0); ++j) {
                    subline = "";
                    std::getline(linestream, subline, delim);
                    if (subline.length() > 0){
                        double value = -1.0;
                        std::stringstream sublinestream(subline);
                        sublinestream >> value;
                        mesh.nods(j, k) = value;
                    }
                }
            }
            nf.close();
        }
        else {
            std::cout << "Can't open the file" << std::endl;
        }
    }

    // overloaded for il::String
    void load_mesh_from_csv
            (const il::String &src_dir,
             const il::String &conn_f_name,
             const il::String &node_f_name,
             // bool is_matlab,
             il::int_t origin,
             il::io_t, Mesh_Geom_T &mesh) {
// This function reads the mesh connectivity matrix (3*N_elements)
// and node coordinates matrix (3*N_nodes) from CSV files.
//
        il::String c_f_path(src_dir); c_f_path.append(conn_f_name);
        il::String n_f_path(src_dir); n_f_path.append(node_f_name);

        std::string line = "", subline = "";
        char delim = ',';

        std::ifstream cf; // connectivity file stream
        cf.open(c_f_path.asCString());
        if(cf.is_open()) {
            // counting the array size
            il::int_t n_row = 0, n_col = 0;
            while(!cf.eof()) {
                std::getline(cf, line);
                if (line.length() == 0) {
                    std::cout << "Empty row" << std::endl;
                } else {
                    std::stringstream linestream(line);
                    ++n_row;
                    il::int_t n_col_t = 0;
                    do {
                        subline = "";
                        std::getline(linestream, subline, delim);
                        if (subline.length() > 0){
                            ++n_col_t;
                        }
                    } while(subline.length() > 0);
                    if (n_col_t != n_col && n_col > 0) {
                        std::cout << "Row size is not constant" << std::endl;
                    }
                    else n_col = n_col_t;
                }
            }
            // resizing the output array (mesh connectivity)
            mesh.conn.resize(n_col, n_row);
            // importing the array
            for (il::int_t k = 0; k < mesh.conn.size(1); ++k) {
                std::getline(cf, line);
                std::stringstream linestream(line);
                for (il::int_t j = 0; j < mesh.conn.size(0); ++j) {
                    subline = "";
                    std::getline(linestream, subline, delim);
                    if (subline.length() > 0){
                        il::int_t value = -1;
                        std::stringstream sublinestream(subline);
                        sublinestream >> value;
                        mesh.conn(j, k) = value;
                    }
                }
            }
            if (origin !=0) {
            // if (is_matlab) {
                // conversion from Matlab (array numbering starting with 1)
                // to C++ standard (array numbering starting with 0)
                il::int_t n_r = (mesh.conn.size(0) >= 3) ? 3 : mesh.conn.size(0);
                for (il::int_t k = 0; k < mesh.conn.size(1); ++k) {
                    for (il::int_t j = 0; j < n_r; ++j) {
                        mesh.conn(j, k) -= origin;
                    }
                }
            }
            cf.close();
        }
        else {
            std::cout << "Can't open the file" << std::endl;
        }

        std::ifstream nf; // node coordinate file stream
        nf.open(n_f_path.asCString());
        if(nf.is_open()) {
            // counting the array size
            il::int_t n_row = 0, n_col = 0;
            while(!nf.eof()) {
                std::getline(nf, line);
                if (line.length() == 0) {
                    std::cout << "Empty row" << std::endl;
                } else {
                    std::stringstream linestream(line);
                    ++n_row;
                    il::int_t n_col_t = 0;
                    do {
                        subline = "";
                        std::getline(linestream, subline, delim);
                        if (subline.length() > 0){
                            ++n_col_t;
                        }
                    } while(subline.length() > 0);
                    if (n_col_t != n_col && n_col > 0) {
                        std::cout << "Row size is not constant" << std::endl;
                    }
                    else n_col = n_col_t;
                }
            }
            // resizing the output array (node coordinates)
            mesh.nods.resize(n_col, n_row);
            // importing the array
            for (il::int_t k = 0; k < mesh.nods.size(1); ++k) {
                std::getline(nf, line);
                std::stringstream linestream(line);
                for (il::int_t j = 0; j < mesh.nods.size(0); ++j) {
                    subline = "";
                    std::getline(linestream, subline, delim);
                    if (subline.length() > 0){
                        double value = -1.0;
                        std::stringstream sublinestream(subline);
                        sublinestream >> value;
                        mesh.nods(j, k) = value;
                    }
                }
            }
            nf.close();
        }
        else {
            std::cout << "Can't open the file" << std::endl;
        }
    }



}

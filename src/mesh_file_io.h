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
#include <complex>
#include <cstdio>
#include <fstream>
#include <iostream>
#include "surface_mesh_utilities.h"

namespace hfp3d {

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

    template <typename T>
    void save_data_to_csv
        (const il::Array<T> &vector,
         const std::string &trg_dir,
         const std::string &of_name) {
        std::string f_path = trg_dir + of_name;
        const char *format;
        if (sizeof(T) == sizeof(std::complex<double>)) {
            format = "%.16g%+.16gi\n";
        } else {
            format = "%.16g\n";
        }
        FILE *of = std::fopen(f_path.c_str(), "w");
        for (int j = 0; j < vector.size(); ++j) {
            std::fprintf(of, format, vector[j]);
        }
        std::fclose(of);
    }

    template <typename T, il::int_t m>
    void save_data_to_csv
        (const il::StaticArray<T, m> &vector,
         const std::string &trg_dir,
         const std::string &of_name) {
        std::string f_path = trg_dir + of_name;
        const char *format;
        if (sizeof(T) == sizeof(std::complex<double>)) {
            format = "%.16g%+.16gi\n";
        } else {
            format = "%.16g\n";
        }
        FILE *of = std::fopen(f_path.c_str(), "w");
        for (int j = 0; j < vector.size(); ++j) {
            T out = vector[j];
            std::fprintf(of, format, out);
        }
        std::fclose(of);
    }

    // overload for il::String
    template <typename T, il::int_t m>
    void save_data_to_csv
        (const il::StaticArray<T, m> &vector,
         const il::String &trg_dir,
         const il::String &of_name) {
        il::String f_path = trg_dir;
        f_path.append(of_name);
        const char *f_p_c = f_path.asCString();
        const char *format;
        if (sizeof(T) == sizeof(std::complex<double>)) {
            format = "%.16g%+.16gi\n";
        } else {
            format = "%.16g\n";
        }
        FILE *of = std::fopen(f_p_c, "w");
        for (int j = 0; j < vector.size(); ++j) {
            T out = vector[j];
            std::fprintf(of, format, out);
        }
        std::fclose(of);
    }

    template <typename T>
    void save_data_to_csv
            (const il::Array2D<T> &matrix,
             const std::string &trg_dir,
             const std::string &of_name) {
        std::string f_path = trg_dir + of_name;
        const char *format;
        if (sizeof(T) == sizeof(std::complex<double>)) {
            format = "%.16g%+.16gi";
        } else {
            format = "%.16g";
        }
        FILE *of = std::fopen(f_path.c_str(), "w");
        for (int j = 0; j < matrix.size(0); ++j) {
            for (int k = 0; k < matrix.size(1); ++k) {
                T out = matrix(j, k);
                std::fprintf(of, format, out);
                if (k < matrix.size(1) - 1) std::fprintf(of, ",");
            }
            std::fprintf(of, "\n");
        }
        std::fclose(of);
    }

// overload for il::String
    template <typename T>
    void save_data_to_csv
            (const il::Array2D<T> &matrix,
             const il::String &trg_dir,
             const il::String &of_name,
             il::io_t, bool &ok) {
    //todo: it will use bool output instead of il::Status for now
        il::String f_path = trg_dir;
        f_path.append(of_name);
        const char *f_p_c = f_path.asCString();
//        const char *format;
//        if (sizeof(T) == sizeof(std::complex<double>)) {
//            format = "%.16g%+.16gi";
//        } else {
//            format = "%.16g";
//        }
        std::ofstream ofs(f_p_c);
        if (ofs.is_open()) {
            for (int j = 0; j < matrix.size(0); ++j) {
                for (int k = 0; k < matrix.size(1); ++k) {
                    T out = matrix(j, k);
                    ofs << out;
                    if (k < matrix.size(1) - 1) {
                        ofs << ",";
                    }
                }
                ofs << std::endl;
            }
            ofs.close();
            ok = true;
        } else {
            std::cout << "Can't open the file ("
                      << f_p_c << ")"<< std::endl;
            ok = false;
        }
    }

    template <class T, il::int_t m, il::int_t n>
    void save_data_to_csv
            (const il::StaticArray2D<T, m, n> &matrix,
             const std::string &trg_dir,
             const std::string &of_name) {
        std::string f_path = trg_dir + of_name;
        const char *format;
        if (sizeof(T) == sizeof(std::complex<double>)) {
            format = "%.16g%+.16gi";
        } else {
            format = "%.16g";
        }
        FILE *of = std::fopen(f_path.c_str(), "w");
        for (int j = 0; j < matrix.size(0); ++j) {
            for (int k = 0; k < matrix.size(1); ++k) {
                T out = matrix(j, k);
                std::fprintf(of, format, out);
                if (k < matrix.size(1) - 1) std::fprintf(of, ",");
            }
            std::fprintf(of, "\n");
        }
        std::fclose(of);
    }

    // overload for il::String
    template <class T, il::int_t m, il::int_t n>
    void save_data_to_csv
            (const il::StaticArray2D<T, m, n> &matrix,
             const il::String &trg_dir,
             const il::String &of_name) {
        il::String f_path = trg_dir;
        f_path.append(of_name);
        const char *f_p_c = f_path.asCString();
        const char *format;
        if (sizeof(T) == sizeof(std::complex<double>)) {
            format = "%.16g%+.16gi";
        } else {
            format = "%.16g";
        }
        FILE *of = std::fopen(f_p_c, "w");
        for (int j = 0; j < matrix.size(0); ++j) {
            for (int k = 0; k < matrix.size(1); ++k) {
                T out = matrix(j, k);
                std::fprintf(of, format, out);
                if (k < matrix.size(1) - 1) std::fprintf(of, ",");
            }
            std::fprintf(of, "\n");
        }
        std::fclose(of);
    }
}

#endif  // INC_HFPX3D_MESH_FILE_IO_H

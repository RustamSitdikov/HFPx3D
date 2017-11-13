//
// This file is part of HFPx3D.
//
// Created by nikolski on 11/13/2017.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//

#include <il/Array2D.h>
#include <il/Status.h>
#include <il/String.h>
#include <fstream>
#include <iostream>

namespace hfp3d {

// loading data

    il::Array2D<double> load_data_from_csv
            (const std::string &src_dir,
             const std::string &crd_f_name,
             il::io_t, il::Status &status) {
        //
        std::string f_path = src_dir + crd_f_name;
        std::string line = "", subline = "";
        char delim = ',';
        il::Array2D<double> data_table;
        std::ifstream nf; // file stream
        nf.open(f_path.c_str());
        if (nf.is_open()) {
            // counting the array size
            il::int_t n_row = 0, n_col = 0;
            while (!nf.eof()) {
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
                        if (subline.length() > 0) {
                            ++n_col_t;
                        }
                    } while (subline.length() > 0);
                    if (n_col_t != n_col && n_col > 0) {
                        std::cout << "Row size is not constant" << std::endl;
                    } else n_col = n_col_t;
                }
            }
            nf.clear();
            nf.seekg(0);
            // resizing the output array
            data_table.resize(n_row, n_col);
            // importing the array
            for (il::int_t k = 0; k < data_table.size(0); ++k) {
                std::getline(nf, line);
                std::stringstream linestream(line);
                for (il::int_t j = 0; j < data_table.size(1); ++j) {
                    subline = "";
                    std::getline(linestream, subline, delim);
                    if (subline.length() > 0) {
                        double value = -1.0;
                        std::stringstream sublinestream(subline);
                        sublinestream >> value;
                        data_table(k, j) = value;
                    }
                }
            }
            nf.close();
        } else {
            std::cout << "Can't open the file" << std::endl;
        }
        return data_table;
    }

//    // overload for il::String
//    il::Array2D<double> load_data_from_csv
//            (const il::String &src_dir,
//             const il::String &crd_f_name,
//             il::io_t, il::Status &status) {
//        //
//        il::String f_path(src_dir);
//        if (f_path.end() != "/") // (!f_path.hasSuffix("/"))
//            f_path.append('/');
//        f_path.append(crd_f_name);
//        std::string line = "", subline = "";
//        char delim = ',';
//        il::Array2D<double> data_table;
//        std::ifstream nf; // file stream
//        nf.open(f_path.asCString());
//        if(nf.is_open()) {
//            // counting the array size
//            il::int_t n_row = 0, n_col = 0;
//            while(!nf.eof()) {
//                std::getline(nf, line);
//                if (line.length() == 0) {
//                    std::cout << "Empty row" << std::endl;
//                } else {
//                    std::stringstream linestream(line);
//                    ++n_row;
//                    il::int_t n_col_t = 0;
//                    do {
//                        subline = "";
//                        std::getline(linestream, subline, delim);
//                        if (subline.length() > 0){
//                            ++n_col_t;
//                        }
//                    } while(subline.length() > 0);
//                    if (n_col_t != n_col && n_col > 0) {
//                        std::cout << "Row size is not constant" << std::endl;
//                    }
//                    else n_col = n_col_t;
//                }
//            }
//            nf.clear();
//            nf.seekg(0);
//            // resizing the output array
//            data_table.resize(n_row, n_col);
//            // importing the array
//            for (il::int_t k = 0; k < data_table.size(0); ++k) {
//                std::getline(nf, line);
//                std::stringstream linestream(line);
//                for (il::int_t j = 0; j < data_table.size(1); ++j) {
//                    subline = "";
//                    std::getline(linestream, subline, delim);
//                    if (subline.length() > 0){
//                        double value = -1.0;
//                        std::stringstream sublinestream(subline);
//                        sublinestream >> value;
//                        data_table(k, j) = value;
//                    }
//                }
//            }
//            nf.close();
//        }
//        else {
//            std::cout << "Can't open the file" << std::endl;
//        }
//        return data_table;
//    }

}
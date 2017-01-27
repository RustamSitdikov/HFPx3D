//
// This file is part of 3d_bem.
//
// Created by D. Nikolski on 1/27/2017.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//

#include <il/Array2D.h>
#include <il/io/numpy.h>

namespace hfp3d {
    il::StaticArray<il::int_t, 2> load_mesh_from_numpy
            (std::string src_dir,
             std::string conn_f_name,
             std::string node_f_name,
             bool is_matlab,
             il::io_t, il::Array2D<il::int_t> &mesh_conn,
             il::Array2D<double> &nodes_crd) {
// This function reads the mesh connectivity matrix (3*N_elements)
// and node coordinates matrix (3*N_nodes) from numpy binary files
        il::Status status{};
        mesh_conn = il::load<il::Array2D<il::int_t>>
                (src_dir + conn_f_name, il::io, status);
        status.abort_on_error();

        nodes_crd = il::load<il::Array2D<double>>
                (src_dir + node_f_name, il::io, status);
        status.abort_on_error();

        il::StaticArray<il::int_t, 2> nums;

        nums[0] = mesh_conn.size(1);
        nums[1] = nodes_crd.size(1);

        if (is_matlab) {
            // conversion from Matlab (array numbering starting with 1)
            // to C++ standard (array numbering starting with 0)
            for (il::int_t k = 0; k < mesh_conn.size(1); ++k) {
                for (il::int_t j = 0; j < mesh_conn.size(0); ++j) {
                    mesh_conn(j, k) -=1;
                }
            }
        }
        return nums;
    }

    void save_matrix_to_csv(
            il::Array2D<double> matrix,
            std::string trg_dir,
            std::string of_name) {
        std::string path = trg_dir + of_name;
        FILE* of=std::fopen(path.c_str(),"w");
        for (int j=0; j < matrix.size(0); ++j){
            for (int k=0; k < matrix.size(1); ++k){
                std::fprintf(of, "%21.16g", matrix(j,k));
                if (k < matrix.size(1)-1) std::fprintf(of, ",");
            }
            std::fprintf(of, "\n");
        }
        std::fclose(of);
    }

}

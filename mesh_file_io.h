//
// This file is part of 3d_bem.
//
// Created by nikolski on 1/27/2017.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//

#ifndef INC_3D_BEM_LOAD_MESH_FROM_FILE_H
#define INC_3D_BEM_LOAD_MESH_FROM_FILE_H

#include <il/Array2D.h>

namespace hfp3d {

    il::StaticArray<il::int_t, 2> load_mesh_from_numpy
            (std::string src_dir,
             std::string conn_f_name,
             std::string node_f_name,
             bool is_matlab,
             il::io_t, il::Array2D <il::int_t> &mesh_conn,
             il::Array2D<double> &nodes_crd);

    void save_matrix_to_csv(
            il::Array2D<double> matrix,
            std::string trg_dir,
            std::string of_name);

}

#endif //INC_3D_BEM_LOAD_MESH_FROM_FILE_H

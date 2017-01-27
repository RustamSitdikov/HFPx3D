//
// This file is part of 3d_bem.
//
// Created by D. Nikolski on 1/10/2017.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//

#ifndef INC_3D_BEM_ELE_BASE_H
#define INC_3D_BEM_ELE_BASE_H

namespace hfp3d {

    struct HZ {
        double h;
        std::complex<double> z;
    };

    // Element's local coordinate system manipulations

    il::StaticArray2D<double, 3, 3> make_el_r_tensor
            (const il::StaticArray2D<double, 3, 3> &el_vert);

    il::StaticArray<std::complex<double>, 3> make_el_tau_crd
            (const il::StaticArray2D<double, 3, 3> &el_vert,
             const il::StaticArray2D<double, 3, 3> &r_tensor);

    HZ make_el_pt_hz
            (const il::StaticArray2D<double, 3, 3> &el_vert,
             const il::StaticArray<double, 3> &X0,
             const il::StaticArray2D<double, 3, 3> &r_tensor);

    il::StaticArray2D<std::complex<double>, 2, 2> make_el_tau_2_mc
            (const il::StaticArray2D<double, 3, 3> &el_vert,
             const il::StaticArray2D<double, 3, 3> &r_tensor);

    // Element's basis (shape) functions

    il::StaticArray2D<std::complex<double>, 6, 6> make_el_sfm_uniform
            (const il::StaticArray2D<double, 3, 3> &el_vert,
             il::io_t, il::StaticArray2D<double, 3, 3> &r_tensor);

    il::StaticArray2D<std::complex<double>, 6, 6> make_el_sfm_nonuniform
            (const il::StaticArray2D<double, 3, 3> &el_vert,
             const il::StaticArray<double, 3> &vertex_wts,
             il::io_t, il::StaticArray2D<double, 3, 3> &r_tensor);

    //il::StaticArray2D<std::complex<double>, 6, 6> make_el_sfm_beta
    // (const il::StaticArray2D<double, 3, 3> &el_vert,
    // const il::StaticArray<double, 3> &vertex_wts,
    // double beta,
    // il::io_t, il::StaticArray2D<double, 3, 3> &r_tensor);

    il::StaticArray2D<std::complex<double>, 6, 6> shift_el_sfm
            (std::complex<double> z);

    // Collocation points

    il::StaticArray<il::StaticArray<double, 3>, 6> el_cp_uniform
            (const il::StaticArray2D<double, 3, 3> &el_vert, double beta);

    il::StaticArray<il::StaticArray<double, 3>, 6> el_cp_nonuniform
            (const il::StaticArray2D<double, 3, 3> &el_vert,
             const il::StaticArray<double, 3> &vertex_wts,
             double beta);

    // auxiliary functions (norm, cross product)

    double l2norm(const il::StaticArray<double, 3> &a);
    il::StaticArray<double, 3> normalize(const il::StaticArray<double, 3> &a);
    il::StaticArray<double, 3> cross(const il::StaticArray<double, 3> &a,
                                     const il::StaticArray<double, 3> &b);

}

#endif //INC_3D_BEM_ELE_BASE_H

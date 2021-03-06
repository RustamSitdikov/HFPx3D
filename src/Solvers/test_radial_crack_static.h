//
// This file is part of HFPx3D.
//
// Created by nikolski on 11/14/2017.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//

#ifndef HFPX3D_SOLVE_RADIAL_CRACK_STATIC_H
#define HFPX3D_SOLVE_RADIAL_CRACK_STATIC_H

#include <cstring>
#include <il/StaticArray.h>

namespace hfp3d {

    double test_dd_radial_crack_static
            (std::string work_directory,
             std::string m_c_f_name,
             std::string m_n_f_name,
             int tt, double tol, bool echo);

    double test_stresses_radial_crack_static
            (std::string work_directory,
             std::string m_c_f_name,
             std::string m_n_f_name,
             il::int_t n_m_p,
             double tol, bool echo);
}

// Reference solutions

namespace ref {

// analytical solution for normally loaded penny-shaped crack
// in infinite elastic solid medium

// DD on the surface
    il::StaticArray<double, 3> get_dd_at_pt
            (double shear_m,
             double poiss_r,
             double a,
             double p,
             il::StaticArray<double, 3> crd);

// Stress field
    il::StaticArray<double, 6> get_stress_at_pt
            (double shear_m,
             double poiss_r,
             double a,
             double p,
             il::StaticArray<double, 3> crd);
}

#endif //HFPX3D_SOLVE_RADIAL_CRACK_STATIC_H

//
// This file is part of 3d_bem.
//
// Created by D. Nikolski on 1/24/2017.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//

//

#ifndef INC_3D_BEM_H_POTENTIAL_H
#define INC_3D_BEM_H_POTENTIAL_H

#include "Elast_Ker_Int.h"

namespace hfp3d {

    class H_Potential : public Kernel_Integration {
    public:
        il::StaticArray3D<std::complex<double>, 6, 3, 9> S11_22
                (double nu, std::complex<double> eix,
                 double h, std::complex<double> d);

        il::StaticArray3D<std::complex<double>, 6, 3, 9> S11_22_12
                (double nu, std::complex<double> eix,
                 double h, std::complex<double> d);

        il::StaticArray3D<std::complex<double>, 6, 3, 9> S13_23
                (double nu, std::complex<double> eix,
                 double h, std::complex<double> d);

        il::StaticArray3D<std::complex<double>, 6, 3, 9> S33
                (double nu, std::complex<double> eix,
                 double h, std::complex<double> d);

        il::StaticArray3D<std::complex<double>, 6, 3, 5> S11_22_red
                (double nu, std::complex<double> eix,
                 double h, std::complex<double> d);

        il::StaticArray3D<std::complex<double>, 6, 3, 5> S11_22_12_red
                (double nu, std::complex<double> eix,
                 double h, std::complex<double> d);

        il::StaticArray3D<std::complex<double>, 6, 3, 5> S13_23_red
                (double nu, std::complex<double> eix,
                 double h, std::complex<double> d);

        il::StaticArray3D<std::complex<double>, 6, 3, 5> S33_red
                (double nu, std::complex<double> eix,
                 double h, std::complex<double> d);

        il::StaticArray3D<std::complex<double>, 6, 4, 3> SijLim
                (double nu, std::complex<double> eix, std::complex<double> d);

    };

}
#endif //INC_3D_BEM_H_POTENTIAL_H

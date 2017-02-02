//
// This file is part of HFPx3D.
//
// Created by Brice Lecampion on 30.01.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//

#ifndef ESSAI_RECTANGULARELEMENTP0_H
#define ESSAI_RECTANGULARELEMENTP0_H

#include <il/StaticArray2D.h>

il::StaticArray2D<double, 3, 6> StressesKernelRectangularP0DD(
    double& x, double& y, double& z, double& a, double& b, double& G,
    double& nu) ;

#endif //ESSAI_RECTANGULARELEMENTP0_H

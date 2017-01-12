//
// This file is part of 3d_bem.
//
// Created by nikolski on 1/10/2017.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//

#ifndef INC_3D_BEM_SUBMATRIX_H
#define INC_3D_BEM_SUBMATRIX_H

#endif //INC_3D_BEM_SUBMATRIX_H

void get_submatrix(il::Array2D<double>& , int, int, int, int, const il::Array2D<double>&);
void get_static_submatrix(il::StaticArray2D& , int, int, const il::Array2D<double>&);
void get_static_sub_from_static(il::StaticArray2D& , int, int, const il::StaticArray2D&);
void set_submatrix(il::Array2D<double>& , int, int, const il::Array2D<double>&);
void set_submatrix_2_static(il::StaticArray2D& , int, int, const il::Array2D<double>&);
void set_static_sub_2_static(il::StaticArray2D& , int, int, const il::StaticArray2D&);
void add_submatrix(il::Array2D<double>& , int, int, const il::Array2D<double>&);
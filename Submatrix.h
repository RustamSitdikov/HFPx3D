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

//#include <il/Array2D.h>
//#include <il/StaticArray2D.h>

template <typename T_sub, typename T_A>
void get_submatrix(T_sub& sub, il::int_t i0, il::int_t i1, il::int_t j0, il::int_t j1, T_A& A);
//template <typename T_A>
//void get_submatrix(il::StaticArray2D& sub, il::int_t i0, il::int_t j0, T_A& A);
template <typename T_sub, typename T_A>
void set_submatrix(T_A& A, il::int_t i0, il::int_t i1, T_sub& B);
//void set_submatrix(il::StaticArray2D& A, il::int_t i0, il::int_t i1, il::Array2D<double>& B);

template <typename T_sub, typename T_A>
void get_submatrix(T_sub& sub, il::int_t i0, il::int_t i1, il::int_t j0, il::int_t j1, T_A& A) {
    IL_ASSERT((i1-i0+1) == sub.size(0));
    IL_ASSERT((j1-j0+1) == sub.size(1));
    IL_ASSERT(i1 <= A.size(0));
    IL_ASSERT(j1 <= A.size(1));

    for(il::int_t i = i0; i <= i1;++i) {
        for (il::int_t j=j0; j<= j1;++j){
            sub(i-i0,j-j0)=A(i,j);
        }
    }
}

//template <typename T_A>
//void get_submatrix(il::StaticArray2D& sub, il::int_t i0, il::int_t j0, T_A& A) {
//    il::int_t i1 = i0 + sub.size(0) - 1; il::int_t j1 = j0 + sub.size(1) - 1;

//    IL_ASSERT(i1 <= A.size(0));
//    IL_ASSERT(j1 <= A.size(1));

//    for(il::int_t i = i0; i <= i1;++i) {
//        for (il::int_t j=j0; j<= j1;++j){
//            sub(i-i0, j-j0) = A(i, j);
//        }
//    }
//}

template <typename T_sub, typename T_A>
void set_submatrix(T_A& A, il::int_t i0, il::int_t i1, T_sub& B) {
    IL_ASSERT(i0 + B.size(0) <= A.size(0));
    IL_ASSERT(i1 + B.size(1) <= A.size(1));

    for (il::int_t j1 = 0; j1 < B.size(1); ++j1) {
        for (il::int_t j0 = 0; j0 < B.size(0); ++j0) {
            A(i0 + j0, i1 + j1) = B(j0, j1);
        }
    }
}

//void set_submatrix(il::StaticArray2D& A, il::int_t i0, il::int_t i1, il::Array2D<double>& B) {
//    IL_ASSERT(i0 + B.size(0) <= A.size(0));
//    IL_ASSERT(i1 + B.size(1) <= A.size(1));

//    for (il::int_t j1 = 0; j1 < B.size(1); ++j1) {
//        for (il::int_t j0 = 0; j0 < B.size(0); ++j0) {
//            A(i0 + j0, i1 + j1) = B(j0, j1);
//        }
//    }
//}


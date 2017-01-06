//
// Created by D. Nikolski on 1/6/2017.
//

#ifndef INC_3D_BEM_LOCAL_IM_H
#define INC_3D_BEM_LOCAL_IM_H

#endif //INC_3D_BEM_LOCAL_IM_H

il::StaticArray<std::complex<double>, 9> ICFns(double, std::complex<double>, double, double, std::complex<double>);
il::StaticArray<std::complex<double>, 5> ICFns_red(double, std::complex<double>, double);

il::StaticArray2D<double, 6, 18> Local_IM_H(double, double, il::StaticArray<std::complex<double>,3>, il::StaticArray<std::complex<double>,3>, il::StaticArray2D<std::complex<double>,6,6>);
//il::StaticArray2D<double, 6, 18> Local_IM_T(double, double, il::StaticArray<std::complex<double>,3>, il::StaticArray<std::complex<double>,3, il::StaticArray2D<std::complex<double>,6,6>>);

//
// This file is part of 3d_bem.
//
// Created by D. Nikolski on 1/24/2017.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//

#include <cmath>
#include <complex>
#include <il/StaticArray.h>
#include <il/StaticArray2D.h>
#include <il/linear_algebra/dense/blas/dot.h>
//#include <il/linear_algebra/dense/norm.h>
#include "Ele_Base.h"

namespace hfp3d {

// auxiliary functions (norm, cross product)

    double L2norm(const il::StaticArray<double, 3> &a) {
        // L2 norm of a vector
        double N = 0.0;
        for (int k = 0; k < a.size(); ++k) {
            N += a[k] * a[k];
        }
        N = std::sqrt(N);
        return N;
    }

    il::StaticArray<double, 3> normalize
            (const il::StaticArray<double, 3> &a) {
        // normalized 3D vector
        il::StaticArray<double, 3> e;
        double N = L2norm(a); // il::norm(a, il::Norm::L2);
        for (int k = 0; k < a.size(); ++k) {
            e[k] = a[k] / N;
        }
        return e;
    }

    il::StaticArray<double, 3> cross(const il::StaticArray<double, 3> &a,
                                     const il::StaticArray<double, 3> &b) {
        // cross product of two 3D vectors
        IL_ASSERT(a.size() == 3);
        IL_ASSERT(b.size() == 3);
        il::StaticArray<double, 3> c;
        c[0] = a[1] * b[2] - a[2] * b[1];
        c[1] = a[2] * b[0] - a[0] * b[2];
        c[2] = a[0] * b[1] - a[1] * b[0];
        return c;
    }

// Element's local coordinate system manipulations

    il::StaticArray2D<double, 3, 3> El_LB_RT
            (const il::StaticArray2D<double, 3, 3> &EV) {
        // This function calculates the rotation tensor -
        // coordinate transform from the element's local Cartesian coordinate system
        // with origin at the first vertex of the element (EV(j, 0))
        // to the "global" (reference) Cartesian coordinate system
        il::StaticArray2D<double, 3, 3> RT;
        il::StaticArray<double, 3> a1{}, a2{}, a3{},
                e1{}, e2{}, e3{};
        for (int j = 0; j < 3; ++j) {
            a1[j] = EV(j, 1) - EV(j, 0);
            a2[j] = EV(j, 2) - EV(j, 0);
        }
        e1 = normalize(a1);
        a3 = cross(e1, a2);
        e3 = normalize(a3);
        e2 = normalize(cross(e3, e1));
        for (int j = 0; j < 3; ++j) {
            RT(j, 0) = e1[j];
            RT(j, 1) = e2[j];
            RT(j, 2) = e3[j];
        }
        return RT;
    }

    il::StaticArray<std::complex<double>, 3> El_RT_Tr
            (const il::StaticArray2D<double, 3, 3> &EV,
             const il::StaticArray2D<double, 3, 3> &RT) {
        // This function calculates the tau-coordinates
        // of the element's vertices
        il::StaticArray<std::complex<double>, 3> tau{0.0};
        il::StaticArray<double, 3> V0;
        il::StaticArray2D<double, 3, 3> BV{0.0};
        for (int k = 0; k < 3; ++k) {
            // Here the 1st vertex is chosen as the origin
            V0[k] = EV(k, 0);
        }
        for (int k = 0; k < 3; ++k) {
            for (int j = 0; j < 3; ++j) {
                // Basis vectors (rows)
                BV(j, k) = EV(k, j) - V0[k];
            }
        }
        // Rotated BV
        il::StaticArray2D<double, 3, 3> BVr = il::dot(BV, RT);
        for (int k = 0; k < 3; ++k) {
            tau[k] = std::complex<double>(BVr(k, 0), BVr(k, 1));
        }
        return tau;
    }

    il::StaticArray2D<std::complex<double>, 2, 2> El_CT
            (const il::StaticArray2D<double, 3, 3> &EV,
             const il::StaticArray2D<double, 3, 3> &RT) {
        // This function calculates the coordinate transform
        // from local Cartesian coordinates to "master element"
        // ([tau, conj(tau)] to [x,y]) with origin at the first vertex
        // of the element (EV(j, 0))
        il::StaticArray2D<std::complex<double>, 2, 2> MI{0.0};
        il::StaticArray<std::complex<double>, 2> z23{0.0};
        il::StaticArray<double, 3> xsi{0.0}, VV{0.0};
        std::complex<double> Dt;
        for (int k = 0; k < 2; ++k) {
            for (int n = 0; n < 3; ++n) {
                VV[n] = EV(n, k + 1) - EV(n, 0);
            }
            xsi = il::dot(VV, RT); // or transposed RT dot VV
            z23[k] = std::complex<double>(xsi[0], xsi[1]);
        }
        // common denominator (determinant)
        Dt = z23[0] * std::conj(z23[1]) - z23[1] * std::conj(z23[0]);
        // inverse transform
        MI(0, 0) = std::conj(z23[1]) / Dt;
        MI(0, 1) = -z23[1] / Dt;
        MI(1, 0) = -std::conj(z23[0]) / Dt;
        MI(1, 1) = z23[0] / Dt;
        return MI;
    }

    hfp3d::el_x_cr El_X_CR
            (const il::StaticArray2D<double, 3, 3> &EV,
             const il::StaticArray<double, 3> &X0,
             const il::StaticArray2D<double, 3, 3> &RT) {
        // This function calculates the h- and tau-coordinates of the point X0
        // with respect to the element's local Cartesian coordinate system
        // with origin at the first vertex of the element (EV(j, 0))
        el_x_cr h_z;
        il::StaticArray<double, 3> V0, X0r;
        //el_x_cr h_z;
        for (int k = 0; k < 3; ++k) {
            V0[k] = EV(k, 0);
            X0r[k] = X0[k] - V0[k];
        }
        X0r = il::dot(X0r, RT); // il::dot(RT_transposed, X0r)
        h_z.h = -X0r[2];
        h_z.z = std::complex<double>(X0r[0], X0r[1]);
        return h_z;
    }

// Element's basis (shape) functions

    il::StaticArray2D<std::complex<double>, 6, 6> El_SFM_S
            (const il::StaticArray2D<double, 3, 3> &EV,
             il::io_t, il::StaticArray2D<double, 3, 3> &RT) {
        // This function calculates the basis (shape) functions' coefficients (rows of SFM)
        // for a triangular boundary element with 2nd order polynomial approximation of unknowns
        // in terms of complex (tau, conj(tau)) representation of local element's coordinates
        // with trivial (middle) edge partitioning;
        // returns the same as El_SFM_N(EV, {1.0, 1.0, 1.0})

        // CT defines inverse coordinate transform [tau, conj(tau)] to [x,y] (see El_CT)
        RT = El_LB_RT(EV);
        il::StaticArray2D<std::complex<double>, 2, 2> CT = El_CT(EV, RT);
        il::StaticArray2D<std::complex<double>, 3, 3> CQ{0.0};
        il::StaticArray2D<std::complex<double>, 6, 6> SFM{0.0}, SFM_M{
                0.0}, CTau{0.0};

        // coefficients of shape functions (rows) for master element (0<=x,y<=1); ~[1, x, y, x^2, y^2, x*y]
        SFM_M(0, 0) = 1.0;
        SFM_M(0, 1) = -3.0;
        SFM_M(0, 2) = -3.0;
        SFM_M(0, 3) = 2.0;
        SFM_M(0, 4) = 2.0;
        SFM_M(0, 5) = 4.0;
        SFM_M(1, 1) = -1.0;
        SFM_M(1, 3) = 2.0;
        SFM_M(2, 2) = -1.0;
        SFM_M(2, 4) = 2.0;
        SFM_M(3, 5) = 4.0;
        SFM_M(4, 2) = 4.0;
        SFM_M(4, 4) = -4.0;
        SFM_M(4, 5) = -4.0;
        SFM_M(5, 1) = 4.0;
        SFM_M(5, 3) = -4.0;
        SFM_M(5, 5) = -4.0;

        // inverse coordinate transform [1, tau, tau_c, tau^2, tau_c^2, tau*tau_c] to [1, x, y, x^2, y^2, x*y]
        for (int j = 0; j <= 1; ++j) {
            for (int k = 0; k <= 1; ++k) {
                CQ(j, k) = CT(j, k) * CT(j, k);
            }
            CQ(j, 2) = 2.0 * CT(j, 0) * CT(j, 1);
            CQ(2, j) = CT(0, j) * CT(1, j);
        }
        CQ(2, 2) = CT(0, 0) * CT(1, 1) + CT(1, 0) * CT(0, 1);
        CTau(0, 0) = 1.0;
        for (int j = 0; j <= 1; ++j) {
            for (int k = 0; k <= 1; ++k) {
                CTau(j + 1, k + 1) = CT(j, k);
            }
        }
        for (int j = 0; j <= 2; ++j) {
            for (int k = 0; k <= 2; ++k) {
                CTau(j + 3, k + 3) = CQ(j, k);
            }
        }
        // assembly of SFM
        SFM = il::dot(SFM_M, CTau);
        return SFM;
    }

    il::StaticArray2D<std::complex<double>, 6, 6> El_SFM_N
            (const il::StaticArray2D<double, 3, 3> &EV,
             const il::StaticArray<double, 3> &VW,
             il::io_t, il::StaticArray2D<double, 3, 3> &RT) {
        // This function calculates the basis (shape) functions' coefficients (rows of SFM)
        // for a triangular boundary element with 2nd order polynomial approximation of unknowns
        // in terms of complex (tau, conj(tau)) representation of local element's coordinates
        // with non-trivial edge partitioning defined by "weights" VW

        double P12 = VW[0] / VW[1], P13 = VW[0] / VW[2], P23 = VW[1] / VW[2],
                C122 = P12 + 1.0,     // (VW[0]+w(2))/VW[1];
                C121 = 1.0 / P12 + 1.0,    // (VW[0]+w(2))/VW[0];
                C12q = C121 + C122,
                C233 = P23 + 1.0,     // (VW[1]+w(3))/VW[2];
                C232 = 1.0 / P23 + 1.0,    // (VW[1]+w(3))/VW[1];
                C23q = C232 + C233,
                C133 = P13 + 1.0,     // (VW[0]+w(3))/VW[2];
                C131 = 1.0 / P13 + 1.0,    // (VW[0]+w(3))/VW[0];
                C13q = C131 + C133;

        RT = El_LB_RT(EV);
        // CT defines inverse coordinate transform [tau, conj(tau)] to [x,y] (see El_CT)
        il::StaticArray2D<std::complex<double>, 2, 2> CT = El_CT(EV, RT);
        il::StaticArray2D<std::complex<double>, 3, 3> CQ{0.0};
        il::StaticArray2D<std::complex<double>, 6, 6> SFM{0.0}, SFM_M{
                0.0}, CTau{0.0};

        // coefficients of shape functions (rows) for master element (0<=x,y<=1); ~[1, x, y, x^2, y^2, x*y]
        SFM_M(0, 0) = 1.0;
        SFM_M(0, 1) = -P12 - 2.0;
        SFM_M(0, 2) = -P13 - 2.0;
        SFM_M(0, 3) = C122;
        SFM_M(0, 4) = C133;
        SFM_M(0, 5) = P13 + P12 + 2.0;
        SFM_M(1, 1) = -1.0 / P12;
        SFM_M(1, 3) = C121;
        SFM_M(1, 5) = 1.0 / P12 - P23;
        SFM_M(2, 2) = -1.0 / P13;
        SFM_M(2, 4) = C131;
        SFM_M(2, 5) = 1.0 / P13 - 1.0 / P23;
        SFM_M(3, 5) = C23q;
        SFM_M(4, 2) = C13q;
        SFM_M(4, 4) = -C13q;
        SFM_M(4, 5) = -C13q;
        SFM_M(5, 1) = C12q;
        SFM_M(5, 3) = -C12q;
        SFM_M(5, 5) = -C12q;

        // inverse coordinate transform [1, tau, tau_c, tau^2, tau_c^2, tau*tau_c] to [1, x, y, x^2, y^2, x*y]
        for (int j = 0; j <= 1; ++j) {
            for (int k = 0; k <= 1; ++k) {
                CQ(j, k) = CT(j, k) * CT(j, k);
            }
            CQ(j, 2) = 2.0 * CT(j, 0) * CT(j, 1);
            CQ(2, j) = CT(0, j) * CT(1, j);
        }
        CQ(2, 2) = CT(0, 0) * CT(1, 1) + CT(1, 0) * CT(0, 1);
        CTau(0, 0) = 1.0;
        for (int j = 0; j <= 1; ++j) {
            for (int k = 0; k <= 1; ++k) {
                CTau(j + 1, k + 1) = CT(j, k);
            }
        }
        for (int j = 0; j <= 2; ++j) {
            for (int k = 0; k <= 2; ++k) {
                CTau(j + 3, k + 3) = CQ(j, k);
            }
        }
        // assembly of SFM
        SFM = il::dot(SFM_M, CTau);
        return SFM;
    };

//il::StaticArray2D<std::complex<double>, 6, 6> El_SFM_C
// (const il::StaticArray2D<double, 3, 3> &EV,
// const il::StaticArray<double,3> &VW
// double beta,
// il::io_t, il::StaticArray2D<double, 3, 3> &RT) {
// This function calculates the basis (shape) functions' coefficients (rows of SFM)
// for a triangular boundary element with 2nd order polynomial approximation of unknowns
// in terms of complex (tau, conj(tau)) representation of local element's coordinates
// with nodes' offset to the centroid (e.g. at collocation points) defined by beta
// and non-trivial edge partitioning
//};

    il::StaticArray2D<std::complex<double>, 6, 6> El_Shift_SFM
            (std::complex<double> z) {
        // "shifted" SFM from z, tau[m], and local SFM
        std::complex<double> zc = std::conj(z);
        il::StaticArray2D<std::complex<double>, 6, 6> ShiftZ{0.0};
        ShiftZ(0, 0) = 1.0;
        ShiftZ(1, 0) = z;
        ShiftZ(1, 1) = 1.0;
        ShiftZ(2, 0) = zc;
        ShiftZ(2, 2) = 1.0;
        ShiftZ(3, 0) = z * z;
        ShiftZ(3, 1) = 2.0 * z;
        ShiftZ(3, 3) = 1.0;
        ShiftZ(4, 0) = zc * zc;
        ShiftZ(4, 2) = 2.0 * zc;
        ShiftZ(4, 4) = 1.0;
        ShiftZ(5, 0) = z * zc;
        ShiftZ(5, 1) = zc;
        ShiftZ(5, 2) = z;
        ShiftZ(5, 5) = 1.0;
        return ShiftZ;
    };

// Collocation points

    il::StaticArray<il::StaticArray<double, 3>, 6> El_CP_S
            (const il::StaticArray2D<double, 3, 3> &EV, double beta) {
        // This function calculates the coordinates
        // of the collocation points on a triangular boundary element
        // with 2nd order polynomial approximation of unknowns
        // and trivial (middle) edge partitioning;
        // offset of the points to the centroid is defined by beta;
        // returns the same as El_CP_N(EV, {1.0, 1.0, 1.0}, beta)
        il::StaticArray<il::StaticArray<double, 3>, 6> CP;
        il::StaticArray<double, 3> EC{0.0};
        int l, m;

        for (int j = 0; j < 3; ++j) {
            for (int k = 0; k < 3; ++k) {
                EC[j] += EV(j, k) / 3.0;
            }
        }
        for (int n = 0; n < 3; ++n) {
            m = (n + 1) % 3;
            l = (m + 1) % 3; // the edge across the (n-3)-th node
            for (int j = 0; j < 3; ++j) {
                (CP[n])[j] = (1.0 - beta) * EV(j, n) + beta * EC[j];
                (CP[n + 3])[j] = 0.5 * (1.0 - beta) * (EV(j, m) + EV(j, l)) +
                                 beta * EC[j];
            }
        }
        return CP;
    };

    il::StaticArray<il::StaticArray<double, 3>, 6> El_CP_N
            (const il::StaticArray2D<double, 3, 3> &EV,
             const il::StaticArray<double, 3> &VW,
             double beta) {
        // This function calculates the coordinates
        // of the collocation points on a triangular boundary element
        // with 2nd order polynomial approximation of unknowns
        // and non-trivial (middle) edge partitioning;
        // offset of the points to the centroid is defined by beta;
        // returns the same as El_CP_N(EV, {1.0, 1.0, 1.0}, beta)
        il::StaticArray<il::StaticArray<double, 3>, 6> CP;
        il::StaticArray<double, 3> EC{0.0};
        int l, m;

        for (int j = 0; j < 3; ++j) {
            for (int k = 0; k < 3; ++k) {
                EC[j] += EV(j, k) / 3.0;
            }
        }
        for (int n = 0; n < 3; ++n) {
            m = (n + 1) % 3;
            l = (m + 1) % 3; // the edge across the (n-3)-th node
            for (int j = 0; j < 3; ++j) {
                (CP[n])[j] = (1.0 - beta) * EV(j, n) + beta * EC[j];
                (CP[n + 3])[j] =
                        (1.0 - beta) * (VW[m] * EV(j, m) + VW[l] * EV(j, l)) /
                        (VW[m] + VW[l]) + beta * EC[j];
            }
        }
        return CP;
    };

}
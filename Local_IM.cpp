//
// This file is part of 3d_bem.
//
// Created by D. Nikolski on 1/6/2017.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//

// Calculation of local element-to-point influence matrices
// for a triangular boundary element with
// 2nd order polynomial approximation of unknowns

#include <il/StaticArray.h>
#include <il/StaticArray2D.h>
#include <il/StaticArray3D.h>
//#include <cmath>
#include <complex>
#include <SijH.h>
#include <ICFns.h>
//#include <il/linear_algebra/dense/blas/blas.h>
#include <il/linear_algebra/dense/blas/dot.h>

il::StaticArray2D<double, 6, 18> Local_IM_H(double mu, double nu, double h, std::complex<double> z, il::StaticArray<std::complex<double>,3> tau, il::StaticArray2D<std::complex<double>,6,6> SFM) {
    // This function assembles a local "stiffness" sub-matrix
    // (influence of DD at the element nodes to stresses at the point z)
    // in terms of a triangular element's local coordinates
    //
    // tau (3) are coordinates of element's vertices and
    // the rows of SFM (6*6) are coefficients of shape functions
    // in terms of the element's own local coordinate system (tau-coordinates);
    // h and z define the position of the (collocation) point x in the same coordinates

    const std::complex<double> I(0.0,1.0);

    // scaling ("-" sign comes from traction Somigliana ID, H-term)
    double scale = -mu/(4.0*M_PI*(1.0-nu));
    // tolerance parameters
    double HTol = 1.0E-9, DTol = 1.0E-8;
    // geometrical
    double an, am;
    std::complex<double> eixm, eixn, eipm, eipn, zc=std::conj(z), ntau2;
    // tz[m] and d[m] can be calculated here
    il::StaticArray<std::complex<double>,3> tz, d, dtau;
    int j, k, l, m, n;
    for (j=0; j<=2; ++j) {
        tz[j] = tau[j] - z;
        l = (j+1)%3;
        dtau[j]=tau[l]-tau[j];
        ntau2=dtau[j]/std::conj(dtau[j]);
        d[j]=0.5*(tz[j]-std::conj(tz[j])*ntau2);
    }
    // also, "shifted" SFM from z, tau[m], and local SFM
    il::StaticArray2D<std::complex<double>, 6, 6> ShiftZ {0.0};
    ShiftZ(0, 0) = 1.0;
    ShiftZ(1, 0) = z; ShiftZ(1, 1) = 1.0;
    ShiftZ(2, 0) = zc; ShiftZ(2, 2) = 1.0;
    ShiftZ(3, 0) = z*z; ShiftZ(3, 1) = 2.0*z; ShiftZ(3, 3) = 1.0;
    ShiftZ(4, 0) = zc*zc; ShiftZ(4, 2) =  2.0*zc; ShiftZ(4, 4) = 1.0;
    ShiftZ(5, 0) = z*zc; ShiftZ(5, 1) = zc; ShiftZ(5, 2) = z; ShiftZ(5, 5) = 1.0;
    il::StaticArray2D<std::complex<double>, 6, 6> SFMz = il::dot(SFM, ShiftZ);

    // constituents of the integrals
    il::StaticArray<std::complex<double>, 9> Fn, Fm;
    il::StaticArray3D<std::complex<double>, 6, 3, 9> Cn, Cm;
    il::StaticArray<std::complex<double>, 5> FnD, FmD;
    il::StaticArray3D<std::complex<double>, 6, 3, 5> CnD, CmD;
    // DD-to-stress influence
    // [(S11+S22)/2; (S11-S22)/2+i*S12; (S13+i*S23)/2; S33]
    // vs SF monomials (Sij_M) and nodal values (Sij_N)
    il::StaticArray3D<std::complex<double>, 6, 4, 3> Sij_M {0.0}, Sij_N {0.0},
            SincLn {0.0}, SincLm {0.0};
    // increment of DD-to-stress influence
    il::StaticArray2D<std::complex<double>, 6, 3> Sincm {0.0}, Sincn {0.0};

    il::StaticArray2D<double, 6, 18> LIM {0.0};

    // searching for "degenerate" edges: point x (collocation pt) projects onto an edge or a vertex
    bool IsDegen = std::abs(d[0])<DTol || std::abs(d[1])<DTol || std::abs(d[2])<DTol; // (d[0]*d[1]*d[2]==0);

    // calculating angles (phi, psi, chi)
    il::StaticArray<double,3> phi {0.0}, psi {0.0};
    il::StaticArray2D<double, 2, 3> chi {0.0};
    for (j=0; j<=2; ++j) {
        phi[j] = std::arg(tz[j]);
        // make sure it's between -pi and pi (add or subtract 2*pi)
        // phi[j] = (std::fmod(phi[j]/M_PI+1.0, 2.0)-1.0)*M_PI;
        psi[j] = std::arg(d[j]);
    }
    for (j=0; j<=2; ++j) {
        for (k=0; k<=1; ++k) {
            l = (j+k)%3;
            chi(k,j) = phi[l]-psi[j];
            // make sure it's between -pi and pi (add or subtract 2*pi)
            // chi(k,j) = (std::fmod(chi(k,j)/M_PI+1.0, 2.0)-1.0)*M_PI;
            chi(k,j) = (chi(k,j)<=-M_PI)? chi(k,j)+M_PI : ((chi(k,j)>M_PI)? chi(k,j)-M_PI : chi(k,j));
            // reprooving for "degenerate" edges
            if (fabs(M_PI_2-std::fabs(chi(k,j)))<DTol) IsDegen = true;
        }
    }

    // summation over edges
    for (m=0; m<=2; ++m) {
        n = (m+1)%3;
        if (std::abs(d[m])>=DTol &&
                std::fabs(M_PI_2-std::fabs(chi(0,m)))>=DTol &&
                std::fabs(M_PI_2-std::fabs(chi(1,m)))>=DTol) {
            eixm = std::exp(I*chi(0,m)); eixn = std::exp(I*chi(1,m));
            if(std::fabs(h)<HTol) { // limit case (point x on the element's plane)
                SincLn = SijLimH(nu, eixn, d[m]);
                SincLm = SijLimH(nu, eixm, d[m]);
                for (j=0; j<=5; ++j) {
                    for (l=0; l<=3; ++l) {
                        for (k=0; k<=2; ++k) {
                            Sij_M(j, l, k) += SincLn(j, l, k) - SincLm(j, l, k);
                        }
                    }
                }
            }
            else { // out-of-plane case
                an = std::abs(tz[n]-d[m])*((chi(1,m)<0)? -1.0 : double((chi(1,m)>0)));
                am = std::abs(tz[m]-d[m])*((chi(0,m)<0)? -1.0 : double((chi(0,m)>0)));
                Fn = ICFns(h, d[m], an, chi(1,m), eixn);
                Fm = ICFns(h, d[m], am, chi(0,m), eixm);
                // S11+S22
                Cn = S11_22H(nu, eixn, h, d[m]);
                Cm = S11_22H(nu, eixm, h, d[m]);
                Sincn = il::dot(Cn, Fn); Sincm = il::dot(Cm, Fm);
                for (j=0; j<=5; ++j) {
                    for (k=0; k<=2; ++k) {
                        Sij_M(j, 0, k) += Sincn(j, k) - Sincm(j, k);
                    }
                }
                if (IsDegen) { // "degenerate" case
                    eipm = std::exp(I*phi[n]); eipn = std::exp(I*phi[m]);
                    FnD = ICFns_red(h, d[m], an);
                    FmD = ICFns_red(h, d[m], am);
                    CnD = S11_22H_red(nu, eipn, h, d[m]);
                    CmD = S11_22H_red(nu, eipm, h, d[m]);
                    Sincn = il::dot(CnD, FnD); Sincm = il::dot(CmD, FmD);
                    for (j=0; j<=5; ++j) {
                        for (k=0; k<=2; ++k) {
                            Sij_M(j, 0, k) += Sincn(j, k) - Sincm(j, k);
                        }
                    }
                }
                // S11-S22+2*I*S12
                Cn = S11_22_12H(nu, eixn, h, d[m]);
                Cm = S11_22_12H(nu, eixm, h, d[m]);
                Sincn = il::dot(Cn, Fn); Sincm = il::dot(Cm, Fm);
                for (j=0; j<=5; ++j) {
                    for (k=0; k<=2; ++k) {
                        Sij_M(j, 1, k) += Sincn(j, k) - Sincm(j, k);
                    }
                }
                if (IsDegen) { // "degenerate" case
                    CnD = S11_22_12H_red(nu, eipn, h, d[m]);
                    CmD = S11_22_12H_red(nu, eipm, h, d[m]);
                    Sincn = il::dot(CnD, FnD); Sincm = il::dot(CmD, FmD);
                    for (j=0; j<=5; ++j) {
                        for (k=0; k<=2; ++k) {
                            Sij_M(j, 1, k) += Sincn(j, k) - Sincm(j, k);
                        }
                    }
                }
                // S13+I*S23
                Cn = S13_23H(nu, eixn, h, d[m]);
                Cm = S13_23H(nu, eixm, h, d[m]);
                Sincn = il::dot(Cn, Fn); Sincm = il::dot(Cm, Fm);
                for (j=0; j<=5; ++j) {
                    for (k=0; k<=2; ++k) {
                        Sij_M(j, 2, k) += Sincn(j, k) - Sincm(j, k);
                    }
                }
                if (IsDegen) { // "degenerate" case
                    CnD = S13_23H_red(nu, eipn, h, d[m]);
                    CmD = S13_23H_red(nu, eipm, h, d[m]);
                    Sincn = il::dot(CnD, FnD); Sincm = il::dot(CmD, FmD);
                    for (j=0; j<=5; ++j) {
                        for (k=0; k<=2; ++k) {
                            Sij_M(j, 2, k) += Sincn(j, k) - Sincm(j, k);
                        }
                    }
                }
                // S33
                Cn = S33H(nu, eixn, h, d[m]);
                Cm = S33H(nu, eixm, h, d[m]);
                Sincn = il::dot(Cn, Fn); Sincm = il::dot(Cm, Fm);
                for (j=0; j<=5; ++j) {
                    for (k=0; k<=2; ++k) {
                        Sij_M(j, 3, k) += Sincn(j, k) - Sincm(j, k);
                    }
                }
                if (IsDegen) { // "degenerate" case
                    CnD = S33H_red(nu, eipn, h, d[m]);
                    CmD = S33H_red(nu, eipm, h, d[m]);
                    Sincn = il::dot(CnD, FnD); Sincm = il::dot(CmD, FmD);
                    for (j=0; j<=5; ++j) {
                        for (k=0; k<=2; ++k) {
                            Sij_M(j, 3, k) += Sincn(j, k) - Sincm(j, k);
                        }
                    }
                }
            }
        }
    }

    // here comes contraction with "shifted" SFM (left)
    Sij_N = il::dot(SFMz, Sij_M);

    // re-shaping of the resulting matrix
    // and scaling (comment out if not necessary)
    for (j=0; j<=5; ++j) {
        int q = j*3;
        for (k=0; k<=2; ++k) {
            // [S11; S22; S12; S13; S23; S33] vs \delta{u}_k at j-th node
            LIM(0,q+k) = scale*(std::real(Sij_N(j,0,k))+std::real(Sij_N(j,1,k)));
            LIM(1,q+k) = scale*(std::real(Sij_N(j,0,k))-std::real(Sij_N(j,1,k)));
            LIM(2,q+k) = scale*std::imag(Sij_N(j,1,k));
            LIM(3,q+k) = scale*2.0*std::real(Sij_N(j,2,k));
            LIM(4,q+k) = scale*2.0*std::imag(Sij_N(j,2,k));
            LIM(5,q+k) = scale*std::real(Sij_N(j,3,k));
        }
    }

    return LIM;
}

// Constituing functions for the integrals
// of any kernel of the elasticity equation
// over a part of a polygonal element.
// Example of usage:
// dot(S11_22H(nu, eix, h, d), ICFns(h, d, a, x, eix))
// dot(S11_22H_red(nu, eip, h, d), ICFns_red(h, d, a))
// dot(S13_23T(nu, eix, h, d), ICFns(h, d, a, x, eix))
// dot(S33T_red(nu, eip, h, d), ICFns_red(h, d, a))
// where eip = std::exp(I*std::arg(t-z));
// eix = std::exp(I*x); x = std::arg(t-z)-std::arg(d);
// a = std::fabs(t-z-d)*sign(x);

// General case (h!=0, collocation point projected into or outside the element)
// powers of r, G0=arctan((ah)/(dr)), H0=arctanh(a/r) and its derivatives w.r. to h

il::StaticArray<std::complex<double>, 9> ICFns(double h, std::complex<double> d, double a, double x, std::complex<double> eix) {

    double D1 = std::abs(d), D2 = D1*D1, a2 = a*a,
            r = std::sqrt(h*h + a2 + D2),
            r2 = r*r, r3 = r2*r, r5 = r3*r2,
            ar = a/r, ar2 = ar*ar,
            hr = std::fabs(h/r),
            B = 1.0/(r2 - a2), B2=B*B, B3=B2*B;
    // evaluation of x and eix can be added here; then change arguments to pointer type
    double TanHi = std::imag(eix)/std::real(eix), tr = hr*TanHi,
            G0 = std::atan(tr), H0 = std::atanh(ar),
            H1 = -0.5*ar*B, H2 = 0.25*(3.0 - ar2)*ar*B2,
            H3 = -0.125*(15.0 - 10.0*ar2 + 3.0*ar2*ar2)*ar*B3;

    il::StaticArray<std::complex<double>, 9> BCE {0.0};
    BCE[0] = r; BCE[1] = 1.0/r; BCE[2] = 1.0/r3; BCE[3] = 1.0/r5;
    BCE[4] = G0-x; BCE[5] = H0; BCE[6] = H1; BCE[7] = H2; BCE[8] = H3;

    return BCE;
}

// Special case (reduced summation, collocation point projected onto the element contour) - additional terms

il::StaticArray<std::complex<double>, 5> ICFns_red(double h, std::complex<double> d, double a) {

    double h2 = h*h, h4 = h2*h2, h6 = h4*h2,
            D1 = std::abs(d), D2 = D1*D1, a2 = a*a,
            ro = std::sqrt(a2 + D2),
            r = std::sqrt(h2 + a2 + D2),
            rr = ro/r, rr2 = rr*rr, rr4 = rr2*rr2,
            L0 = std::atanh(rr), L1 = -0.5*rr/h2, L2 = 0.25*(3.0 - rr2)*rr/h4,
            L3 = -0.125*(15.0 - 10.0*rr2 + 3.0*rr4)*rr/h6;

    il::StaticArray<std::complex<double>, 5> BCE {0.0};
    BCE[0] = 1.0; BCE[1] = L0; BCE[2] = L1; BCE[3] = L2; BCE[4] = L3;

    return BCE;
}

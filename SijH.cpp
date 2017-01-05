//
// Created by D. Nikolski on 1/5/2017.
//
// Integration of the hypersingular kernel of the elasticity equation
// over a part of a polygonal element (a sector associated with one edge)
// with 2nd order polynomial approximating (shape) functions
//
// Stress components (vs local Cartesian coordinate system of the element)
// combined as S11+S22, S11-S22+2*I*S12, S13+I*S23, S33

#include <il/StaticArray.h>
#include <il/StaticArray3D.h>
#include <complex>
#include <cmath>

// General case (h!=0, collocation point projected into or outside the element)

il::StaticArray3D<std::complex<double>, 6, 3, 9> S11_22H(double nu, std::complex<double> eix, double h, std::complex<double> d) {
    const std::complex<double> I(0.0,1.0);

    double C1_1 = 1.0+nu;
    double C1_2 = 1.0+2.0*nu;
    double C2_1 = 2.0+nu;
    double C3_2 = 3.0+2.0*nu;
    double C4_1 = 4.0+nu;
    double C5_4 = 5.0+4.0*nu;
    double C7_2 = 7.0+2.0*nu;
    double C7_5 = 7.0+5.0*nu;
    double C7_6 = 7.0+6.0*nu;
    double C11_4 = 11.0+4.0*nu;
    double C11_5 = 11.0+5.0*nu;
    double C13_2 = 13.0+2.0*nu;
    double C13_10 = 13.0+10.0*nu;

    double CosX = std::real(eix);
    double TanX = std::imag(eix)/CosX;
    std::complex<double> TCosX = CosX*eix;

    double h2 = h*h; double h4 = h2*h2;
    double sgh = (( h < 0 )? -1.0 : double(( h > 0 ))); // sign(h)

    double D1 = std::abs(d); double D2 = D1*D1; double D4 = D2*D2;
    // std::complex<double> dC = std::conj(d);
    std::complex<double> D0 = std::polar(1.0,std::arg(d)); //  = d/D1
    std::complex<double> D02 = D0*D0; //  = d^2/D1^2

    std::complex<double> D1_h1 = D2+h2;
    std::complex<double> D1_h3 = D2+3.0*h2;
    std::complex<double> D1mh3 = D2-3.0*h2;

    std::complex<double> P0, P1, P2;

    il::StaticArray3D<std::complex<double>, 6, 3, 9> C{0.0};

    C(0, 0, 2) = h*std::imag(d);
    C(0, 1, 2) = -h*std::real(d);
    P0 = 0.1875*h; P1 = 3.0*h2; P2 = D2*TanX;
    C(0, 0, 3) = -P0*(P1*std::imag(d)+P2*std::real(d));
    C(0, 1, 3) = P0*(P1*std::real(d)-P2*std::imag(d));
    P0 = C7_2*h;
    C(0, 0, 6) = P0*std::real(D0);
    C(0, 1, 6) = P0*std::imag(D0);
    C(0, 2, 6) = -C1_2*D1;
    P0 = 3.0*h*D1_h3;
    C(0, 0, 7) = P0*std::real(D0);
    C(0, 1, 7) = P0*std::imag(D0);
    C(0, 2, 7) = -2.0*h2*D1;
    P0 = -0.5*h*D1mh3*D1_h1;
    C(0, 0, 8) = P0*std::real(D0);
    C(0, 1, 8) = P0*std::imag(D0);

    C(1, 1, 1) = 0.2*C11_5*h*D02*TCosX;
    C(1, 0, 1) = I*C(1, 1, 1);
    P1 = 0.1*(D2*(7.0+2.0*I*TanX)+16.0*h2*TCosX)*D02; P2 = C7_5/60.0*D2*TanX;
    C(1, 0, 2) = -(I*P1+P2)*h;
    C(1, 1, 2) = -(P1+I*P2)*h;
    P1 = (D2*h2*(0.4+0.11875*I*TanX)+0.09375*I*D4*TanX+0.4*h4*TCosX)*D02;
    P2 = (-0.05625*D2+0.11875*h2)*D2*TanX;
    C(1, 0, 3) = (I*P1+P2)*h;
    C(1, 1, 3) = (P1+I*P2)*h;
    C(1, 0, 4) = -C1_1*sgh;
    C(1, 1, 4) = -I*C1_1*sgh;
    C(1, 2, 5) = 0.5*C1_2*D0;
    P1 = 0.5*C7_2*d*D0; P2 = D1*(0.3+4.0/3.0*C2_1);
    C(1, 0, 6) = (P1+P2)*h;
    C(1, 1, 6) = I*(-P1+P2)*h;
    C(1, 2, 6) = 2.0*C2_1*h2*D0;
    P1 = 1.5*d*D0*D1_h3; P2 = D1*(1.0/6.0*C1_2*D2+h2*(43.0/30.0+C2_1/3.0));
    C(1, 0, 7) = (P1+P2)*h;
    C(1, 1, 7) = I*(-P1+P2)*h;
    C(1, 2, 7) = 2.0*h4*D0;
    P0 = h*D1_h1; P1 = 0.15*D1*D2-1.9/6.0*D1*h2; P2 = 0.25*d*D0*D1mh3;
    C(1, 0, 8) = -P0*(P1+P2);
    C(1, 1, 8) = I*P0*(-P1+P2);

    for (int k = 0; k < C.size(2); ++k) {
        for (int j = 0; j < C.size(1); ++j) {
            C(2, j, k) = std::conj(C(1, j, k));
        }
    }

    C(3, 2, 0) = I*C1_2*D02*TCosX;
    P0 = d*h; P1 = 0.0625*C13_10; P2 = D02*(0.0625*C13_10+0.5*C3_2*TCosX);
    C(3, 0, 1) = -I*P0*(P1-P2);
    C(3, 1, 1) = P0*(P1+P2);
    C(3, 2, 1) = 2.0*I*C2_1*D02*h2*TCosX;
    //P1 = ; P2 = ;
    C(3, 0, 2) = d*h*(0.09375*I*C7_2*h2-0.03125*C3_2*D2*TanX+D02*(D2*(-1.0/12.0*I*C7_6+0.09375*C3_2*TanX)-I*h2*(0.09375*C7_2+0.25*C4_1*TCosX)));
    C(3, 1, 2) = d*h*(-0.09375*C7_2*h2-I*0.03125*C3_2*D2*TanX-D02*(D2*(1.0/12.0*C7_6+0.09375*I*C3_2*TanX)+h2*(0.09375*C7_2+0.25*C4_1*TCosX)));
    C(3, 2, 2) = -I*D02*h4*TCosX;
    //P0 = ; P1 = ; P2 = ;
    C(3, 0, 3) = d*h*(D2*TanX*((0.09375-0.28125*D02)*h2-(0.046875+0.109375*D02)*D2)+I*h2*(0.625*D02*D2-0.234375*h2+D02*h2*(0.234375+0.3125*TCosX)));
    C(3, 1, 3) = d*h*(I*D2*TanX*((0.09375+0.28125*D02)*h2+(-0.046875+0.109375*D02)*D2)+h2*(0.625*D02*D2+0.234375*h2+D02*h2*(0.234375+0.3125*TCosX)));
    C(3, 0, 5) = C3_2*h*D0*(0.25*D02-0.75);
    C(3, 1, 5) = -I*C3_2*h*D0*(0.25*D02+0.75);
    C(3, 2, 5) = 0.5*C1_2*d*D0;
    P0 = h*D0; P1 = D02*(0.75*C5_4*D2+0.25*C11_4*h2); P2 = 0.25*C5_4*D2+0.75*C11_4*h2;
    C(3, 0, 6) = h*D0*(P1-P2);
    C(3, 1, 6) = -I*h*D0*(P1+P2);
    C(3, 2, 6) = 2.0*C2_1*h2*d*D0;
    //P1 = ; P2 = ;
    C(3, 0, 7) = h*D0*(0.125*C1_2*D4-0.25*C7_2*D2*h2-0.375*C13_2*h4+D02*(0.625*C1_2*D4+0.75*C7_2*D2*h2+0.125*C13_2*h4));
    C(3, 1, 7) = I*h*D0*(0.125*C1_2*D4-0.25*C7_2*D2*h2-0.375*C13_2*h4-D02*(0.625*C1_2*D4+0.75*C7_2*D2*h2+0.125*C13_2*h4));
    C(3, 2, 7) = 2.0*h4*d*D0;
    P0 = h*D0*D1_h1; P1 = D02*(7.0/24.0*D4-11.0/12.0*D2*h2-5.0/24.0*h4); P2 = 0.125*D4-0.25*D2*h2+0.625*h4;
    C(3, 0, 8) = -P0*(P1+P2);
    C(3, 1, 8) = I*P0*(P1-P2);

    for (int k = 0; k < C.size(2); ++k) {
        for (int j = 0; j < C.size(1); ++j) {
            C(4, j, k) = std::conj(C(3, j, k));
        }
    }

    P0 = 0.125*C13_10*h;
    C(5, 0, 1) = P0*std::imag(d);
    C(5, 1, 1) = -P0*std::real(d);
    C(5, 2, 1) = 1.0/12.0*C1_2*D2*TanX;
    P1 = 0.0625*C3_2*D2*TanX; P2 = 0.1875*C7_2*h2;
    C(5, 0, 2) = -h*(P1*std::real(d)+P2*std::imag(d));
    C(5, 1, 2) = -h*(P1*std::imag(d)-P2*std::real(d));
    C(5, 2, 2) = -1.0/6.0*C2_1*D2*h2*TanX;
    P1 = 0.09375*D2*TanX*(2*h2-D2); P2 = 0.46875*h4;
    C(5, 0, 3) = h*(P1*std::real(d)+P2*std::imag(d));
    C(5, 1, 3) = h*(P1*std::imag(d)-P2*std::real(d));
    C(5, 2, 3) = 0.25*D2*h4*TanX;
    C(5, 2, 4) = -4.0*C1_1*std::fabs(h);
    P0 = -1.5*C3_2*h;
    C(5, 0, 5) = P0*std::real(D0);
    C(5, 1, 5) = P0*std::imag(D0);
    C(5, 2, 5) = -0.5*C1_2*D1;
    P0 = -h*(0.5*C5_4*D2+1.5*C11_4*h2);
    C(5, 0, 6) = P0*std::real(D0);
    C(5, 1, 6) = P0*std::imag(D0);
    C(5, 2, 6) = (1.0/6.0*C1_2*D2+1.5*h2*(5+2*nu))*D1;
    P0 = h*(0.25*C1_2*D4-0.5*C7_2*D2*h2-0.75*C13_2*h4);
    C(5, 0, 7) = P0*std::real(D0);
    C(5, 1, 7) = P0*std::imag(D0);
    C(5, 2, 7) = 2.0/3.0*h2*(C2_1*D2+h2*(7+nu))*D1;
    P0 = -h*D1_h1*(0.25*D4-0.5*D2*h2+1.25*h4);
    C(5, 0, 8) = P0*std::real(D0);
    C(5, 1, 8) = P0*std::imag(D0);
    C(5, 2, 8) = 2.0/3.0*h4*D1_h1*D1;

    return C;
}

il::StaticArray3D<std::complex<double>, 6, 3, 9> S11_22_12H(double nu, std::complex<double> eix, double h, std::complex<double> d) {
    const std::complex<double> I(0.0,1.0);

    double C1m1 = 1.0-nu;
    double C1m2 = 1.0-2.0*nu;
    double C2m1 = 2.0-nu;
    double C3m1 = 3.0-nu;
    double C3m4 = 3.0-4.0*nu;
    double C5m2 = 5.0-2.0*nu;
    double C5m4 = 5.0-4.0*nu;
    double C6m5 = 6.0-5.0*nu;
    double C7m2 = 7.0-2.0*nu;
    double C8m5 = 8.0-5.0*nu;
    double C9m2 = 9.0-2.0*nu;
    double C9m4 = 9.0-4.0*nu;
    double C9m8 = 9.0-8.0*nu;
    double C13m2 = 13.0-2.0*nu;
    double C15m4 = 15.0-4.0*nu;
    double C15m8 = 15.0-8.0*nu;
    double C115m38d80 = 1.4375-0.475*nu;

    double CosX = std::real(eix);
    double TanX = std::imag(eix)/CosX;
    std::complex<double> TCosX = CosX*eix;
    std::complex<double> TCosX1m1 = TCosX-1.0;
    std::complex<double> CTCosX4_3 = 3.0+4.0*TCosX;
    std::complex<double> e2X = eix*eix;
    std::complex<double> CTCosX1 = TCosX1m1*TCosX;
    std::complex<double> WCTCosX1 = (13.0+e2X-10.0*TCosX)*TCosX;

    double h2 = h*h; double h4 = h2*h2;
    double sgh = (( h < 0 )? -1.0 : double(( h > 0 ))); // sign(h)

    double D1 = std::abs(d); double D2 = D1*D1; double D4 = D2*D2;
    std::complex<double> dC = std::conj(d);
    std::complex<double> D0 = std::polar(1.0,std::arg(d)); //  = d/D1
    std::complex<double> D0C = std::conj(D0);
    std::complex<double> D02 = D0*D0; //  = d^2/D1^2
    std::complex<double> D03 = D0*D02; //  = d^3/D1^3
    std::complex<double> D04 = D02*D02; //  = d^4/D1^4

    std::complex<double> D2h2 = D2*h2;

    std::complex<double> D1_h1 = D2+h2;
    std::complex<double> D1_h3 = D2+3.0*h2;
    std::complex<double> D3_h1 = 3.0*D2+h2;
    std::complex<double> D1mh3 = D2-3.0*h2;

    std::complex<double> P0, P1, P2;

    il::StaticArray3D<std::complex<double>, 6, 3, 9> C{0.0};

    C(0, 2, 1) = -I*C1m2*D02*TCosX;
    C(0, 0, 2) = I*d*h*(-0.5+D02*(0.5+0.75*TCosX));
    C(0, 1, 2) = d*h*(0.5+D02*(0.5+0.75*TCosX));
    C(0, 2, 2) = I*h2*D02*TCosX;
    //P0 = d*h; P1 = ; P2 = ;
    C(0, 0, 3) = d*h*(0.28125*I*h2-0.09375*D2*TanX+D02*(D2*(-0.75*I+0.28125*TanX)-I*h2*(0.28125+0.375*TCosX)));
    C(0, 1, 3) = -d*h*(0.28125*h2+0.09375*I*D2*TanX+D02*(D2*(0.75+0.28125*I*TanX)+h2*(0.28125+0.375*TCosX)));
    P0 = 0.5*D0*h; P1 = 3.0*D02;
    C(0, 0, 6) = P0*(-P1+C9m4);
    C(0, 1, 6) = I*P0*(P1+C9m4);
    C(0, 2, 6) = -C1m2*D0*d;
    P0 = 1.5*h*D0; P2 = D3_h1*D02;
    C(0, 0, 7) = P0*(D1_h3-P2);
    C(0, 1, 7) = I*P0*(D1_h3+P2);
    C(0, 2, 7) = -2.0*h2*D0*d;
    P0 = 0.25*h*D0*D1_h1; P2 = D02*(5.0*D2+h2);
    C(0, 0, 8) = -P0*(D1mh3+P2);
    C(0, 1, 8) = I*P0*(-D1mh3+P2);

    P0 = 0.4*D02*h*TCosX; P2 = 8.0*D02*(-1.0+TCosX);
    C(1, 0, 1) = I*P0*(C8m5+P2);
    C(1, 1, 1) = P0*(-C8m5+P2);
    C(1, 2, 1) = -I*C1m2*d*D02*(0.625+TCosX);
    P0 = D02*h; //P1 = ; P2 = ;
    C(1, 0, 2) = P0*(D2*(-0.7*I+0.2*TanX)-1.6*I*h2*TCosX-2.0*D02*(0.8*I*h2*CTCosX1+D2*(0.1*TanX-I*(0.7+0.4*TCosX))));
    C(1, 1, 2) = P0*(D2*(0.7+0.2*I*TanX)+1.6*h2*TCosX+2.0*D02*(-0.8*h2*CTCosX1+D2*(0.7+0.1*I*TanX+0.4*TCosX)));
    C(1, 2, 2) = d*D02*(-D2*C1m2*(-0.5*I+0.1875*TanX)+I*h2*(1.1875+1.75*TCosX-0.125*nu*CTCosX4_3));
    //P1 = ; P2 = ;
    C(1, 0, 3) = P0*(D2h2*(0.4*I-0.11875*TanX)-0.09375*D4*TanX+0.4*I*h4*TCosX+
                     D02*(3.0*D4*(-0.4*I+0.18125*TanX)+0.4*I*h4*CTCosX1+D2h2*(0.11875*TanX-0.4*I*(2.0+TCosX))));
    C(1, 1, 3) = -P0*(D2h2*(0.4+0.11875*I*TanX)+0.09375*I*D4*TanX+0.4*h4*TCosX+
                      D02*(3.0*D4*(0.4+0.18125*I*TanX)-0.4*h4*CTCosX1+D2h2*(0.11875*I*TanX+0.4*(2.0+TCosX))));
    C(1, 2, 3) = -I*d*D02*h2*(D2*(1.5+0.5625*I*TanX)+0.1875*h2*CTCosX4_3);
    C(1, 2, 5) = -0.5*C1m2*D03;
    P0 = d*D0*h;
    C(1, 0, 6) = -P0*(4.5*D02-0.5*C9m4);
    C(1, 1, 6) = I*P0*(4.5*D02+0.5*C9m4);
    C(1, 2, 6) = -D03*(2.0*C2m1*h2+3.0*C1m2*D2);
    P1 = 1.5*D2+4.5*h2; P2 = D02*(7.5*D2+4.5*h2);
    C(1, 0, 7) = P0*(P1-P2);
    C(1, 1, 7) = I*P0*(P1+P2);
    C(1, 2, 7) = -0.25*D03*(5.0*D4*C1m2+6.0*D2h2*C7m2+h4*C13m2);
    P0 = 0.25*P0*D1_h1; P1 = D2-3.0*h2; P2 = D02*(7.0*D2+3.0*h2);
    C(1, 0, 8) = -P0*(P1+P2);
    C(1, 1, 8) = I*P0*(-P1+P2);
    C(1, 2, 8) = -0.5*D03*h2*D1_h1*(5.0*D2+h2);

    C(2, 0, 1) = 3.2*I*h*D02*TCosX;
    C(2, 1, 1) = 3.2*h*D02*TCosX;
    C(2, 2, 1) = 0.625*I*C1m2*d;
    P1 = 0.1*D02*(D2*(7.0+2.0*I*TanX)+16.0*h2*TCosX); P2 = I*C6m5/30.0*D2*TanX;
    C(2, 0, 2) = -I*h*(P1-P2);
    C(2, 1, 2) = -h*(P1+P2);
    C(2, 2, 2) = d*(0.0625*C1m2*D2*TanX-I*h2*(1.1875-0.375*nu));
    P1 = D2*TanX*(-0.05625*D2+0.11875*h2); P2 = D02*(D2*TanX*(0.11875*h2+0.09375*D2)-0.4*I*h2*(D2+h2*TCosX));
    C(2, 0, 3) = h*(P1-P2);
    C(2, 1, 3) = I*h*(P1+P2);
    C(2, 2, 3) = d*h2*(-0.1875*D2*TanX+0.5625*I*h2);
    C(2, 0, 4) = -2.0*C1m1*sgh;
    C(2, 1, 4) = I*C(2, 0, 4);
    C(2, 2, 5) = 1.5*C1m2*D0;
    P1 = 4.5*d*D0*h; P2 = D1*h*(4.3-8.0/3.0*nu);
    C(2, 0, 6) = P1+P2;
    C(2, 1, 6) = I*(-P1+P2);
    C(2, 2, 6) = D0*(6.0*C2m1*h2+C1m2*D2);
    P1 = 1.5*d*D0*h*D1_h3; P2 = 1.0/3.0*D1*h*((2.0*C2m1+3.3)*h2+0.5*C3m4*D2);
    C(2, 0, 7) = P1+P2;
    C(2, 1, 7) = I*(-P1+P2);
    C(2, 2, 7) = D0*(3.0*h2*D1_h3-0.25*C1m2*D1mh3*D1_h1);
    P0 = h*D1*D1_h1; P1 = 0.15*D2-0.95/3.0*h2; P2 = 0.25*D02*D1mh3;
    C(2, 0, 8) = -P0*(P1+P2);
    C(2, 1, 8) = I*P0*(-P1+P2);
    C(2, 2, 8) = -0.5*D0*h2*D1mh3*D1_h1;

    C(3, 2, 0) = 6.4/3.0*I*D04*C1m2*CTCosX1;
    P0 = D02*d*h; P1 = D02*(1.4625+0.375*WCTCosX1); P2 = 1.25*(0.15+C1m1)+0.5*C5m4*TCosX;
    C(3, 0, 1) = -I*P0*(P1-P2);
    C(3, 1, 1) = -P0*(P1+P2);
    C(3, 2, 1) = I*D04*(12.8/3.0*C2m1*CTCosX1*h2-C1m2/3.0*(5.2+0.725*I*TanX+3.2*TCosX)*D2);
    //P1 = ; P2 = ;
    C(3, 0, 2) = I*P0*(D02*(D2*(3.275+TanX*(0.78125*I+0.025*TanX)+TCosX)+h2*(0.86875+0.1875*WCTCosX1))-
                       D2*(C1m1+0.25/3.0+0.09375*I*C5m4*TanX)-h2*(0.0625*C5m2*CTCosX4_3-0.09375));
    C(3, 1, 2) = P0*(D02*(D2*(3.275+TanX*(0.78125*I+0.025*TanX)+TCosX)+h2*(0.86875+0.1875*WCTCosX1))+
                     D2*(C1m1+0.25/3.0+0.09375*I*C5m4*TanX)+h2*(0.0625*C5m2*CTCosX4_3-0.09375));
    C(3, 2, 2) = D04*(-D4*C1m2*(-0.8*I+0.3625*TanX)-0.8/3.0*I*h4*C13m2*CTCosX1+
                      D2h2/3.0*(-C115m38d80*TanX+0.4*I*(1.0+8.0*C3m1+2.0*C7m2*TCosX)));
    //P1 = ; P2 = ;
    C(3, 0, 3) = P0*((D2h2*(0.625*I-0.28125*TanX)-0.109375*D4*TanX+0.078125*I*h4*CTCosX4_3)+
                     D02*(D4*(-2.0*I+1.015625*TanX)-I*h4*(0.234375+0.046875*WCTCosX1)+D2h2*(0.46875*TanX-0.125*I*(15.0+4.0*TCosX))));
    C(3, 1, 3) = -P0*((D2h2*(0.625+0.28125*I*TanX)+0.109375*I*D4*TanX+0.078125*h4*CTCosX4_3)+
                      D02*(D4*(2.0+1.015625*I*TanX)+h4*(0.234375+0.046875*WCTCosX1)+D2h2*(0.46875*I*TanX+0.125*(15.0+4.0*TCosX))));
    C(3, 2, 3) = D04*h2*(3.0*D4*(-0.8*I+0.3625*TanX)+0.8*I*h4*CTCosX1+D2h2*(0.2375*TanX-0.8*I*(2.0+TCosX)));
    P0 = 0.25*h*D03;
    C(3, 0, 5) = P0*(-3.0*D02+C5m4);
    C(3, 1, 5) = I*P0*(3.0*D02+C5m4);
    C(3, 2, 5) = -1.5*d*D03*C1m2;
    P1 = 3.0*D2*C9m8+h2*C15m8; P2 = 9.0*D02*(5.0*D2+h2);
    C(3, 0, 6) = P0*(P1-P2);
    C(3, 1, 6) = I*P0*(P1+P2);
    C(3, 2, 6) = -d*D03*(6.0*h2*C2m1+5.0*D2*C1m2);
    P0 = 0.5*P0; P1 = 5.0*D4*C3m4+6.0*D2h2*C9m4+h4*C15m4; P2 = 3.0*D02*(35.0*D4+30.0*D2h2+3.0*h4);
    C(3, 0, 7) = P0*(P1-P2);
    C(3, 1, 7) = I*P0*(P1+P2);
    C(3, 2, 7) = -d*D03*(0.75*h4*C13m2+2.5*D2h2*C7m2+1.75*D4*C1m2);
    P0 = P0*D1_h1; P1 = (-7.0*D4+22.0*D2h2+5.0*h4)/3.0; P2 = D02*(21.0*D4+14.0*D2h2+h4);
    C(3, 0, 8) = P0*(P1-P2);
    C(3, 1, 8) = I*P0*(P1+P2);
    C(3, 2, 8) = -d*D03*h2*D1_h1*(3.5*D2+1.5*h2);

    P1 = 1.25*nu*dC;
    C(4, 0, 1) = h*(2.875*std::imag(d)-I*P1);
    C(4, 1, 1) = h*(P1-2.875*std::real(d));
    C(4, 2, 1) = 0.725/3.0*C1m2*D2*TanX;
    P0 = 0.0625*h; P1 = 6.0*nu*I*h2-C5m2*D2*TanX; P2 = I*(3.0*C9m2*I*h2-2.0*nu*D2*TanX);
    C(4, 0, 2) = P0*(P1*std::real(d)+P2*std::imag(d));
    C(4, 1, 2) = P0*(P1*std::imag(d)-P2*std::real(d));
    C(4, 2, 2) = D2*TanX*(0.0375*C1m2*D2-C115m38d80/3.0*h2);
    P0 = 0.09375*h; P1 = 5.0*h4; P2 = (-D4+2.0*D2h2)*TanX;
    C(4, 0, 3) = P0*(P1*std::imag(d)+P2*std::real(d));
    C(4, 1, 3) = P0*(-P1*std::real(d)+P2*std::imag(d));
    C(4, 2, 3) = D2*h2*TanX*(0.2375*h2-0.1125*D2);
    C(4, 2, 4) = -8.0*C1m1*std::fabs(h);
    P0 = 1.5*h; P2 = 2.0*I*nu;
    C(4, 0, 5) = P0*(-C5m2*std::real(D0)-P2*std::imag(D0));
    C(4, 1, 5) = P0*(-C5m2*std::imag(D0)+P2*std::real(D0));
    C(4, 2, 5) = -1.5*C1m2*D1;
    P1 = 2.0*h*D1_h3*nu*D0C; P2 = 4.5*h*(D2+5.0*h2);
    C(4, 0, 6) = P1-P2*std::real(D0);
    C(4, 1, 6) = I*P1-P2*std::imag(D0);
    C(4, 2, 6) = D1*(1.0/3.0*C1m2*D2+(3.6*C3m1-0.4)*h2);
    P1 = 0.5*nu*h*D1_h1*D1mh3*D0C; P2 = 0.75*h*(D4-6.0*D2h2-15.0*h4);
    C(4, 0, 7) = -P1+P2*std::real(D0);
    C(4, 1, 7) = -I*P1+P2*std::imag(D0);
    C(4, 2, 7) = D1*(-0.15*C1m2*D4+1.0/6.0*C7m2*D2h2+(0.35+1.9*(8-nu))/3.0*h4);
    P2 = -0.25*h*D1_h1*(D4-2.0*D2h2+5.0*h4);
    C(4, 0, 8) = P2*std::real(D0);
    C(4, 1, 8) = P2*std::imag(D0);
    C(4, 2, 8) = D1*h2*D1_h1*(1.9/3.0*h2-0.3*D2);

    C(5, 2, 0) = 6.4/3.0*I*C1m2*D02*TCosX;
    P0 = d*h; P1 = 0.1875+1.25*C1m1; P2 = D02*(1.4375+2.5*TCosX);
    C(5, 0, 1) = I*P0*(-P1+P2);
    C(5, 1, 1) = P0*(P1+P2);
    C(5, 2, 1) = D02/3.0*(C1m2*D2*(2.6*I-0.725*TanX)+12.8*C2m1*I*h2*TCosX);
    //P1 = ; P2 = ;
    C(5, 0, 2) = P0*(0.09375*I*h2*C9m4-0.03125*D2*C5m4*TanX+
                     D02*(D2*(-3.25/3.0*I+0.46875*TanX)-0.3125*I*h2*(2.7+4.0*TCosX)));
    C(5, 1, 2) = -P0*(0.09375*h2*C9m4+0.03125*I*D2*C5m4*TanX+
                      D02*(D2*(3.25/3.0+0.46875*I*TanX)+0.3125*h2*(2.7+4.0*TCosX)));
    C(5, 2, 2) = D02*(0.0625*C1m2*TanX*D4+(I*(-5.0+1.6*nu)+C115m38d80*TanX)/3.0*D2h2-0.8/3.0*I*C13m2*TCosX*h4);
    //P1 = ; P2 = ;
    C(5, 0, 3) = P0*(-(0.234375*I*h4-0.09375*D2h2*TanX+0.046875*D4*TanX)+
                     D02*(0.078125*I*h4*CTCosX4_3+(0.625*I-0.28125*TanX)*D2h2-0.109375*D4*TanX));
    C(5, 1, 3) = P0*((0.234375*h4+0.09375*I*D2h2*TanX-0.046875*I*D4*TanX)+
                     D02*(0.078125*h4*CTCosX4_3+(0.625+0.28125*I*TanX)*D2h2+0.109375*I*D4*TanX));
    C(5, 2, 3) = D02*h2*(-0.1875*D4*TanX+(0.8*I-0.2375*TanX)*D2h2+0.8*I*h4*TCosX);
    P0 = D0*h; P1 = 1.25*D02; P2 = 3.75-3.0*nu;
    C(5, 0, 5) = P0*(P1-P2);
    C(5, 1, 5) = -I*P0*(P1+P2);
    C(5, 2, 5) = 1.5*d*D0*C1m2;
    P0 = 0.25*P0; P1 = D02*(15.0*h2+27.0*D2); P2 = 3.0*C15m8*h2+C9m8*D2;
    C(5, 0, 6) = P0*(P1-P2);
    C(5, 1, 6) = -I*P0*(P1+P2);
    C(5, 2, 6) = d*D0*(C1m2*D2+6.0*C2m1*h2);
    P0 = 0.5*P0; P1 = C3m4*D4-2.0*C9m4*D2h2-3.0*C15m4*h4; P2 = 3.0*D02*(5.0*D4+18.0*D2h2+5.0*h4);
    C(5, 0, 7) = P0*(P1+P2);
    C(5, 1, 7) = -I*P0*(-P1+P2);
    C(5, 2, 7) = 0.25*d*D0*(-C1m2*D4+2.0*C7m2*D2h2+3.0*C13m2*h4);
    P0 = P0*D1_h1; P1 = (-D4+2.0*D2h2-5.0*h4); P2 = D02/3.0*(-7.0*D4+22.0*D2h2+5.0*h4);
    C(5, 0, 8) = P0*(P1+P2);
    C(5, 1, 8) = I*P0*(P1-P2);
    C(5, 2, 8) = -0.5*d*D0*h2*D1_h1*D1mh3;

    return C;
}

il::StaticArray3D<std::complex<double>, 6, 3, 9> S13_23H(double nu, std::complex<double> eix, double h, std::complex<double> d) {
    const std::complex<double> I(0.0,1.0);

    double C1_1 = 1.0+nu;
    double C2m1 = 2.0-nu;
    double C3_1 = 3.0+nu;
    double C3m1 = 3.0-nu;
    double C5m1 = 5.0-nu;
    double C6_1 = 6.0+nu;
    double C12_1 = 12.0+nu;

    double CosX = std::real(eix);
    double TanX = std::imag(eix)/CosX;
    std::complex<double> CTanX8_3i = 8.0+3.0*I*TanX;
    std::complex<double> TCosX = CosX*eix;
    std::complex<double> TCosXC = std::conj(TCosX);
    std::complex<double> TCosX1m1 = TCosX-1.0;
    std::complex<double> CTCosX3_4 = 3.0+4.0*TCosX;
    std::complex<double> CTCosX5_8 = 5.0+8.0*TCosX;
    std::complex<double> CTCosX1 = TCosX1m1*TCosX;
//  std::complex<double> e2X = eix*eix;
//  std::complex<double> WCTCosX1 = (13.0+e2X-10.0*TCosX)*TCosX;

    double h2 = h*h; double h3 = h2*h; double h4 = h2*h2;
//  double sgh = (( h < 0 )? -1.0 : double(( h > 0 ))); // sign(h)

    double D1 = std::abs(d); double D2 = D1*D1; double D4 = D2*D2;
    std::complex<double> dC = std::conj(d);
    std::complex<double> D0 = std::polar(1.0,std::arg(d)); // =d/D1
    std::complex<double> D0C = std::conj(D0);
    std::complex<double> D02 = D0*D0; // =d^2/D1^2
    std::complex<double> D02C = std::conj(D02);
    std::complex<double> D03 = D0*D02; // =d^3/D1^3
    std::complex<double> D04 = D02*D02; // =d^4/D1^4

    std::complex<double> D2h2 = D2*h2;

    std::complex<double> D1_h1 = D2+h2;
//  std::complex<double> D1_h3 = D2+3.0*h2;
//  std::complex<double> D3_h1 = 3.0*D2+h2;
    std::complex<double> D1mh3 = D2-3.0*h2;
//  std::complex<double> D3mh1 = 3.0*D2-h2;

    std::complex<double> P0, P1, P2;

    il::StaticArray3D<std::complex<double>, 6, 3, 9> C{0.0};

    C(0, 1, 1) = -0.5*nu*D02*TCosX;
    C(0, 0, 1) = I*C(0, 1, 1);
    C(0, 1, 2) = 0.5*D02*h2*TCosX;
    C(0, 0, 2) = I*C(0, 1, 2);
    C(0, 0, 6) = -0.5*(C2m1*D1+nu*d*D0);
    C(0, 1, 6) = 0.5*I*(-C2m1*D1+nu*d*D0);
    C(0, 2, 6) = -D0*h;
    C(0, 0, 7) = -h2*D1*(D02+1.0);
    C(0, 1, 7) = I*h2*D1*(D02-1.0);
    C(0, 2, 7) = -2.0*h3*D0;

    C(1, 1, 1) = -0.0625*d*D02*nu*CTCosX5_8;
    C(1, 0, 1) = I*C(1, 1, 1);
    C(1, 2, 1) = -I*D02*h*TCosX;
    C(1, 1, 2) = d*D02*(0.03125*nu*D2*CTanX8_3i+h2*(0.5+0.09375*nu+0.125*(6.0+nu)*TCosX));
    C(1, 0, 2) = I*C(1, 1, 2);
    C(1, 2, 2) = I*D02*h3*TCosX;
    C(1, 1, 3) = -0.09375*d*D02*h2*(D2*CTanX8_3i+h2*CTCosX3_4);
    C(1, 0, 3) = I*C(1, 1, 3);
    C(1, 0, 5) = 0.25*D0*(C2m1-nu*D02);
    C(1, 1, 5) = 0.25*I*D0*(C2m1+nu*D02);
    P2 = 3.0*nu*D2+C3_1*h2;
    C(1, 0, 6) = 0.5*D0*(C5m1*h2-D02*P2);
    C(1, 1, 6) = 0.5*I*D0*(C5m1*h2+D02*P2);
    C(1, 2, 6) = -d*D0*h;
    P2 = D02*(0.625*nu*D4+0.75*C6_1*D2h2+0.125*C12_1*h4);
    C(1, 0, 7) = D0*(h4-P2);
    C(1, 1, 7) = I*D0*(h4+P2);
    C(1, 2, 7) = -2.0*d*D0*h3;
    C(1, 0, 8) = -0.25*D03*h2*D1_h1*(5.0*D2+h2);
    C(1, 1, 8) = -I*C(1, 0, 8);

    C(2, 1, 1) = 0.3125*nu*d;
    C(2, 0, 1) = I*C(2, 1, 1);
    C(2, 1, 2) = -0.03125*d*(h2*(16+3.0*nu)+I*nu*D2*TanX);
    C(2, 0, 2) = I*C(2, 1, 2);
    C(2, 2, 2) = D2/12.0*h*TanX;
    C(2, 1, 3) = 0.09375*d*h2*(3.0*h2+I*D2*TanX);
    C(2, 0, 3) = I*C(2, 1, 3);
    C(2, 2, 3) = -0.25*D2*h3*TanX;
    P1 = 3.0*nu*D0; P2 = C2m1*D0C;
    C(2, 0, 5) = 0.25*(P1+P2);
    C(2, 1, 5) = -0.25*I*(P1-P2);
    P1 = nu*D2+3.0*C3_1*h2; P2 = C5m1*h2;
    C(2, 0, 6) = 0.5*(P1*D0+P2*D0C);
    C(2, 1, 6) = -0.5*I*(P1*D0-P2*D0C);
    C(2, 2, 6) = -10.0/3.0*D1*h;
    P1 = 0.125*(nu*D4-2.0*C6_1*D2h2-3.0*C12_1*h4);
    C(2, 0, 7) = -P1*D0+h4*D0C;
    C(2, 1, 7) = I*(P1*D0+h4*D0C);
    C(2, 2, 7) = -D1*h*(D2+11.0*h2)/3.0;
    C(2, 0, 8) = -0.25*D0*h2*D1mh3*D1_h1;
    C(2, 1, 8) = -I*C(2, 0, 8);
    C(2, 2, 8) = -2.0/3.0*D1*h3*D1_h1;

    P1 = 0.5*C2m1*TCosX; P2 = 3.2/3.0*nu*D02*CTCosX1;
    C(3, 0, 0) = I*D02*(P1+P2);
    C(3, 1, 0) = D02*(-P1+P2);
    P1 = 0.5*C5m1*h2*TCosX; P2 = D02*(-3.2/3.0*C3_1*h2*CTCosX1+nu*D2/3.0*(2.6+0.3625*I*TanX+1.6*TCosX));
    C(3, 0, 1) = I*D02*(P1-P2);
    C(3, 1, 1) = -D02*(P1+P2);
    C(3, 2, 1) = -0.125*I*d*D02*h*CTCosX5_8;
    P1 = D02*(D4*nu*(0.4+0.18125*I*TanX)-h4*0.4/3.0*C12_1*CTCosX1+
              D2h2*(1.4+0.8/3.0*nu+0.4/3.0*C6_1*TCosX+I*(0.2+0.11875/3.0*nu)*TanX));
    P2 = 0.5*TCosX*h4;
    C(3, 0, 2) = I*D02*(P1-P2);
    C(3, 1, 2) = D02*(P1+P2);
    C(3, 2, 2) = 0.0625*I*d*D02*h*(D2*CTanX8_3i+h2*(19.0+28.0*TCosX));
    P0 = D04*h2; //P1 = ; P2 = ;
    C(3, 0, 3) = P0*(D4*(-1.2*I+0.54375*TanX)+0.4*I*h4*CTCosX1+D2h2*(0.11875*TanX-0.4*I*(2.0+TCosX)));
    C(3, 1, 3) = P0*(D4*(-1.2-0.54375*I*TanX)+0.4*h4*CTCosX1+D2h2*(-0.11875*I*TanX-0.4*(2.0+TCosX)));
    C(3, 2, 3) = -0.1875*I*d*D02*h3*(D2*CTanX8_3i+h2*CTCosX3_4);
    P1 = 0.25*C2m1*d*D0; P2 = 0.75*nu*d*D03;
    C(3, 0, 5) = P1-P2;
    C(3, 1, 5) = I*(P1+P2);
    C(3, 2, 5) = -0.5*D03*h;
    P1 = 0.5*C5m1*d*D0*h2; P2 = 0.5*d*D03*(5.0*nu*D2+3.0*C3_1*h2);
    C(3, 0, 6) = P1-P2;
    C(3, 1, 6) = I*(P1+P2);
    C(3, 2, 6) = -D03*h*(3.0*D2+4.0*h2);
    P1 = d*D0*h4; P2 = 0.125*d*D03*(7.0*nu*D4+10.0*C6_1*D2h2+3.0*C12_1*h4);
    C(3, 0, 7) = P1-P2;
    C(3, 1, 7) = I*(P1+P2);
    C(3, 2, 7) = -0.25*D03*h*(5.0*D4+42.0*D2h2+13.0*h4);
    C(3, 0, 8) = -0.25*d*D03*h2*D1_h1*(7.0*D2+3.0*h2);
    C(3, 1, 8) = -I*C(3, 0, 8);
    C(3, 2, 8) = -0.5*D03*h3*D1_h1*(5.0*D2+h2);

    C(4, 1, 0) = 0.5*C2m1*D02C*TCosXC;
    C(4, 0, 0) = -I*C(4, 1, 0);
    P1 = 0.3625/3.0*nu*D2*TanX; P2 = 0.5*C5m1*h2*D02C*TCosXC;
    C(4, 0, 1) = P1-I*P2;
    C(4, 1, 1) = -I*P1+P2;
    P1 = (0.01875*nu*D2-h2*(0.2+0.11875/3.0*nu))*D2*TanX; P2 = 0.5*h4*D02C*TCosXC;
    C(4, 0, 2) = P1+I*P2;
    C(4, 1, 2) = -I*P1-P2;
    C(4, 0, 3) = D2*h2*(-0.05625*D2+0.11875*h2)*TanX;
    C(4, 1, 3) = -I*C(4, 0, 3);
    C(4, 0, 4) = -2.0*C1_1*std::fabs(h);
    C(4, 1, 4) = -I*C(4, 0, 4);
    P1 = 0.75*nu*D1;  P2 = 0.25*C2m1*dC*D0C;
    C(4, 0, 5) = -P1+P2;
    C(4, 1, 5) = I*(P1+P2);
    P1 = D1*(0.5/3.0*nu*D2+h2*(4.3+0.9*nu)); P2 = 0.5*C5m1*h2*dC*D0C;
    C(4, 0, 6) = P1+P2;
    C(4, 1, 6) = -I*(P1-P2);
    P1 = D1*(-0.075*nu*D4+0.25/3.0*C6_1*D2h2+h4/3.0*(7.3+0.475*nu)); P2 = h4*dC*D0C;
    C(4, 0, 7) = P1+P2;
    C(4, 1, 7) = -I*(P1-P2);
    C(4, 0, 8) = D1*h2*D1_h1*(-0.15*D2+0.95/3.0*h2);
    C(4, 1, 8) = -I*C(4, 0, 8);

    C(5, 1, 0) = 3.2/3.0*nu*D02*TCosX;
    C(5, 0, 0) = I*C(5, 1, 0);
    P1 = 0.125/3.0*C2m1*D2*TanX; P2 = D02/3.0*(nu*D2*(1.3+0.3625*I*TanX)+3.2*C3_1*h2*TCosX);
    C(5, 0, 1) = P1+I*P2;
    C(5, 1, 1) = I*P1+P2;
    P1 = 0.125/3.0*C5m1*D2*h2*TanX;
    P2 = D02*(0.03125*nu*D4*TanX+D2h2/3.0*(-I*(2.1+0.4*nu)+(0.6+0.11875*nu)*TanX)-0.4/3.0*I*h4*C12_1*TCosX);
    C(5, 0, 2) = -P1+P2;
    C(5, 1, 2) = -I*(P1+P2);
    P1 = 0.125*D2*h4*TanX; P2 = D02*h2*(0.4*TCosX*h4+(0.4+0.11875*I*TanX)*D2h2+0.09375*I*TanX*D4);
    C(5, 0, 3) = P1+I*P2;
    C(5, 1, 3) = I*P1+P2;
    C(5, 0, 4) = -C3m1*std::fabs(h);
    C(5, 1, 4) = I*C(5, 0, 4);
    P1 = 0.25*C2m1*D1; P2 = 0.75*nu*d*D0;
    C(5, 0, 5) = -P1+P2;
    C(5, 1, 5) = -I*(P1+P2);
    P1 = D1*(0.75*(6.0-nu)*h2+0.25/3.0*C2m1*D2); P2 = 0.5*d*D0*(nu*D2+3.0*C3_1*h2);
    C(5, 0, 6) = P1+P2;
    C(5, 1, 6) = I*(P1-P2);
    P1 = 1.0/6.0*D1*h2*((15.0-nu)*h2+C5m1*D2); P2 = 0.125*d*D0*(3.0*C12_1*h4+2.0*C6_1*D2h2-nu*D4);
    C(5, 0, 7) = P1+P2;
    C(5, 1, 7) = I*(P1-P2);
    P0 = h2*D1_h1; P1 = 1.0/3.0*h2*D1; P2 = 0.25*d*D0*D1mh3;
    C(5, 0, 8) = P0*(P1-P2);
    C(5, 1, 8) = I*P0*(P1+P2);

    C(5, 2, 1) = 0.625*I*h*d;
    C(5, 2, 2) = 0.0625*h*d*(-19.0*I*h2+D2*TanX);
    C(5, 2, 3) = 0.1875*h3*d*(3.0*I*h2-D2*TanX);
    C(5, 2, 5) = 1.5*h*D0;
    C(5, 2, 6) = (D2+12.0*h2)*h*D0;
    C(5, 2, 7) = 0.25*(-D4+14.0*D2h2+39.0*h4)*h*D0;
    C(5, 2, 8) = -0.5*D0*h3*D1_h1*D1mh3;

    for (int j = 0; j < C.size(2); ++j) {
        C(4, 2, j) = std::conj(C(5, 2, j));
    }

    return C;
}

il::StaticArray3D<std::complex<double>, 6, 3, 9> S33H(double nu, std::complex<double> eix, double h, std::complex<double> d) {
    const std::complex<double> I(0.0,1.0);

    double CosX = std::real(eix);
    double TanX = std::imag(eix)/CosX;
    std::complex<double> TCosX = CosX*eix;
    std::complex<double> CTCosX4_3 = 3.0+4.0*TCosX;
    std::complex<double> CTanX8_3i = 8.0+3.0*I*TanX;

    double h2 = h*h; double h3 = h*h2; double h4 = h2*h2;
//  double sgh = (( h < 0 )? -1.0 : double(( h > 0 ))); // sign(h)

    double D1 = std::abs(d); double D2 = D1*D1; double D4 = D2*D2;
    // std::complex<double> dC = std::conj(d);
    std::complex<double> D0 = std::polar(1.0,std::arg(d)); //  = d/D1
    std::complex<double> D02 = D0*D0; //  = d^2/D1^2

    std::complex<double> D1_h1 = D2+h2;
//  std::complex<double> D1_h3 = D2+3.0*h2;
    std::complex<double> D1mh3 = D2-3.0*h2;

    std::complex<double> P0, P1, P2;

    il::StaticArray3D<std::complex<double>, 6, 3, 9> C{0.0};

    C(0, 0, 6) = -2.0*h*std::real(D0);
    C(0, 1, 6) = -2.0*h*std::imag(D0);
    C(0, 2, 6) = -2.0*D1;
    C(0, 0, 7) = -4.0*h3*std::real(D0);
    C(0, 1, 7) = -4.0*h3*std::imag(D0);
    C(0, 2, 7) = 4.0*h2*D1;

    C(1, 1, 1) = -D02*h*TCosX;
    C(1, 0, 1) = I*C(1, 1, 1);
    P1 = 1.0/12.0*TanX*D2; P2 = D02*h2*TCosX;
    C(1, 0, 2) = (P1+I*P2)*h;
    C(1, 1, 2) = (I*P1+P2)*h;
    C(1, 0, 3) = -0.25*D2*h3*TanX;
    C(1, 1, 3) = I*C(1, 0, 3);
    C(1, 2, 5) = D0;
    P2 = 10.0/3.0*D1;
    C(1, 0, 6) = (-d*D0-P2)*h;
    C(1, 1, 6) = I*(d*D0-P2)*h;
    C(1, 2, 6) = -4.0*h2*D0;
    P0 = D1*h; P1 = (D2+11.0*h2)/3.0; P2 = 2.0*D02*h2;
    C(1, 0, 7) = -(P1+P2)*P0;
    C(1, 1, 7) = -I*(P1-P2)*P0;
    C(1, 2, 7) = -4.0*h4*D0;
    C(1, 0, 8) = -2.0/3.0*h3*D1*D1_h1;
    C(1, 1, 8) = I*C(1, 0, 8);

    for (int k = 0; k < C.size(2); ++k) {
        for (int j = 0; j < C.size(1); ++j) {
            C(2, j, k) = std::conj(C(1, j, k));
        }
    }

    C(3, 2, 0) = 2.0*I*D02*TCosX;
    P0 = d*h; P2 = D02*(0.625+TCosX);
    C(3, 0, 1) = I*P0*(0.625-P2);
    C(3, 1, 1) = -P0*(0.625+P2);
    C(3, 2, 1) = -4.0*I*D02*h2*TCosX;
    P0 = 0.0625*d*h; P1 = -19.0*I*h2+D2*TanX; P2 = D02*(D2*CTanX8_3i+h2*(19.0+28.0*TCosX));
    C(3, 0, 2) = P0*(P1+I*P2);
    C(3, 1, 2) = P0*(I*P1+P2);
    C(3, 2, 2) = 2.0*I*D02*h4*TCosX;
    P0 = 0.1875*d*h3; P1 = 3.0*I*h2-D2*TanX; P2 = D02*(D2*CTanX8_3i+h2*CTCosX4_3);
    C(3, 0, 3) = P0*(P1-I*P2);
    C(3, 1, 3) = P0*(I*P1-P2);
    P0 = h*D0;
    C(3, 0, 5) = P0*(-0.5*D02+1.5);
    C(3, 1, 5) = I*P0*(0.5*D02+1.5);
    C(3, 2, 5) = d*D0;
    P1 = D2+12.0*h2; P2 = D02*(3.0*D2+4.0*h2);
    C(3, 0, 6) = P0*(P1-P2);
    C(3, 1, 6) = I*P0*(P1+P2);
    C(3, 2, 6) = -4.0*h2*d*D0;
    P1 = -0.25*D4+3.5*D2*h2+9.75*h4; P2 = D02*(1.25*D4+10.5*D2*h2+3.25*h4);
    C(3, 0, 7) = P0*(P1-P2);
    C(3, 1, 7) = I*P0*(P1+P2);
    C(3, 2, 7) = -4.0*h4*d*D0;
    P0 = h3*D0*D1_h1; P1 = 0.5*D1mh3; P2 = D02*(2.5*D2+0.5*h2);
    C(3, 0, 8) = -P0*(P1+P2);
    C(3, 1, 8) = I*P0*(-P1+P2);

    for (int k = 0; k < C.size(2); ++k) {
        for (int j = 0; j < C.size(1); ++j) {
            C(4, j, k) = std::conj(C(3, j, k));
        }
    }

    C(5, 0, 1) = -1.25*h*std::imag(d);
    C(5, 1, 1) = 1.25*h*std::real(d);
    C(5, 2, 1) = 1.0/6.0*D2*TanX;
    P1 = 0.125*D2*TanX; P2 = 2.375*h2;
    C(5, 0, 2) = h*(P1*std::real(d)+P2*std::imag(d));
    C(5, 1, 2) = h*(P1*std::imag(d)-P2*std::real(d));
    C(5, 2, 2) = 1.0/3.0*D2*h2*TanX;
    P1 = 3.0*P1; P2 = 1.125*h2;
    C(5, 0, 3) = -h3*(P1*std::real(d)+P2*std::imag(d));
    C(5, 1, 3) = -h3*(P1*std::imag(d)-P2*std::real(d));
    C(5, 2, 3) = -0.5*D2*h4*TanX;
    C(5, 0, 5) = 3.0*h*std::real(D0);
    C(5, 1, 5) = 3.0*h*std::imag(D0);
    C(5, 2, 5) = -D1;
    P0 = 2.0*h*(D2+12.0*h2);
    C(5, 0, 6) = P0*std::real(D0);
    C(5, 1, 6) = P0*std::imag(D0);
    C(5, 2, 6) = (1.0/3.0*D2-9.0*h2)*D1;
    P0 = h*(-0.5*D4+7.0*D2*h2+19.5*h4);
    C(5, 0, 7) = P0*std::real(D0);
    C(5, 1, 7) = P0*std::imag(D0);
    C(5, 2, 7) = -h2*(4.0/3.0*D2+8.0*h2)*D1;
    P0 = h3*D1_h1*D1mh3;
    C(5, 0, 8) = -P0*std::real(D0);
    C(5, 1, 8) = -P0*std::imag(D0);
    C(5, 2, 8) = -4.0/3.0*h4*D1_h1*D1;

    return C;
}

// Special case (reduced summation, collocation point projected onto the element contour) - additional terms

il::StaticArray3D<std::complex<double>, 6, 3, 5> S11_22H_red(double nu, std::complex<double> eix, double h, std::complex<double> d) {
    const std::complex<double> I(0.0,1.0);

    double C1_1 = 1.0+nu;
    double C1_2 = 1.0+2.0*nu;
    double C2_1 = 2.0+nu;
    double C3_2 = 3.0+2.0*nu;
    double C7_2 = 7.0+2.0*nu;
    double C11_4 = 11.0+4.0*nu;
    double C13_2 = 13.0+2.0*nu;

    double CosX = std::real(eix);
    double SinX = std::imag(eix);
    std::complex<double> e2x = eix*eix;
    std::complex<double> Ce1x3_1 = eix*(3.0+e2x);
    std::complex<double> Ce1x3m1 = eix*(3.0-e2x);

    double h2 = h*h; double h3 = h2*h; double h4 = h2*h2; double h5 = h4*h; double h7 = h5*h2;
    double sgh = (( h < 0 )? -1.0 : static_cast<double>(( h > 0 ))); // sign(h)

    std::complex<double> P1, P2, P3, P4;

    il::StaticArray3D<std::complex<double>, 6, 3, 5> C{0.0};

    C(0, 0, 2) = -C7_2*h*SinX; C(0, 1, 2) = C7_2*h*CosX;
    C(0, 0, 3) = -9.0*h3*SinX; C(0, 1, 3) = 9.0*h3*CosX;
    C(0, 0, 4) = -1.5*h5*SinX; C(0, 1, 4) = 1.5*h5*CosX;

    C(1, 0, 0) = -0.5*I*C1_1*e2x*sgh; C(1, 1, 0) = -0.5*C1_1*e2x*sgh;
    C(1, 2, 1) = 0.5*I*C1_2*eix;
    C(1, 2, 2) = 2.0*I*C2_1*h2*eix;
    C(1, 2, 3) = 2.0*I*h4*eix;

    for (int k = 0; k < C.size(2); ++k) {
        for (int j = 0; j < C.size(1); ++j) {
            C(2, j, k) = std::conj(C(1, j, k));
        }
    }

    P1 = 0.25*C3_2*h; P2 = 0.25*C11_4*h3; P3 = 0.125*C13_2*h5; P4 = 0.625*h7;

    C(3, 2, 0) = -2.0*I*C1_1*e2x*std::fabs(h);
    C(3, 0, 1) = -I*P1*Ce1x3_1; C(3, 1, 1) = P1*Ce1x3m1;
    C(3, 0, 2) = -I*P2*Ce1x3_1; C(3, 1, 2) = P2*Ce1x3m1;
    C(3, 0, 3) = -I*P3*Ce1x3_1; C(3, 1, 3) = P3*Ce1x3m1;
    C(3, 0, 4) = -I*P4/3.0*Ce1x3_1; C(3, 1, 4) = P4/3.0*Ce1x3m1;

    for (int k = 0; k < C.size(2); ++k) {
        for (int j = 0; j < C.size(1); ++j) {
            C(4, j, k) = std::conj(C(3, j, k));
        }
    }

    C(5, 0, 1) = 6.0*P1*SinX; C(5, 1, 1) = -6.0*P1*CosX;
    C(5, 0, 2) = 6.0*P2*SinX; C(5, 1, 2) = -6.0*P2*CosX;
    C(5, 0, 3) = 6.0*P3*SinX; C(5, 1, 3) = -6.0*P3*CosX;
    C(5, 0, 4) = 2.0*P4*SinX; C(5, 1, 4) = -2.0*P4*CosX;

    return C;
}

il::StaticArray3D<std::complex<double>, 6, 3, 5> S11_22_12H_red(double nu, std::complex<double> eix, double h, std::complex<double> d) {
    const std::complex<double> I(0.0,1.0);

    double C1m1 = 1.0-nu;
    double C1m2 = 1.0-2.0*nu;
    double C2m1 = 2.0-nu;
    double C5m4 = 5.0-4.0*nu;
    double C9m4 = 9.0-4.0*nu;
    double C13m2 = 13.0-2.0*nu;
    double C15m4 = 15.0-4.0*nu;
    double C15m8 = 15.0-8.0*nu;

    double CosX = std::real(eix);
    double SinX = std::imag(eix);
    std::complex<double> e2x = eix*eix;
    std::complex<double> emx = std::conj(eix);
    //std::complex<double> Ce2x3_1 = 3.0+e2x;
    //std::complex<double> Ce2x3m1 = 3.0-e2x;
    std::complex<double> Ce1x3_1 = eix*(3.0+e2x);
    std::complex<double> Ce1x3m1 = eix*(3.0-e2x);

    double h2 = h*h; double h3 = h2*h; double h4 = h2*h2;
    double h5 = h4*h; double h7 = h5*h2;
    double sgh = (( h < 0 )? -1.0 : static_cast<double>(( h > 0 ))); // sign(h)
    double abh = std::fabs(h);

    std::complex<double> P1, P2, P3, P4;

    il::StaticArray3D<std::complex<double>, 6, 3, 5> C{0.0};

    C(0, 2, 0) = -I*nu*e2x/abh;
    C(0, 0, 2) = 0.5*I*h*(3.0*Ce1x3_1-4.0*nu*eix); C(0, 1, 2) = -0.5*h*(3.0*Ce1x3m1-4.0*nu*eix);
    C(0, 0, 3) = 1.5*I*h3*Ce1x3_1; C(0, 1, 3) = -1.5*h3*Ce1x3m1;
    C(0, 0, 4) = 0.25*I*h5*Ce1x3_1; C(0, 1, 4) = -0.25*h5*Ce1x3m1;

    P1 = 0.5*I*C1m2*eix;
    P2 = 2.0*I*C2m1*h2*eix;
    P3 = 0.25*I*C13m2*h4*eix;
    P4 = 0.5*I*h2*h4*eix;

    C(1, 0, 0) = -I*sgh*e2x*(C1m1+0.5*e2x); C(1, 1, 0) = sgh*e2x*(C1m1-0.5*e2x);
    C(1, 2, 1) = P1*e2x;
    C(1, 2, 2) = P2*e2x;
    C(1, 2, 3) = P3*e2x;
    C(1, 2, 4) = P4*e2x;

    C(2, 1, 0) = -sgh*e2x; C(2, 0, 0) = I*C(2, 1, 0);
    C(2, 2, 1) = 3.0*P1;
    C(2, 2, 2) = 3.0*P2;
    C(2, 2, 3) = 3.0*P3;
    C(2, 2, 4) = 3.0*P4;

    P1 = 0.25*h*eix;
    P2 = 0.25*h3*eix;
    P3 = 0.125*h5*eix;
    P4 = 0.125*h7*eix;

    C(3, 2, 0) = -2.0*I*C1m1*abh*e2x*e2x;
    C(3, 0, 1) = -I*e2x*P1*(C5m4+3.0*e2x);
    C(3, 0, 2) = -I*e2x*P2*(C15m8+9.0*e2x);
    C(3, 0, 3) = -I*e2x*P3*(C15m4+9.0*e2x);
    C(3, 0, 4) = -I*e2x*P4*(5.0/3.0+e2x);

    C(3, 1, 1) = e2x*P1*(C5m4-3.0*e2x);
    C(3, 1, 2) = e2x*P2*(C15m8-9.0*e2x);
    C(3, 1, 3) = e2x*P3*(C15m4-9.0*e2x);
    C(3, 1, 4) = e2x*P4*(5.0/3.0-e2x);

    C(5, 2, 0) = -4.0*I*C1m1*abh*e2x;
    C(5, 0, 1) = -I*P1*(3.0*C5m4+5.0*e2x);
    C(5, 0, 2) = -3.0*I*P2*(C15m8+5.0*e2x);
    C(5, 0, 3) = -3.0*I*P3*(C15m4+5.0*e2x);
    C(5, 0, 4) = -5.0*I*P4*(1.0+e2x/3.0);

    C(5, 1, 1) = P1*(3.0*C5m4-5.0*e2x);
    C(5, 1, 2) = 3.0*P2*(C15m8-5.0*e2x);
    C(5, 1, 3) = 3.0*P3*(C15m4-5.0*e2x);
    C(5, 1, 4) = 5.0*P4*(1.0-e2x/3.0);

    C(4, 0, 1) = 3.0*I*std::conj(P1)*(C5m4-5.0*e2x);
    //C(4, 0, 1) = -0.75*I*h*(5.0*eix-C5m4*emx);
    C(4, 1, 1) = -3.0*std::conj(P1)*(C5m4+5.0*e2x);
    //C(4, 1, 1) = -0.75*h*(5.0*eix+C5m4*emx);
    C(4, 0, 2) = 3.0*I*std::conj(P2)*(C15m8-15.0*e2x);
    //C(4, 0, 2) = -0.75*I*h3*(15.0*eix-C15m8*emx);
    C(4, 1, 2) = -3.0*std::conj(P2)*(C15m8+15.0*e2x);
    //C(4, 1, 2) = -0.75*h3*(15.0*eix+C15m8*emx);
    C(4, 0, 3) = 3.0*I*std::conj(P3)*(C15m4-15.0*e2x);
    //C(4, 0, 3) = -0.375*I*h5*(15.0*eix-C15m4*emx);
    C(4, 1, 3) = -3.0*std::conj(P3)*(C15m4+15.0*e2x);
    //C(4, 1, 3) = -0.375*h5*(15.0*eix+C15m4*emx);
    C(4, 0, 4) = 1.25*h7*SinX; C(4, 1, 4) = -1.25*h7*CosX;

    return C;
}

il::StaticArray3D<std::complex<double>, 6, 3, 5> S13_23H_red(double nu, std::complex<double> eix, double h, std::complex<double> d) {
    const std::complex<double> I(0.0,1.0);

    double C1_1 = 1.0+nu;
    double C1m1 = 1.0-nu;
    double C2m1 = 2.0-nu;
    double C3_1 = 3.0+nu;
    double C3m1 = 3.0-nu;
    double C5m1 = 5.0-nu;
    double C12_1 = 12.0+nu;

    // double CosX = std::real(eix);
    // double SinX = std::imag(eix);
    std::complex<double> e2x = eix*eix;
    std::complex<double> e3x = e2x*eix;
    std::complex<double> emx = std::conj(eix);
    std::complex<double> em2 = std::conj(e2x);
    std::complex<double> Ce2x3 = 3.0*emx+e3x;
    std::complex<double> Ce1x3_1 = eix*(3.0+e2x);
    std::complex<double> Ce1x3m1 = eix*(3.0-e2x);

    double h2 = h*h; double h3 = h2*h; double h4 = h2*h2;
    double h5 = h4*h; double h6 = h4*h2; double h7 = h5*h2;
    // double sgh = (( h < 0 )? -1.0 : static_cast<double>(( h > 0 ))); // sign(h)
    double abh = std::fabs(h);

    std::complex<double> P1, P2, P3, P4;

    il::StaticArray3D<std::complex<double>, 6, 3, 5> C{0.0};

    C(0, 1, 0) = -0.25*C1m1*e2x/abh; C(0, 0, 0) = I*C(0, 1, 0);
    C(0, 2, 2) = -I*h*eix; C(0, 2, 3) = -2.0*I*h3*eix;

    P1 = 0.25*nu*Ce2x3;
    P2 = 0.5*h2*C3_1*Ce2x3;
    P3 = 0.125*h4*C12_1*Ce2x3;
    // P4 = 0.25*h6*Ce2x3;

    C(1, 0, 1) = 0.25*I*eix*(C2m1+nu*e2x); C(1, 1, 1) = -0.25*eix*(C2m1-nu*e2x);
    C(1, 0, 2) = 0.5*I*h2*eix*(C5m1+C3_1*e2x); C(1, 1, 2) = -0.5*h2*eix*(C5m1-C3_1*e2x);
    C(1, 0, 3) = 0.125*I*h4*eix*(8.0+C12_1*e2x); C(1 ,1 ,3) = -0.125*h4*eix*(8.0-C12_1*e2x);
    C(1, 1, 4) = 0.25*h6*e3x; C(1, 0, 4) = I*C(1 ,1 ,4);

    C(2, 0, 1) = std::conj(C(1, 0, 1)-I*P1); C(2, 1, 1) = std::conj(P1-C(1, 1, 1));
    C(2, 0, 2) = std::conj(C(1, 0, 2)-I*P2); C(2, 1, 2) = std::conj(P2-C(1, 1, 2));
    C(2, 0, 3) = std::conj(C(1, 0, 3)-I*P3); C(2, 1, 3) = std::conj(P3-C(1, 1, 3));
    C(2, 1, 4) = 0.75*eix*h6; C(2, 0, 4) = I*C(2, 1, 4);

    P1 = 0.5*I*h*eix;
    P2 = 4.0*I*h3*eix;
    P3 = 3.25*I*h5*eix;
    P4 = 0.5*I*h7*eix;

    C(5, 1, 0) = -C1_1*abh*e2x; C(5, 0, 0) = I*C(5, 1, 0);
    C(5, 2, 1) = 3.0*P1;
    C(5, 2, 2) = 3.0*P2;
    C(5, 2, 3) = 3.0*P3;
    C(5, 2, 4) = 3.0*P4;

    C(3, 1, 0) = 0.5*(C3m1-C1_1*e2x)*abh*e2x; C(3, 0, 0) = -0.5*I*(C3m1+C1_1*e2x)*abh*e2x;
    C(3, 2, 1) = e2x*P1;
    C(3, 2, 2) = e2x*P2;
    C(3, 2, 3) = e2x*P3;
    C(3, 2, 4) = e2x*P4;

    C(4, 1, 0) = -0.5*C3m1*abh*em2; C(4, 0, 0) = -I*C(4, 1, 0);
    C(4, 2, 1) = std::conj(C(5, 2, 1));
    C(4, 2, 2) = std::conj(C(5, 2, 2));
    C(4, 2, 3) = std::conj(C(5, 2, 3));
    C(4, 2, 4) = std::conj(C(5, 2, 4));

    return C;
}

il::StaticArray3D<std::complex<double>, 6, 3, 5> S33H_red(double nu, std::complex<double> eix, double h, std::complex<double> d) {
    const std::complex<double> I(0.0,1.0);

    double CosX = std::real(eix);
    double SinX = std::imag(eix);
    std::complex<double> e2x = eix*eix;
    std::complex<double> Ce1x3_1 = eix*(3.0+e2x);
    std::complex<double> Ce1x3m1 = eix*(3.0-e2x);

    double h2 = h*h; double h3 = h2*h; double h4 = h2*h2;
    double h5 = h4*h; double h7 = h5*h2;
    // double sgh = (( h < 0 )? -1.0 : double(( h > 0 ))); // sign(h)

    il::StaticArray3D<std::complex<double>, 6, 3, 5> C{0.0};

    C(0, 0, 2) = 2.0*h*SinX; C(0, 1, 2) = -2.0*h*CosX;
    C(0, 0, 3) = 4.0*h3*SinX; C(0, 1, 3) = -4.0*h3*CosX;

    C(1, 2, 1) = I*eix;
    C(1, 2, 2) = -4.0*I*h2*eix;
    C(1, 2, 3) = -4.0*I*h4*eix;

    for (int j = 1; j < C.size(2)-1; ++j) {
        C(2, 2, j) = std::conj(C(1, 2, j));
    }

    C(3, 0, 1) = 0.5*I*h*Ce1x3_1; C(3, 1, 1) = -0.5*h*Ce1x3m1;
    C(3, 0, 2) = 4.0*I*h3*Ce1x3_1; C(3, 1, 2) = -4.0*h3*Ce1x3m1;
    C(3, 0, 3) = 3.25*I*h5*Ce1x3_1; C(3, 1, 3) = -3.25*h5*Ce1x3m1;
    C(3, 0, 4) = 0.5*I*h7*Ce1x3_1; C(3, 1, 4) = -0.5*h7*Ce1x3m1;

    for (int k = 1; k < C.size(2); ++k) {
        for (int j = 0; j < C.size(1)-1; ++j) {
            C(4, j, k) = std::conj(C(3, j, k));
        }
    }

    C(5, 0, 1) = -3.0*h*SinX; C(5, 1, 1) = 3.0*h*CosX;
    C(5, 0, 2) = -24.0*h3*SinX; C(5, 1, 2) = 24.0*h3*CosX;
    C(5, 0, 3) = -19.5*h5*SinX; C(5, 1, 3) = 19.5*h5*CosX;
    C(5, 0, 4) = -3.0*h7*SinX; C(5, 1, 4) = 3.0*h7*CosX;

    return C;
}

// Limit case (h==0, plane) - all stress components

il::StaticArray3D<std::complex<double>, 6, 4, 3> SijLimH(double nu, std::complex<double> eix, std::complex<double> d) {
    const std::complex<double> I(0.0,1.0);

    double C1_2 = 1.0+2.0*nu;
    double C1m2 = 1.0-2.0*nu;
    double C2m1 = 2.0-nu;

    double CosX = std::real(eix);
    double SinX = std::imag(eix);
    // double TanX = SinX/CosX;
    std::complex<double> e2x = eix*eix;
    double H0Lim = std::atanh(SinX); // std::complex<double> ?
    // double H0Lim = 0.5*(std::log(1.0+SinX)-std::log(1.0-SinX)); // std::complex<double> ?

    double D1 = std::abs(d); // double D2 = D1*D1; double D4 = D2*D2;
    std::complex<double> D0 = std::polar(1.0,std::arg(d)); // = d/D1
    std::complex<double> D02 = D0*D0; //  = d^2/D1^2
    std::complex<double> D03 = D0*D02; //  = d^3/D1^3
    std::complex<double> D04 = D02*D02; //  = d^4/D1^4

    il::StaticArray3D<std::complex<double>, 6, 4, 3> C{0.0};

    il::StaticArray<std::complex<double>, 6> V1{0.0}, V2{0.0};

    V1[0] = SinX/D1;
    V1[1] = H0Lim*D0;
    V1[2] = std::conj(V1[1]);
    V1[3] = d*D0*(H0Lim+2.0*I*eix);
    V1[4] = std::conj(V1[3]);
    V1[5] = -H0Lim*D1;

    V2[0] = 0.5*D02/D1*(SinX-2.0*I*eix*CosX*CosX);
    V2[1] = -0.125*D03*(4.0*H0Lim+I*eix*(e2x+8.0));
    V2[2] = 0.125*D03*(12.0*H0Lim+5.0*I*eix);
    V2[3] = -0.5*D04*D1*(3.0*H0Lim-2.0*I*eix*(e2x-3.0));
    V2[4] = -1.5*D1*H0Lim;
    V2[5] = 1.5*D02*D1*(H0Lim+2.0*I*eix);

    for (int j = 0; j < C.size(0); ++j) {
        C(j, 0, 2) = 0.5*C1_2*V1[j];
        C(j, 1, 2) = C1m2*V2[j];
        C(j, 2, 0) = 0.25*C2m1*V1[j];
        C(j, 2, 1) = 0.5*nu*V2[j];
        C(j, 3, 2) = V1[j];
    }

    return C;
}

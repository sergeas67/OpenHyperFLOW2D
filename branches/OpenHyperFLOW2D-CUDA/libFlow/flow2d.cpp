/*************************************************
*  Libflow v2.2                                  *
*  Copyright (C)  1995-2014 by Serge A. Suchkov  *
*  Copyright policy: LGPL V3                     *
*************************************************/
#include "libFlow/flow2d.hpp"

Flow2D::Flow2D(double u, double v):Flow() {
    double Tmp;
    UU = u;
    VV = v;
    Tmp =sqrt(UU*UU+VV*VV+1.e-12);
    Wg(Tmp);
}

Flow2D::Flow2D(Flow& f,double u, double v):Flow(f) {
    double Tmp;

    UU = u;
    VV = v;
    Tmp = sqrt(UU*UU+VV*VV+1.e-12);
    Wg(Tmp);
}

Flow2D::Flow2D(double _mu,
               double _lam,
               double Cp,
               double T,
               double P,
               double R,
               double u,
               double v):Flow(Cp,T,P,R,_lam,_mu) {
    double Tmp;
    UU = u;
    VV = v;
    Tmp = sqrt(UU*UU+VV*VV+1.e-12);
    Wg(Tmp);
}

double Flow2D::U(double u) {
    double Tmp;
    UU = u;
    Tmp = sqrt(UU*UU+VV*VV+1.e-12);
    Wg(Tmp);
    return UU;
}

double Flow2D::V(double v) {
    double Tmp;
    VV = v;
    Tmp = sqrt(UU*UU+VV*VV+1.e-12);
    Wg(Tmp);
    return VV;
}

double Flow2D::Wg(double u, double v) {
    double Tmp;
    UU=u;
    VV=v;
    Tmp = sqrt(UU*UU+VV*VV+1.e-12);
    return Wg(Tmp);
}

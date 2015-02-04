/*************************************************
*  Libflow v2.2                                  *
*  Copyright (C)  1995-2015 by Serge A. Suchkov  *
*  Copyright policy: LGPL V3                     *
*   http://openhyperflow2d.googlecode.com        *
*************************************************/
#include "libFlow/flow2d.hpp"

Flow2D::Flow2D(FP u, FP v):Flow() {
    FP Tmp;
    UU = u;
    VV = v;
    Tmp =sqrt(UU*UU+VV*VV+1.e-12);
    Wg(Tmp);
}

Flow2D::Flow2D(Flow& f,FP u, FP v):Flow(f) {
    FP Tmp;

    UU = u;
    VV = v;
    Tmp = sqrt(UU*UU+VV*VV+1.e-12);
    Wg(Tmp);
}

Flow2D::Flow2D(FP _mu,
               FP _lam,
               FP Cp,
               FP T,
               FP P,
               FP R,
               FP u,
               FP v):Flow(Cp,T,P,R,_lam,_mu) {
    FP Tmp;
    UU = u;
    VV = v;
    Tmp = sqrt(UU*UU+VV*VV+1.e-12);
    Wg(Tmp);
}

FP Flow2D::U(FP u) {
    FP Tmp;
    UU = u;
    Tmp = sqrt(UU*UU+VV*VV+1.e-12);
    Wg(Tmp);
    return UU;
}

FP Flow2D::V(FP v) {
    FP Tmp;
    VV = v;
    Tmp = sqrt(UU*UU+VV*VV+1.e-12);
    Wg(Tmp);
    return VV;
}

FP Flow2D::Wg(FP u, FP v) {
    FP Tmp;
    UU=u;
    VV=v;
    Tmp = sqrt(UU*UU+VV*VV+1.e-12);
    return Wg(Tmp);
}

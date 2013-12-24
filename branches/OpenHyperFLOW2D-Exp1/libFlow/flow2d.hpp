/*************************************************
*  Libflow v2.2                                  *
*  Copyright (C)  1995-2013 by Serge A. Suchkov  *
*  Copyright policy: LGPL V3                     *
*************************************************/
#ifndef _flow2d_hpp
#define _flow2d_hpp

#include "libFlow/flow.hpp"
#include <string.h>

class   Flow2D: public Flow {

    double UU;
    double VV;

    double Wg(double w) {
        return Flow::Wg(w);
    }
    double LAM(double l) {
        return Flow::LAM(l);
    }

public:

    double MACH(double m) {
        return Flow::MACH(m);
    }
    double Wg() {
        return /*Flow::Wg()*/sqrt(UU*UU+VV*VV+1.e-5);
    }
    double LAM() {
        return Flow::LAM();
    }
    double MACH() {
        return Flow::MACH();
    }


    Flow2D() {
        UU = Wg();VV = 0;
    }
    Flow2D(Flow& f):Flow(f) {
        UU = Wg();VV = 0;
    }
    Flow2D(double u, double v);
    Flow2D(Flow& f,double u, double v);
    Flow2D(double _mu,double _lam,double Cp,double T,double P,double R,double u,double v);

    virtual ~Flow2D() {
        ;
    }

    double U() {
        return UU;
    }
    double V() {
        return VV;
    }
    double U(double u);
    double V(double v);
    double Wg(double u, double v);

    Flow2D& operator = (Flow2D& NewFlow2D) {
        memcpy(this, &NewFlow2D, sizeof(Flow2D));return *this;
    }
};

#endif //_flow2d_hpp


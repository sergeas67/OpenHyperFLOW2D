/*************************************************
*  Libflow v2.2                                  *
*  Copyright (C)  1995-2016 by Serge A. Suchkov  *
*  Copyright policy: LGPL V3                     *
*  http://github.com/sergeas67/openhyperflow2d   *
*************************************************/
#ifndef _flow2d_hpp
#define _flow2d_hpp

#include "libFlow/flow.hpp"
#include <string.h>

class   Flow2D: public Flow {

    FP UU;
    FP VV;

    FP Wg(FP w) {
        return Flow::Wg(w);
    }

    FP LAM(FP l) {

        if(VV != 0.0 && VV !=0.0) {
          FP angle = atan(VV/UU);

          Flow::LAM(l);

          UU = Flow::Wg()*cos(angle);
          VV = Flow::Wg()*sin(angle);
        } else {

          Flow::LAM(l);

          if(VV == 0.0) {
             UU = Flow::Wg();
          } else if (UU == 0.0) {
             VV = Flow::Wg();
          }
       }
        return Flow::LAM();
    }

public:

    FP MACH(FP m) {

        if(VV != 0.0 && VV !=0.0) {
          FP angle = atan(VV/UU);
          Flow::MACH(m);

          UU = Flow::Wg()*cos(angle);
          VV = Flow::Wg()*sin(angle);

        } else {

        Flow::MACH(m);

          if(VV == 0.0) {
             UU = Flow::Wg();
          } else if (UU == 0.0) {
             VV = Flow::Wg();
          }
       }
        return Flow::MACH();
    }

    FP Wg() {
        return /*Flow::Wg()*/sqrt(UU*UU+VV*VV+1.e-5);
    }

    FP LAM() {
        return Flow::LAM();
    }

    FP MACH() {
        return Flow::MACH();
    }


    Flow2D() {
        UU = Wg();VV = 0;
    }

    Flow2D(Flow& f):Flow(f) {
        UU = Wg();VV = 0;
    }

    Flow2D(FP u, FP v);
    Flow2D(Flow& f,FP u, FP v);
    Flow2D(FP _mu,FP _lam,FP Cp,FP T,FP P,FP R,FP u,FP v);

    virtual ~Flow2D() {
        ;
    }

    FP U() {
        return UU;
    }
    FP V() {
        return VV;
    }
    FP U(FP u);
    FP V(FP v);
    FP Wg(FP u, FP v);

    Flow2D& operator = (Flow2D& NewFlow2D) {
        memcpy(this, &NewFlow2D, sizeof(Flow2D));return *this;
    }
};

#endif //_flow2d_hpp


/*******************************************************************************
*   OpenHyperFLOW2D                                                            *
*                                                                              *
*   Version  1.0.3                                                             *
*   Copyright (C)  1995-2016 by Serge A. Suchkov                               *
*   Copyright policy: LGPL V3                                                  *
*   http://github.com/sergeas67/openhyperflow2d                                *
*                                                                              *
*   last update: 07/04/2016                                                    *
*******************************************************************************/
#include "libExcept/except.hpp"
#include "libOpenHyperFLOW2D/hyper_flow_bound_circle.hpp"
//  Bound Circle
BoundCircle2D::BoundCircle2D(char* name,
                             UMatrix2D< FlowNode2D< FP, NUM_COMPONENTS> >* JM,  // Computation area 
                             FP  x,                      // Start circle 
                             FP  y,                      // coordinates (x,y)
                             FP  x1,                     // Center circle
                             FP  y1,                     // coordinates (x,y) 
                             FP  dx,                     // dx - step
                             FP  dy,                     // dy - step
                             ulong    ct,                // condition type
                             int      MaterialID,        // Material ID 0- gas !0 - solid
                             Flow2D*  pInFlow2D,         // init flow2d object on circle bound
                             FP*      Y,                 // component matrix
                             ulong    bctt,              // Bound contour turbulence type
                             ostream* dbg_output):BoundContour2D(name,JM,
                                                                 (int)(x/dx+0.4999),
                                                                 (int)(y/dy+0.4999)
                                                                 ),Area2D(name,JM) {
    int    k,i;
    FP xx1,yy1,xx2,yy2,fi0,r = sqrt((x-x1)*(x-x1)+(y-y1)*(y-y1)+1.e-30);
    xx1 = x;
    yy1 = y;
    fi0 = atan2((y1-y),(x1-x));
    if (pInFlow2D) {
        pInFlow2D->U(0.);
        pInFlow2D->V(0.);
    }
    int    ix,iy,ix0,iy0;

    k   = max(1,(int)(2*pi*r/sqrt(dx*dx+dy*dy)));
    
    ix0 = (unsigned int)(x1/dx+0.499999);
    iy0 = (unsigned int)(y1/dx+0.499999);
    
    /*
    ix0 = max(0,ix0);
    iy0 = max(0,iy0);
    ix0 = min(ix0,JM->GetX() - 1);
    iy0 = min(iy0,JM->GetY() - 1); 
    */
    for (i=0; i<k; i++) {
        xx2 = x1+(r*sin(fi0+(2.*pi*i)/k-pi/2.)); 
        yy2 = y1+(r*cos(fi0+(2.*pi*i)/k-pi/2.));
        ix = (unsigned int)(xx2/dx+0.499999);
        iy = (unsigned int)(yy2/dy+0.499999);
        if (dbg_output)
            *dbg_output << "\ndebug: Add circle node ("<< xx2 <<","<< yy2 << ")"
            "->[" << ix <<","<< iy << "]"
            << "\n"  << flush;
     if( ix >= 0 && iy >= 0 &&
         ix <= ((int)JM->GetX() - 1) && 
         iy <= ((int)JM->GetY() - 1)) {
        AddBound2D(name,ix,iy,ct, NULL, pInFlow2D,Y,bctt); 
      } /* else {
        ix = max(0,ix);
        iy = max(0,iy);
        ix = min(ix,JM->GetX() - 1);
        iy = min(iy,JM->GetY() - 1); 
      }*/
      xx1 = xx2;
      yy1 = yy2;
    }
    CloseContour2D(name,ct, NULL, pInFlow2D,Y,bctt);
    SetBounds();
    if (dbg_output)
        *dbg_output << "\n \ndebug:Start point for fill circle ("<< x1 <<","<< y1 << ")"
        "->[" << x1/dx <<","<< y1/dy << "]"
        << "\n" << flush;
    if(MaterialID) {
        FillArea2D(x1,y1,NT_S_2D,MaterialID);
    } else {
        FillArea2D(x1,y1,NT_F_2D,pInFlow2D,Y,bctt,MaterialID);
    }
}

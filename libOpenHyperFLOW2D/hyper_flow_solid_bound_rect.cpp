/*******************************************************************************
*   OpenHyperFLOW2D                                                            *
*                                                                              *
*   Version  1.0.3                                                             *
*   Copyright (C)  1995-2016 by Serge A. Suchkov                               *
*   Copyright policy: LGPL V3                                                  *
*   http://github.com/sergeas67/openhyperflow2d                                *
*                                                                              *
*   last update: 14/04/2016                                                    *
*******************************************************************************/
#include "libExcept/except.hpp"
#include "libOpenHyperFLOW2D/hyper_flow_solid_bound_rect.hpp"
// SolidBoundRect constructor
SolidBoundRect2D::SolidBoundRect2D(char* name,
                                   UMatrix2D< FlowNode2D< FP, NUM_COMPONENTS> >* JM,// Computation area reference
                                   FP  x,                                           // Start rect 
                                   FP  y,                                           // coordinates (x,y) left-down corner
                                   FP  DX,                                          // Rect width
                                   FP  DY,                                          // Rect hight
                                   FP  dx,                                          // dx - step
                                   FP  dy,                                          // dy - step
                                   ulong   ct,                                      // condition type
                                   Flow2D* pInFlow2D,                               // init flow2d object on circle bound
                                   FP* _Y,                                          // component matrix
                                   ulong   bctt,                                    // bound contour turbulence type
                                   ostream* dbg_output):BoundContour2D(name,JM,     // embedded BoundContour object
                                                                       (unsigned int)(x/dx+0.4999),
                                                                       (unsigned int)(y/dy+0.4999)),Area2D(name,JM) {
    FP xx1,yy1,xx2,yy2;
    unsigned int    ix,iy;
    int             tt;   
    dX = DX;
    dY = DY;
    X  = x;
    Y  = y;
    if (pInFlow2D) {          // ??
        pInFlow2D->U(0.);
        pInFlow2D->V(0.);
    }
    xx1 = x;
    yy1 = y;
//      1        
    xx2 = xx1+DX;
    yy2 = yy1;
    ix = (unsigned int)(xx2/dx+0.4999);
    iy = (unsigned int)(yy2/dy+0.4999);
    
    if (isTurbulenceCond2D(bctt, TCT_k_eps_Model_2D)) {
        tt = bctt | TCT_dkdy_NULL_2D| TCT_k_CONST_2D | TCT_eps_mud2kdy2_WALL_2D;
    } else {
        tt = bctt;
    }
    AddBound2D(name,ix,iy,ct, NULL, pInFlow2D,_Y,tt); 
    xx1 = xx2;
    yy1 = yy2;
//      2        
    xx2 = xx1;
    yy2 = yy1+DY;
    ix = (unsigned int)(xx2/dx+0.4999);
    iy = (unsigned int)(yy2/dy+0.4999);

    if (isTurbulenceCond2D(bctt, TCT_k_eps_Model_2D)) {
        tt = bctt | TCT_dkdx_NULL_2D | TCT_k_CONST_2D | TCT_eps_mud2kdx2_WALL_2D;
    } else {
        tt = bctt;
    }
    AddBound2D(name,ix,iy,ct, NULL, pInFlow2D,_Y,tt); 
    xx1 = xx2;
    yy1 = yy2;
//      3        
    xx2 = xx1-DX;
    yy2 = yy1;
    
    ix = (unsigned int)(xx2/dx+0.4999);
    iy = (unsigned int)(yy2/dy+0.4999);

    if (isTurbulenceCond2D(bctt, TCT_k_eps_Model_2D)) {
        tt = bctt | TCT_dkdy_NULL_2D | TCT_k_CONST_2D | TCT_eps_mud2kdy2_WALL_2D; 
    } else {
        tt = bctt;
    }
    AddBound2D(name,ix,iy,ct, NULL, pInFlow2D,_Y,tt); 

    if (isTurbulenceCond2D(bctt, TCT_k_eps_Model_2D)) {
        tt = bctt | TCT_dkdx_NULL_2D  | TCT_k_CONST_2D | TCT_eps_mud2kdx2_WALL_2D; 
    } else {
        tt = bctt;
    }
    CloseContour2D(name, ct, NULL, pInFlow2D,_Y,tt);
    SetBounds();
    ix = (int)((x+DX/2)/dx+0.4999);
    iy = (int)((y+DY/2)/dy+0.4999);
    FillArea2D(ix,iy,NT_S_2D);
}


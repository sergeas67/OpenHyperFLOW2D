/*******************************************************************************
*   OpenHyperFLOW2D                                                            *
*                                                                              *
*   Version  1.0.1                                                             *
*   Copyright (C)  1995-2014 by Serge A. Suchkov                               *
*   Copyright policy: LGPL V3                                                  *
*                                                                              *
*   last update: 01/07/2014                                                    *
*******************************************************************************/
#include "libExcept/except.hpp"
#include "libOpenHyperFLOW2D/hyper_flow_solid_bound_circle.hpp"
// Solid Bound Circle
SolidBoundCircle2D::SolidBoundCircle2D(char* name,
                                      UMatrix2D< FlowNode2D< double, NUM_COMPONENTS> >* JM,  // Computation area 
                                      double  x,                      // Start circle 
                                      double  y,                      // coordinates (x,y)
                                      double  x1,                     // Center circle
                                      double  y1,                     // coordinates (x,y) 
    #ifdef _UNIFORM_MESH_
                                      double  dx,                     // dx - step
                                      double  dy,                     // dy - step
    #else 
                                      unsigned int     num_segments,
                                      Mesh2D* p_mesh,
    #endif // _UNIFORM_MESH_
                                      ulong   ct,                     // condition type
                                      Flow2D* pInFlow2D,              // init flow2d object on circle bound
                                      double* Y    ,                  // component matrix
                                      ulong   bctt,                   // Bound contour turbulence type
                                      ostream* dbg_output):BoundContour2D(name,JM,
    #ifndef _UNIFORM_MESH_
                                                                          p_mesh,x,y

    #else
                                                                         (int)(x/dx+0.4999),
                                                                         (int)(y/dy+0.4999)
    #endif // _UNIFORM_MESH_
                                                                         ),Area2D(name,JM
    #ifndef _UNIFORM_MESH_
                                                                                ,p_mesh
    #endif // _UNIFORM_MESH_
                                                                               ) {
    int    k,i;
    double xx1,yy1,xx2,yy2,fi0,r = sqrt((x-x1)*(x-x1)+(y-y1)*(y-y1)+1.e-30);
    xx1 = x;
    yy1 = y;
    fi0 = atan2((y1-y),(x1-x));
    if (pInFlow2D) {
        pInFlow2D->U(0.);
        pInFlow2D->V(0.);
    }
#ifdef _UNIFORM_MESH_
    int    ix,iy;
    k = max(1,(int)(2*pi*r/sqrt(dx*dx+dy*dy)));
#else
    k = num_segments;
#endif // _UNIFORM_MESH_
    for (i=0; i<k; i++) {
        xx2 = x1+(r*sin(fi0+(2.*pi*i)/k-pi/2.)); 
        yy2 = y1+(r*cos(fi0+(2.*pi*i)/k-pi/2.));
#ifdef _UNIFORM_MESH_
        ix = (unsigned int)(xx2/dx+0.499999);
        iy = (unsigned int)(yy2/dy+0.499999);
#endif // _UNIFORM_MESH_
        if (dbg_output)
            *dbg_output << "\ndebug: Add circle node ("<< xx2 <<","<< yy2 << ")"
#ifdef _UNIFORM_MESH_
            "->[" << ix <<","<< iy << "]"
#endif // _UNIFORM_MESH_
            << "\n"  << flush;
     if( ix >= 0 && iy >= 0 &&
         ix <= ((int)JM->GetX() - 1) && 
         iy <= ((int)JM->GetY() - 1)) {
#ifdef _UNIFORM_MESH_
        AddBound2D(name,ix,iy,ct, NULL, pInFlow2D,Y,bctt); 
#else
        AddBound2D(name,xx2,yy2,ct, NULL, pInFlow2D,Y,bctt); 
#endif // _UNIFORM_MESH_
      }
        xx1 = xx2;
        yy1 = yy2;
    }
    CloseContour2D(name,ct, NULL, pInFlow2D,Y,bctt);
    SetBounds();
    if (dbg_output)
        *dbg_output << "\n \ndebug:Start point for fill circle ("<< x1 <<","<< y1 << ")"
#ifdef _UNIFORM_MESH_
        "->[" << x1/dx <<","<< y1/dy << "]"
#endif // _UNIFORM_MESH_
        << "\n" << flush;
#ifdef _UNIFORM_MESH_
    //FillArea2D((int)(x1/dx+0.4999),(int)(y1/dy+0.4999),NT_S_2D);
    FillArea2D(x1,y1,NT_S_2D);
#else
    FillArea2D(x1,y1,NT_S_2D);
#endif // _UNIFORM_MESH_

}

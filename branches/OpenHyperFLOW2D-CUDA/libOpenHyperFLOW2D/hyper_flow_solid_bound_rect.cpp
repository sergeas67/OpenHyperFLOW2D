/*******************************************************************************
*   OpenHyperFLOW2D                                                            *
*                                                                              *
*   Version  1.0.0                                                             *
*   Copyright (C)  1995-2013 by Serge A. Suchkov                               *
*   Copyright policy: LGPL V3                                                  *
*                                                                              *
*   last update: 15/12/2013                                                    *
*******************************************************************************/
#include "libExcept/except.hpp"
#include "libOpenHyperFLOW2D/hyper_flow_solid_bound_rect.hpp"
// SolidBoundRect constructor
SolidBoundRect2D::SolidBoundRect2D(char* name,
                                   UMatrix2D< FlowNode2D< double, NUM_COMPONENTS> >* JM,// Computation area reference
                                   double  x,                                           // Start rect 
                                   double  y,                                           // coordinates (x,y) left-down corner
                                   double  DX,                                          // Rect width
                                   double  DY,                                          // Rect hight
    #ifdef _UNIFORM_MESH_
                                   double  dx,                                          // dx - step
                                   double  dy,                                          // dy - step
    #else 
                                   Mesh2D* p_mesh,
    #endif // _UNIFORM_MESH_
                                   ulong   ct,                                          // condition type
                                   Flow2D* pInFlow2D,                                   // init flow2d object on circle bound
                                   double* _Y,                                          // component matrix
                                   ulong   bctt,                                        // bound contour turbulence type
                                   ostream* dbg_output):BoundContour2D(name,JM,         // embedded BoundContour object
    #ifndef _UNIFORM_MESH_
                                                                       p_mesh,x,y
    #else
                                                                       (unsigned int)(x/dx+0.4999),
                                                                       (unsigned int)(y/dy+0.4999)
    #endif // _UNIFORM_MESH_
                                                                       ),Area2D(name,JM   // embedded Area object
    #ifndef _UNIFORM_MESH_
                                                                        ,p_mesh
    #endif // _UNIFORM_MESH_
                                                                       ) {
    double xx1,yy1,xx2,yy2;
#ifdef _UNIFORM_MESH_
    unsigned int    ix,iy;
#else
    pMesh = p_mesh;
#endif // _UNIFORM_MESH_
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
#ifdef _UNIFORM_MESH_
    ix = (unsigned int)(xx2/dx+0.4999);
    iy = (unsigned int)(yy2/dy+0.4999);
#endif // _UNIFORM_MESH_
    
    if (isTurbulenceCond2D(bctt, TCT_k_eps_Model_2D)) {
        tt = bctt | TCT_dkdy_NULL_2D| TCT_k_CONST_2D | TCT_eps_mud2kdy2_WALL_2D;
    } else {
        tt = bctt;
    }
#ifdef _UNIFORM_MESH_
    AddBound2D(name,ix,iy,ct, NULL, pInFlow2D,_Y,tt); 
#else
    AddBound2D(name,xx2,yy2,ct, NULL, pInFlow2D,_Y,tt); 
#endif // _UNIFORM_MESH_
    xx1 = xx2;
    yy1 = yy2;
//      2        
    xx2 = xx1;
    yy2 = yy1+DY;
#ifdef _UNIFORM_MESH_
    ix = (unsigned int)(xx2/dx+0.4999);
    iy = (unsigned int)(yy2/dy+0.4999);
#endif // _UNIFORM_MESH_

    if (isTurbulenceCond2D(bctt, TCT_k_eps_Model_2D)) {
        tt = bctt | TCT_dkdx_NULL_2D | TCT_k_CONST_2D | TCT_eps_mud2kdx2_WALL_2D;
    } else {
        tt = bctt;
    }
#ifdef _UNIFORM_MESH_
    AddBound2D(name,ix,iy,ct, NULL, pInFlow2D,_Y,tt); 
#else
    AddBound2D(name,xx2,yy2,ct, NULL, pInFlow2D,_Y,tt); 
#endif // _UNIFORM_MESH_
    xx1 = xx2;
    yy1 = yy2;
//      3        
    xx2 = xx1-DX;
    yy2 = yy1;
#ifdef _UNIFORM_MESH_
    ix = (unsigned int)(xx2/dx+0.4999);
    iy = (unsigned int)(yy2/dy+0.4999);
#endif // _UNIFORM_MESH_

    if (isTurbulenceCond2D(bctt, TCT_k_eps_Model_2D)) {
        tt = bctt | TCT_dkdy_NULL_2D | TCT_k_CONST_2D | TCT_eps_mud2kdy2_WALL_2D; 
    } else {
        tt = bctt;
    }
#ifdef _UNIFORM_MESH_
    AddBound2D(name,ix,iy,ct, NULL, pInFlow2D,_Y,tt); 
#else
    AddBound2D(name,xx2,yy2,ct, NULL, pInFlow2D,_Y,tt); 
#endif // _UNIFORM_MESH_

    if (isTurbulenceCond2D(bctt, TCT_k_eps_Model_2D)) {
        tt = bctt | TCT_dkdx_NULL_2D  | TCT_k_CONST_2D | TCT_eps_mud2kdx2_WALL_2D; 
    } else {
        tt = bctt;
    }
    CloseContour2D(name, ct, NULL, pInFlow2D,_Y,tt);
    SetBounds();
#ifdef _UNIFORM_MESH_
    ix = (int)((x+DX/2)/dx+0.4999);
    iy = (int)((y+DY/2)/dy+0.4999);
    FillArea2D(ix,iy,NT_S_2D);
#else
    FillArea2D(x(x+DX/2),(y+DY/2),NT_S_2D);
#endif // _UNIFORM_MESH_
}


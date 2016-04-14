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
#ifndef _hyper_flow_solid_bound_rect_hpp_
#define _hyper_flow_solid_bound_rect_hpp_

#include <iostream>
#include <fstream>
using namespace std;
#include "utl/umatrix2d.hpp"

//#include "utl/umatrix3d.hpp"

#define _SAFE_ACCESS_
#include "utl/uarray.hpp"
#undef   _SAFE_ACCESS_

#include "obj_data/obj_data.hpp"
#include "libOpenHyperFLOW2D/hyper_flow_node.hpp"
#include "libOpenHyperFLOW2D/hyper_flow_bound_contour.hpp"
#include "libOpenHyperFLOW2D/hyper_flow_area.hpp"

#ifndef _UNIFORM_MESH_
#include "libMesh/mesh.hpp"
#endif  // _UNIFORM_MESH_

#ifndef  max
    #define max(a,b) (((a)>(b))?(a):(b))
#endif 

#ifndef  min
    #define min(a,b) (((a)<(b))?(a):(b))
#endif 
//////////////////////////////////////////////////
//  SolidBoundRect2D object                     //
//////////////////////////////////////////////////
class SolidBoundRect2D : protected BoundContour2D,
                         protected Area2D {
    FP X, Y, dX,dY;
    char* SolidBoundRectName;

public:

  SolidBoundRect2D(char* name,
                   UMatrix2D< FlowNode2D< FP, NUM_COMPONENTS> >* JM, // Computation area reference
                   FP  x,                      // Start rect 
                   FP  y,                      // coordinates (x,y) left-down coner
                   FP  DX,                     // Rect width
                   FP  DY,                     // Rect hight
#ifdef _UNIFORM_MESH_
                   FP  dx,                     // dx - step
                   FP  dy,                     // dy - step
#else 
                   Mesh2D* p_mesh,
#endif // _UNIFORM_MESH_
                   ulong   ct,                     // condition type
                   Flow2D* pInFlow2D=0,            // init flow2d object on circle bound
                   FP* Y=0,                    // component matrix
                   ulong   bctt=TCT_No_Turbulence_2D,// Bound contour turbulence type
                   ostream* dbg_output=NULL);                  

    //void ReplaceSolidBoundRect(FP  x, FP  y);

    FP GetX() {
        return X;
    }
    FP GetY() {
        return Y;
    }

    FP GetDX() {
        return dX;
    }
    FP GetDY() {
        return dY;
    }
 char* GetSolidBoundRectName2D() {
      return SolidBoundRectName;
 }


};
#endif // _hyper_flow_solid_bound_rect_hpp_



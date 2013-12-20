/*******************************************************************************
*   OpenHyperFLOW2D                                                            *
*                                                                              *
*   Version  1.0.0                                                             *
*   Copyright (C)  1995-2013 by Serge A. Suchkov                               *
*   Copyright policy: LGPL V3                                                  *
*                                                                              *
*   last update: 15/12/2013                                                    *
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
    double X, Y, dX,dY;
    char* SolidBoundRectName;
#ifndef _UNIFORM_MESH_
    Mesh2D*  pMesh;       //  anisotropic mesh
#else

#endif // _UNIFORM_MESH_
public:

  SolidBoundRect2D(char* name,
                   UMatrix2D< FlowNode2D< double, NUM_COMPONENTS> >* JM, // Computation area reference
                   double  x,                      // Start rect 
                   double  y,                      // coordinates (x,y) left-down coner
                   double  DX,                     // Rect width
                   double  DY,                     // Rect hight
#ifdef _UNIFORM_MESH_
                   double  dx,                     // dx - step
                   double  dy,                     // dy - step
#else 
                   Mesh2D* p_mesh,
#endif // _UNIFORM_MESH_
                   ulong   ct,                     // condition type
                   Flow2D* pInFlow2D=0,            // init flow2d object on circle bound
                   double* Y=0,                    // component matrix
                   ulong   bctt=TCT_No_Turbulence_2D,// Bound contour turbulence type
                   ostream* dbg_output=NULL);                  

    //void ReplaceSolidBoundRect(double  x, double  y);

    double GetX() {
        return X;
    }
    double GetY() {
        return Y;
    }

    double GetDX() {
        return dX;
    }
    double GetDY() {
        return dY;
    }
 char* GetSolidBoundRectName2D() {
      return SolidBoundRectName;
 }


};
#endif // _hyper_flow_solid_bound_rect_hpp_



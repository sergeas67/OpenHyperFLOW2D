/*******************************************************************************
*   OpenHyperFLOW2D                                                            *
*                                                                              *
*   Version  1.0.2                                                             *
*   Copyright (C)  1995-2014 by Serge A. Suchkov                               *
*   Copyright policy: LGPL V3                                                  *
*                                                                              *
*   last update: 25/02/2014                                                    *
*******************************************************************************/
#ifndef _hyper_flow_solid_bound_circle_hpp_
#define _hyper_flow_solid_bound_circle_hpp_

#include <iostream>
#include <fstream>
using namespace std;
#include "utl/umatrix2d.hpp"

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
//  SolidBoundCircle2D object                   //
//////////////////////////////////////////////////
class SolidBoundCircle2D : protected BoundContour2D,
                           protected Area2D  {
char* SolidBoundCircleName;

#ifndef _UNIFORM_MESH_
    Mesh2D*  pMesh;       //  anisotropic mesh
#endif // _UNIFORM_MESH_

public:

  SolidBoundCircle2D(char* name,
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
                     Flow2D* pInFlow2D=NULL,         // init flow2d object on circle bound
                     double* Y=NULL,                 // component matrix
                     ulong   bctt=TCT_No_Turbulence_2D,// Bound contour turbulence type
                     ostream* dbg_output=NULL);                  
  SolidBoundCircle2D(char* name,
                     UMatrix2D< FlowNode2D< double, NUM_COMPONENTS> >* JM,  // Computation area 
                     double  x,                     // Center circle
                     double  y,                     // coordinates (x,y) 
                     double  r,
#ifdef _UNIFORM_MESH_
                     double  dx,                     // dx - step
                     double  dy,                     // dy - step
#else 
                     unsigned int     num_segments,
                     Mesh2D* p_mesh,
#endif // _UNIFORM_MESH_
                     ulong   ct,                     // condition type
                     Flow2D* pInFlow2D=NULL,         // init flow2d object on circle bound
                     double* Y=NULL,                 // component matrix
                     ulong   bctt=TCT_No_Turbulence_2D,// Bound contour turbulence type
                     ostream* dbg_output=NULL);                  
    
    char* GetSolidBoundCircleName2D() {
         return SolidBoundCircleName;
    }
#ifndef _UNIFORM_MESH_
    Mesh2D*    GetMesh() {
        return pMesh;
    }
#endif // _UNIFORM_MESH_
};
#endif // _hyper_flow_solid_bound_circle_hpp_



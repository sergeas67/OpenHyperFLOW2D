/*******************************************************************************
*   OpenHyperFLOW2D                                                            *
*                                                                              *
*   Version  1.0.2                                                             *
*   Copyright (C)  1995-2014 by Serge A. Suchkov                               *
*   Copyright policy: LGPL V3                                                  *
*                                                                              *
*   last update: 06/12/2014                                                    *
*******************************************************************************/
#ifndef _hyper_flow_bound_circle_hpp_
#define _hyper_flow_bound_circle_hpp_

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

#ifndef  max
    #define max(a,b) (((a)>(b))?(a):(b))
#endif 

#ifndef  min
    #define min(a,b) (((a)<(b))?(a):(b))
#endif
 
//////////////////////////////////////////////////
//  BoundCircle2D object                        //
//////////////////////////////////////////////////
class BoundCircle2D : protected BoundContour2D,
                      protected Area2D  {
char* BoundCircleName;

public:

  BoundCircle2D(char* name,
                UMatrix2D< FlowNode2D< FP, NUM_COMPONENTS> >* JM,  // Computation area 
                FP  x,                            // Start circle 
                FP  y,                            // coordinates (x,y)
                FP  x1,                           // Center circle
                FP  y1,                           // coordinates (x,y) 
                FP  dx,                           // dx - step
                FP  dy,                           // dy - step
                ulong    ct,                       // condition type
                int      MaterialID=GAS_ID,            // Material ID 0- gas !0 - solid
                Flow2D*  pInFlow2D=NULL,           // init flow2d object on circle bound
                FP*      Y=NULL,                       // component matrix
                ulong    bctt=TCT_No_Turbulence_2D,// Bound contour turbulence type
                ostream* dbg_output=NULL);                  
  
  BoundCircle2D(char* name,
                UMatrix2D< FlowNode2D< FP, NUM_COMPONENTS> >* JM,  // Computation area 
                FP  x,                            // Center circle
                FP  y,                            // coordinates (x,y) 
                FP  r,
                FP  dx,                           // dx - step
                FP  dy,                           // dy - step
                ulong    ct,                       // condition type
                int      MaterialID=GAS_ID,            // Material ID 0- gas !0 - solid
                Flow2D*  pInFlow2D=NULL,           // init flow2d object on circle bound
                FP*      Y=NULL,                   // component matrix
                ulong    bctt=TCT_No_Turbulence_2D,// Bound contour turbulence type
                ostream* dbg_output=NULL);                  
    
    char* GetBoundCircleName2D() {
         return BoundCircleName;
    }
};
#endif // _hyper_flow_bound_circle_hpp_



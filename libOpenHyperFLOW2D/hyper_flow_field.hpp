/*******************************************************************************
*   OpenHyperFLOW2D                                                            *
*                                                                              *
*   Version  1.0.1                                                             *
*   Copyright (C)  1995-2014 by Serge A. Suchkov                               *
*   Copyright policy: LGPL V3                                                  *
*                                                                              *
*   last update: 01/07/2014                                                    *
*******************************************************************************/
#ifndef _hyper_flow_field_hpp_
#define _hyper_flow_field_hpp_

#include <iostream>
#include <fstream>
using namespace std;
#include "utl/umatrix2d.hpp"

#define _SAFE_ACCESS_
#include "utl/uarray.hpp"
#undef   _SAFE_ACCESS_

#include "obj_data/obj_data.hpp"
#include "libOpenHyperFLOW2D/hyper_flow_node.hpp"

#ifndef _UNIFORM_MESH_
#include "libMesh/mesh.hpp"
#endif  // _UNIFORM_MESH_

#ifndef  max
    #define max(a,b) (((a)>(b))?(a):(b))
#endif 

#ifndef  min
    #define min(a,b) (((a)<(b))?(a):(b))
#endif 

//////////////////////////////////////////////////////
//  FlowField2D object                              //
//////////////////////////////////////////////////////
class FlowField2D : public UMatrix2D< FlowNode2D< FP, NUM_COMPONENTS> >
{
     char*    FlowFieldName;
     int      fd;
     int      p;
 public:
     //FlowField2D(InputData*  data);
     FlowField2D(char* name, InputData*  data);
     FlowField2D(char* name, int x, int y);
     FlowField2D(FlowField2D*);
     FlowField2D(char* name, UMatrix2D< FlowNode2D< FP, NUM_COMPONENTS> >*);
     ~FlowField2D();
     char* GetFlowFieldName() {return FlowFieldName;}
     int SaveFlowField2D(char* filename);
};
#endif // _hyper_flow_field_hpp_





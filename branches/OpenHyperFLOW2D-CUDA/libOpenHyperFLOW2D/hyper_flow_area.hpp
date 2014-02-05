/*******************************************************************************
*   OpenHyperFLOW2D                                                            *
*                                                                              *
*   Version  1.0.1                                                             *
*   Copyright (C)  1995-2014 by Serge A. Suchkov                               *
*   Copyright policy: LGPL V3                                                  *
*                                                                              *
*   last update: 16/01/2014                                                    *
*******************************************************************************/
#ifndef _hyper_flow_area_hpp_
#define _hyper_flow_area_hpp_
#include <iostream>
#include <fstream>
using namespace std;
#include "utl/umatrix2d.hpp"

#define GAS_ID  0

#define _SAFE_ACCESS_
#include "utl/uarray.hpp"
#undef   _SAFE_ACCESS_

#include "obj_data/obj_data.hpp"
#include "libOpenHyperFLOW2D/hyper_flow_node.hpp"
#include "libOpenHyperFLOW2D/hyper_flow_field.hpp"

#ifndef _UNIFORM_MESH_
#include "libMesh/mesh.hpp"
#endif  // _UNIFORM_MESH_

#ifndef  max
    #define max(a,b) (((a)>(b))?(a):(b))
#endif

#ifndef  min
    #define min(a,b) (((a)<(b))?(a):(b))
#endif

enum FieldSetMode {
  FSM_REWRITE,
  FSM_WRITE_IF_EMPTY,
  FSM_WRITE_IF_GAS,
  FSM_WRITE_IF_SOLID
};

enum AreaState {
    AS_OK,
    AS_ERR_OUT_OF_RANGE,
    AS_ERR_INIT_POINT
};

//////////////////////////////////////////////////
//  Area2D object                               //
//////////////////////////////////////////////////
class 
Area2D  {
    char*                         AreaName;                // Name of area
    UMatrix2D<FlowNode2D<double,NUM_COMPONENTS> >* pMFN;     // Reference to computational area matrix
#ifndef _UNIFORM_MESH_
    Mesh2D*                       pMesh;                   //  anisotropic mesh
#else
    unsigned int                  StartX,StartY;           // Initialisation start point coordinates (nodes)
#endif // _UNIFORM_MESH_
    double                        fStartX,fStartY;         // Initialisation start point coordinates (m) 
    AreaState                     as;                      // Area state
    double*                       pY;                      // Y[a] - Y1,Y2,Y3,...,Ya
    ulong                         ANT;                     // Area nodes type
    ulong                         ATT;                     // Area turbulence type

public:
#ifdef _UNIFORM_MESH_
    void FillArea2D(unsigned int x,
                    unsigned int y,
                    ulong,ulong att=TCT_No_Turbulence_2D,int MaterialID=GAS_ID);
    void FillArea2D(unsigned int x,
                    unsigned int y,
                    ulong, Flow2D* pf2d, double* Y=NULL,ulong att=TCT_No_Turbulence_2D,int MaterialID=GAS_ID);
    void FillArea2D(unsigned int,
                    unsigned int,
                    ulong, Flow* pf, double* Y=NULL,ulong att=TCT_No_Turbulence_2D,int MaterialID=GAS_ID);
    void FillArea2D(double x,
                    double y,
                    ulong, ulong att=TCT_No_Turbulence_2D,int MaterialID=GAS_ID);
    void FillArea2D(double x,
                    double y,
                    ulong, Flow2D* pf2d, double* Y=NULL,ulong att=TCT_No_Turbulence_2D,int MaterialID=GAS_ID);
    void FillArea2D(double x,
                    double y,
                    ulong, Flow* pf, double* Y=NULL,ulong att=TCT_No_Turbulence_2D,int MaterialID=GAS_ID);
#else
    void FillArea2D(double x,
                    double y,
                    ulong, int att=TCT_No_Turbulence_2D,int MaterialID=GAS_ID);
    void FillArea2D(double x,
                    double y,
                    ulong, Flow2D* pf2d, double* Y=NULL,ulong att=TCT_No_Turbulence_2D,int MaterialID=GAS_ID);
    void FillArea2D(double x,
                    double y,
                    ulong, Flow* pf, double* Y=NULL,ulong att=TCT_No_Turbulence_2D,int MaterialID=GAS_ID);
#endif // _UNIFORM_MESH_
    Area2D(char* name, UMatrix2D<FlowNode2D< double,NUM_COMPONENTS > >* p_j
#ifndef _UNIFORM_MESH_
          ,Mesh2D* p_mesh       //  anisotropic mesh
#endif // _UNIFORM_MESH_
        );
    ~Area2D();
    
AreaState GetAreaState() {
        return as;
    }

#ifdef _UNIFORM_MESH_
    unsigned int GetStartX() {
        return StartX;
    }
    unsigned int GetStartY() {
        return StartY;
    }
#else
    double GetStartFX() {
        return fStartX;
    }
    double GetStartFY() {
        return fStartY;
    }
    Mesh2D*    GetMesh() {
        return pMesh;
    }
#endif // _UNIFORM_MESH_

 char* GetAreaName() {return AreaName;}
};

void Abort_OpenHyperFLOW2D(int num_cuda_dev=1);
void Exit_OpenHyperFLOW2D(int num_cuda_dev=1);

// <------------- 2D --------------->           
#endif // _hyper_flow_area_hpp_



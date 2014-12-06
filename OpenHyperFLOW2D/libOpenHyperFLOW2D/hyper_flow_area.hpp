/*******************************************************************************
*   OpenHyperFLOW2D                                                            *
*                                                                              *
*   Version  1.0.2                                                             *
*   Copyright (C)  1995-2014 by Serge A. Suchkov                               *
*   Copyright policy: LGPL V3                                                  *
*                                                                              *
*   last update: 06/12/2014                                                    *
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
    char*                         AreaName;               // Name of area
    UMatrix2D<FlowNode2D<FP,NUM_COMPONENTS> >* pMFN;      // Reference to computational area matrix
    unsigned int                  StartX,StartY;          // Initialisation start point coordinates (nodes)
    FP                            fStartX,fStartY;        // Initialisation start point coordinates (m) 
    AreaState                     as;                     // Area state
    FP*                           pY;                     // Y[a] - Y1,Y2,Y3,...,Ya
    ulong                         ANT;                    // Area nodes type
    ulong                         ATT;                    // Area turbulence type

public:
    void FillArea2D(unsigned int x,
                    unsigned int y,
                    ulong,ulong att=TCT_No_Turbulence_2D,int MaterialID=GAS_ID);
    void FillArea2D(unsigned int x,
                    unsigned int y,
                    ulong, Flow2D* pf2d, FP* Y=NULL,ulong att=TCT_No_Turbulence_2D,int MaterialID=GAS_ID);
    void FillArea2D(unsigned int,
                    unsigned int,
                    ulong, Flow* pf, FP* Y=NULL,ulong att=TCT_No_Turbulence_2D,int MaterialID=GAS_ID);
    void FillArea2D(FP x,
                    FP y,
                    ulong, ulong att=TCT_No_Turbulence_2D,int MaterialID=GAS_ID);
    void FillArea2D(FP x,
                    FP y,
                    ulong, Flow2D* pf2d, FP* Y=NULL,ulong att=TCT_No_Turbulence_2D,int MaterialID=GAS_ID);
    void FillArea2D(FP x,
                    FP y,
                    ulong, Flow* pf, FP* Y=NULL,ulong att=TCT_No_Turbulence_2D,int MaterialID=GAS_ID);
    Area2D(char* name, UMatrix2D<FlowNode2D< FP,NUM_COMPONENTS > >* p_j);
    ~Area2D();
    
AreaState GetAreaState() {
        return as;
    }

    unsigned int GetStartX() {
        return StartX;
    }
    unsigned int GetStartY() {
        return StartY;
    }

 char* GetAreaName() {return AreaName;}
};

void Abort_OpenHyperFLOW2D();
void Exit_OpenHyperFLOW2D();

// <------------- 2D --------------->           
#endif // _hyper_flow_area_hpp_



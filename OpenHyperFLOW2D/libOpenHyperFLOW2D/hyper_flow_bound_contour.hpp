/*******************************************************************************
*   OpenHyperFLOW2D                                                            *
*                                                                              *
*   Version  1.0.2                                                             *
*   Copyright (C)  1995-2014 by Serge A. Suchkov                               *
*   Copyright policy: LGPL V3                                                  *
*                                                                              *
*   last update: 06/12/2014                                                    *
*******************************************************************************/
#ifndef _hyper_flow_bound_contour_hpp_
#define _hyper_flow_bound_contour_hpp_

#include <iostream>
#include <fstream>
using namespace std;
#include "utl/umatrix2d.hpp"

#define _SAFE_ACCESS_
#include "utl/uarray.hpp"
#undef   _SAFE_ACCESS_

#include "obj_data/obj_data.hpp"
#include "libOpenHyperFLOW2D/hyper_flow_node.hpp"
#include "libOpenHyperFLOW2D/hyper_flow_bound.hpp"

#ifndef _UNIFORM_MESH_
#include "libMesh/mesh.hpp"
#endif  // _UNIFORM_MESH_

#ifndef  max
    #define max(a,b) (((a)>(b))?(a):(b))
#endif 

#ifndef  min
    #define min(a,b) (((a)<(b))?(a):(b))
#endif 
/*
enum BoundState {
    BND_INACTIVE, BND_OK, BND_ERR
};

enum BoundMoveType {
     BMT_TAIL,
     BMT_NOTAIL
};
*/
//////////////////////////////////////////////////
//  Bound 2D Contour object                     //
//////////////////////////////////////////////////
class BoundContour2D : private UArray<Bound2D*> {
    FP   f_current_x,f_first_x;
    FP   f_current_y,f_first_y;
    Table*   BoundContourTable;
    char*    BoundContourName;
    UMatrix2D< FlowNode2D< FP, NUM_COMPONENTS> >* FlowNodeMatrixPtr;
    
#ifndef _UNIFORM_MESH_
    Mesh2D*  pMesh;       //  anisotropic mesh
#else
    unsigned int current_x,first_x;
    unsigned int current_y,first_y;
#endif // _UNIFORM_MESH_
    int      isContourClosed;
    int      isActivateContour;

public:
    BoundContour2D(char* Name, UMatrix2D< FlowNode2D< FP, NUM_COMPONENTS> >* fnm,
#ifndef _UNIFORM_MESH_
                 Mesh2D* p_mesh,       //  anisotropic mesh
                 FP x=0.,
                 FP y=0.
#else
                 unsigned int x=0,
                 unsigned int y=0
#endif // _UNIFORM_MESH_
                 );

    ~BoundContour2D();

    int LoadBoundContourTable(Table* t);
    
    int AddBound2D(char* name,
                   unsigned int x,
                   unsigned int y,
                   ulong        bt,
                   Flow*        pInFlow=0,
                   Flow2D*      pInFlow2D=0,
                   FP*      y1=0,
                   ulong        btt=TCT_No_Turbulence_2D);
    
#ifndef _UNIFORM_MESH_    
    int AddBound2D(char* name,
                   Mesh2D*      p_mesh,       //  anisotropic mesh
                   FP       f_x,
                   FP       f_y,
                   ulong        bt,
                   Flow*        pInFlow=0,
                   Flow2D*      pInFlow2D=0,
                   FP*      y1=0,
                   int          btt=TCT_No_Turbulence_2D);
#endif // _UNIFORM_MESH_

    int CloseContour2D(char*    name,
                       ulong    bt,
                       Flow*    pInFlow=0,
                       Flow2D*  pInFlow2D=0,
                       FP*  y=0,
                       ulong    btt=TCT_No_Turbulence_2D);

    int InsertBound2D(char* name,
                      unsigned int nb,
                      unsigned int x,
                      unsigned int y,
                      ulong        bt,
                      Flow*        pInFlow=0,
                      Flow2D*      pInFlow2D=0,
                      FP*      y1=0,
                      ulong        btt=TCT_No_Turbulence_2D);

#ifndef _UNIFORM_MESH_
    int InsertBound2D(char* name,
                      unsigned int nb,
                      Mesh2D* p_mesh,       //  anisotropic mesh
                      FP       f_x,
                      FP       f_y,
                      ulong        bt,
                      Flow*        pInFlow=0,
                      Flow2D*      pInFlow2D=0,
                      FP*      y1=0,
                      int          btt=TCT_No_Turbulence_2D);
#endif // _UNIFORM_MESH_
    
    int    DelBound2D(int nb);
    int    DelBound2D(char name);
    int    SetBounds(int MaterialID=0);
    int    SetBounds(UArray<FlowNode2D< FP, NUM_COMPONENTS>* >* node_array, int MaterialID=0);
#ifdef _UNIFORM_MESH_
    int    GetCurrentX();
    int    GetCurrentY();
    void   SetCurrentX(int x);
    void   SetCurrentY(int y);
    void   SetFirstX(int x);
    void   SetFirstY(int y);
#endif // _UNIFORM_MESH_
    FP GetCurrentFX();
    FP GetCurrentFY();
    void   SetCurrentFX(FP x);
    void   SetCurrentFY(FP y);
    void   SetFirstFX(FP x);
    void   SetFirstFY(FP y);
    int    GetNumBounds();
    Bound2D* GetBound(int nb);
    int    IsContourClosed();
    int    IsContourActivate();
    void   ClearBoundContour();
    int    RotateBoundContour2D(FP x0, 
                                FP y0, 
                                FP angle);
    char*  GetBoundContourName() {return BoundContourName;}

#ifndef _UNIFORM_MESH_
Mesh2D* GetMesh() {
        return pMesh;
    }
#endif // _UNIFORM_MESH_
};

// <------------- 2D --------------->           
#endif // _hyper_flow_bound_contour_hpp_



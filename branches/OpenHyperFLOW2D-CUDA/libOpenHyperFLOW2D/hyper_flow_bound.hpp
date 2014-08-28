/*******************************************************************************
*   OpenHyperFLOW2D                                                            *
*                                                                              *
*   Version  1.0.1                                                             *
*   Copyright (C)  1995-2014 by Serge A. Suchkov                               *
*   Copyright policy: LGPL V3                                                  *
*                                                                              *
*   last update: 01/07/2014                                                    *
*******************************************************************************/
#ifndef _hyper_flow_bound_hpp_
#define _hyper_flow_bound_hpp_

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

const double   pi  = M_PI;
const double   _2pi=2*pi;

enum ChemicalReactionsModel {
     NO_REACTIONS,
     CRM_ZELDOVICH,
     CRM_ARRENIUS,
     CRM_EDM,
};

struct ChemicalReactionsModelData2D {
 double K0;         // Stechiometric ratio (OX/Fuel)
 double gamma;      // Chemical reaction "completation" coefficient
 double Tf;         // Ignition temperature 
// Fuel              
 double R_Fuel;
 double H_Fuel;
 Table* Entalpy_Fuel;
 Table* lam_Fuel; 
 Table* mu_Fuel;  
 Table* Cp_Fuel;  
// Oxidizer 
 double R_OX;
 double H_OX;
 Table* Entalpy_OX;
 Table* lam_OX;
 Table* mu_OX; 
 Table* Cp_OX; 
 // Air (inert gas)
 double R_air; 
 double H_air; 
 Table* Entalpy_air;
 Table* lam_air;
 Table* mu_air; 
 Table* Cp_air; 
 // Combustion products
 double R_cp; 
 double H_cp; 
 Table* Entalpy_cp;
 Table* lam_cp;
 Table* mu_cp; 
 Table* Cp_cp; 
};

struct cudaChemicalReactionsModelData2D {
 double  K0;         // Stechiometric ratio (OX/Fuel)
 double  gamma;      // Chemical reaction "completation" coefficient
 double  Tf;         // Ignition temperature 
 int     TablesSize;
 double* Tg; 
// Fuel              
 double  R_Fuel;
 double  H_Fuel;
 double* Entalpy_Fuel;
 double* lam_Fuel; 
 double* mu_Fuel;  
 double* Cp_Fuel;  
// Oxidizer 
 double R_OX;
 double H_OX;
 double* Entalpy_OX;
 double* lam_OX;
 double* mu_OX; 
 double* Cp_OX; 
 // Air (inert gas)
 double R_air; 
 double H_air; 
 double* Entalpy_air;
 double* lam_air;
 double* mu_air; 
 double* Cp_air; 
 // Combustion products
 double R_cp; 
 double H_cp; 
 double* Entalpy_cp;
 double* lam_cp;
 double* mu_cp; 
 double* Cp_cp; 
};

enum BoundState {
    BND_INACTIVE, BND_OK, BND_ERR
};

enum BoundMoveType {
     BMT_TAIL,
     BMT_NOTAIL
};

//////////////////////////////////////////////////
//  Bound2D object                              //
//////////////////////////////////////////////////
class Bound2D {
    char*                           BoundName;            // Name of bound
    UMatrix2D< FlowNode2D< double, NUM_COMPONENTS> >* pMFN; //  N-components 
    XY<double>                      fStart;               //  X-Y coordinate  start  point of bound (meters)
    XY<double>                      fEnd;                 //  X-Y coordinate  end  point of bound   (meters)
#ifdef _UNIFORM_MESH_
    XY<unsigned int>                Start;                //  X-Y coordinate  start  point of bound (nodes)
    XY<unsigned int>                End;                  //  X-Y coordinate  end  point of bound   (nodes)
    double                          dx,dy;                //  dx and dy step...  
#else
    Mesh2D*                         pMesh;                //  anisotropic mesh
#endif // _UNIFORM_MESH_
    Flow2D*                         pBoundFlow2D;         //  Flow2D object (if it need) for this bound
    Flow*                           pBoundFlow;           //  Flow  object (if it need) for this bound
    double*                         pY;                   //  Y[a] - Y1,Y2,Y3,...,Ya
    BoundState                      bs;                   //  Bound state
    int                             BNT;                  //  Bound condition bit flag set 
    int                             BTC;                  //  Bound turbulence condition bit flag set 

    void  InitBound(char* name, 
                    UMatrix2D< FlowNode2D< double, NUM_COMPONENTS> >*,
#ifndef _UNIFORM_MESH_
                    Mesh2D* p_mesh,                       //  anisotropic mesh
#endif // _UNIFORM_MESH_
                    int     bt,double* y=0,int  btc=TCT_No_Turbulence_2D);
public:
    Bound2D(char* name, 
          UMatrix2D< FlowNode2D< double, NUM_COMPONENTS> >*,
#ifndef _UNIFORM_MESH_
          Mesh2D* p_mesh,                                 //  anisotropic mesh
#endif // _UNIFORM_MESH_
          double  x1,
          double  y1,
          double  x2,
          double  y2,
          int     bt,
          Flow*   pInFlow,
          double* Y=0,
          int     btc=TCT_No_Turbulence_2D);

    Bound2D(char* name, 
          UMatrix2D< FlowNode2D< double, NUM_COMPONENTS> >*,
#ifndef _UNIFORM_MESH_
          Mesh2D* p_mesh,                                 //  anisotropic mesh
#endif // _UNIFORM_MESH_
          double  x1,
          double  y1,
          double  x2,
          double  y2,
          int     bt,
          Flow2D* pInFlow2D,
          double* Y=0,
          int     btc=TCT_No_Turbulence_2D);

#ifdef _UNIFORM_MESH_
    Bound2D(char* name, 
          UMatrix2D< FlowNode2D< double, NUM_COMPONENTS> >*,
          unsigned int x1,
          unsigned int y1,
          unsigned int x2,
          unsigned int y2,
          int          bt,
          Flow*        pInFlow,
          double*      Y=0,
          int          btc=TCT_No_Turbulence_2D);

    Bound2D(char* name, 
          UMatrix2D< FlowNode2D< double, NUM_COMPONENTS> >*,
          unsigned int x1,
          unsigned int y1,
          unsigned int x2,
          unsigned int y2,
          int          bt,
          Flow2D*      pInFlow2D,
          double*      Y=0,
          int          btc=TCT_No_Turbulence_2D);
#endif // _UNIFORM_MESH_

    Bound2D(char* name, 
          UMatrix2D< FlowNode2D< double, NUM_COMPONENTS> >*,
#ifndef _UNIFORM_MESH_
          Mesh2D* p_mesh,                                 //  anisotropic mesh
#endif // _UNIFORM_MESH_
          XY<double>*  p_start,
          XY<double>*  p_end,
          int          bt,
          Flow*        pInFlow,
          double*      Y=0,
          int          btc=TCT_No_Turbulence_2D);

    Bound2D(char* name, 
          UMatrix2D< FlowNode2D< double, NUM_COMPONENTS> >*,
#ifndef _UNIFORM_MESH_
          Mesh2D*      p_mesh,                          //  anisotropic mesh
#endif // _UNIFORM_MESH_
          XY<double>*  p_start,
          XY<double>*  p_end,
          int          bt,
          Flow2D*      pInFlow2D,
          double*      Y=0,
          int          btc=TCT_No_Turbulence_2D);

#ifdef _UNIFORM_MESH_
    Bound2D(char* name, 
          UMatrix2D< FlowNode2D< double, NUM_COMPONENTS> >*,
          XY<unsigned int>*  p_start,
          XY<unsigned int>*  p_end,
          int                bt,
          Flow*              pInFlow=0,
          double*            Y=0,
          int                btc=TCT_No_Turbulence_2D);

    Bound2D(char* name, 
          UMatrix2D< FlowNode2D< double, NUM_COMPONENTS> >*,
          XY<unsigned int>*  p_start,
          XY<unsigned int>*  p_end,
          int                bt,
          Flow2D*            pInFlow2D,
          double*            Y=0,
          int                btc=TCT_No_Turbulence_2D);
#endif // _UNIFORM_MESH_

    virtual     ~Bound2D();

    BoundState   SetBound(int MaterialID=0);
    BoundState   SetBound(UArray< FlowNode2D<double, NUM_COMPONENTS>* >* node_array, int MaterialID=0);
    BoundState   GetBoundState();
    int          GetBoundCond();
    double*      GetYArray();
    Flow*        GetBoundFlow();
    Flow2D*      GetBoundFlow2D();

    char*        GetBoundName() {return BoundName;}
    unsigned int GetNumComponents();
#ifdef _UNIFORM_MESH_
    unsigned int GetStartX();       
    unsigned int GetStartY();
    unsigned int GetEndX();      
    unsigned int GetEndY();
#endif // _UNIFORM_MESH_
    double       GetStartFX();
    double       GetStartFY();
    
    double       GetEndFX();      
    double       GetEndFY();
#ifdef _UNIFORM_MESH_
    void         SetStartX(unsigned int x);
    void         SetStartY(unsigned int y);
    
    void         SetEndX(unsigned int x);
    void         SetEndY(unsigned int y);
    
    void         SetStartXY(XY<unsigned int>* xy);
    void         SetEndXY(XY<unsigned int>* xy);
#endif // _UNIFORM_MESH_
    void         SetStartFX(double x);
    void         SetStartFY(double y);
    
    void         SetEndFX(double x);
    void         SetEndFY(double y);
    
    void         SetStartFXY(XY<double>* xy);
    void         SetEndFXY(XY<double>* xy);
    
    int          RotateBound2D(double x0, 
                               double y0, 
                               double angle);
    
    int          TestRotateBound2D(double x0,
                                   double y0,
                                   double angle);
#ifndef _UNIFORM_MESH_
    Mesh2D*
    GetMesh() {
        return pMesh;
    }
#endif // _UNIFORM_MESH_
};
// <------------- 2D --------------->           

#endif // _hyper_flow_bound_hpp_



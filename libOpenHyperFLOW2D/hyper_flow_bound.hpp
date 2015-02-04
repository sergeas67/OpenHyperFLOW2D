/*******************************************************************************
*   OpenHyperFLOW2D                                                            *
*                                                                              *
*   Version  1.0.2                                                             *
*   Copyright (C)  1995-2015 by Serge A. Suchkov                               *
*   Copyright policy: LGPL V3                                                  *
*   http://openhyperflow2d.googlecode.com                                      *
*                                                                              *
*   last update: 01/02/2015                                                    *
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

#ifndef  max
    #define max(a,b) (((a)>(b))?(a):(b))
#endif 

#ifndef  min
    #define min(a,b) (((a)<(b))?(a):(b))
#endif 

const FP   pi  = M_PI;
const FP   _2pi=2*pi;

enum ChemicalReactionsModel {
     NO_REACTIONS,
     CRM_ZELDOVICH,
     CRM_ARRENIUS,
     CRM_EDM,
};

struct ChemicalReactionsModelData2D {
 FP K0;         // Stechiometric ratio (OX/Fuel)
 FP gamma;      // Chemical reaction "completation" coefficient
 FP Tf;         // Ignition temperature 
// Fuel              
 FP R_Fuel;
 FP H_Fuel;
 Table* Entalpy_Fuel;
 Table* lam_Fuel; 
 Table* mu_Fuel;  
 Table* Cp_Fuel;  
// Oxidizer 
 FP R_OX;
 FP H_OX;
 Table* Entalpy_OX;
 Table* lam_OX;
 Table* mu_OX; 
 Table* Cp_OX; 
 // Air (inert gas)
 FP R_air; 
 FP H_air; 
 Table* Entalpy_air;
 Table* lam_air;
 Table* mu_air; 
 Table* Cp_air; 
 // Combustion products
 FP R_cp; 
 FP H_cp; 
 Table* Entalpy_cp;
 Table* lam_cp;
 Table* mu_cp; 
 Table* Cp_cp; 
};

struct cudaChemicalReactionsModelData2D {
 FP  K0;         // Stechiometric ratio (OX/Fuel)
 FP  gamma;      // Chemical reaction "completation" coefficient
 FP  Tf;         // Ignition temperature 
 int     TablesSize;
 FP* Tg; 
// Fuel              
 FP  R_Fuel;
 FP  H_Fuel;
 FP* Entalpy_Fuel;
 FP* lam_Fuel; 
 FP* mu_Fuel;  
 FP* Cp_Fuel;  
// Oxidizer 
 FP R_OX;
 FP H_OX;
 FP* Entalpy_OX;
 FP* lam_OX;
 FP* mu_OX; 
 FP* Cp_OX; 
 // Air (inert gas)
 FP R_air; 
 FP H_air; 
 FP* Entalpy_air;
 FP* lam_air;
 FP* mu_air; 
 FP* Cp_air; 
 // Combustion products
 FP R_cp; 
 FP H_cp; 
 FP* Entalpy_cp;
 FP* lam_cp;
 FP* mu_cp; 
 FP* Cp_cp; 
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
    char*                           BoundName;        // Name of bound
    UMatrix2D< FlowNode2D< FP, NUM_COMPONENTS> >* pMFN; //  N-components 
    XY<FP>                      fStart;               //  X-Y coordinate  start  point of bound (meters)
    XY<FP>                      fEnd;                 //  X-Y coordinate  end  point of bound   (meters)
    XY<unsigned int>            Start;                //  X-Y coordinate  start  point of bound (nodes)
    XY<unsigned int>            End;                  //  X-Y coordinate  end  point of bound   (nodes)
    FP                          dx,dy;                //  dx and dy step...  
    Flow2D*                     pBoundFlow2D;         //  Flow2D object (if it need) for this bound
    Flow*                       pBoundFlow;           //  Flow  object (if it need) for this bound
    FP*                         pY;                   //  Y[a] - Y1,Y2,Y3,...,Ya
    BoundState                  bs;                   //  Bound state
    int                         BNT;                  //  Bound condition bit flag set 
    int                         BTC;                  //  Bound turbulence condition bit flag set 

    void  InitBound(char* name, 
                    UMatrix2D< FlowNode2D< FP, NUM_COMPONENTS> >*,
                    int     bt,FP* y=0,int  btc=TCT_No_Turbulence_2D);
public:
    Bound2D(char* name, 
          UMatrix2D< FlowNode2D< FP, NUM_COMPONENTS> >*,
          FP  x1,
          FP  y1,
          FP  x2,
          FP  y2,
          int     bt,
          Flow*   pInFlow,
          FP* Y=0,
          int     btc=TCT_No_Turbulence_2D);

    Bound2D(char* name, 
          UMatrix2D< FlowNode2D< FP, NUM_COMPONENTS> >*,
          FP  x1,
          FP  y1,
          FP  x2,
          FP  y2,
          int     bt,
          Flow2D* pInFlow2D,
          FP* Y=0,
          int     btc=TCT_No_Turbulence_2D);

    Bound2D(char* name, 
          UMatrix2D< FlowNode2D< FP, NUM_COMPONENTS> >*,
          unsigned int x1,
          unsigned int y1,
          unsigned int x2,
          unsigned int y2,
          int          bt,
          Flow*        pInFlow,
          FP*      Y=0,
          int          btc=TCT_No_Turbulence_2D);

    Bound2D(char* name, 
          UMatrix2D< FlowNode2D< FP, NUM_COMPONENTS> >*,
          unsigned int x1,
          unsigned int y1,
          unsigned int x2,
          unsigned int y2,
          int          bt,
          Flow2D*      pInFlow2D,
          FP*      Y=0,
          int          btc=TCT_No_Turbulence_2D);

    Bound2D(char* name, 
          UMatrix2D< FlowNode2D< FP, NUM_COMPONENTS> >*,
          XY<FP>*  p_start,
          XY<FP>*  p_end,
          int          bt,
          Flow*        pInFlow,
          FP*      Y=0,
          int          btc=TCT_No_Turbulence_2D);

    Bound2D(char* name, 
          UMatrix2D< FlowNode2D< FP, NUM_COMPONENTS> >*,
          XY<FP>*  p_start,
          XY<FP>*  p_end,
          int          bt,
          Flow2D*      pInFlow2D,
          FP*      Y=0,
          int          btc=TCT_No_Turbulence_2D);

    Bound2D(char* name, 
          UMatrix2D< FlowNode2D< FP, NUM_COMPONENTS> >*,
          XY<unsigned int>*  p_start,
          XY<unsigned int>*  p_end,
          int                bt,
          Flow*              pInFlow=0,
          FP*            Y=0,
          int                btc=TCT_No_Turbulence_2D);

    Bound2D(char* name, 
          UMatrix2D< FlowNode2D< FP, NUM_COMPONENTS> >*,
          XY<unsigned int>*  p_start,
          XY<unsigned int>*  p_end,
          int                bt,
          Flow2D*            pInFlow2D,
          FP*            Y=0,
          int                btc=TCT_No_Turbulence_2D);

    virtual     ~Bound2D();

    BoundState   SetBound(int MaterialID=0);
    BoundState   SetBound(UArray< FlowNode2D<FP, NUM_COMPONENTS>* >* node_array, int MaterialID=0);
    BoundState   GetBoundState();
    int          GetBoundCond();
    FP*      GetYArray();
    Flow*        GetBoundFlow();
    Flow2D*      GetBoundFlow2D();

    char*        GetBoundName() {return BoundName;}
    unsigned int GetNumComponents();
    unsigned int GetStartX();       
    unsigned int GetStartY();
    unsigned int GetEndX();      
    unsigned int GetEndY();
    FP       GetStartFX();
    FP       GetStartFY();
    
    FP       GetEndFX();      
    FP       GetEndFY();
    void         SetStartX(unsigned int x);
    void         SetStartY(unsigned int y);
    
    void         SetEndX(unsigned int x);
    void         SetEndY(unsigned int y);
    
    void         SetStartXY(XY<unsigned int>* xy);
    void         SetEndXY(XY<unsigned int>* xy);
    void         SetStartFX(FP x);
    void         SetStartFY(FP y);
    
    void         SetEndFX(FP x);
    void         SetEndFY(FP y);
    
    void         SetStartFXY(XY<FP>* xy);
    void         SetEndFXY(XY<FP>* xy);
    
    int          RotateBound2D(FP x0, 
                               FP y0, 
                               FP angle);
    
    int          TestRotateBound2D(FP x0,
                                   FP y0,
                                   FP angle);
};
// <------------- 2D --------------->           

#endif // _hyper_flow_bound_hpp_



/*******************************************************************************
*   OpenHyperFLOW2D                                                            *
*                                                                              *
*   Version  1.0.3                                                             *
*   Copyright (C)  1995-2016 by Serge A. Suchkov                               *
*   Copyright policy: LGPL V3                                                  *
*   http://openhyperflow2d.googlecode.com                                      *
*                                                                              *
*   last update: 04/07/2016                                                    *
*******************************************************************************/
#ifndef _hyper_flow_airfoil_hpp_
#define _hyper_flow_airfoil_hpp_

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

enum AirfoilType {
 AFT_NACA=3700,
 AFT_TsAGI,
};

//////////////////////////////////////////////////////////
// Airfoil                                              //
//////////////////////////////////////////////////////////
/*
# "Bezier presentation of airfoils", by Wolfgang Boehm,
#  Computer Aided Geometric Design 4 (1987) pp 17-22.
#
#				Gershon Elber, Nov. 1992
#
*/
// Macros set for Airfoil primitive
#define bez_d4_i0(x) pow((1.-x),4)
#define bez_d4_i1(x) 4.*pow((1.-x),3)*x
#define bez_d4_i2(x) 6.*pow((1.-x),2)*x*x
#define bez_d4_i3(x) 4.*(1.-x)*pow(x,3)
#define bez_d4_i4(x) pow(x,4)

#define bez_d8_i0(x) pow((1.-x),8)
#define bez_d8_i1(x) 8.*pow((1.-x),7)*x
#define bez_d8_i2(x) 28.*pow((1.-x),6)*x*x
#define bez_d8_i3(x) 56.*pow((1.-x),5)*pow(x,3)
#define bez_d8_i4(x) 70.*pow((1.-x),4)*pow(x,4)
#define bez_d8_i5(x) 56.*pow((1.-x),3)*pow(x,5)
#define bez_d8_i6(x) 28.*pow((1.-x),2)*pow(x,6)
#define bez_d8_i7(x) 8.*(1.-x)*pow(x,7)
#define bez_d8_i8(x) pow(x,8)
//////////////////////////////////////////////////
//  SolidBoundAirfoil2D object                  //
//////////////////////////////////////////////////
class SolidBoundAirfoil2D : protected BoundContour2D,
                            protected Area2D {
    FP pp,mm;
    FP mean_y(FP t);
    FP mean_x(FP t);
    FP z_x(FP x);
    FP z_y(FP x, FP tk);
    FP airfoil_y1(FP t,FP  thick);
    FP airfoil_y2(FP t,FP  thick);
    FP airfoil_y(FP t);
    FP airfoil_x(FP t);
    AirfoilType  airfoil_type;
    char*  SolidBoundAirfoilName;

public:

  SolidBoundAirfoil2D(char* name,
                      UMatrix2D< FlowNode2D< FP, NUM_COMPONENTS> >* JM, // Computation area reference
                      FP  x,                      // Start profile 
                      FP  y,                      // coordinates (x,y)
                      FP  m_m,                    //
                      FP  p_p,                    //
                      FP  thick,                  // airfoil thick
                      FP  dx,                     // dx - step
                      FP  dy,                     // dy - step
                      ulong   ct,                 // condition type
                      Flow2D* pInFlow2D,          // init flow2d object on circle bound
                      FP* Y=NULL,                 // component matrix
                      ulong   bctt=TCT_No_Turbulence_2D, // Bound contour turbulence type
                      FP  scale=1.,               // airfoil scale
                      FP  attack_angle=0,         // Angle of attack
                      ostream* dbg_output=NULL);
  
  SolidBoundAirfoil2D(char* name,
                      UMatrix2D< FlowNode2D< FP, NUM_COMPONENTS> >* JM, // Computation area reference
                      FP  x,                      // Start profile 
                      FP  y,                      // coordinates (x,y)
                      UArray< XY<FP> >* airfoil_points,
                      FP  dx,                     // dx - step
                      FP  dy,                     // dy - step
                      ulong   ct,                 // condition type
                      Flow2D* pInFlow2D,          // init flow2d object on circle bound
                      FP* Y=NULL,                 // component matrix
                      ulong   bctt=TCT_No_Turbulence_2D, // Bound contour turbulence type
                      FP  scale=1.,               // airfoil scale
                      FP  attack_angle=0,         // Angle of attack
                      ostream* dbg_output=NULL);
  
  SolidBoundAirfoil2D(char* name,
                      UMatrix2D< FlowNode2D< FP, NUM_COMPONENTS> >* JM, // Computation area reference
                      FP  x,                      // Start profile 
                      FP  y,                      // coordinates (x,y)
                      InputData* airfoil_input_data, // Airfoil input data
                      FP  dx,                     // dx - step
                      FP  dy,                     // dy - step
                      ulong   ct,                 // condition type
                      Flow2D* pInFlow2D,          // init flow2d object on circle bound
                      FP* Y=NULL,                 // component matrix
                      ulong   bctt=TCT_No_Turbulence_2D, // Bound contour turbulence type
                      FP  scale=1.,               // airfoil scale
                      FP  attack_angle=0,         // Angle of attack
                      ostream* dbg_output=NULL);


  char* GetSolidBoundAirfoilName2D() {
         return SolidBoundAirfoilName;
    }

  AirfoilType GetAirfoilType() { return airfoil_type; }
};

#endif // _hyper_flow_airfoil_hpp_



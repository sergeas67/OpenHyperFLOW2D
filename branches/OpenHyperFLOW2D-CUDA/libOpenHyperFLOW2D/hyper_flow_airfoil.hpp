/*******************************************************************************
*   OpenHyperFLOW2D                                                            *
*                                                                              *
*   Version  1.0.1                                                             *
*   Copyright (C)  1995-2014 by Serge A. Suchkov                               *
*   Copyright policy: LGPL V3                                                  *
*                                                                              *
*   last update: 16/01/2014                                                    *
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

#ifndef _UNIFORM_MESH_
#include "libMesh/mesh.hpp"
#endif  // _UNIFORM_MESH_

#ifndef  max
    #define max(a,b) (((a)>(b))?(a):(b))
#endif 

#ifndef  min
    #define min(a,b) (((a)<(b))?(a):(b))
#endif 

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
    double pp,mm;
    double mean_y(double t);
    double mean_x(double t);
    double z_x(double x);
    double z_y(double x, double tk);
    double airfoil_y1(double t,double  thick);
    double airfoil_y2(double t,double  thick);
    double airfoil_y(double t);
    double airfoil_x(double t);
    char*  SolidBoundAirfoilName;

public:

  SolidBoundAirfoil2D(char* name,
                      UMatrix2D< FlowNode2D< double, NUM_COMPONENTS> >* JM, // Computation area reference
                      double  x,                      // Start profile 
                      double  y,                      // coordinates (x,y)
                      double  m_m,                    //
                      double  p_p,                    //
                      double  thick,                  // airfoil thick
#ifdef _UNIFORM_MESH_
                      double  dx,                     // dx - step
                      double  dy,                     // dy - step
#else 
                      unsigned int num_segments,
                      Mesh2D* p_mesh,
#endif // _UNIFORM_MESH_
                      ulong   ct,                     // condition type
                      Flow2D* pInFlow2D,              // init flow2d object on circle bound
                      double* Y=NULL,                 // component matrix
                      ulong   bctt=TCT_No_Turbulence_2D, // Bound contour turbulence type
                      double  scale=1.,               // airfoil scale
                      double  attack_angle=0,         // Angle of attack
                      ostream* dbg_output=NULL);

  SolidBoundAirfoil2D(char* name,
                    UMatrix2D< FlowNode2D< double, NUM_COMPONENTS> >* JM, // Computation area reference
                    double  x,                      // Start profile 
                    double  y,                      // coordinates (x,y)
                    void*   airfoil_points,
                    //UArray< XY < double > >* airfoil_points,
#ifdef _UNIFORM_MESH_
                    double  dx,                     // dx - step
                    double  dy,                     // dy - step
#else 
                    unsigned int num_segments,
                    Mesh2D* p_mesh,
#endif // _UNIFORM_MESH_
                    ulong   ct,                     // condition type
                    Flow2D* pInFlow2D,              // init flow2d object on circle bound
                    double* Y=NULL,                 // component matrix
                    ulong   bctt=TCT_No_Turbulence_2D, // Bound contour turbulence type
                    double  scale=1.,               // airfoil scale
                    double  attack_angle=0,         // Angle of attack
                    ostream* dbg_output=NULL);


  char* GetSolidBoundAirfoilName2D() {
         return SolidBoundAirfoilName;
    }
#ifndef _UNIFORM_MESH_
    Mesh2D*    GetMesh() {
        return pMesh;
    }
#endif // _UNIFORM_MESH_
};

#endif // _hyper_flow_airfoil_hpp_



/*******************************************************************************
*   OpenHyperFLOW2D                                                            *
*                                                                              *
*   Version  1.0.1                                                             *
*   Copyright (C)  1995-2014 by Serge A. Suchkov                               *
*   Copyright policy: LGPL V3                                                  *
*                                                                              *
*   last update: 01/07/2014                                                    *
*******************************************************************************/
#ifndef _out_cfd_param_hpp_
#define _out_cfd_param_hpp_

#ifdef _OPENMP
#include <omp.h>
#endif //OPENMP

#include "utl/uarray.hpp"
#include "utl/umatrix2d.hpp"

#include "obj_data/obj_data.hpp"
#include "libOpenHyperFLOW2D/hyper_flow_node.hpp"

#ifdef _DEBUG_0
#include "libExcept/except.hpp"
#endif // _DEBUG_0

struct XCut {
FP x0;  // X start point
FP y0;  // Y start point
FP dy;  // cut width
};

FP T_asterisk(FlowNode2D<FP,NUM_COMPONENTS>* node );
FP p_asterisk(FlowNode2D<FP,NUM_COMPONENTS>* node );
FP Schliren(FlowNode2D<FP,NUM_COMPONENTS>* node );

FP CalcaveragePressure2D(UMatrix2D< FlowNode2D<FP,NUM_COMPONENTS> >* pJ,
                             FP x0, // initial point of probed area
                             FP l,  // length of probed area
                             FP d   // diameter (or cross-section size) of probed area
                            );

FP CalcaverageTemperature2D(UMatrix2D< FlowNode2D<FP,NUM_COMPONENTS> >* pJ,
                                FP x0, // initial point of probed area
                                FP l,  // length of probed area
                                FP d,  // diameter (or cross-section size) of probed area
                                int is_mid_ethalpy_T = 0
                                );

FP CalcMassFlowRateX2D(UMatrix2D< FlowNode2D<FP,NUM_COMPONENTS> >* pJ,
                           FP x0, // initial point of probed area
                           FP y0,
                           FP dy   // diameter (or cross-section size) of probed area)
                          );

FP CalcXForceYSym2D(UMatrix2D< FlowNode2D<FP,NUM_COMPONENTS> >* pJ,
                         FP x0, // initial point of probed area
                         FP l,  // length of probed area
                         FP d   // diameter (or cross-section size) of probed area
                        );

FP CalcXForce2D(UMatrix2D< FlowNode2D<FP,NUM_COMPONENTS> >* pJ,
                    FP x0, // initial point of probed area
                    FP l,  // length of probed area
                    FP d  // diameter (or cross-section size) of probed area
                   );

FP CalcYForce2D(UMatrix2D< FlowNode2D<FP,NUM_COMPONENTS> >* pJ,
                    FP x0, // initial point of probed area
                    FP l,  // length of probed area
                    FP d  // diameter (or cross-section size) of probed area
                   );
FP CalcXForce2D(UMatrix2D< FlowNode2D<FP,NUM_COMPONENTS> >* pJ,
                    FP x0, // X initial point of probed area
                    FP y0, // Y initial point of probed area
                    FP dx,  // X size of probed area
                    FP dy  // Y of probed area
                    );

FP CalcYForce2D(UMatrix2D< FlowNode2D<FP,NUM_COMPONENTS> >* pJ,
                    FP x0, // X initial point of probed area
                    FP y0, // Y initial point of probed area
                    FP dx, // X size of probed area
                    FP dy  // Y of probed area
                    );

void SaveXHeatFlux2D(ofstream* OutputData,
                     UMatrix2D< FlowNode2D<FP,NUM_COMPONENTS> >* pJ,
                     Flow2D* TestFlow, 
                     FP  Ts,
                     int     y_max,
                     int     y_min);

void SaveYHeatFlux2D(ofstream* OutputData,
                     UMatrix2D< FlowNode2D<FP,NUM_COMPONENTS> >* pJ,
                     FP  Ts);

FP   Calc_XSigmaFi(UMatrix2D< FlowNode2D<FP,NUM_COMPONENTS> >* pJ,
                      FP x0,
                      FP y0,
                      FP dy,
                      Flow2D* pF);

FP   Calc_Cp(FlowNode2D<FP,NUM_COMPONENTS>* CurrentNode, Flow2D* pF);
FP   Calc_Cx_2D(UMatrix2D< FlowNode2D<FP,NUM_COMPONENTS> >* pJ, 
                    FP x, 
                    FP y,
                    FP d_x, 
                    FP d_y,
                    Flow2D* pF);

FP   Calc_Cy_2D(UMatrix2D< FlowNode2D<FP,NUM_COMPONENTS> >* pJ, 
                    FP x, 
                    FP y,
                    FP d_x, 
                    FP d_y,
                    Flow2D* pF);

FP Calc_Cv(UMatrix2D< FlowNode2D<FP,NUM_COMPONENTS> >* pJ,               
               FP x0, // initial point of probed area                        
               FP y0,                                                        
               FP dy, // diameter (or cross-section size) of probed area)
               FP p_amb,
               Flow2D* pF);
                                                                                 

FP Calc_Cd(UMatrix2D< FlowNode2D<FP,NUM_COMPONENTS> >* pJ,               
               FP x0, // initial point of probed area                        
               FP y0,                                                        
               FP dy,   // diameter (or cross-section size) of probed area)   
               Flow2D* pF);

FP CalcArea2D(UMatrix2D< FlowNode2D<FP,NUM_COMPONENTS> >* pJ,             
                  FP x0,  // initial X  point of probed area                  
                  FP y0,  // initial Y  point of probed area                  
                  FP dy   // diameter (or cross-section size) of probed area) 
                 );                                                             
#endif // _out_cfd_param_hpp_




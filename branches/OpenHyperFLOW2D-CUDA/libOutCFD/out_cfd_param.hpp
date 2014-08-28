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
double x0;  // X start point
double y0;  // Y start point
double dy;  // cut width
};

double T_asterisk(FlowNode2D<double,NUM_COMPONENTS>* node );
double p_asterisk(FlowNode2D<double,NUM_COMPONENTS>* node );
double Schliren(FlowNode2D<double,NUM_COMPONENTS>* node );

double CalcaveragePressure2D(UMatrix2D< FlowNode2D<double,NUM_COMPONENTS> >* pJ,
                             double x0, // initial point of probed area
                             double l,  // length of probed area
                             double d   // diameter (or cross-section size) of probed area
                            );

double CalcaverageTemperature2D(UMatrix2D< FlowNode2D<double,NUM_COMPONENTS> >* pJ,
                                double x0, // initial point of probed area
                                double l,  // length of probed area
                                double d,  // diameter (or cross-section size) of probed area
                                int is_mid_ethalpy_T = 0
                                );

double CalcMassFlowRateX2D(UMatrix2D< FlowNode2D<double,NUM_COMPONENTS> >* pJ,
                           double x0, // initial point of probed area
                           double y0,
                           double dy   // diameter (or cross-section size) of probed area)
                          );

double CalcXForceYSym2D(UMatrix2D< FlowNode2D<double,NUM_COMPONENTS> >* pJ,
                         double x0, // initial point of probed area
                         double l,  // length of probed area
                         double d   // diameter (or cross-section size) of probed area
                        );

double CalcXForce2D(UMatrix2D< FlowNode2D<double,NUM_COMPONENTS> >* pJ,
                    double x0, // initial point of probed area
                    double l,  // length of probed area
                    double d  // diameter (or cross-section size) of probed area
                   );

double CalcYForce2D(UMatrix2D< FlowNode2D<double,NUM_COMPONENTS> >* pJ,
                    double x0, // initial point of probed area
                    double l,  // length of probed area
                    double d  // diameter (or cross-section size) of probed area
                   );
double CalcXForce2D(UMatrix2D< FlowNode2D<double,NUM_COMPONENTS> >* pJ,
                    double x0, // X initial point of probed area
                    double y0, // Y initial point of probed area
                    double dx,  // X size of probed area
                    double dy  // Y of probed area
                    );

double CalcYForce2D(UMatrix2D< FlowNode2D<double,NUM_COMPONENTS> >* pJ,
                    double x0, // X initial point of probed area
                    double y0, // Y initial point of probed area
                    double dx, // X size of probed area
                    double dy  // Y of probed area
                    );

void SaveXHeatFlux2D(ofstream* OutputData,
                     UMatrix2D< FlowNode2D<double,NUM_COMPONENTS> >* pJ,
                     double  Ts);

void SaveYHeatFlux2D(ofstream* OutputData,
                     UMatrix2D< FlowNode2D<double,NUM_COMPONENTS> >* pJ,
                     double  Ts);

double   Calc_XSigmaFi(UMatrix2D< FlowNode2D<double,NUM_COMPONENTS> >* pJ,
                      double x0,
                      double y0,
                      double dy,
                      Flow2D* pF);

double   Calc_Cp(UMatrix2D< FlowNode2D<double,NUM_COMPONENTS> >* pJ, int i, int j, Flow2D* pF);
double   Calc_Cx_2D(UMatrix2D< FlowNode2D<double,NUM_COMPONENTS> >* pJ, 
                    double x, 
                    double y,
                    double d_x, 
                    double d_y,
                    Flow2D* pF);

double   Calc_Cy_2D(UMatrix2D< FlowNode2D<double,NUM_COMPONENTS> >* pJ, 
                    double x, 
                    double y,
                    double d_x, 
                    double d_y,
                    Flow2D* pF);
#endif //_out_cfd_param_hpp_




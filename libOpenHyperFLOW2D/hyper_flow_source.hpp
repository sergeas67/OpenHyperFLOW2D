/*******************************************************************************
*   OpenHyperFLOW2D                                                            *
*                                                                              *
*   Version  2.0.1                                                             *
*   Copyright (C)  1995-2016 by Serge A. Suchkov                               *
*   Copyright policy: LGPL V3                                                  *
*   http://github.com/sergeas67/openhyperflow2d                                *
*                                                                              *
*   last update: 14/04/2016                                                    *
*******************************************************************************/
#ifndef _hyper_flow_source_hpp_
#define _hyper_flow_source_hpp_

#include "utl/uarray.hpp"
#include "utl/umatrix2d.hpp"
#include "obj_data/obj_data.hpp"
#include "libOpenHyperFLOW2D/hyper_flow_node.hpp"

#ifndef ComputationalMatrix2D
#define ComputationalMatrix2D  UMatrix2D< FlowNode2D<FP,NUM_COMPONENTS> >
#endif //ComputationalMatrix
// <------------- 2D --------------->           

class Source2D
{
 
 ComputationalMatrix2D* F;                //   Reference to computational matrix
 int                    sx,sy,ex,ey;      //   Start and end points (in nodes)
 int                    c_index;          //   Component index
 FP                     Cp, M_s0, T, T_f; //   Cp, Ms, T, Tf of source
 int                    StartSrcIter;     //   Start source iteration

public:
  
  Source2D(ComputationalMatrix2D* f,
           int s_x, int s_y, 
           int e_x,int e_y, 
           int c_idx,
           FP cp, 
           FP ms, 
           FP t, 
           FP t_f,
           int StartIter=0);
 
 ~Source2D();
 
 void   SetSource2D( int i = 0 );
 void   ClearSource2D();
 
 int    GetCompIndex2D() { return  c_index; }

 int    GetSX()          { return  sx; }
 int    GetEX()          { return  ex; }
 int    GetSY()          { return  sy; }
 int    GetEY()          { return  ey; }
 int    GetStartIter()   { return StartSrcIter; }
 
 void   SetSX(int s_x)   {sx=s_x;}
 void   SetSY(int s_y)   {sy=s_y;}
 
 void   SetEX(int e_x)   {ex=e_x;}
 void   SetEY(int e_y)   {ey=e_y;}
 
 void   SetStartIter(int s) { StartSrcIter = s;}

 FP GetCp()          { return  Cp;  } 
 FP GetMs()          { return  M_s0;} 
 FP GetT()           { return  T;   } 
 FP GetTf()          { return  T_f; } 

 void   SetCp(FP cp) {Cp=cp;} 
 void   SetMs(FP ms) {M_s0=ms;} 
 void   SetT(FP t)   {T=t;}
 void   SetTf(FP t_f){T_f=t_f;} 

};

class SourceList2D : public UArray<Source2D*> {
    ComputationalMatrix2D* F;                //   Reference to computational matrix
    InputData*      data;                    //   Reference to input data
public:
    SourceList2D(ComputationalMatrix2D*, InputData*);
    ~SourceList2D() {};

void  SetSources2D(int iter = 0);
void  ClearSources2D();

};

// <------------- 2D --------------->           
#endif //_hyper_flow_source_hpp_


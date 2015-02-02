/*******************************************************************************
*   OpenHyperFLOW2D                                                            *
*                                                                              *
*   Version  1.0.2                                                             *
*   Copyright (C)  1995-2015 by Serge A. Suchkov                               *
*   Copyright policy: LGPL V3                                                  *
*                                                                              *
*   last update: 01/02/2015                                                    *
*******************************************************************************/
#include "libOpenHyperFLOW2D/hyper_flow_source.hpp"
// <------------- 2D --------------->           
Source2D::Source2D(ComputationalMatrix2D* f,
                   int s_x, int s_y, 
                   int e_x,int e_y, 
                   int c_idx,
                   FP cp, 
                   FP ms, 
                   FP t, 
                   FP t_f,
                   int si):F(f),sx(s_x),sy(s_y),ex(e_x),ey(e_y),c_index(c_idx),Cp(cp),M_s0(ms),T(t),T_f(t_f),StartSrcIter(si) {}

Source2D::~Source2D() {
    ClearSource2D();
}

void   Source2D::SetSource2D(int start_iter) {

if(start_iter < StartSrcIter) return;

int DX = sx-ex;
int DY = sy-ey;
int SKX, SKY;
FP dF, DR, DR2;
int i;
unsigned int x,y;


 // Source in single point
    if ( DX == 0 && DY == 0 ) {
        
        if(FlowNode2D<FP,NUM_COMPONENTS>::FT==FT_AXISYMMETRIC) {
            if(sy == 0 || ey ==0) {
              F->GetValue(sx,sy).Src[i2d_Ro]      =  M_s0/(M_PI*F->GetValue(sx,sy).dx*F->GetValue(sx,sy).dy*F->GetValue(sx,sy).dy);
            } else {
              F->GetValue(sx,sy).Src[i2d_Ro]      =  M_s0/(2*M_PI*F->GetValue(sx,sy).dx*F->GetValue(sx,sy).dy*F->GetValue(sx,sy).y);
            }
        } else {
              F->GetValue(sx,sy).Src[i2d_Ro]      =  M_s0/(F->GetValue(sx,sy).dx*F->GetValue(sx,sy).dy); 
        }
           
        F->GetValue(sx,sy).SrcAdd[i2d_Ro] =  0.;
        F->GetValue(sx,sy).Src[i2d_RoU]   =  0;
        F->GetValue(sx,sy).Tf             =  T_f;
        if(c_index < 4)
          F->GetValue(sx,sy).Src[c_index+4] =  F->GetValue(sx,sy).Src[i2d_Ro];
        F->GetValue(sx,sy).Src[i2d_RoE]     =  Cp*T*F->GetValue(sx,sy).Src[i2d_Ro];
        return;
    }

    
// Source on line
    if ( fabs(DX) > fabs(DY) ) {
        if ( DX>0 )SKX=1;
        else    SKX=-1;
        if ( DY>0 )SKY=1;
        else    SKY=-1;
        
        dF = fabs(DY)/fabs(DX);
        
        for ( i=0;i!=DX+SKX;i+=SKX ) {
            x = (unsigned int)(sx+i*SKX);
            y = (unsigned int)(sy+i*dF*SKY);
            if(FlowNode2D<FP,NUM_COMPONENTS>::FT==FT_AXISYMMETRIC) {
             if(sy == 0 || ey ==0) {
                DR = DY*F->GetValue(sx,sy).dy;
                F->GetValue(x,y).Src[i2d_Ro]      =  M_s0/(M_PI*(F->GetValue(sx,sy).dx*DR*DR));
             } else {
                DR2 = M_PI*fabs(sy*sy*F->GetValue(sx,sy).dy*F->GetValue(sx,sy).dy-
                           ey*ey*F->GetValue(sx,sy).dy*F->GetValue(sx,sy).dy);
                F->GetValue(x,y).Src[i2d_Ro]      =  M_s0/(F->GetValue(sx,sy).dx*DR2);
             }
                
            } else { 
                F->GetValue(x,y).Src[i2d_Ro]      =  M_s0/(F->GetValue(sx,sy).dx*F->GetValue(sx,sy).dy);
            }
               
            F->GetValue(x,y).SrcAdd[i2d_Ro] =  0.;
            F->GetValue(x,y).Tf             =  T_f;
            F->GetValue(x,y).Src[i2d_RoU]   =  0;
            F->GetValue(x,y).Src[i2d_RoV]   =  0;
            F->GetValue(x,y).Src[c_index+4] =  F->GetValue(x,y).Src[i2d_Ro];
            F->GetValue(x,y).Src[i2d_RoE]   =  Cp*T*F->GetValue(x,y).Src[i2d_Ro];
        }
    } else {
        if ( DY>0 )SKY=1;
        else    SKY=-1;
        if ( DX>0 )SKX=1;
        else    SKX=-1;
        
        dF = fabs(DX)/fabs(DY);
        
        for ( i=0;i!=DY+SKY;i+=SKY ) {
            x = (unsigned int)(sx+i*dF*SKX);
            y = (unsigned int)(sy+i*SKY);
            if(FlowNode2D<FP,NUM_COMPONENTS>::FT==FT_AXISYMMETRIC) {
                if(sy == 0 || ey ==0) {
                   DR = DY*F->GetValue(sx,sy).dy;
                   F->GetValue(x,y).Src[i2d_Ro]      =  M_s0/(M_PI*(F->GetValue(sx,sy).dx*DR*DR));
                } else   {
                    DR2 = M_PI*fabs(sy*sy*F->GetValue(sx,sy).dy*F->GetValue(sx,sy).dy-
                                    ey*ey*F->GetValue(sx,sy).dy*F->GetValue(sx,sy).dy);
                    F->GetValue(x,y).Src[i2d_Ro]  =  M_s0/(F->GetValue(sx,sy).dx*DR2);
                }
              }

            F->GetValue(x,y).SrcAdd[i2d_Ro] =  0.;
            F->GetValue(x,y).Tf             =  T_f;
            F->GetValue(x,y).Src[i2d_RoU]   =  0;
            F->GetValue(x,y).Src[i2d_RoV]   =  0;
            F->GetValue(x,y).Src[c_index+4] =  F->GetValue(x,y).Src[i2d_Ro];
            F->GetValue(x,y).Src[i2d_RoE]   =  Cp*T*F->GetValue(x,y).Src[i2d_Ro];
        }
    }
    return;
}


void   Source2D::ClearSource2D() {
int DX = sx-ex;
int DY = sy-ey;
int SKX, SKY;
int i;
unsigned int x,y;
FP dF;

 // Source in single point
    if ( DX == 0 && DY == 0 ) {
        
        for (int k=0;k<FlowNode2D<FP, NUM_COMPONENTS>::NumEq;k++ )
            F->GetValue(sx,sy).Src[k] = F->GetValue(sx,sy).SrcAdd[k] = 0.;
        return;
    }

    
// Source on line
    if ( fabs(DX) > fabs(DY) ) {
        if ( DX>0 )SKX=1;
        else    SKX=-1;
        if ( DY>0 )SKY=1;
        else    SKY=-1;
        
        dF = fabs(DY)/fabs(DX);
        
        for ( i=0;i!=DX+SKX;i+=SKX ) {
            x = (unsigned int)(sx+i*SKX);
            y = (unsigned int)(sy+i*dF*SKY);
            for (int k=0;k<FlowNode2D<FP, NUM_COMPONENTS>::NumEq;k++ )
                F->GetValue(x,y).Src[k] = F->GetValue(x,y).SrcAdd[k] = 0.;
        }
    } else {
        if ( DY>0 )SKY=1;
        else    SKY=-1;
        if ( DX>0 )SKX=1;
        else    SKX=-1;
        
        dF = fabs(DX)/fabs(DY);
        
        for ( i=0;i!=DY+SKY;i+=SKY ) {
            x = (unsigned int)(sx+i*dF*SKX);
            y = (unsigned int)(sy+i*SKY);
            for (int k=0;k<FlowNode2D<FP, NUM_COMPONENTS>::NumEq;k++ )
                F->GetValue(x,y).Src[k] = F->GetValue(x,y).SrcAdd[k] = 0.;
        }
    }
    return;
}

SourceList2D::SourceList2D(ComputationalMatrix2D* f, InputData* d) {

    int     NumSrc;
    int     GasSource_SX;
    int     GasSource_SY;  
    int     GasSource_EX;  
    int     GasSource_EY;  
    int     GasSourceIndex;
    int     StartSrcIter;
    
    FP  Msrc;          
    FP  Tsrc;          
    FP  Tf_src;
    FP  Cp=0.;
    FP  Y_mix[4];

    char    Str[255];
    
    Source2D* TmpSrc;
    
    Table*  Cp_cp;
    Table*  Cp_Ox;
    Table*  Cp_Fuel;
    Table*  Cp_air;
    
    data =  d;
    F    =  f;

    if(data) {
        Cp_cp   = data->GetTable((char*)"Cp_cp");
        Cp_Ox   = data->GetTable((char*)"Cp_OX");
        Cp_Fuel = data->GetTable((char*)"Cp_Fuel");
        Cp_air  = data->GetTable((char*)"Cp_air");

   NumSrc  = data->GetIntVal((char*)"NumSrc");
    
    for (int i=0; i< NumSrc;i++ ) {
        snprintf(Str,255,"Src%i.GasSrcSX",i+1);
        GasSource_SX    = data->GetIntVal(Str);
        snprintf(Str,255,"Src%i.GasSrcSY",i+1);
        GasSource_SY    = data->GetIntVal(Str);
        snprintf(Str,255,"Src%i.GasSrcEX",i+1);
        GasSource_EX    = data->GetIntVal(Str);
        snprintf(Str,255,"Src%i.GasSrcEY",i+1);
        GasSource_EY    = data->GetIntVal(Str);
        snprintf(Str,255,"Src%i.GasSrcIndex",i+1);
        GasSourceIndex  = data->GetIntVal(Str);
        snprintf(Str,255,"Src%i.Msrc",i+1);
        Msrc           =  data->GetFloatVal(Str);
        snprintf(Str,255,"Src%i.Tsrc",i+1);
        Tsrc          =   data->GetFloatVal(Str);
        snprintf(Str,255,"Src%i.Tf_src",i+1);
        Tf_src        =   data->GetFloatVal(Str);
        snprintf(Str,255,"Src%i.StartIter",i+1);
        StartSrcIter  =   data->GetFloatVal(Str);
        
        if ( GasSourceIndex==0 )        // Fuel
            Cp = Cp_Fuel->GetVal(Tsrc); 
        else if ( GasSourceIndex==1 )   // Ox
            Cp = Cp_Ox->GetVal(Tsrc);   
        else if ( GasSourceIndex==2 )   // Combustion products 
            Cp = Cp_cp->GetVal(Tsrc);   
        else if ( GasSourceIndex==3 )   // Air
            Cp = Cp_air->GetVal(Tsrc);  
        else if ( GasSourceIndex==4 ) { // Mixture
            snprintf(Str,255,"Src%i.Y_fuel",i+1);
            Y_mix[0] =  data->GetFloatVal(Str);
            snprintf(Str,255,"Src%i.Y_ox",i+1);
            Y_mix[1] =  data->GetFloatVal(Str);
            snprintf(Str,255,"Src%i.Y_cp",i+1);
            Y_mix[2] =  data->GetFloatVal(Str);
            snprintf(Str,255,"Src%i.Y_air",i+1);
            Y_mix[3] = 1 - Y_mix[0] + Y_mix[1] + Y_mix[2];
            Cp = Y_mix[0]*Cp_Fuel->GetVal(Tsrc) + 
                 Y_mix[1]*Cp_Ox->GetVal(Tsrc) + 
                 Y_mix[2]*Cp_cp->GetVal(Tsrc)+
                 Y_mix[3]*Cp_air->GetVal(Tsrc);
        }  if ( GasSourceIndex > 4 ) {
               *data->GetMessageStream() << "\nERROR: ";
               *data->GetMessageStream() << "Error component index" << GasSourceIndex << " in point gas source.\n" << flush;
               data->GetMessageStream()->flush();
       }
       
       TmpSrc = new Source2D(F,GasSource_SX,GasSource_SY,GasSource_EX,GasSource_EY,GasSourceIndex,Cp,Msrc,Tsrc,Tf_src,StartSrcIter);
       AddElement(&TmpSrc);
     }
   }
}

void SourceList2D::SetSources2D(int ii) {
    for(unsigned int i=0;i<GetNumElements();i++) {
        GetElement(i)->SetSource2D(ii);
    }
  return;
}

void SourceList2D::ClearSources2D() {
    for(unsigned int i=0;i<GetNumElements();i++) {
        GetElement(i)->ClearSource2D();
    }
  return;
}
// <------------- 2D --------------->           


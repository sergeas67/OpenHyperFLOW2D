/*******************************************************************************
*   OpenHyperFLOW2D-CUDA                                                       *
*                                                                              *
*   Transient, Density based Effective Explicit Parallel Hybrid Solver         *
*   TDEEPHS (CUDA+MPI)                                                         *
*   Version  1.0.2                                                             *
*   Copyright (C)  1995-2015 by Serge A. Suchkov                               *
*   Copyright policy: LGPL V3                                                  *
*                                                                              *
*   deeps2d_core.cpp: OpenHyperFLOW2D solver core code....                     *
*                                                                              *
*  last update: 11/01/2015                                                     *
********************************************************************************/
#include "deeps2d_core.hpp"

#include <sys/time.h>
#include <sys/timeb.h>
#include <sys/file.h>

int NumContour;
int start_iter = 5;
int n_s;
int num_threads = 1;
int num_blocks  = 1;

SourceList2D*  SrcList = NULL;
int            isGasSource=0;
int            isRecalcYplus;
int            isHighOrder;
int            TurbStartIter;
int            TurbExtModel;
int            err_i, err_j;
int            turb_mod_name_index = 0;
FP             Ts0,A,W,Mach;

unsigned int*  dt_min_host;
unsigned int*  dt_min_device;

UArray< unsigned int* >* dt_min_host_Array;
UArray< unsigned int* >* dt_min_device_Array;

UArray< XY<int> >* GlobalSubmatrix;
UArray< XY<int> >* WallNodes;

UArray<UMatrix2D< FlowNode2D<FP,NUM_COMPONENTS> >*>*     SubmatrixArray;
UArray<UMatrix2D< FlowNodeCore2D<FP,NUM_COMPONENTS> >*>* CoreSubmatrixArray;
UArray<XCut>*                                            XCutArray;

FP  makeZero;

int     Nstep;
FP  ExitMonitorValue;
int     MonitorNumber;
int     MonitorCondition; // 0 - equal
                          // 1 - less than
                          // 2 - great than

unsigned int     AddSrcStartIter = 0;
FP  beta[6+NUM_COMPONENTS];
FP  beta0;
FP  CFL;
Table*  CFL_Scenario;
Table*  beta_Scenario;
#ifdef _OPEN_MP
FP  DD_max_var;
#endif // OPEN_MP

//------------------------------------------
// Cx,Cy,Cd,Cv,Cp,St
//------------------------------------------
int     is_Cx_calc,is_Cd_calc;
FP  p_ambient;
FP  x0_body,y0_body,dx_body,dy_body;
FP  x0_nozzle,y0_nozzle,dy_nozzle;

int     Cx_Flow_index;
int     Cd_Flow_index;
int     Cp_Flow_index;
int     SigmaFi_Flow_index;
int     y_max,y_min;
//------------------------------------------
int     NOutStep;
int     isFirstStart=0;
int     ScaleFLAG=1;
int     isAdiabaticWall;
int     isOutHeatFluxX;
int     isOutHeatFluxY;
int     is_p_asterisk_out;

ChemicalReactionsModelData2D chemical_reactions;
BlendingFactorFunction       bFF;

FP* Y=NULL;
FP  Cp=0.;
FP  lam=0.;
FP  mu=0.;
FP  Tg=0.;
FP  Rg=0.;
FP  Pg=0.;
FP  Wg=0.;
FP  Ug,Vg;
int     CompIndex;

FP Y_fuel[4]={1.,0.,0.,0.};  /* fuel */
FP Y_ox[4]  ={0.,1.,0.,0.};  /* OX */
FP Y_cp[4]  ={0.,0.,1.,0.};  /* cp */
FP Y_air[4] ={0.,0.,0.,1.};  /* air */
FP Y_mix[4] ={0.,0.,0.,0.};  /* mixture */
#ifdef _RMS_
const  char*  RMS_Name[11] = {"Rho","Rho*U","Rho*V","Rho*E","Rho*Yfu","Rho*Yox","Rho*Ycp","Rho*k","Rho*eps","Rho*omega","nu_t"};
char          RMSFileName[255];
ofstream*     pRMS_OutFile;      // Output RMS stream
#endif //_RMS_

char          MonitorsFileName[255];
ofstream*     pMonitors_OutFile; // Output Monitors stream

int    useSwapFile=0;
char   OldSwapFileName[255];
void*  OldSwapData;
u_long OldFileSizeGas;
int    Old_fd;
int    isScan;
int    isPiston;

int MemLen;
int isRun=0;
int isDraw=0;
int isSnapshot=0;

char*                                        ProjectName;
char                                         GasSwapFileName[255];
char                                         ErrFileName[255];
char                                         OutFileName[255];
char                                         TecPlotFileName[255];
int                                          fd_g;
void*                                        GasSwapData;
int                                          TurbMod    = 0;
int                                          EndFLAG    = 1;
int                                          PrintFLAG;
ofstream*                                    pInputData;        // Output data stream
ofstream*                                    pHeatFlux_OutFile; // Output HeatFlux stream
unsigned int                                 iter = 0;          // iteration number
unsigned int                                 last_iter=0;       // Global iteration number
int                                          isStop=0;          // Stop flag
InputData*                                   Data=NULL;         // Object data loader
UMatrix2D< FlowNode2D<FP,NUM_COMPONENTS> >*  J=NULL;        // Main computation area
UMatrix2D< FlowNodeCore2D<FP,NUM_COMPONENTS> >*  C=NULL;    // Main computation area
UArray<Flow*>*                               FlowList=NULL;     // List of 'Flow' objects
UArray<Flow2D*>*                             Flow2DList=NULL;   // List of 'Flow2D' objects
UArray<Bound2D*>                             SingleBoundsList;  // Single Bounds List;

FP                                           dt;                // time step
FP*                                          RoUx=NULL;
FP*                                          RoVy=NULL;
FP                                           GlobalTime=0.;
FP                                           CurrentTimePart=0;

ofstream*                                    pOutputData;     // output data stream (file)

int                                          I,NSaveStep;
unsigned int                                 MaxX=0;          // X dimension of computation area
unsigned int                                 MaxY=0;          // Y dimension of computation area

FP                                           dxdy,dx2,dy2;
FP                                           SigW,SigF,delta_bl;

unsigned long FileSizeGas      = 0;
int      isVerboseOutput       = 0;
int      isTurbulenceReset     = 0;

void ChkCudaState(cudaError_t cudaState, char* name) {
    if(cudaState != cudaSuccess) {
     printf("\nError in %s\n",name);
     if(cudaState == cudaErrorMemoryAllocation) {
        printf("Memory allocation error.\n");
     } else if(cudaState == cudaErrorLaunchTimeout ) {
       printf("Timeout.\n");
     }else if(cudaState == cudaErrorLaunchOutOfResources) {
       printf("Resources temporary insufficient.\n");
     }else if(cudaState ==  cudaErrorInvalidConfiguration ) {
       printf("Resources insufficient for this device\n");
     }else if(cudaState == cudaErrorInvalidValue) {
       printf("Invalid value.\n");
     }else if(cudaState == cudaErrorInvalidHostPointer ) {
       printf("Invalid host pointer.\n");
     }else if(cudaState == cudaErrorInvalidDevicePointer) {
       printf("Invalid device pointer.\n");
     }else if(cudaState == cudaErrorNotReady) {
       printf("Device not ready.\n");
     } else {
      printf("Unknown error.\n");
     }
       Exit_OpenHyperFLOW2D(1);
    }
}

int CalibrateThreadBlockSize(int  cur_block_size,
                             int* opt_block_size,
                             FP*  opt_round_trip,
                             FP   round_trip){
   if(opt_block_size[1] != 1) {
       if(opt_round_trip[0] == 0.0) {
         opt_block_size[0] = opt_block_size[1] =  cur_block_size;
      } else if (opt_round_trip[0] > round_trip) {
        opt_block_size[0] = cur_block_size;
        opt_block_size[1] = cur_block_size/2;
      } else {
        opt_block_size[1] = cur_block_size/2; 
      }
   }

   opt_round_trip[0] = round_trip;
   
   if(opt_block_size[1] == 1)
     return opt_block_size[0];
   else
     return max(1,opt_block_size[1]);
};

inline void  DEEPS2D_Stage1(UMatrix2D< FlowNode2D<FP,NUM_COMPONENTS> >*     pLJ,
                     UMatrix2D< FlowNodeCore2D<FP,NUM_COMPONENTS> >* pLC,
                     int MIN_X, int MAX_X, int MAX_Y,
                     FP dxx, FP dyy,
                     FP dtdx, FP dtdy) {

    for (int i = MIN_X;i<MAX_X;i++ ) {
       for (int j = 0;j<MAX_Y;j++ ) {

          FlowNode2D< FP,NUM_COMPONENTS >* CurrentNode=NULL;

          CurrentNode = &(pLJ->GetValue(i,j)); 

          if ( CurrentNode->isCond2D(CT_NODE_IS_SET_2D) &&
                 !CurrentNode->isCond2D(CT_SOLID_2D)
              && !CurrentNode->isCond2D(NT_FC_2D)
              ) {

              FlowNodeCore2D< FP,NUM_COMPONENTS >* NextNode=NULL;

              FlowNode2D< FP,NUM_COMPONENTS >* UpNode=NULL;          // near
              FlowNode2D< FP,NUM_COMPONENTS >* DownNode=NULL;        // nodes
              FlowNode2D< FP,NUM_COMPONENTS >* LeftNode=NULL;        // references
              FlowNode2D< FP,NUM_COMPONENTS >* RightNode=NULL;

              FP beta[NUM_COMPONENTS+6];  
              FP _beta[NUM_COMPONENTS+6]; 
              FP dXX,dYY;

              int Num_Eq = FlowNode2D<FP,NUM_COMPONENTS>::NumEq-SetTurbulenceModel(CurrentNode);

              NextNode    = &(pLC->GetValue(i,j)); 

              int  n1=CurrentNode->idXl;
              int  n2=CurrentNode->idXr;
              int  n3=CurrentNode->idYu;
              int  n4=CurrentNode->idYd;

              int  N1 = i - n1;
              int  N2 = i + n2;
              int  N3 = j + n3;
              int  N4 = j - n4;

              int  n_n = max(n1+n2,1);
              int  m_m = max(n3+n4,1);

              FP  n_n_1 = 1./n_n;
              FP  m_m_1 = 1./m_m;

              UpNode    = &(pLJ->GetValue(i,N3));
              DownNode  = &(pLJ->GetValue(i,N4));
              RightNode = &(pLJ->GetValue(N2,j));
              LeftNode  = &(pLJ->GetValue(N1,j));

              // Scan equation system ... k - number of equation
              for (int k=0;k<Num_Eq;k++ ) {
                    int      c_flag = 0;
                    int      dx_flag, dx2_flag;
                    int      dy_flag, dy2_flag;

                    beta[k]  = CurrentNode->beta[k];
                    _beta[k] = 1. - beta[k];

                // Precomputed variables for current node ...
                    c_flag  = dx_flag = dy_flag = dx2_flag = dy2_flag = 0;
                    if ( k < 4 ) { // Make bit flags for future test for current equation
                        c_flag   = CT_Ro_CONST_2D     << k; 
                        dx_flag  = CT_dRodx_NULL_2D   << k;
                        dy_flag  = CT_dRody_NULL_2D   << k;
                        dx2_flag = CT_d2Rodx2_NULL_2D << k;
                        dy2_flag = CT_d2Rody2_NULL_2D << k;
                    } else if (k < (4+NUM_COMPONENTS)) {
                        c_flag   = CT_Y_CONST_2D;
                        dx_flag  = CT_dYdx_NULL_2D;
                        dy_flag  = CT_dYdy_NULL_2D;
                        dx2_flag = CT_d2Ydx2_NULL_2D;
                        dy2_flag = CT_d2Ydy2_NULL_2D;
                    } else if ((CurrentNode->isTurbulenceCond2D(TCT_k_eps_Model_2D) ||
                                CurrentNode->isTurbulenceCond2D(TCT_Spalart_Allmaras_Model_2D))) { //
                      if( k == i2d_k) {
                          c_flag   = TCT_k_CONST_2D     << (k-4-NUM_COMPONENTS); 
                          dx_flag  = TCT_dkdx_NULL_2D   << (k-4-NUM_COMPONENTS);
                          dy_flag  = TCT_dkdy_NULL_2D   << (k-4-NUM_COMPONENTS);
                          dx2_flag = TCT_d2kdx2_NULL_2D << (k-4-NUM_COMPONENTS);
                          dy2_flag = TCT_d2kdy2_NULL_2D << (k-4-NUM_COMPONENTS);
                      } else if (k == i2d_eps) {
                          c_flag   = TCT_eps_CONST_2D     << (k-4-NUM_COMPONENTS); 
                          dx_flag  = TCT_depsdx_NULL_2D   << (k-4-NUM_COMPONENTS);
                          dy_flag  = TCT_depsdy_NULL_2D   << (k-4-NUM_COMPONENTS);
                          dx2_flag = TCT_d2epsdx2_NULL_2D << (k-4-NUM_COMPONENTS);
                          dy2_flag = TCT_d2epsdy2_NULL_2D << (k-4-NUM_COMPONENTS);
                      }
                    }
                    // Check BC for current equation
                    if (k<(4+NUM_COMPONENTS)) {

                        if ( CurrentNode->isCond2D((CondType2D)c_flag) )
                            c_flag  = 0;
                        else
                            c_flag  = 1;

                        if ( CurrentNode->isCond2D((CondType2D)dx_flag) ) {
                            dx_flag = 0;
                        } else {
                            dx_flag = 1;
                        }

                        if ( CurrentNode->isCond2D((CondType2D)dy_flag) ) {
                            dy_flag = 0;
                        } else {
                            dy_flag = 1;
                        }

                        if ( CurrentNode->isCond2D((CondType2D)dx2_flag) ) {
                            dx2_flag = 1;
                        } else {
                            dx2_flag = 0;
                        }

                        if ( CurrentNode->isCond2D((CondType2D)dy2_flag) ) {
                            dy2_flag = 1;
                        } else {
                            dy2_flag = 0;
                        }
                    } else if((CurrentNode->isTurbulenceCond2D(TCT_k_eps_Model_2D) ||
                               CurrentNode->isTurbulenceCond2D(TCT_Spalart_Allmaras_Model_2D)) ) {
                        if ( CurrentNode->isTurbulenceCond2D((TurbulenceCondType2D)c_flag) )
                            c_flag  = 0;
                        else
                            c_flag  = 1;

                        if ( CurrentNode->isTurbulenceCond2D((TurbulenceCondType2D)dx_flag) ) {
                            dx_flag = 0;
                        } else {
                            dx_flag = 1;
                        }
                        if ( CurrentNode->isTurbulenceCond2D((TurbulenceCondType2D)dy_flag) ) {
                            dy_flag = 0;
                        } else {
                            dy_flag = 1;
                        }
                        if ( CurrentNode->isTurbulenceCond2D((TurbulenceCondType2D)dx2_flag) ) {
                            dx2_flag = 1;
                        } else {
                            dx2_flag = 0;
                        }
                        if ( CurrentNode->isTurbulenceCond2D((TurbulenceCondType2D)dy2_flag) ) {
                            dy2_flag = 1;
                        } else {
                            dy2_flag = 0;
                        }
                    }
                    if ( c_flag ) {
                        if ( dx_flag ) {
                            dXX = CurrentNode->dSdx[k] = (RightNode->A[k]-LeftNode->A[k])*n_n_1;
                        } else {
                            CurrentNode->S[k] = (LeftNode->S[k]*n2+RightNode->S[k]*n1)*n_n_1;
                            dXX = CurrentNode->dSdx[k] = 0.;
                        }
                        if ( dy_flag ) {
                            dYY = CurrentNode->dSdy[k] = (UpNode->B[k]-DownNode->B[k])*m_m_1;

                        } else {
                            CurrentNode->S[k] =  (UpNode->S[k]*n3+DownNode->S[k]*n4)*m_m_1;
                            dYY = CurrentNode->dSdy[k] = 0;
                        }
                        if ( dx2_flag ) {
                            dXX = (LeftNode->dSdx[k]+RightNode->dSdx[k])*0.5;
                        }
                        if ( dy2_flag ) {
                            dYY = (UpNode->dSdy[k]+DownNode->dSdy[k])*0.5;
                        }

                        if ( CurrentNode->FT ) {
                            NextNode->S[k] = CurrentNode->S[k]*beta[k]+_beta[k]*(dxx*(LeftNode->S[k]+RightNode->S[k])+dyy*(UpNode->S[k]+DownNode->S[k]))*0.5
                                          - (dtdx*dXX+dtdy*(dYY+CurrentNode->F[k]/(j+1))) + (CurrentNode->Src[k])*dt+CurrentNode->SrcAdd[k];
                        } else {
                            NextNode->S[k] = CurrentNode->S[k]*beta[k]+_beta[k]*(dxx*(LeftNode->S[k]+RightNode->S[k])+dyy*(UpNode->S[k]+DownNode->S[k]))*0.5
                                          - (dtdx*dXX+dtdy*dYY) + (CurrentNode->Src[k])*dt+CurrentNode->SrcAdd[k];
                        }
                }
            }
        }
     }
 }
}


inline int ConvertTurbMod(int input_tm) {
int      AddEq = 2;
       if (input_tm == 2) {          /*(TCT_Prandtl_Model_2D)*/
           AddEq = 2;
       } else if ( input_tm == 4) {  /*(TCT_k_eps_Model_2D)*/
           AddEq = 0;
       } else if ( input_tm == 3) {  /*(TCT_Spalart_Allmaras_Model_2D)*/
           AddEq = 1;
       } else {
           AddEq = 2;
       }
 return AddEq;

}
#ifdef _MPI_
inline FP DEEPS2D_Stage2(UMatrix2D< FlowNode2D<FP,NUM_COMPONENTS> >*     pLJ,
                             UMatrix2D< FlowNodeCore2D<FP,NUM_COMPONENTS> >* pLC,
                            int MIN_X, int MAX_X, int MAX_Y,
                            FP  beta0, int b_FF, ChemicalReactionsModelData2D* pCRMD,
                            int iter, int last_iter, int TurbStartIter,
                            FP SigW, FP SigF, FP dx_1, FP dy_1, FP delta_bl,
                            FP CFL0,FP beta_init,
#ifdef _RMS_
#ifdef _MPI
                            Var_pack* DD_max, int rank,
#else
                            UMatrix2D<FP>& RMS, 
                            UMatrix2D<int>&    iRMS,
                            UMatrix2D<FP>& DD_max,
                            int* i_c, int* j_c,
                            int ii,
#endif //_MPI 
#endif // _RMS_
                            TurbulenceExtendedModel TurbExtModel ) {

    FP dt_min_local;
    FP beta_min;

    beta_min = min(beta0,beta_init);

    for (int i = MIN_X;i<MAX_X;i++ ) {
       for (int j = 0;j<MAX_Y;j++ ) {

           FlowNode2D< FP,NUM_COMPONENTS >* CurrentNode=NULL;

           CurrentNode = &(pLJ->GetValue(i,j)); 

           if (CurrentNode->isCond2D(CT_NODE_IS_SET_2D) &&
               !CurrentNode->isCond2D(CT_SOLID_2D) &&
               !CurrentNode->isCond2D(NT_FC_2D)) {

               FlowNodeCore2D< FP,NUM_COMPONENTS >* NextNode=NULL;

               FlowNode2D< FP,NUM_COMPONENTS >* UpNode=NULL;          // near
               FlowNode2D< FP,NUM_COMPONENTS >* DownNode=NULL;        // nodes
               FlowNode2D< FP,NUM_COMPONENTS >* LeftNode=NULL;        // references
               FlowNode2D< FP,NUM_COMPONENTS >* RightNode=NULL;
               
               FP DD_local[NUM_COMPONENTS+6];

               int Num_Eq = FlowNode2D<FP,NUM_COMPONENTS>::NumEq-SetTurbulenceModel(CurrentNode);

               NextNode    = &(pLC->GetValue(i,j)); 

               int  n1=CurrentNode->idXl;
               int  n2=CurrentNode->idXr;
               int  n3=CurrentNode->idYu;
               int  n4=CurrentNode->idYd;

               int  N1 = i - n1;
               int  N2 = i + n2;
               int  N3 = j + n3;
               int  N4 = j - n4;

               int  n_n = max(n1+n2,1);
               int  m_m = max(n3+n4,1);

               FP  n_n_1 = 1./n_n;
               FP  m_m_1 = 1./m_m;

               UpNode    = &(pLJ->GetValue(i,N3));
               DownNode  = &(pLJ->GetValue(i,N4));
               RightNode = &(pLJ->GetValue(N2,j));
               LeftNode  = &(pLJ->GetValue(N1,j));

               // Scan equation system ... k - number of equation
               for (int k=0;k<Num_Eq;k++ ) {

                   int c_flag = 0;

                   if ( k < 4 ) // Make bit flags for future test for current equation // FlowNode2D<FP,NUM_COMPONENTS>::NumEq-AddEq-NUM_COMPONENTS ?
                       c_flag  = CT_Ro_CONST_2D   << k;
                   else if (k<(4+NUM_COMPONENTS))  // 7 ?
                       c_flag  = CT_Y_CONST_2D;
                   else if((CurrentNode->isTurbulenceCond2D(TCT_k_eps_Model_2D) ||
                            CurrentNode->isTurbulenceCond2D(TCT_Spalart_Allmaras_Model_2D) )) 
                       c_flag  = TCT_k_CONST_2D << (k-4-NUM_COMPONENTS); 

                   DD_local[k] = 0;

                   if ( !CurrentNode->isCond2D((CondType2D)c_flag) && 
                         CurrentNode->S[k] != 0. ) {

                         FP Tmp;

                         if(k == i2d_RoU && k == i2d_RoV ) {
                             Tmp = sqrt(CurrentNode->S[i2d_RoU]*CurrentNode->S[i2d_RoU]+
                                        CurrentNode->S[i2d_RoV]*CurrentNode->S[i2d_RoV]+1.e-30); // Flux
                         } else {
                             Tmp = CurrentNode->S[k];
                         }

                         if(fabs(Tmp) > 1.e-15)
                            DD_local[k] = fabs((NextNode->S[k]-CurrentNode->S[k])/Tmp);
                         else
                            DD_local[k] = 0.0;

                         if( b_FF == BFF_L) {
                         //LINEAR locally adopted blending factor function  (LLABFF)
                           CurrentNode->beta[k] = min(beta_min,(beta_min*beta_min)/(beta_min+DD_local[k]));
                         } else if( b_FF == BFF_LR) {
                         //LINEAR locally adopted blending factor function with relaxation (LLABFFR)
                           CurrentNode->beta[k] = min((beta_min+CurrentNode->beta[k])*0.5,(beta_min*beta_min)/(beta_min+DD_local[k]));
                         } else if( b_FF == BFF_S) {
                         //SQUARE locally adopted blending factor function (SLABF)
                           CurrentNode->beta[k] = min(beta_min,(beta_min*beta_min)/(beta_min+DD_local[k]*DD_local[k]));
                         } else if (b_FF == BFF_SR) {
                         //SQUARE locally adopted blending factor function with relaxation (SLABFFR)
                           CurrentNode->beta[k] = min((beta_min+CurrentNode->beta[k])*0.5,(beta_min*beta_min)/(beta_min+DD_local[k]*DD_local[k]));
                         } else if( b_FF == BFF_SQR) {
                         //SQRT() locally adopted blending factor function (SQRLABF)
                           CurrentNode->beta[k] = min(beta_min,(beta_min*beta_min)/(beta_min+sqrt(DD_local[k])));
                         } else if( b_FF == BFF_SQRR) {
                         //SQRT() locally adopted blending factor function with relaxation (SQRLABFFR)
                           CurrentNode->beta[k] = min((beta_min+CurrentNode->beta[k])*0.5,(beta_min*beta_min)/(beta_min+sqrt(DD_local[k]))); 
                         } else if( b_FF == BFF_LG ) {
                           FP LGAF = sqrt(CurrentNode->dSdx[k]*CurrentNode->dSdx[k] +
                                              CurrentNode->dSdy[k]*CurrentNode->dSdy[k] + 1.0e-30) * (dx+dy)/CurrentNode->S[k];
                           CurrentNode->beta[k] = min(beta_min,(beta_min*beta_min)/(beta_min+LGAF));
                         } else if (b_FF == BFF_MACH) {
                           //Mach number depended locally adapted blending factor function  (MLABFF)
                           CurrentNode->beta[k] = min(beta0,(beta_min*beta_min)/(beta_min+pow(DD_local[k],1.0/(1.0+Mach)))); // 1/(1+M) or 1/M ???
                         } else if (b_FF == BFF_HYBRID) {
                           // Hybrid BFF
                           FP LGAF = sqrt(CurrentNode->dSdx[k]*CurrentNode->dSdx[k] +
                                              CurrentNode->dSdy[k]*CurrentNode->dSdy[k] + 1.0e-30) * (dx+dy)/CurrentNode->S[k];
                           CurrentNode->beta[k] = min(beta_min,(beta_min*beta_min)/(beta_min+(LGAF+
                                                                               DD_local[k]*DD_local[k]+
                                                                               pow(DD_local[k],1.0/(1.0+Mach)))/3.));
                         } else if( b_FF == BFF_MIXED ) {
                           FP LGAF = sqrt(CurrentNode->dSdx[k]*CurrentNode->dSdx[k] +
                                              CurrentNode->dSdy[k]*CurrentNode->dSdy[k] + 1.0e-30) * (dx+dy)/CurrentNode->S[k];
                           CurrentNode->beta[k] = min(beta_min,(beta_min*beta_min)/(beta_min+(LGAF+DD_local[k]*DD_local[k])*0.5));
                         } else if (b_FF == BFF_SR_LIMITED) {
                           //LIMITED locally adopted blending factor function with relaxation (SLLABFFR)
                           CurrentNode->beta[k] = min((beta_min+CurrentNode->beta[k])*0.5,beta_min/(1.0+DD_local[k]));
                         } else {
                           // Default->SQRLABF
                           CurrentNode->beta[k] = min(beta_min,(beta_min*beta_min)/(beta_min+sqrt(DD_local[k])));
                         }
#ifdef _RMS_
                         RMS(k,ii) += DD_local[k];
                         iRMS(k,ii)++;
                         DD_max(k,ii) = max(DD_max(k,ii),DD_local[k]);

                         if ( DD_max(k,ii) == DD_local[k] ) {
                              i_c[ii] = i;
                              j_c[ii] = j;
                         }
#endif // RMS
                   }
                   if (k<(4+NUM_COMPONENTS)) {
                       if ( !CurrentNode->isCond2D((CondType2D)c_flag) )
                             CurrentNode->S[k]   = NextNode->S[k];
                   } else if ((CurrentNode->isTurbulenceCond2D(TCT_k_eps_Model_2D) ||
                               CurrentNode->isTurbulenceCond2D(TCT_Spalart_Allmaras_Model_2D)) ){
                       if ( !CurrentNode->isTurbulenceCond2D((TurbulenceCondType2D)c_flag) )
                             CurrentNode->S[k]   =  NextNode->S[k];
                   }
               }

               CurrentNode->droYdx[NUM_COMPONENTS]=CurrentNode->droYdy[NUM_COMPONENTS]=0.;

               for (int k=4;k<FlowNode2D<FP,NUM_COMPONENTS>::NumEq-2;k++ ) {
                   if ( !CurrentNode->isCond2D(CT_dYdx_NULL_2D) ) {
                       CurrentNode->droYdx[k-4]=(RightNode->S[k]-LeftNode->S[k])*dx_1*0.5;
                       CurrentNode->droYdx[NUM_COMPONENTS]+=(RightNode->S[k]-LeftNode->S[k])*dx_1*0.5;
                   }
                   if ( !CurrentNode->isCond2D(CT_dYdy_NULL_2D) ) {
                         CurrentNode->droYdy[k-4]=(UpNode->S[k]-DownNode->S[k])*dy_1*0.5;
                         CurrentNode->droYdy[NUM_COMPONENTS]+=(DownNode->S[k]-UpNode->S[k])*dy_1*0.5;
                   }
               }

               if (CurrentNode->isCond2D(CT_WALL_NO_SLIP_2D) || CurrentNode->isCond2D(CT_WALL_LAW_2D) )  {
                   CurrentNode->dUdx=(RightNode->U*n1-LeftNode->U*n2)*dx_1*n_n_1;
                   CurrentNode->dVdx=(RightNode->V*n1-LeftNode->V*n2)*dx_1*n_n_1;

                   CurrentNode->dUdy=(UpNode->U*n3-DownNode->U*n4)*dy_1*m_m_1;
                   CurrentNode->dVdy=(UpNode->V*n3-DownNode->V*n4)*dy_1*m_m_1;

                   if(CurrentNode->isTurbulenceCond2D(TCT_k_eps_Model_2D)){
                     CurrentNode->dkdx   =(RightNode->S[i2d_k]*n1-LeftNode->S[i2d_k]*n2)*dx_1/CurrentNode->S[i2d_Ro]*n_n_1;
                     CurrentNode->depsdx =(RightNode->S[i2d_eps]*n1-LeftNode->S[i2d_eps]*n2)*dx_1/CurrentNode->S[i2d_Ro]*n_n_1;

                     CurrentNode->dkdy   =(UpNode->S[i2d_k]*n3-DownNode->S[i2d_k]*n4)*dy_1/CurrentNode->S[i2d_Ro]*m_m_1;
                     CurrentNode->depsdy =(UpNode->S[i2d_eps]*n3-DownNode->S[i2d_eps]*n4)*dy_1/CurrentNode->S[i2d_Ro]*m_m_1;
                   } else if (CurrentNode->isTurbulenceCond2D(TCT_Spalart_Allmaras_Model_2D)) {
                              CurrentNode->dkdx   =(RightNode->S[i2d_k]*n1-LeftNode->S[i2d_k]*n2)*dx_1/CurrentNode->S[i2d_Ro]*n_n_1;
                              CurrentNode->dkdy   =(UpNode->S[i2d_k]*n3-DownNode->S[i2d_k]*n4)*dy_1/CurrentNode->S[i2d_Ro]*m_m_1;
                   }
               } else {
                   CurrentNode->dUdx   =(RightNode->U-LeftNode->U)*dx_1*n_n_1;
                   CurrentNode->dVdx   =(RightNode->V-LeftNode->V)*dx_1*n_n_1;

                   CurrentNode->dUdy   =(UpNode->U-DownNode->U)*dy_1*m_m_1;
                   CurrentNode->dVdy   =(UpNode->V-DownNode->V)*dy_1*m_m_1;
                   if(CurrentNode->isTurbulenceCond2D(TCT_k_eps_Model_2D)){
                     CurrentNode->dkdx   =(RightNode->S[i2d_k]-LeftNode->S[i2d_k])*dx_1/CurrentNode->S[i2d_Ro]*n_n_1;
                     CurrentNode->depsdx =(RightNode->S[i2d_eps]-LeftNode->S[i2d_eps])*dx_1/CurrentNode->S[i2d_Ro]*n_n_1;

                     CurrentNode->dkdy   =(UpNode->S[i2d_k]-DownNode->S[i2d_k])*dy_1/CurrentNode->S[i2d_Ro]*m_m_1;
                     CurrentNode->depsdy =(UpNode->S[i2d_eps]-DownNode->S[i2d_eps])*dy_1/CurrentNode->S[i2d_Ro]*m_m_1;
                   } else if (CurrentNode->isTurbulenceCond2D(TCT_Spalart_Allmaras_Model_2D)) {
                              CurrentNode->dkdx   =(RightNode->S[i2d_k]-LeftNode->S[i2d_k])*dx_1/CurrentNode->S[i2d_Ro]*n_n_1;
                              CurrentNode->dkdy   =(UpNode->S[i2d_k]-DownNode->S[i2d_k])*dy_1/CurrentNode->S[i2d_Ro]*m_m_1;
                   }
               }

               CurrentNode->dTdx=(RightNode->Tg-LeftNode->Tg)*dx_1*n_n_1;
               CurrentNode->dTdy=(UpNode->Tg-DownNode->Tg)*dy_1*m_m_1;

               CalcChemicalReactions(CurrentNode,CRM_ZELDOVICH, (void*)(pCRMD));

               if((int)(iter+last_iter) < TurbStartIter) {
                  CurrentNode->FillNode2D(0,1,SigW,SigF,(TurbulenceExtendedModel)TurbExtModel,delta_bl);
               } else {
                  CurrentNode->FillNode2D(1,0,SigW,SigF,(TurbulenceExtendedModel)TurbExtModel,delta_bl);
               }

               if( CurrentNode->Tg < 0. ) {
                   *f_stream << "\nTg=" << CurrentNode->Tg << " K. p=" << CurrentNode->p <<" Pa dt=" << dt << " sec.\n" << flush;
                   *f_stream << "\nERROR: Computational unstability in UMatrix2D< FlowNode2D<FP,NUM_COMPONENTS> >(" << i <<","<< j <<") on iteration " << iter+last_iter<< "...\n";
                   return 0.0;
               }  else {
                   FP AAA          = sqrt(CurrentNode->k*CurrentNode->R*CurrentNode->Tg); 
                   dt_min_local        = min(CFL,CFL0)*min(dx/(AAA+fabs(CurrentNode->U)),dy/(AAA+fabs(CurrentNode->V)));
               }
        }
      }
    }
 return dt_min_local;
}
#endif // _MPI_
// CFL_Scenario->GetVal(iter+last_iter)

#ifdef _CUDA_
void DEEPS2D_Run(ofstream* f_stream, 
                 UMatrix2D<FlowNode2D<FP,NUM_COMPONENTS> >*     pJ,
                 UMatrix2D<FlowNodeCore2D<FP,NUM_COMPONENTS> >* pC,
                 UArray< FlowNode2D<FP,NUM_COMPONENTS>* >*      cudaSubmatrixArray,
                 UArray< FlowNodeCore2D<FP,NUM_COMPONENTS>* >*  cudaCoreSubmatrixArray,
                 UArray< XY<int> >*                             cudaDimArray,
                 UArray< XY<int>* >*                            cudaWallNodesArray,
                 UArray< ChemicalReactionsModelData2D* >*       cudaCRM2D,
                 int num_mp,
                 cudaStream_t *cuda_streams,
                 cudaEvent_t  *cuda_events) {

   // local variables
    FlowNode2D<FP,NUM_COMPONENTS>* TmpMatrixPtr;
    //FlowNodeCore2D<FP,NUM_COMPONENTS>* TmpCoreMatrixPtr;
    //int    SubMaxX, SubStartIndex; 
    int      TmpMaxX;
    int      current_div;
    int      opt_thread_block_size[2];
    FP       opt_round_trip[1];

    int      num_cuda_threads;
    int      num_cuda_blocks; 

    FP       dtdx;
    FP       dtdy;
    FP       dyy;
    FP       dxx;
    FP       dx_1,dy_1; // 1/dx, 1/dy
    FP       d_time;
    FP       VCOMP;
    timeval  start, stop, mark1, mark2;
    FP       int2float_scale;
#ifdef _RMS_
    FP       max_RMS;
    int      k_max_RMS;
#endif //_RMS_
    FP*       dt_min;
    int*      i_c;
    int*      j_c;
    FP        dtmin=1.0;
    //FP      DD[FlowNode2D<FP,NUM_COMPONENTS>::NumEq];
#ifdef _RMS_
    FP                sum_RMS[FlowNode2D<FP,NUM_COMPONENTS>::NumEq];
    unsigned long     sum_iRMS[FlowNode2D<FP,NUM_COMPONENTS>::NumEq];
    UMatrix2D<FP>     RMS(FlowNode2D<FP,NUM_COMPONENTS>::NumEq,cudaArraySubmatrix->GetNumElements());
    UMatrix2D<int>    iRMS(FlowNode2D<FP,NUM_COMPONENTS>::NumEq,cudaArraySubmatrix->GetNumElements());
#endif //_RMS_
    UMatrix2D<FP>     DD_max(FlowNode2D<FP,NUM_COMPONENTS>::NumEq,cudaDimArray->GetNumElements());

    isScan = 0;
    dyy    = dx/(dx+dy);
    dxx    = dy/(dx+dy);
    dxdy   = dx*dy;
    dx2    = dx*dx;
    dy2    = dy*dy;
    d_time = 0.;
    
    dt_min = new FP[cudaArraySubmatrix->GetNumElements()];
    i_c    = new int[cudaArraySubmatrix->GetNumElements()];
    j_c    = new int[cudaArraySubmatrix->GetNumElements()];

    for(int ii=0;ii<(int)cudaArraySubmatrix->GetNumElements();ii++) {
        dt_min[ii] = dtmin = dt;
        i_c[ii] = j_c[ii] = 0;
     }
#ifdef _RMS_    
    snprintf(RMSFileName,255,"RMS-%s",OutFileName);
    CutFile(RMSFileName);
    pRMS_OutFile = OpenData(RMSFileName);
    SaveRMSHeader(pRMS_OutFile);
#endif // _RMS_
    if(MonitorPointsArray && MonitorPointsArray->GetNumElements() > 0) {
        snprintf(MonitorsFileName,255,"Monitors-%s",OutFileName);
        CutFile(MonitorsFileName);
        pMonitors_OutFile = OpenData(MonitorsFileName);
        SaveMonitorsHeader(pMonitors_OutFile, MonitorPointsArray);
    }

#ifdef _DEBUG_0
         ___try {
#endif  // _DEBUG_0
                    I = 0;
                    isRun = 1;
                    FlowNode2D<FP,NUM_COMPONENTS>*      cudaJ;
                    FlowNodeCore2D<FP,NUM_COMPONENTS>*  cudaC;
                    dtmin = dt;
                    current_div =  4*warp_size;
                    int2float_scale  = (FP)(INT_MAX)/(256*dt);

             do {

                  gettimeofday(&mark2,NULL);
                  gettimeofday(&start,NULL);

                  if( AddSrcStartIter < iter + last_iter) {
                   FlowNode2D<FP,NUM_COMPONENTS>::isSrcAdd = 1;
                  } else {
                   FlowNode2D<FP,NUM_COMPONENTS>::isSrcAdd = 0;
                  }

                  n_s = cudaDimArray->GetNumElements();

                  iter = 0;

                  do {

#ifdef _RMS_
                       max_RMS = 0.5 * ExitMonitorValue;
                       k_max_RMS = -1;

#endif // _RMS_
#ifdef _RMS_
                       for (int k=0;k<(int)FlowNode2D<FP,NUM_COMPONENTS>::NumEq;k++ ) {

                           for(int ii=0;ii<n_s;ii++) {
                               sum_RMS[k]  = RMS(k,ii)  = 0.;    // Clean sum residual
                               sum_iRMS[k] = iRMS(k,ii) = 0;     // num involved nodes
                           }
                           DD[k]      = 0.;
                       }
#endif // _RMS_

                       unsigned int dtest_int = (unsigned int)(int2float_scale*10);

                       int iX0=0;

                       dtmin = 10*dt;
#pragma unroll
                       for(int ii=0;ii<n_s;ii++) {  // CUDA version

                       int max_X = cudaDimArray->GetElement(ii).GetX();
                       int max_Y = cudaDimArray->GetElement(ii).GetY();

                       int r_Overlap;
                       int l_Overlap;

                       if(cudaSetDevice(ii) != cudaSuccess ) {
                          *f_stream << "\nError set CUDA device no: "<< ii << endl;
                          Exit_OpenHyperFLOW2D(n_s);
                       }

#ifdef _DEVICE_MMAP_
                       dt_min_host = dt_min_host_Array->GetElement(ii);
                       *dt_min_host = dtest_int;
#else
                       dt_min_device = dt_min_device_Array->GetElement(ii);
                       CopyHostToDevice(&dtest_int,dt_min_device,sizeof(unsigned int));
#endif //_DEVICE_MMAP_
                       cudaJ = cudaSubmatrixArray->GetElement(ii);
                       cudaC = cudaCoreSubmatrixArray->GetElement(ii);

                       num_cuda_threads =  4*warp_size/max(1,current_div);
                       num_cuda_blocks  = (max_X*max_Y)/num_cuda_threads;

                       if(num_cuda_blocks*num_cuda_threads != max_X*max_Y)
                          num_cuda_blocks++;

                       if(ii == n_s-1)
                         r_Overlap = 0;
                       else
                         r_Overlap = 1;

                       if(ii == 0)
                         l_Overlap = 0;
                       else
                         l_Overlap = 1;

                       for (int k=0;k<(int)FlowNode2D<FP,NUM_COMPONENTS>::NumEq;k++ ) {
                          DD_max(k,ii) =  0.;
                         }

                       dx_1 = 1.0/dx;
                       dy_1 = 1.0/dy;

                       // dt = ??? 

                       dtdx = dt/dx;
                       dtdy = dt/dy;

                       cuda_DEEPS2D_Stage1<<<num_cuda_blocks,num_cuda_threads, 0, cuda_streams[ii]>>>(cudaJ,cudaC,
                                                                                                      //cudaCRM2D->GetElement(ii),
                                                                                                      max_X*max_Y,
                                                                                                      max_X,max_Y,
                                                                                                      iX0 + max_X - r_Overlap,
                                                                                                      iX0 + l_Overlap,
                                                                                                      dxx,dyy,dtdx,dtdy,dt,
                                                                                                      FlowNode2D<FP,NUM_COMPONENTS>::FT,
                                                                                                      FlowNode2D<FP,NUM_COMPONENTS>::NumEq-ConvertTurbMod(TurbMod),
                                                                                                      ProblemType);
                       iX0 += max_X;
                      }
                   CUDA_BARRIER((char*)"cuda_DEEPS2D_Stage1");
                   iX0 = 0;
 #pragma unroll
                   for(int ii=0;ii<n_s;ii++) {  // CUDA version

                       if(cudaSetDevice(ii) != cudaSuccess ) {
                          *f_stream << "\nError set CUDA device no: "<< ii << endl;
                          Exit_OpenHyperFLOW2D(n_s);
                       }

                       int max_X = cudaDimArray->GetElement(ii).GetX();
                       int max_Y = cudaDimArray->GetElement(ii).GetY();

                       int r_Overlap;
                       int l_Overlap;

                       cudaJ = cudaSubmatrixArray->GetElement(ii);
                       cudaC = cudaCoreSubmatrixArray->GetElement(ii);

                       num_cuda_threads = min(max_num_threads,4*warp_size/max(1,current_div));
                       num_cuda_blocks  = (max_X*max_Y)/num_cuda_threads;

                       if(num_cuda_blocks*num_cuda_threads != max_X*max_Y)
                          num_cuda_blocks++;

                       if(ii == n_s-1)
                         r_Overlap = 0;
                       else
                         r_Overlap = 1;
                       if(ii == 0)
                         l_Overlap = 0;
                       else
                         l_Overlap = 1;

                        int noTurbCond;
                        if(iter+last_iter<TurbStartIter) 
                           noTurbCond = 1;
                        else 
                           noTurbCond = 0;

                        cudaHu =  cudaHuArray->GetElement(ii);
                        dt_min_device = dt_min_device_Array->GetElement(ii);

                        cuda_DEEPS2D_Stage2<<<num_cuda_blocks,num_cuda_threads, 0, cuda_streams[ii]>>>(cudaJ,cudaC,
                                                                                                       max_X*max_Y,
                                                                                                       max_X,max_Y,
                                                                                                       iX0 + max_X - r_Overlap,
                                                                                                       iX0 + l_Overlap,
                                                                                                       beta_Scenario->GetVal(iter+last_iter),beta0,
                                                                                                       bFF, CFL_Scenario->GetVal(iter+last_iter),
                                                                                                       cudaCRM2D->GetElement(ii),
                                                                                                       noTurbCond,
                                                                                                       SigW,SigF,dx_1,dy_1,delta_bl,
                                                                                                       FlowNode2D<FP,NUM_COMPONENTS>::FT,
                                                                                                       FlowNode2D<FP,NUM_COMPONENTS>::NumEq-ConvertTurbMod(TurbMod),
                                                                                                       cudaHu,
                                                                                                       FlowNode2D<FP,NUM_COMPONENTS>::isSrcAdd,
                                                                                                       dt_min_device, int2float_scale,
                                                                                                       (TurbulenceExtendedModel)TurbExtModel,
                                                                                                       ProblemType);
                        if (MonitorPointsArray &&
                            iter/NOutStep*NOutStep == iter ) {
                            for(int ii_monitor=0;ii_monitor<(int)MonitorPointsArray->GetNumElements();ii_monitor++) {
                                if(ii == MonitorPointsArray->GetElement(ii_monitor).rank) {
                                    int i_i = MonitorPointsArray->GetElement(ii_monitor).MonitorXY.GetX()/FlowNode2D<FP,NUM_COMPONENTS>::dx - iX0;
                                    int j_j = MonitorPointsArray->GetElement(ii_monitor).MonitorXY.GetY()/FlowNode2D<FP,NUM_COMPONENTS>::dy;
                                    CopyDeviceToHost(&cudaJ[i_i*max_Y + j_j].p,&MonitorPointsArray->GetElement(ii_monitor).p,sizeof(FP),cuda_streams[ii]);
                                    CopyDeviceToHost(&cudaJ[i_i*max_Y + j_j].Tg,&MonitorPointsArray->GetElement(ii_monitor).T,sizeof(FP),cuda_streams[ii]);
                                }
                            }                                       
                        }
                        
                        iX0 += max_X;
                   }

                  CUDA_BARRIER((char*)"cuda_DEEPS2D_Stage2");

                   for(int ii=0;ii<n_s;ii++) {

                       int r_Overlap;
                       int l_Overlap;
                       
                       int max_X;
                       int max_Y;

                       if(cudaSetDevice(ii) != cudaSuccess ) {
                          *f_stream << "\nError set CUDA device no: "<< ii << endl;
                          Exit_OpenHyperFLOW2D(n_s);
                       }

#ifdef _DEVICE_MMAP_
                       dt_min_host = dt_min_host_Array->GetElement(ii);
                       unsigned int  int2float = *dt_min_host;
#else
                       unsigned int  int2float;
                       dt_min_device = dt_min_device_Array->GetElement(ii);
                       CopyDeviceToHost(dt_min_device,&int2float,sizeof(unsigned int),cuda_streams[ii]);
#endif //_DEVICE_MMAP_
                       dtmin  =  min(dtmin,(FP)(int2float)/int2float_scale);

                       if(dtmin == 0) {
                          *f_stream << "\nERROR: Computational unstability  on iteration " << iter+last_iter<< endl;
                          Abort_OpenHyperFLOW2D(n_s);
                        }

                       if(ii == n_s-1)
                         r_Overlap = 0;
                       else
                         r_Overlap = 1;
                       if(ii == 0)
                         l_Overlap = 0;
                       else
                         l_Overlap = 1;
// --- Halo exchange ---
                        if(r_Overlap > 0) {

                            max_X = cudaDimArray->GetElement(ii).GetX();
                            max_Y = cudaDimArray->GetElement(ii).GetY();

                            cudaJ = cudaSubmatrixArray->GetElement(ii);
                            size_t cuda_HalloSize = pJ->GetColSize();

                            void*  cuda_Src  = (void*)((ulong)cudaJ+(max_X*max_Y*sizeof(FlowNode2D<FP,NUM_COMPONENTS>) - 2*cuda_HalloSize));
                            void*  cuda_Dst  = (void*)(cudaArraySubmatrix->GetElement(ii+1));

                            CopyDeviceToDeviceP2P(cuda_Src,
                                                  ii,
                                                  cuda_Dst,
                                                  ii+1,
                                                  cuda_HalloSize,
                                                  cuda_streams[ii]);

                            cuda_Dst  = (void*)((ulong)cuda_Src + cuda_HalloSize);
                            cuda_Src  = (void*)((ulong)cudaArraySubmatrix->GetElement(ii+1)+cuda_HalloSize);

                            if(cudaSetDevice(ii+1) != cudaSuccess ) {
                               *f_stream << "\nError set CUDA device no: "<< ii << endl;
                               Exit_OpenHyperFLOW2D(n_s);
                            }

                            CopyDeviceToDeviceP2P(cuda_Src,
                                                  ii+1,
                                                  cuda_Dst,
                                                  ii,
                                                  cuda_HalloSize,
                                                  cuda_streams[ii+1]);
                        }
// --- Halo exchange ---
                   }

                   if(n_s > 1) {
                       CUDA_BARRIER((char*)"Halo exchange");
                   }
/*
             if(!isAdiabaticWall)
                CalcHeatOnWallSources(pJ,dx,dy,dt);
*/


        CurrentTimePart += dt;
        
        dt = dtmin;
         
        if ( isVerboseOutput && iter/NOutStep*NOutStep == iter ) {
             
             gettimeofday(&mark1,NULL);
             d_time = (FP)(mark1.tv_sec-mark2.tv_sec)+(FP)(mark1.tv_usec-mark2.tv_usec)*1.e-6; 

             if(d_time > 0.)
                 VCOMP = (FP)(NOutStep)/d_time;
             else
                 VCOMP = 0.;
             
             if(iter > 2)
               current_div = max(1,CalibrateThreadBlockSize(current_div,
                                                            opt_thread_block_size,
                                                            opt_round_trip,
                                                            d_time));
             
             memcpy(&mark2,&mark1,sizeof(mark1));
#ifdef _RMS_
             SaveRMS(pRMS_OutFile,last_iter+iter, sum_RMS);

             if(k_max_RMS == i2d_nu_t)
                k_max_RMS +=turb_mod_name_index;

             if(k_max_RMS != -1 && (MonitorNumber == 0 || MonitorNumber == 5))
             *f_stream << "Step No " << iter+last_iter << " maxRMS["<< RMS_Name[k_max_RMS] << "]="<< (FP)(max_RMS*100.) \
                        <<  " % step_time=" << (FP)d_time << " sec (" << (FP)VCOMP <<" step/sec) dt="<< dt <<"\n" << flush;
             else if(MonitorNumber > 0 &&  MonitorNumber < 5 )
                 *f_stream << "Step No " << iter+last_iter << " maxRMS["<< RMS_Name[MonitorNumber-1] << "]="<< (FP)(max_RMS*100.) \
                  <<  " % step_time=" << (FP)d_time << " sec (" << (FP)VCOMP <<" step/sec) dt="<< dt <<"\n" << flush;
             else
             *f_stream << "Step No " << iter+last_iter << " maxRMS["<< k_max_RMS << "]="<< (FP)(max_RMS*100.) \
                        <<  " % step_time=" << (FP)d_time << " sec (" << (FP)VCOMP <<" step/sec) dt="<< dt <<"\n" << flush;
#else
             *f_stream << "Step No " << iter+last_iter <<  "/" << Nstep <<" step_time=" << (FP)d_time << " sec (" << (FP)VCOMP <<" step/sec) dt="<< dt <<
              " blocks=" <<  num_cuda_blocks << " threads="  << num_cuda_threads << " div=" << current_div <<" \n"
                 
                 << flush;
#endif // _RMS_
             if(MonitorPointsArray && MonitorPointsArray->GetNumElements() > 0) {
                SaveMonitors(pMonitors_OutFile,GlobalTime+CurrentTimePart,MonitorPointsArray);
             }
        }
  iter++;
} while((int)iter < Nstep);


  for (unsigned int i=0;i<GlobalSubmatrix->GetNumElements();i++) {

      int r_Overlap, l_Overlap;
      int SubStartIndex = GlobalSubmatrix->GetElementPtr(i)->GetX();
      int SubMaxX = GlobalSubmatrix->GetElementPtr(i)->GetY();

      if(cudaSetDevice(i) != cudaSuccess ) {
          *f_stream << "\nError set CUDA device no: "<< i << endl;
          Exit_OpenHyperFLOW2D(n_s);
       }

      if(i == GlobalSubmatrix->GetNumElements()-1)
        r_Overlap = 0;
      else
        r_Overlap = 1;
      if(i == 0)
        l_Overlap = 0;
      else
        l_Overlap = 1;

       TmpMaxX = (SubMaxX-SubStartIndex) - r_Overlap;
       TmpMatrixPtr = (FlowNode2D<FP,NUM_COMPONENTS>*)((ulong)J->GetMatrixPtr()+(ulong)(sizeof(FlowNode2D<FP,NUM_COMPONENTS>)*(SubStartIndex)*MaxY));

#ifdef _PARALLEL_RECALC_Y_PLUS_

       if(ProblemType == SM_NS &&
          (TurbExtModel == TEM_Spalart_Allmaras ||
           TurbExtModel == TEM_vanDriest ||
           TurbExtModel == TEM_k_eps_Chien)) {

           FP  x0 = SubStartIndex*FlowNode2D<FP,NUM_COMPONENTS>::dx;

           *f_stream << "Parallel recalc y+ on CUDA device No  " << i  << endl;

           cuda_Recalc_y_plus<<<num_cuda_blocks,num_cuda_threads, 0, cuda_streams[i]>>>(cudaSubmatrixArray->GetElement(i),
                                                                                        TmpMaxX*MaxY,
                                                                                        cudaWallNodesArray->GetElement(i),
                                                                                        NumWallNodes,
                                                                                        min(dx,dy),
                                                                                        max((x0+FlowNode2D<FP,NUM_COMPONENTS>::dx*TmpMaxX),
                                                                                                       (FlowNode2D<FP,NUM_COMPONENTS>::dy*MaxY)),
                                                                                        FlowNode2D<FP,NUM_COMPONENTS>::dx,
                                                                                        FlowNode2D<FP,NUM_COMPONENTS>::dy,
                                                                                        MaxY);
           CUDA_BARRIER((char*)"y+ recalc");
       }

#endif // _PARALLEL_RECALC_Y_PLUS_

       CopyDeviceToHost(cudaArraySubmatrix->GetElement(i),TmpMatrixPtr,(sizeof(FlowNode2D<FP,NUM_COMPONENTS>))*(TmpMaxX*MaxY),cuda_streams[i]);
  }
       CUDA_BARRIER((char*)"Data collection");

       if ( isGasSource && SrcList) {
             *f_stream << "\nSet gas sources...";
             SrcList->SetSources2D();
        }
       
       *f_stream << "OK" << endl;

       if ( XCutArray->GetNumElements() > 0 ) {
           for(int i_xcut = 0; i_xcut<(int)XCutArray->GetNumElements(); i_xcut++) {
               XCut* TmpXCut = XCutArray->GetElementPtr(i_xcut);
               *f_stream << "Cut(" << i_xcut + 1  <<") X=" << TmpXCut->x0 << " Y=" << TmpXCut->y0 <<
                  " dY=" << TmpXCut->dy << " MassFlow="<< CalcMassFlowRateX2D(J,TmpXCut->x0,TmpXCut->y0,TmpXCut->dy) << "  (kg/sec*m)" << endl;
           }
       }

#ifdef  _GNUPLOT_
        *f_stream << "\nSave current results in file " << OutFileName << "...\n" << flush; 
        DataSnapshot(OutFileName,WM_REWRITE);
#endif // _GNUPLOT_
        if ( (I/NSaveStep)*NSaveStep == I ) {
#ifdef  _TECPLOT_
             *f_stream << "Add current results to transient solution file " << TecPlotFileName << "...\n" << flush; 
             DataSnapshot(TecPlotFileName,WM_APPEND); 
#endif // _TECPLOT_
         }
         I++;
         gettimeofday(&stop,NULL);
         d_time = (FP)(stop.tv_sec-start.tv_sec)+(FP)(stop.tv_usec-start.tv_usec)*1.e-6; 
         *f_stream << "HyperFLOW/DEEPS computation cycle time=" << (FP)d_time << " sec ( average  speed " << (FP)(Nstep/d_time) <<" step/sec).       \n" << flush;
         f_stream->flush();
         last_iter  += iter;
         iter = 0;
         GlobalTime += CurrentTimePart;
         CurrentTimePart  = 0.;

         if(isOutHeatFluxX) {
          char HeatFluxXFileName[255];
          snprintf(HeatFluxXFileName,255,"HeatFlux-X-%s",OutFileName);
          CutFile(HeatFluxXFileName);
          pHeatFlux_OutFile = OpenData(HeatFluxXFileName);
          SaveXHeatFlux2D(pHeatFlux_OutFile,J,Flow2DList->GetElement(Cp_Flow_index-1),Ts0,y_max,y_min);
          pHeatFlux_OutFile->close();
         }

         if(isOutHeatFluxY) {
          char HeatFluxYFileName[255];
          snprintf(HeatFluxYFileName,255,"HeatFlux-Y-%s",OutFileName);
          CutFile(HeatFluxYFileName);
          pHeatFlux_OutFile = OpenData(HeatFluxYFileName);
          SaveYHeatFlux2D(pHeatFlux_OutFile,J,Ts0);
          pHeatFlux_OutFile->close();
         }
// Sync swap files...
//  Gas area
                     if ( GasSwapData ) {
                       if(isVerboseOutput)
                        *f_stream << "\nSync swap file for gas..." << flush;
#ifdef  _NO_MMAP_
#ifdef _WRITE_LARGE_FILE_
//#warning "Sync file >2Gb..."
                        lseek(fd_g,0,SEEK_SET);
                        ssize_t max_write = 1024L*1024L*1024L;
                        ssize_t one_write = 0L;
                        ssize_t len  = 0L;
                        off_t  off = 0L;
                        char* TmpPtr=(char*)J->GetMatrixPtr();
                        if(J->GetMatrixSize() > max_write) {
                           for(off = 0L,one_write = max_write; len < J->GetMatrixSize(); off += max_write) {
                           len += pwrite64(fd_g,TmpPtr+off,one_write,off);
                           if(J->GetMatrixSize() - len < max_write)
                              one_write = J->GetMatrixSize() - len;
                            }
                        if(len != J->GetMatrixSize())
                           *f_stream << "Error: len(" << len << ") != FileSize(" << J->GetMatrixSize() << ") " << endl << flush;
                        } else {
                           len = pwrite64(fd_g,J->GetMatrixPtr(),J->GetMatrixSize(),0L);
                        }
#else
                        lseek(fd_g,0,SEEK_SET);
                        write(fd_g,J->GetMatrixPtr(),J->GetMatrixSize());
#endif // _WRITE_LARGE_FILE_
#else
                        msync(J->GetMatrixPtr(),J->GetMatrixSize(),MS_SYNC);
#endif // _NO_MMAP_ 
                        *f_stream << "OK" << endl;
                     }

if(MonitorNumber < 5) {
#ifdef _RMS_
    if( max_RMS > ExitMonitorValue )
     MonitorCondition = 1;
   else
     MonitorCondition = 0;
#endif // _RMS_
} else {
   if( GlobalTime + CurrentTimePart <  ExitMonitorValue )
     MonitorCondition = 1;
   else
     MonitorCondition = 0;
}
}while( MonitorCondition );
//---   Save  results ---

#ifdef _DEBUG_0
                ___try {
#endif  // _DEBUG_0
#ifdef  _GNUPLOT_
                    DataSnapshot(OutFileName,WM_REWRITE);
#endif //  _GNUPLOT_
#ifdef _DEBUG_0
                } __except(SysException e) {
                    if ( e == SIGINT ) {
                        *f_stream << "Interrupted by user\n" <<  flush;
                    } else {
                        *f_stream << "Save error data !\n";
                        f_stream->flush();
                        if ( GasSwapData!=0 ) {
#ifdef _REMOVE_SWAPFILE_
                            if ( e != SIGINT ) isDel=1;
#endif //_REMOVE_SWAPFILE_
#ifdef _REMOVE_SWAPFILE_
                                 CloseSwapFile(GasSwapFileName,
                                               GasSwapData,
                                               FileSizeGas,
                                               fd_g,
                                               isDel);
#endif // _REMOVE_SWAPFILE_
                        }
                          Abort_OpenHyperFLOW2D(num_mp);
                    }
                }
                __end_except;
#endif  // _DEBUG_0
//--- End Save results ---
                *f_stream << "\nResults saved in file \"" << OutFileName << "\".\n" << flush;
                f_stream->flush();
#ifdef _DEBUG_0
            }__except( UMatrix2D<FP>*  m) {
                *f_stream << "\n";
                *f_stream << "Error in UMatrix2D<FP> ("<< err_i << "," << err_j <<")." << "\n" << flush;
                f_stream->flush();

                if ( GasSwapData!=0 ) {
#ifdef _REMOVE_SWAPFILE_
                    CloseSwapFile(GasSwapFileName,
                                  GasSwapData,
                                  FileSizeGas,
                                  fd_g,
                                  1);
#endif // _REMOVE_SWAPFILE_
                    useSwapFile = 0;
                } else {
                    delete J;
                }
                Abort_OpenHyperFLOW2D(num_mp);
            } __except( ComputationalMatrix2D*  m) {
                *f_stream << "\n";
                *f_stream << "Error in UMatrix2D< FlowNode2D<FP,NUM_COMPONENTS> ("<< err_i << "," << err_j <<")->";
                if(m->GetMatrixState() == MXS_ERR_OUT_OF_INDEX){
                *f_stream << "MXS_ERR_OUT_OF_INDEX"  << "\n" << flush; 
                } else if(m->GetMatrixState() == MXS_ERR_MEM) {
                *f_stream << "MXS_ERR_MEM"  << "\n" << flush;
                }else if(m->GetMatrixState() == MXS_ERR_MAP) {
                *f_stream << "MXS_ERR_MAP"  << "\n" << flush;
                }
                f_stream->flush();
                if ( GasSwapData!=0 ) {
#ifdef _REMOVE_SWAPFILE_
                    CloseSwapFile(GasSwapFileName,
                                  GasSwapData,
                                  FileSizeGas,
                                  fd_g,
                                  1);
#endif // _REMOVE_SWAPFILE_
                    useSwapFile=0;
                } else {
                    delete J;
                }
                Abort_OpenHyperFLOW2D(num_mp);
            }__except(SysException e) {
                //int isDel=0;
                if ( e == SIGINT ) {
                    *f_stream << "\n";
                    *f_stream << "Interrupted by user\n" <<  flush;
                } else   *f_stream << "Handled system signal: " << (int)e << "\n" << flush;
                f_stream->flush();
                if ( GasSwapData!=0 ) {
#ifdef _REMOVE_SWAPFILE_
                    if ( e != SIGINT ) isDel=1;
#endif //_REMOVE_SWAPFILE_
#ifdef _REMOVE_SWAPFILE_
             CloseSwapFile(GasSwapFileName,
                           GasSwapData,
                           FileSizeGas,
                           fd_g,
                           isDel);
#endif // _REMOVE_SWAPFILE_
                } else {
                  delete J;
                }
                  Abort_OpenHyperFLOW2D(num_mp);
            }__except(void*) {
                *f_stream << "\n";
                *f_stream << "Unknown error  in ("<< err_i << "," << err_j <<")." << "\n" << flush;
                f_stream->flush();
                if ( GasSwapData!=0 ) {
#ifdef _REMOVE_SWAPFILE_
                    CloseSwapFile(GasSwapFileName,
                                  GasSwapData,
                                  FileSizeGas,
                                  fd_g,
                                  1);
#endif // _REMOVE_SWAPFILE_
                    useSwapFile=0;
                } else {
                  delete J;
                }
                  Abort_OpenHyperFLOW2D(num_mp);
            }
            __end_except;
#endif // _DEBUG_0

           *f_stream << "\nReady. Computation finished.\n"  << flush;
};
#endif // _CUDA_

#ifdef _MPI_
void DEEPS2D_Run(ofstream* f_stream
#ifdef _MPI
                ,UMatrix2D< FlowNode2D<FP,NUM_COMPONENTS> >*     pJ,
                 UMatrix2D< FlowNodeCore2D<FP,NUM_COMPONENTS> >* pC,
                 int rank, int last_rank
#endif //_MPI
                 ) {

   // local variables
    FP   n_n,m_m,n_n_1,m_m_1;
    FP   dXX,dYY,AAA;
    int      n1,n2,n3,n4;
    unsigned int j,k;
    unsigned int StartXLocal,MaxXLocal;
    FP   dtdx;
    FP   dtdy;
    FP   dyy;
    FP   dxx;
    FP   dx_1,dy_1; // 1/dx, 1/dy
    FP   d_time;
    FP   t,VCOMP;
    timeval  start, stop, mark1, mark2;
    int      N1,N2,N3,N4;
#ifdef _OPEN_MP
    unsigned int i_max,j_max,k_max;
#endif //_OPEN_MP
#ifdef _MPI
    unsigned int i_max,j_max,k_max;
#endif //_MPI
    FP   dt_min_local;
    FP   _beta[6+NUM_COMPONENTS];
    FP   DD_local[FlowNode2D<FP,NUM_COMPONENTS>::NumEq];
#ifdef _RMS_
    FP   max_RMS;
    int      k_max_RMS;
#endif //_RMS_
#ifndef _MPI
    UMatrix2D< FlowNode2D<FP,NUM_COMPONENTS> >*     pJ=NULL;
    UMatrix2D< FlowNodeCore2D<FP,NUM_COMPONENTS> >* pC=NULL;
    FP*   dt_min;
    int*      i_c;
    int*      j_c;
    FP    dtmin=1.0;
    FP    DD[FlowNode2D<FP,NUM_COMPONENTS>::NumEq];
#ifdef _RMS_
    FP    sum_RMS[FlowNode2D<FP,NUM_COMPONENTS>::NumEq];
    unsigned long sum_iRMS[FlowNode2D<FP,NUM_COMPONENTS>::NumEq];
    UMatrix2D<FP> RMS(FlowNode2D<FP,NUM_COMPONENTS>::NumEq,SubmatrixArray->GetNumElements());
    UMatrix2D<int>    iRMS(FlowNode2D<FP,NUM_COMPONENTS>::NumEq,SubmatrixArray->GetNumElements());
#endif //_RMS_
    UMatrix2D<FP> DD_max(FlowNode2D<FP,NUM_COMPONENTS>::NumEq,SubmatrixArray->GetNumElements());
#else
#ifdef _RMS_
    FP   RMS[FlowNode2D<FP,NUM_COMPONENTS>::NumEq];
    unsigned long iRMS[FlowNode2D<FP,NUM_COMPONENTS>::NumEq];
#endif //_RMS_
    Var_pack DD_max[last_rank+1];
#ifdef _MPI_NB
    MPI::Request  HaloExchange[4];
    MPI::Request  DD_Exchange[2*last_rank];
#endif //_MPI_NB
    unsigned int r_Overlap, l_Overlap;
#endif // _MPI
    int      AddEq=0;

    isScan = 0;
    dyy    = dx/(dx+dy);
    dxx    = dy/(dx+dy);
    dxdy   = dx*dy;
    dx2    = dx*dx;
    dy2    = dy*dy;
    t      = 0;
    n1 = n2 = n3 = n4 = 1;
    d_time = 0.;
#ifndef _MPI
    dt_min = new FP[SubmatrixArray->GetNumElements()];
    i_c    = new int[SubmatrixArray->GetNumElements()];
    j_c    = new int[SubmatrixArray->GetNumElements()];
#endif // _MPI

#ifndef _MPI
    for(int ii=0;ii<(int)SubmatrixArray->GetNumElements();ii++) {
        dt_min[ii] = dtmin = dt;
        i_c[ii] = j_c[ii] = 0;
     }
#ifdef _RMS_
    snprintf(RMSFileName,255,"RMS-%s",OutFileName);
    CutFile(RMSFileName);
    pRMS_OutFile = OpenData(RMSFileName);
    SaveRMSHeader(pRMS_OutFile);
#endif // _RMS_

#else
    MPI::COMM_WORLD.Barrier();
    MPI::COMM_WORLD.Bcast(&dt,1,MPI::DOUBLE,0);
#ifdef _RMS_
    if( rank == 0 ) {
        char    RMSFileName[255];
        snprintf(RMSFileName,255,"RMS-%s",OutFileName);
        CutFile(RMSFileName);
        pRMS_OutFile = OpenData(RMSFileName);
        SaveRMSHeader(pRMS_OutFile);
    }
#endif // _RMS_

    for (int ii=0;(int)ii<last_rank+1;ii++)
         DD_max[ii].dt_min=dt;

#endif // _MPI
#ifdef _DEBUG_0
          ___try {
#endif  // _DEBUG_0

                    I = 0;
                    isRun = 1;
#ifdef _MPI
                    if( rank == 0) {
                        StartXLocal = l_Overlap = 0;
                    } else {
                        StartXLocal = l_Overlap = 1;
                    }

                    if( rank == last_rank) {
                        MaxXLocal=pJ->GetX();
                        r_Overlap = 0;
                    } else {
                        MaxXLocal=pJ->GetX()-1;
                        r_Overlap = 1;
                    }
#else
                    pJ    = J;

                    StartXLocal = 0;
                    MaxXLocal=pJ->GetX();
#endif // _MPI
             do {

                  gettimeofday(&mark2,NULL);
                  gettimeofday(&start,NULL);
                  
                  if( AddSrcStartIter < iter + last_iter){
                   FlowNode2D<FP,NUM_COMPONENTS>::isSrcAdd = 1;
                  } else {
                   FlowNode2D<FP,NUM_COMPONENTS>::isSrcAdd = 0;
                  }
#ifdef _OPENMP
                  n_s = (int)SubmatrixArray->GetNumElements();

#pragma omp parallel shared(f_stream,CoreSubmatrixArray, SubmatrixArray, chemical_reactions,Y_mix,\
                            Cp,i_max,j_max,k_max,Tg,beta0,CurrentTimePart,DD,dx,dy,MaxX,MaxY,dt_min, dtmin,DD_max,i_c,j_c,n_s) \
                     private(iter,j,k,n1,n2,n3,n4,N1,N2,N3,N4,n_n,m_m,pC,pJ,err_i,err_j,\
                             beta,_beta,AddEq,dXX,dYY,DD_local,makeZero,AAA,StartXLocal,MaxXLocal,\
                             dtdx,dtdy,dt)
#endif //_OPENMP

#ifdef _RMS_                           
                            //sum_RMS, sum_iRMS, RMS, iRMS,
#endif // _RMS_
                  iter = 0;
                  do {
#ifdef _MPI
// MPI version      
#ifdef _RMS_
                  if(rank == 0 ) {
                    if(MonitorNumber < 5)
                      max_RMS =  0.5*ExitMonitorValue;
                    else
                      max_RMS = 0.;
                  } else {
                    if(MonitorNumber < 5)
                       max_RMS =  2*ExitMonitorValue;
                    else
                       max_RMS = 0.;
                  }

                  k_max_RMS = -1;
#endif // _RMS_
                  for (int kk=0;kk<FlowNode2D<FP,NUM_COMPONENTS>::NumEq;kk++ ) {
#ifdef _RMS_
                      RMS[kk] = 0.;                                           // Clean sum residual
                      iRMS[kk]= 0;                                            // num involved nodes
                      
                      DD_max[rank].DD[kk].RMS  = 0.;                          // sum residual per rank
                      DD_max[rank].DD[kk].iRMS = 0;                           // num involved nodes per rank
#endif // _RMS_
                      
                      DD_max[rank].DD[kk].DD   = 0.;                          // max residual per rank
                      DD_max[rank].DD[kk].i    = 0;                           // max residual x-coord
                      DD_max[rank].DD[kk].j    = 0;                           // max residual y-coord
                      DD_local[kk]             = 0.;                          // local residual
                    }
#else
#ifdef _OPENMP
//#pragma omp single
#endif //_OPENMP
                   {
#ifdef _RMS_
                       max_RMS = 0.5 * ExitMonitorValue;
                       k_max_RMS = -1;
#endif // _RMS_
                       for ( k=0;k<(int)FlowNode2D<FP,NUM_COMPONENTS>::NumEq;k++ ) {
#ifdef _RMS_
                           for(int ii=0;ii<(int)SubmatrixArray->GetNumElements();ii++) {
                               sum_RMS[k]  = RMS(k,ii)  = 0.;    // Clean sum residual
                               sum_iRMS[k] = iRMS(k,ii) = 0;     // num involved nodes
                           }
#endif // _RMS_
                           DD[k]      = 0.;
                       }
                   }
#endif //_MPI

#ifdef _OPENMP
#pragma omp barrier
#pragma omp for private(CurrentNode,NextNode,UpNode,DownNode,LeftNode,RightNode)  schedule(dynamic) ordered nowait
                   for(int ii=0;ii<n_s;ii++) {  // OpenMP version
#endif //_OPENMP

#ifndef _MPI
#ifndef _OPENMP
                    int ii = 0;                                                    // Single thread version
                    dt = dt_min[0];
#endif //_OPENMP
#else
                     dt  = 1.;
                     for (int ii=0;(int)ii<last_rank+1;ii++) {                     // Choose
                         dt  = min(DD_max[ii].dt_min,dt);                          // minimal
                    }

#endif // _MPI

#ifdef _OPENMP
//#pragma omp critical
                   dtmin = min(dt_min[ii],dtmin);
                   
                   dt    = dtmin;
                   i_c[ii] = j_c[ii] = 0;
#endif //_OPENMP

#ifndef _MPI
                    pJ = SubmatrixArray->GetElement(ii);
                    pC = CoreSubmatrixArray->GetElement(ii);
#endif // _MPI

#ifdef _MPI

#else
                    for ( k=0;k<(int)FlowNode2D<FP,NUM_COMPONENTS>::NumEq;k++ ) {
                       DD_max(k,ii) =  0.;
                      }

                    if( ii == 0)
                       StartXLocal=0;
                    else
                       StartXLocal=1;
                    if( ii == (int)SubmatrixArray->GetNumElements()-1) {
                        MaxXLocal=pJ->GetX();
                    } else {
                        MaxXLocal=pJ->GetX()-1;
                    }
#endif // _MPI
                    dx_1 = 1.0/dx;
                    dy_1 = 1.0/dy;

                    dtdx = dt/dx;
                    dtdy = dt/dy;

#ifdef _OPENMP
                    dt_min[ii]   = 1.;
#else
#ifdef _MPI
                    DD_max[rank].dt_min = dt_min_local = 1.;
#endif // _MPI
#endif // _OPENMP
                    DEEPS2D_Stage1(pJ,pC,StartXLocal,MaxXLocal,MaxY,dxx,dyy,dtdx,dtdy);
#ifdef _OPENMP
           }
#endif // _OPENMP
        
#ifdef _OPENMP
#pragma omp barrier
#pragma omp for private(CurrentNode,NextNode,UpNode,DownNode,LeftNode,RightNode) schedule(dynamic) ordered nowait
             for(int ii=0;ii<n_s;ii++) {

                    pJ = SubmatrixArray->GetElement(ii);
                    pC = CoreSubmatrixArray->GetElement(ii);

                    if( ii == 0)
                      StartXLocal=0;
                    else
                      StartXLocal=1;
                    if( ii == (int)SubmatrixArray->GetNumElements()-1) {
                        MaxXLocal=pJ->GetX();
                    } else {
                        MaxXLocal=pJ->GetX()-1;
                    }
#endif //_OPENMP
#ifdef _MPI                    
                    DD_max[rank].dt_min = min(DD_max[rank].dt_min,
#else
                    dt_min[ii] = min(dt_min[ii],

#endif //_MPI 
                    DEEPS2D_Stage2(pJ,pC,StartXLocal,
                                   MaxXLocal,MaxY,beta0,
                                   bFF,&chemical_reactions,
                                   iter,last_iter,TurbStartIter,
                                   SigW,SigF,dx_1,dy_1,delta_bl,
                                   CFL_Scenario->GetVal(iter+last_iter),
                                   beta_Scenario->GetVal(iter+last_iter),
#ifdef _RMS_        
#ifdef _MPI
                                   DD_max, rank,
#else
                                   RMS, iRMS, DD_max, i_c, j_c, ii,
#endif //_MPI 
#endif // _RMS_
                                   (TurbulenceExtendedModel)TurbExtModel));
                    
#ifdef _MPI
                    if(DD_max[rank].dt_min == 0.0) {
                        *f_stream << "\nERROR: Computational unstability  on iteration " << iter+last_iter<< " in subdomain no " << rank << "\n";
#else
                    if(dt_min[ii] == 0.0) {
                        *f_stream << "\nERROR: Computational unstability  on iteration " << iter+last_iter<< " in subdomain no " << ii << "\n";
#endif //_MPI 
                        Abort_OpenHyperFLOW2D();
                    }
#ifdef _MPI
// --- Halo exchange ---
                 if(rank < last_rank) {
// Send Tail
                 void*  tmp_SendPtr  = (void*)((ulong)pJ->GetMatrixPtr()+
                                                     (pJ->GetMatrixSize()-2*pJ->GetColSize()));
                 u_long tmp_SendSize = pJ->GetColSize();
#ifdef _MPI_NB
                 HaloExchange[0] = MPI::COMM_WORLD.Isend(tmp_SendPtr,
                                                         tmp_SendSize,
                                                         MPI::BYTE,rank+1,
                                                         tag_MatrixTail);
#else
                 MPI::COMM_WORLD.Send(tmp_SendPtr,
                                      tmp_SendSize,
                                      MPI::BYTE,rank+1,tag_MatrixTail);
#endif //_MPI_NB

// Recive Head
                 void*  tmp_RecvPtr  = (void*)((ulong)pJ->GetMatrixPtr()+
                                                     (pJ->GetMatrixSize()-
                                                      pJ->GetColSize()));
                 u_long tmp_RecvSize = pJ->GetColSize();
#ifdef _MPI_NB
                 HaloExchange[1] = MPI::COMM_WORLD.Irecv(tmp_RecvPtr,
                                                         tmp_RecvSize,
                                                         MPI::BYTE,rank+1,
                                                         tag_MatrixHead);
#else
                 MPI::COMM_WORLD.Recv(tmp_RecvPtr,
                                      tmp_RecvSize,
                                      MPI::BYTE,rank+1,tag_MatrixHead);
#endif //_MPI_NB

             }
               if(rank > 0) {
// Recive Tail
                 void*  tmp_RecvPtr  = (void*)(pJ->GetMatrixPtr());
                 u_long tmp_RecvSize = pJ->GetColSize();
#ifdef _MPI_NB
                 HaloExchange[3] = MPI::COMM_WORLD.Irecv(tmp_RecvPtr,
                                                         tmp_RecvSize,
                                                         MPI::BYTE,rank-1,
                                                         tag_MatrixTail);
#else
                 MPI::COMM_WORLD.Recv(tmp_RecvPtr,
                                      tmp_RecvSize,
                                      MPI::BYTE,rank-1,tag_MatrixTail);
#endif //_MPI_NB

//Send  Head
                 void*  tmp_SendPtr  = (void*)((ulong)pJ->GetMatrixPtr()+
                                                      pJ->GetColSize());
                 u_long tmp_SendSize = pJ->GetColSize();
#ifdef _MPI_NB
                 HaloExchange[2] = MPI::COMM_WORLD.Isend(tmp_SendPtr,
                                                         tmp_SendSize,
                                                         MPI::BYTE,rank-1,tag_MatrixHead);
#else
                 MPI::COMM_WORLD.Send(tmp_SendPtr,
                                      tmp_SendSize,
                                      MPI::BYTE,rank-1,tag_MatrixHead);
#endif //_MPI_NB

             }
#endif // _MPI

             if(!isAdiabaticWall)
                CalcHeatOnWallSources(pJ,dx,dy,dt
#ifdef _MPI
                                      ,rank,last_rank
#else
                                      ,ii,(int)SubmatrixArray->GetNumElements()-1 
#endif // _MPI
                                   );
#ifdef _OPENMP
 for (DD_max_var = 1.,k=0;k<(int)FlowNode2D<FP,NUM_COMPONENTS>::NumEq;k++ ) {
    DD_max_var = DD[k] = max(DD[k],DD_max(k,ii));
    if(DD[k] == DD_max(k,ii)) {
       i_max =  i_c[ii];
       j_max =  j_c[ii];
       k_max =  k;
     }
  }
#pragma omp critical
     dtmin = min(dt_min[ii],dtmin); 
     dt    = dtmin;
#endif // _OPENMP

#ifdef _MPI
#ifdef _MPI_NB
        if(rank < last_rank) {
           HaloExchange[0].Wait();
           HaloExchange[1].Wait();
        }
        if(rank > 0) {
           HaloExchange[3].Wait();
           HaloExchange[2].Wait();
        }
       // MPI::Request::Waitall(4, HaloExchange);
#endif //_MPI_NB

       if(rank > 0) {
#ifdef _MPI_NB
           DD_Exchange[rank-1] = MPI::COMM_WORLD.Isend(&DD_max[rank],sizeof(Var_pack),MPI::BYTE,0,tag_DD);
#else
           MPI::COMM_WORLD.Send(&DD_max[rank],sizeof(Var_pack),MPI::BYTE,0,tag_DD);
#endif //_MPI_NB

       } else {
         for(int ii=1;ii<last_rank+1;ii++) {
#ifdef _MPI_NB
           DD_Exchange[last_rank+ii] = MPI::COMM_WORLD.Irecv(&DD_max[ii],sizeof(Var_pack),MPI::BYTE,ii,tag_DD);
#else
           MPI::COMM_WORLD.Recv(&DD_max[ii],sizeof(Var_pack),MPI::BYTE,ii,tag_DD);
#endif //_MPI_NB

         }
       }

#ifdef _MPI_NB
       if(rank > 0) {
          DD_Exchange[rank-1].Wait();
       } else {
         for(int ii=1;ii<last_rank+1;ii++) {
          DD_Exchange[last_rank+ii].Wait();  
         }
       }
       //MPI::Request::Waitall(2*(last_rank), DD_Exchange);
#endif //_MPI_NB

       if(rank == 0) {
           for(int ii=0,DD_max_var=0.;ii<last_rank+1;ii++) {
               for (k=0;k<(int)(FlowNode2D<FP,NUM_COMPONENTS>::NumEq);k++ ) {
                  DD_max_var = max(DD_max[ii].DD[k].DD,DD_max_var);
#ifdef _RMS_
                  RMS[k] += DD_max[ii].DD[k].RMS;
                  iRMS[k]+= DD_max[ii].DD[k].iRMS;
#endif // _RMS_
                  if(DD_max[ii].DD[k].DD == DD_max_var) {
                   i_max = DD_max[ii].DD[k].i;
                   j_max = DD_max[ii].DD[k].j;
                   k_max = k;
                 }
             }
           }
#ifdef _RMS_               
           for (k=0;k<(int)(FlowNode2D<FP,NUM_COMPONENTS>::NumEq);k++ ) {
                   if(iRMS[k] > 0) {
                     RMS[k] = sqrt(RMS[k]/iRMS[k]);
                     if(MonitorNumber == 0 || MonitorNumber > 4) {
                        max_RMS = max(RMS[k],max_RMS);
                        if(max_RMS == RMS[k])
                           k_max_RMS = k;
                     } else {
                        max_RMS = max(RMS[MonitorNumber-1],max_RMS);
                        if(max_RMS == RMS[MonitorNumber-1])
                           k_max_RMS = k;
                     }

                    }
            }
#endif // _RMS_
#else
#ifdef _OPENMP
 }

#pragma omp single
#endif // _OPENMP
    {
#ifdef _RMS_
        for(k=0;k<(int)(FlowNode2D<FP,NUM_COMPONENTS>::NumEq);k++ ) {
         for(int ii=0;ii<(int)SubmatrixArray->GetNumElements();ii++) {
              if(iRMS(k,ii) > 0) {
                 sum_RMS[k]  += RMS(k,ii);
                 sum_iRMS[k] += iRMS(k,ii);
             }
           }

           if(sum_iRMS[k] != 0)
              sum_RMS[k] = sum_RMS[k]/sum_iRMS[k];
           else
              sum_RMS[k] = 0;

           if( MonitorNumber == 0 || MonitorNumber > 4) {
               max_RMS = max(sum_RMS[k],max_RMS);
               if(max_RMS == sum_RMS[k])
                  k_max_RMS = k;
           } else {
               max_RMS = max(sum_RMS[MonitorNumber-1],max_RMS);
               if(max_RMS == sum_RMS[MonitorNumber-1])
                  k_max_RMS = MonitorNumber-1;
           }

         }
#endif // _RMS_

#endif // _MPI
        CurrentTimePart += dt;
         if ( isVerboseOutput && iter/NOutStep*NOutStep == iter ) {
             gettimeofday(&mark1,NULL);
             d_time = (FP)(mark1.tv_sec-mark2.tv_sec)+(FP)(mark1.tv_usec-mark2.tv_usec)*1.e-6; 

             if(d_time > 0.)
                 VCOMP = (FP)(NOutStep)/d_time;
             else
                 VCOMP = 0.;
             memcpy(&mark2,&mark1,sizeof(mark1));
#ifdef _RMS_
             SaveRMS(pRMS_OutFile,last_iter+iter,
#ifdef  _MPI
             RMS);
#else
             sum_RMS);
#endif // _MPI

             if(k_max_RMS == i2d_nu_t)
                k_max_RMS +=turb_mod_name_index;

             if(k_max_RMS != -1 && (MonitorNumber == 0 || MonitorNumber == 5))
             *f_stream << "Step No " << iter+last_iter << " maxRMS["<< RMS_Name[k_max_RMS] << "]="<< (FP)(max_RMS*100.) \
                        <<  " % step_time=" << (FP)d_time << " sec (" << (FP)VCOMP <<" step/sec) dt="<< dt <<"\n" << flush;
             else if(MonitorNumber > 0 &&  MonitorNumber < 5 )
                 *f_stream << "Step No " << iter+last_iter << " maxRMS["<< RMS_Name[MonitorNumber-1] << "]="<< (FP)(max_RMS*100.) \
                  <<  " % step_time=" << (FP)d_time << " sec (" << (FP)VCOMP <<" step/sec) dt="<< dt <<"\n" << flush;
             else
             *f_stream << "Step No " << iter+last_iter << " maxRMS["<< k_max_RMS << "]="<< (FP)(max_RMS*100.) \
                        <<  " % step_time=" << (FP)d_time << " sec (" << (FP)VCOMP <<" step/sec) dt="<< dt <<"\n" << flush;
#else
             *f_stream << "Step No " << iter+last_iter <<  " step_time=" << (FP)d_time << " sec (" << (FP)VCOMP <<" step/sec) dt="<< dt <<"\n" << flush;
#endif // _RMS_
             f_stream->flush();
            }
#ifndef _MPI
#endif // _MPI
        }
     iter++;
   } while((int)iter < Nstep);
#ifdef _OPENMP
#pragma omp single 
          {
#endif //  _OPENMP

#ifdef _MPI
     // Collect all submatrix
        if(rank>0) {
        void*  tmp_SendPtr=(void*)((u_long)(pJ->GetMatrixPtr())+
                                            pJ->GetColSize()*l_Overlap);
        u_long tmp_SendSize=pJ->GetMatrixSize()-
                            pJ->GetColSize()*(l_Overlap);
#ifdef _IMPI_
        LongMatrixSend(0, tmp_SendPtr, tmp_SendSize);  // Low Mem Send subdomain
#else
        MPI::COMM_WORLD.Send(tmp_SendPtr,
                             tmp_SendSize,
                             MPI::BYTE,0,tag_Matrix); 
#endif // _IMPI_
        } else {
        void*  tmp_RecvPtr;
        u_long tmp_RecvSize;
        for(int ii=1;ii<last_rank+1;ii++) {
           tmp_RecvPtr=(void*)((u_long)(ArraySubmatrix->GetElement(ii)->GetMatrixPtr())+
                                        ArraySubmatrix->GetElement(ii)->GetColSize()*(r_Overlap+
                                                                                      l_Overlap));
           tmp_RecvSize=ArraySubmatrix->GetElement(ii)->GetMatrixSize()-
                        ArraySubmatrix->GetElement(ii)->GetColSize()*(r_Overlap-l_Overlap);
#ifdef _IMPI_
           LongMatrixRecv(ii, tmp_RecvPtr, tmp_RecvSize);  // Low Mem Send subdomain
#else
           MPI::COMM_WORLD.Recv(tmp_RecvPtr,
                                tmp_RecvSize,
                                MPI::BYTE,ii,tag_Matrix); 
#endif // _IMPI_
        }
#ifdef _PARALLEL_RECALC_Y_PLUS_
        if(WallNodes && !WallNodesUw_2D && J) 
             WallNodesUw_2D = GetWallFrictionVelocityArray2D(J,WallNodes);
         else if(WallNodes && J)
             RecalcWallFrictionVelocityArray2D(J,WallNodesUw_2D,WallNodes);
#endif // _PARALLEL_RECALC_Y_PLUS_
    }

#ifdef _PARALLEL_RECALC_Y_PLUS_
    if( rank == 0 ) {
        *f_stream << "Parallel recalc y+...";
        for(int ii = 1; (int)ii < last_rank + 1; ii++ ) {
            MPI::COMM_WORLD.Send(WallNodesUw_2D->GetArrayPtr(),
                                 WallNodesUw_2D->GetNumElements()*WallNodesUw_2D->GetElementSize(),
                                 MPI::BYTE,ii,tag_WallFrictionVelocity);
        }
    } else {
            MPI::COMM_WORLD.Recv(WallNodesUw_2D->GetArrayPtr(),
                                 WallNodesUw_2D->GetNumElements()*WallNodesUw_2D->GetElementSize(),
                                 MPI::BYTE,0,tag_WallFrictionVelocity);
    }
    ParallelRecalc_y_plus(pJ,WallNodes,WallNodesUw_2D,x0);
#endif // _PARALLEL_RECALC_Y_PLUS_

    if( rank == 0) {
#ifndef _PARALLEL_RECALC_Y_PLUS_
        *f_stream << "Recalc y+...";
         Recalc_y_plus(J,WallNodes);
#endif // _PARALLEL_RECALC_Y_PLUS_
#else
#ifdef _PARALLEL_RECALC_Y_PLUS_
         *f_stream << "Parallel recalc y+...";
         ParallelRecalc_y_plus(pJ,WallNodes,WallNodesUw_2D,x0);
#else
   *f_stream << "Recalc y+...";
   Recalc_y_plus(J,WallNodes);
#endif // _PARALLEL_RECALC_Y_PLUS_
#endif //  _MPI
         
        if ( isGasSource && SrcList) {
              *f_stream << "\nSet gas sources..." << endl;
              SrcList->SetSources2D();
         }
        if ( XCutArray->GetNumElements() > 0 ) {
            for(int i_xcut = 0; i_xcut<(int)XCutArray->GetNumElements(); i_xcut++) {
                XCut* TmpXCut = XCutArray->GetElementPtr(i_xcut);
                *f_stream << "Cut(" << i_xcut + 1  <<") X=" << TmpXCut->x0 << " Y=" << TmpXCut->y0 <<
                   " dY=" << TmpXCut->dy << " MassFlow="<< CalcMassFlowRateX2D(J,TmpXCut->x0,TmpXCut->y0,TmpXCut->dy) << "  (kg/sec*m)" << endl;
            }
        }
        gettimeofday(&stop,NULL);
#ifdef  _GNUPLOT_
        *f_stream << "\nSave current results in file " << OutFileName << "...\n" << flush; 
        DataSnapshot(OutFileName,WM_REWRITE);
#endif // _GNUPLOT_
        if ( (I/NSaveStep)*NSaveStep == I ) {
#ifdef  _TECPLOT_
             *f_stream << "Add current results to transient solution file " << TecPlotFileName << "...\n" << flush; 
             DataSnapshot(TecPlotFileName,WM_APPEND); 
#endif // _TECPLOT_
         }
         I++;
         d_time = (FP)(stop.tv_sec-start.tv_sec)+(FP)(stop.tv_usec-start.tv_usec)*1.e-6; 
         *f_stream << "HyperFLOW/DEEPS computation cycle time=" << (FP)d_time << " sec ( average  speed " << (FP)(Nstep/d_time) <<" step/sec).       \n" << flush;
         f_stream->flush();
         last_iter  += iter;
         iter = 0;
         GlobalTime += CurrentTimePart;
         CurrentTimePart  = 0.;

         if(isOutHeatFluxX) {
          char HeatFluxXFileName[255];
          snprintf(HeatFluxXFileName,255,"HeatFlux-X-%s",OutFileName);
          CutFile(HeatFluxXFileName);
          pHeatFlux_OutFile = OpenData(HeatFluxXFileName);
          SaveXHeatFlux2D(pHeatFlux_OutFile,J,Ts0);
          pHeatFlux_OutFile->close();
         }

         if(isOutHeatFluxY) {
          char HeatFluxYFileName[255];
          snprintf(HeatFluxYFileName,255,"HeatFlux-Y-%s",OutFileName);
          CutFile(HeatFluxYFileName);
          pHeatFlux_OutFile = OpenData(HeatFluxYFileName);
          SaveYHeatFlux2D(pHeatFlux_OutFile,J,Ts0);
          pHeatFlux_OutFile->close();
         }
// Sync swap files...
//  Gas area
                     if ( GasSwapData ) {
                       if(isVerboseOutput)
                        *f_stream << "\nSync swap file for gas..." << flush;
#ifdef  _NO_MMAP_
#ifdef _WRITE_LARGE_FILE_
//#warning "Sync file >2Gb..."
                        lseek(fd_g,0,SEEK_SET);
                        ssize_t max_write = 1024L*1024L*1024L;
                        ssize_t one_write = 0L;
                        ssize_t len  = 0L;
                        off_t  off = 0L;
                        char* TmpPtr=(char*)J->GetMatrixPtr();
                        if(J->GetMatrixSize() > max_write) {
                           for(off = 0L,one_write = max_write; len < J->GetMatrixSize(); off += max_write) {
                           len += pwrite64(fd_g,TmpPtr+off,one_write,off);
                           if(J->GetMatrixSize() - len < max_write)
                              one_write = J->GetMatrixSize() - len;
                            }
                        if(len != J->GetMatrixSize())
                           *f_stream << "Error: len(" << len << ") != FileSize(" << J->GetMatrixSize() << ") " << endl << flush;
                        } else {
                           len = pwrite64(fd_g,J->GetMatrixPtr(),J->GetMatrixSize(),0L);
                        }
#else
                        lseek(fd_g,0,SEEK_SET);
                        write(fd_g,J->GetMatrixPtr(),J->GetMatrixSize());
#endif // _WRITE_LARGE_FILE_
#else
                        msync(J->GetMatrixPtr(),J->GetMatrixSize(),MS_SYNC);
#endif //  _GPFS
                        *f_stream << "OK" << endl;
                     }
#ifdef _MPI
      }
// Get Data back
    MPI::COMM_WORLD.Bcast(&GlobalTime,1,MPI::DOUBLE,0);
    if(rank>0) {
        void*  tmp_RecvPtr=(void*)((u_long)(pJ->GetMatrixPtr())+pJ->GetColSize()*l_Overlap);
        u_long tmp_RecvSize=pJ->GetMatrixSize()-pJ->GetColSize()*(l_Overlap);
#ifdef _IMPI_
        LongMatrixRecv(0, tmp_RecvPtr, tmp_RecvSize);  // Low Mem Send subdomain
#else
        MPI::COMM_WORLD.Recv(tmp_RecvPtr,
                             tmp_RecvSize,
                             MPI::BYTE,0,tag_Matrix);
#endif //_IMPI_
#ifdef _PARALLEL_RECALC_Y_PLUS_
        if(!WallNodesUw_2D)
            WallNodesUw_2D = new UArray<FP>(NumWallNodes,-1);

        MPI::COMM_WORLD.Recv(WallNodesUw_2D->GetArrayPtr(),
                             NumWallNodes*sizeof(FP),
                             MPI::BYTE,0,tag_WallFrictionVelocity);
#endif // _PARALLEL_RECALC_Y_PLUS_
    } else {
      void*  tmp_SendPtr;
      u_long tmp_SendSize;
      for(int ii=1;ii<last_rank+1;ii++) {
          tmp_SendPtr=(void*)((u_long)(ArraySubmatrix->GetElement(ii)->GetMatrixPtr())+
          ArraySubmatrix->GetElement(ii)->GetColSize()*(r_Overlap+l_Overlap));
          tmp_SendSize=ArraySubmatrix->GetElement(ii)->GetMatrixSize()-
          ArraySubmatrix->GetElement(ii)->GetColSize()*(r_Overlap-l_Overlap);
#ifdef _IMPI_
          LongMatrixSend(ii, tmp_SendPtr, tmp_SendSize);  // Low Mem Send subdomain
#else
          MPI::COMM_WORLD.Send(tmp_SendPtr,
                               tmp_SendSize,
                               MPI::BYTE,ii,tag_Matrix);
#endif //_IMPI_ 
#ifdef _PARALLEL_RECALC_Y_PLUS_
          MPI::COMM_WORLD.Send(WallNodesUw_2D->GetArrayPtr(),
                               NumWallNodes*sizeof(FP),
                               MPI::BYTE,ii,tag_WallFrictionVelocity);
#endif // _PARALLEL_RECALC_Y_PLUS_
      }
    }
     MPI::COMM_WORLD.Bcast(&last_iter,1,MPI::INT,0);
     MPI::COMM_WORLD.Bcast(&isRun,1,MPI::INT,0);
     if(!isRun)
         Abort_OpenHyperFLOW2D();
#endif //  _MPI
#ifdef _PROFILE_
 Abort_OpenHyperFLOW2D();
#endif // _PROFILE_
#ifdef _MPI
 MPI::COMM_WORLD.Bcast(&isRun,1,MPI::INT,0);
#endif //  _MPI
#ifdef _OPENMP
          }
#pragma omp barrier
#endif //  _OPENMP

if(MonitorNumber < 5) {
#ifdef _RMS_
    if( max_RMS > ExitMonitorValue )
     MonitorCondition = 1;
   else
     MonitorCondition = 0;
#endif // _RMS_
} else {
   if( GlobalTime <  ExitMonitorValue )
     MonitorCondition = 1;
   else
     MonitorCondition = 0;
}
}while( MonitorCondition );
//---   Save  results ---
#ifdef _MPI
if (rank == 0) {
#endif //  _MPI

#ifdef _DEBUG_0
                ___try {
#endif  // _DEBUG_0
#ifdef  _GNUPLOT_
                    DataSnapshot(OutFileName,WM_REWRITE);
#endif //  _GNUPLOT_
#ifdef _DEBUG_0
                } __except(SysException e) {
                    if ( e == SIGINT ) {
                        *f_stream << "Interrupted by user\n" <<  flush;
                    } else {
                        *f_stream << "Save error data !\n";
                        f_stream->flush();
                        if ( GasSwapData!=0 ) {
#ifdef _REMOVE_SWAPFILE_
                            if ( e != SIGINT ) isDel=1;
#endif //_REMOVE_SWAPFILE_
#ifdef _REMOVE_SWAPFILE_
                                 CloseSwapFile(GasSwapFileName,
                                               GasSwapData,
                                               FileSizeGas,
                                               fd_g,
                                               isDel);
#endif // _REMOVE_SWAPFILE_
                        }
                          Abort_OpenHyperFLOW2D();
                    }
                }
                __end_except;
#endif  // _DEBUG_0
//--- End Save results ---
                *f_stream << "\nResults saved in file \"" << OutFileName << "\".\n" << flush;
                f_stream->flush();
#ifdef _MPI
                isRun = 0;
                Exit_OpenHyperFLOW2D();
              }
#endif //  _MPI
#ifdef _DEBUG_0
            }__except( UMatrix2D<FP>*  m) {
                *f_stream << "\n";
                *f_stream << "Error in UMatrix2D<FP> ("<< err_i << "," << err_j <<")." << "\n" << flush;
                f_stream->flush();

                if ( GasSwapData!=0 ) {
#ifdef _REMOVE_SWAPFILE_
                    CloseSwapFile(GasSwapFileName,
                                  GasSwapData,
                                  FileSizeGas,
                                  fd_g,
                                  1);
#endif // _REMOVE_SWAPFILE_
                    useSwapFile = 0;
                } else {
                    delete J;
                }
                Abort_OpenHyperFLOW2D();
            } __except( ComputationalMatrix2D*  m) {
                *f_stream << "\n";
                *f_stream << "Error in UMatrix2D< FlowNode2D<FP,NUM_COMPONENTS> ("<< err_i << "," << err_j <<")->";
                if(m->GetMatrixState() == MXS_ERR_OUT_OF_INDEX){
                *f_stream << "MXS_ERR_OUT_OF_INDEX"  << "\n" << flush; 
                } else if(m->GetMatrixState() == MXS_ERR_MEM) {
                *f_stream << "MXS_ERR_MEM"  << "\n" << flush;
                }else if(m->GetMatrixState() == MXS_ERR_MAP) {
                *f_stream << "MXS_ERR_MAP"  << "\n" << flush;
                }
                f_stream->flush();
                if ( GasSwapData!=0 ) {
#ifdef _REMOVE_SWAPFILE_
                    CloseSwapFile(GasSwapFileName,
                                  GasSwapData,
                                  FileSizeGas,
                                  fd_g,
                                  1);
#endif // _REMOVE_SWAPFILE_
                    useSwapFile=0;
                } else {
                    delete J;
                }
                Abort_OpenHyperFLOW2D();
            }__except(SysException e) {
                //int isDel=0;
                if ( e == SIGINT ) {
                    *f_stream << "\n";
                    *f_stream << "Interrupted by user\n" <<  flush;
                } else   *f_stream << "Handled system signal: " << (int)e << "\n" << flush;
                f_stream->flush();
                if ( GasSwapData!=0 ) {
#ifdef _REMOVE_SWAPFILE_
                    if ( e != SIGINT ) isDel=1;
#endif //_REMOVE_SWAPFILE_
#ifdef _REMOVE_SWAPFILE_
             CloseSwapFile(GasSwapFileName,
                           GasSwapData,
                           FileSizeGas,
                           fd_g,
                           isDel);
#endif // _REMOVE_SWAPFILE_
                } else {
                  delete J;
                }
                  Abort_OpenHyperFLOW2D();
            }__except(void*) {
                *f_stream << "\n";
                *f_stream << "Unknown error  in ("<< err_i << "," << err_j <<")." << "\n" << flush;
                f_stream->flush();
                if ( GasSwapData!=0 ) {
#ifdef _REMOVE_SWAPFILE_
                    CloseSwapFile(GasSwapFileName,
                                  GasSwapData,
                                  FileSizeGas,
                                  fd_g,
                                  1);
#endif // _REMOVE_SWAPFILE_
                    useSwapFile=0;
                } else {
                  delete J;
                }
                  Abort_OpenHyperFLOW2D();
            }
            __end_except;
#endif // _DEBUG_0

#ifdef _MPI
           if(rank == 0)
#endif // _MPI
             *f_stream << "\nReady. Computation finished.\n"  << flush;

           Exit_OpenHyperFLOW2D();
};
#endif // _MPI_



#ifdef _PARALLEL_RECALC_Y_PLUS_
void ParallelRecalc_y_plus(ComputationalMatrix2D* pJ, 
                           UArray< XY<int> >* WallNodes,
                           UArray<FP>* WallFrictionVelocity2D,
                           FP x0) {
#ifdef _OPEN_MP
#pragma omp parallel for
#endif //_OPEN_MP
        for (int i=0;i<(int)pJ->GetX();i++ ) {
            for (int j=0;j<(int)pJ->GetY();j++ ) {
                 if ( pJ->GetValue(i,j).isCond2D(CT_NODE_IS_SET_2D) &&
                      !pJ->GetValue(i,j).isCond2D(CT_SOLID_2D)) {

                        for (int ii=0;ii<(int)WallNodes->GetNumElements();ii++) {

                             unsigned int iw,jw;
                             FP U_w,x,y,wx,wy;

                             iw = WallNodes->GetElementPtr(ii)->GetX();
                             jw = WallNodes->GetElementPtr(ii)->GetY();

                             U_w   =  WallFrictionVelocity2D->GetElement(ii);

                             x = x0 + i * FlowNode2D<FP,NUM_COMPONENTS>::dx;
                             y = j * FlowNode2D<FP,NUM_COMPONENTS>::dy;

                             wx = iw * FlowNode2D<FP,NUM_COMPONENTS>::dx;
                             wy = jw * FlowNode2D<FP,NUM_COMPONENTS>::dy;

                             if(x == wx && y == wy ) {
                                 pJ->GetValue(i,j).y_plus = U_w*min((FlowNode2D<FP,NUM_COMPONENTS>::dx),
                                                                    (FlowNode2D<FP,NUM_COMPONENTS>::dy))*pJ->GetValue(i,j).S[i2d_Ro]/pJ->GetValue(i,j).mu;
                             } else {
                                 if(pJ->GetValue(i,j).l_min == min(pJ->GetValue(i,j).l_min,
                                                                   sqrt((x-wx)*(x-wx) + (y-wy)*(y-wy))))
                                 pJ->GetValue(i,j).y_plus = U_w*pJ->GetValue(i,j).l_min*pJ->GetValue(i,j).S[i2d_Ro]/pJ->GetValue(i,j).mu;
                             }
                        }
                 }
            }
        }
}
#else
void Recalc_y_plus(ComputationalMatrix2D* pJ, UArray< XY<int> >* WallNodes) {
    unsigned int iw,jw;
    FP tau_w, U_w;

    for (int ii=0;ii<(int)WallNodes->GetNumElements();ii++) { 

        iw = WallNodes->GetElementPtr(ii)->GetX();
        jw = WallNodes->GetElementPtr(ii)->GetY();

        tau_w = (fabs(pJ->GetValue(iw,jw).dUdy) + 
                 fabs(pJ->GetValue(iw,jw).dVdx)) * pJ->GetValue(iw,jw).mu;

        U_w   = sqrt(tau_w/pJ->GetValue(iw,jw).S[i2d_Ro]+1e-30);

        for (int i=0;i<(int)pJ->GetX();i++ ) {
               for (int j=0;j<(int)pJ->GetY();j++ ) {

                       if (pJ->GetValue(i,j).isCond2D(CT_NODE_IS_SET_2D) &&
                           !pJ->GetValue(i,j).isCond2D(CT_SOLID_2D)) {
                               if ( i == (int)iw && j == (int)jw) {
                                     pJ->GetValue(i,j).y_plus = U_w*min(dx,dy)/pJ->GetValue(i,j).mu;
                               } else { 
                                   if(pJ->GetValue(i,j).l_min == max(min(dx,dy),sqrt((i-iw)*(i-iw)*dx*dx + 
                                                                                     (j-jw)*(j-jw)*dy*dy)))
                                      pJ->GetValue(i,j).y_plus = U_w*pJ->GetValue(i,j).l_min*pJ->GetValue(i,j).S[i2d_Ro]/pJ->GetValue(i,j).mu;
                               }
                         }
                       }
               }
    }
}
#endif // _PARALLEL_RECALC_Y_PLUS_


inline  void CalcHeatOnWallSources(UMatrix2D< FlowNode2D<FP,NUM_COMPONENTS> >* F, FP dx, FP dy, FP dt, int rank, int last_rank) {

        unsigned int StartXLocal,MaxXLocal;
        FP dx_local, dy_local;

        if( rank == 0) {
            StartXLocal = 0;
        } else {
            StartXLocal = 1;
        }

        if( rank == last_rank) {
            MaxXLocal=F->GetX();
        } else {
            MaxXLocal=F->GetX()-1;
        }
        // Clean Q
        for (unsigned int i=StartXLocal;i<MaxXLocal;i++ )
            for ( unsigned int j=0;j<F->GetY();j++ ) {
                 F->GetValue(i,j).Q_conv = 0.;
            }

        for (unsigned int i=StartXLocal;i<MaxXLocal;i++ )
            for ( unsigned int j=0;j<F->GetY();j++ ) {

                FlowNode2D< FP,NUM_COMPONENTS >* CurrentNode=NULL;
                FlowNode2D< FP,NUM_COMPONENTS >* UpNode=NULL;
                FlowNode2D< FP,NUM_COMPONENTS >* DownNode=NULL;
                FlowNode2D< FP,NUM_COMPONENTS >* LeftNode=NULL;
                FlowNode2D< FP,NUM_COMPONENTS >* RightNode=NULL;


                CurrentNode = &F->GetValue(i,j);

                if(!CurrentNode->isCond2D(CT_SOLID_2D) && CurrentNode->isCleanSources == 1) {
                    
                    dx_local = dx;
                    dy_local = dy;

                    if(j < F->GetY()-1)
                      UpNode    = &(F->GetValue(i,j+1));
                    else
                      UpNode    = NULL;

                    if(j>0)
                      DownNode  = &(F->GetValue(i,j-1));
                    else
                      DownNode  = NULL;

                    if(i>0)
                      LeftNode  = &(F->GetValue(i-1,j));
                    else
                      LeftNode  = NULL;

                    if(i < F->GetX()-1)
                      RightNode = &(F->GetValue(i+1,j));
                    else
                      RightNode = NULL;

                if ( CurrentNode->isCond2D(CT_WALL_LAW_2D) || 
                     CurrentNode->isCond2D(CT_WALL_NO_SLIP_2D)) { 
                    FP lam_eff = 0;
                    int num_near_nodes = 0;
                    
                    if( CurrentNode->isCleanSources ) 
                        CurrentNode->SrcAdd[i2d_RoE] = 0.; 
                    
                    
                    if ( DownNode && DownNode->isCond2D(CT_SOLID_2D) ) {
                        
                        num_near_nodes = 1;
                        lam_eff = CurrentNode->lam+CurrentNode->lam_t; 
                        
                        if (CurrentNode->UpNode) {
                            lam_eff +=CurrentNode->UpNode->lam + CurrentNode->UpNode->lam_t;
                            num_near_nodes++;
                        }
                        
                        lam_eff = lam_eff/num_near_nodes;
                        
                        if(DownNode->Q_conv > 0.)
                           DownNode->Q_conv = (DownNode->Q_conv - lam_eff*(DownNode->Tg - CurrentNode->Tg)/dy_local)*0.5;
                        else
                           DownNode->Q_conv = -lam_eff*(DownNode->Tg - CurrentNode->Tg)/dy_local;
                        
                        CurrentNode->SrcAdd[i2d_RoE] = -dt*DownNode->Q_conv/dy;
                        //CurrentNode->SrcAdd[i2d_RoE] += dt*lam_eff*(DownNode->Tg - CurrentNode->Tg)/dy2;
                    }
                    
                    if ( UpNode && UpNode->isCond2D(CT_SOLID_2D) ) {
                        
                        num_near_nodes = 1;
                        lam_eff = CurrentNode->lam+CurrentNode->lam_t; 
                        
                        if (CurrentNode->DownNode) {
                            lam_eff +=CurrentNode->DownNode->lam + CurrentNode->DownNode->lam_t;
                            num_near_nodes++;
                        }
                        
                        lam_eff = lam_eff/num_near_nodes;
                        
                        if(UpNode->Q_conv > 0.)
                           UpNode->Q_conv = (UpNode->Q_conv - lam_eff*(UpNode->Tg - CurrentNode->Tg)/dy_local)*0.5;
                        else
                           UpNode->Q_conv = -lam_eff*(UpNode->Tg - CurrentNode->Tg)/dy_local;
                        
                        CurrentNode->SrcAdd[i2d_RoE] = -dt*UpNode->Q_conv/dy;
                        //CurrentNode->SrcAdd[i2d_RoE] += dt*lam_eff*(UpNode->Tg - CurrentNode->Tg)/dy2;
                    }
                    
                    if ( LeftNode && LeftNode->isCond2D(CT_SOLID_2D) ) {
                        
                        num_near_nodes = 1;
                        lam_eff = CurrentNode->lam+CurrentNode->lam_t; 
                        
                        if (CurrentNode->RightNode) {
                            lam_eff +=CurrentNode->RightNode->lam + CurrentNode->RightNode->lam_t;
                            num_near_nodes++;
                        }
                        
                        lam_eff = lam_eff/num_near_nodes;
                        
                        if(LeftNode->Q_conv > 0.)
                           LeftNode->Q_conv = (LeftNode->Q_conv - lam_eff*(LeftNode->Tg - CurrentNode->Tg)/dx_local)*0.5;
                        else
                           LeftNode->Q_conv = -lam_eff*(LeftNode->Tg - CurrentNode->Tg)/dx_local;
                        
                        CurrentNode->SrcAdd[i2d_RoE] = -dt*LeftNode->Q_conv/dx;
                        //CurrentNode->SrcAdd[i2d_RoE] += dt*lam_eff*(LeftNode->Tg - CurrentNode->Tg)/dx2;
                    }
                    
                    if ( RightNode && RightNode->isCond2D(CT_SOLID_2D) ) {
                        
                        num_near_nodes = 1;
                        lam_eff = CurrentNode->lam+CurrentNode->lam_t; 
                        
                        if (CurrentNode->LeftNode) {
                            lam_eff +=CurrentNode->LeftNode->lam + CurrentNode->LeftNode->lam_t;
                            num_near_nodes++;
                        }
                        
                        lam_eff = lam_eff/num_near_nodes;
                        
                        if(RightNode->Q_conv > 0.)
                           RightNode->Q_conv = (RightNode->Q_conv - lam_eff*(RightNode->Tg - CurrentNode->Tg)/dx_local)*0.5;
                        else
                           RightNode->Q_conv = -lam_eff*(RightNode->Tg - CurrentNode->Tg)/dx_local;
                        
                        CurrentNode->SrcAdd[i2d_RoE] = -dt*RightNode->Q_conv/dx;
                        //CurrentNode->SrcAdd[i2d_RoE] += dt*lam_eff*(RightNode->Tg - CurrentNode->Tg)/dx2;
                    }
                }
            }
         }
}


#ifdef _CUDA_
 __host__ __device__ 
#endif //_CUDA_ 
inline int SetTurbulenceModel(FlowNode2D<FP,NUM_COMPONENTS>* pJ) {
int      AddEq = 2;
       if (pJ->isTurbulenceCond2D(TCT_Prandtl_Model_2D)) {
           AddEq = 2;
       } else if ( pJ->isTurbulenceCond2D(TCT_k_eps_Model_2D)) {
           AddEq = 0;
       } else if ( pJ->isTurbulenceCond2D(TCT_Spalart_Allmaras_Model_2D)) {
           AddEq = 1;
       } else {
           AddEq = 2;
       }
 return AddEq;
}



void SetMinDistanceToWall2D(ComputationalMatrix2D* pJ2D,
                            UArray< XY<int> >* WallNodes2D
                            ,FP x0
                            ) {

FP  min_l_min = min((FlowNode2D<FP,NUM_COMPONENTS>::dx),
                        (FlowNode2D<FP,NUM_COMPONENTS>::dy));
#ifdef _OPEN_MP
#pragma omp parallel  for
#endif //_OPEN_MP

for (int i=0;i<(int)pJ2D->GetX();i++ ) {
    for (int j=0;j<(int)pJ2D->GetY();j++ ) {
           if (pJ2D->GetValue(i,j).isCond2D(CT_NODE_IS_SET_2D) &&
              !pJ2D->GetValue(i,j).isCond2D(CT_SOLID_2D)) {
              if(pJ2D->GetValue(i,j).Tg != 0 && pJ2D->GetValue(i,j).p == 0.) {
                 pJ2D->GetValue(i,j).SetCond2D(CT_SOLID_2D);
             } else {   
                        unsigned int iw,jw;
                        FP wx, wy;
                        FP x, y;

                        pJ2D->GetValue(i,j).l_min = max((x0+FlowNode2D<FP,NUM_COMPONENTS>::dx*pJ2D->GetX()),
                                                           (FlowNode2D<FP,NUM_COMPONENTS>::dy*pJ2D->GetY()));
                        
                        x = x0 + i * FlowNode2D<FP,NUM_COMPONENTS>::dx;
                        y = j * FlowNode2D<FP,NUM_COMPONENTS>::dy;

                        for (int ii=0;ii<(int)WallNodes2D->GetNumElements();ii++) {

                             iw = WallNodes2D->GetElementPtr(ii)->GetX();
                             jw = WallNodes2D->GetElementPtr(ii)->GetY();

                             wx = iw * FlowNode2D<FP,NUM_COMPONENTS>::dx;
                             wy = jw * FlowNode2D<FP,NUM_COMPONENTS>::dy;

                             pJ2D->GetValue(i,j).l_min = min(pJ2D->GetValue(i,j).l_min,
                                                             sqrt((x-wx)*(x-wx) + (y-wy)*(y-wy)));
                             pJ2D->GetValue(i,j).l_min =  max(min_l_min,pJ2D->GetValue(i,j).l_min);
                        }
             }
         }
      }
   }
}

#ifdef _IMPI_
void LongSend(int rank, void* src,  size_t len, int data_tag) {
const size_t fix_buf_size=1024*1024; // 1M
char  TmpBuff[fix_buf_size+1];
int   num;
off_t offset=0L;

      if(len <= fix_buf_size) {
                MPI::COMM_WORLD.Send(src,
                                     len,
                                     MPI::BYTE,rank,data_tag);
      } else {
       num = len/fix_buf_size;
       for(int i=0;i<num;i++) {
                memcpy(TmpBuff,src+offset,fix_buf_size);
                MPI::COMM_WORLD.Send(TmpBuff,
                                     fix_buf_size,
                                     MPI::BYTE,rank,data_tag);
                offset += fix_buf_size;
       }

       if(len - num*fix_buf_size > 0L) {
          memcpy(TmpBuff,src+offset,len - num*fix_buf_size);
          MPI::COMM_WORLD.Send(TmpBuff,
                              len - num*fix_buf_size,
                              MPI::BYTE,rank,data_tag);
       }
     }
}

void LongRecv(int rank, void* dst,  size_t len, int data_tag) {
const size_t fix_buf_size=1024*1024; // 1M
char  TmpBuff[fix_buf_size+1];
int   num;
off_t offset=0L;
      if(len <= fix_buf_size) {
            MPI::COMM_WORLD.Recv(dst,
                                 len,
                                 MPI::BYTE,rank,data_tag);
      } else {
       num = len/fix_buf_size;
       for(int i=0;i<num;i++) {
                MPI::COMM_WORLD.Recv(TmpBuff,
                                     fix_buf_size,
                                     MPI::BYTE,rank,data_tag);
                memcpy(dst+offset,TmpBuff,fix_buf_size);
                offset += fix_buf_size;
       }

       if(len - num*fix_buf_size > 0L) {
          MPI::COMM_WORLD.Recv(TmpBuff,
                              len - num*fix_buf_size,
                              MPI::BYTE,rank,data_tag);
          memcpy(dst+offset, TmpBuff, len - num*fix_buf_size);
       }
      }
}

void LongMatrixSend(int rank, void* src,  size_t len) {
const size_t fix_buf_size=1024*1024; // 1M
char  TmpBuff[fix_buf_size+1];
int   num;
off_t offset=0L;

      if(len <= fix_buf_size) {
                MPI::COMM_WORLD.Send(src,
                                     len,
                                     MPI::BYTE,rank,tag_Matrix);
      } else {
       num = len/fix_buf_size;
       for(int i=0;i<num;i++) {
                memcpy(TmpBuff,src+offset,fix_buf_size);
                MPI::COMM_WORLD.Send(TmpBuff,
                                     fix_buf_size,
                                     MPI::BYTE,rank,tag_Matrix);
                offset += fix_buf_size;
       }

       if(len - num*fix_buf_size > 0L) {
          memcpy(TmpBuff,src+offset,len - num*fix_buf_size);
          MPI::COMM_WORLD.Send(TmpBuff,
                              len - num*fix_buf_size,
                              MPI::BYTE,rank,tag_Matrix);
       }
     }
}

void LongMatrixRecv(int rank, void* dst,  size_t len) {
const size_t fix_buf_size=1024*1024; // 1M
char  TmpBuff[fix_buf_size+1];
int   num;
off_t offset=0L;
      if(len <= fix_buf_size) {
            MPI::COMM_WORLD.Recv(dst,
                                 len,
                                 MPI::BYTE,rank,tag_Matrix);
      } else {
       num = len/fix_buf_size;
       for(int i=0;i<num;i++) {
                MPI::COMM_WORLD.Recv(TmpBuff,
                                     fix_buf_size,
                                     MPI::BYTE,rank,tag_Matrix);
                memcpy(dst+offset,TmpBuff,fix_buf_size);
                offset += fix_buf_size;
       }

       if(len - num*fix_buf_size > 0L) {
          MPI::COMM_WORLD.Recv(TmpBuff,
                              len - num*fix_buf_size,
                              MPI::BYTE,rank,tag_Matrix);
          memcpy(dst+offset, TmpBuff, len - num*fix_buf_size);
       }
      }
}
#endif // _IMPI_



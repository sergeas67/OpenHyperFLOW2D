/*******************************************************************************
*   OpenHyperFLOW2D                                                            *
*                                                                              *
*   Transient, Density based Effective Explicit Parallel Solver (TDEEPS2D)     *
*                                                                              *
*   Version  1.0.3                                                             *
*   Copyright (C)  1995-2016 by Serge A. Suchkov                               *
*   Copyright policy: LGPL V3                                                  *
*   http://github.com/sergeas67/openhyperflow2d                                *
*   deeps2d_core.cpp: OpenHyperFLOW2D solver core code....                     *
*                                                                              *
*  last update: 07/04/2016                                                     *
********************************************************************************/
#include "deeps2d_core.hpp"

#include <sys/time.h>
#include <sys/timeb.h>
#include <sys/file.h>

int NumContour;
int start_iter = 5;
int n_s;
int PreloadFlag = 0,p_g=0;

SourceList2D*  SrcList = NULL;
int            isGasSource=0;
int            TurbStartIter;
int            TurbExtModel;
int            err_i, err_j;
int            turb_mod_name_index = 0;
int            isAlternateRMS;
int            isIgnoreUnsetNodes;
FP             Ts0,A,W,Mach;

UArray< XY<int> >* GlobalSubDomain;
UArray< XY<int> >* WallNodes;

UArray<UMatrix2D< FlowNode2D<FP,NUM_COMPONENTS> >*>*     SubDomainArray;
UArray<UMatrix2D< FlowNodeCore2D<FP,NUM_COMPONENTS> >*>* CoreSubDomainArray;
UArray<XCut>*                                            XCutArray;

int             Nstep;
FP              ExitMonitorValue;
int             MonitorIndex;
int             MonitorCondition; // 0 - equal
                                  // 1 - less than
                                  // 2 - great than
UArray< MonitorPoint >*   MonitorPointsArray = NULL;

unsigned int     AddSrcStartIter = 0;
FP      beta0;
FP      CFL;
Table*  CFL_Scenario;
Table*  beta_Scenario;
FP      chord = 0;
#ifdef _OPEN_MP
FP      DD_max_var;
#endif // OPEN_MP

//------------------------------------------
// Cx,Cy,Cd,Cv,Cp,St
//------------------------------------------
int     is_Cx_calc,is_Cd_calc;
FP      p_ambient;
FP      x0_body,y0_body,dx_body,dy_body;
FP      x0_nozzle,y0_nozzle,dy_nozzle;

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
BlendingFactorFunction     bFF;

FP*     Y=NULL;
FP      Cp=0.;
FP      lam=0.;
FP      mu=0.;
FP      Tg=0.;
FP      Rg=0.;
FP      Pg=0.;
FP      Wg=0.;
FP      Ug,Vg;
int     CompIndex;

FP      Y_fuel[4]={1.,0.,0.,0.};  /* fuel */
FP      Y_ox[4]  ={0.,1.,0.,0.};  /* OX */
FP      Y_cp[4]  ={0.,0.,1.,0.};  /* cp */
FP      Y_air[4] ={0.,0.,0.,1.};  /* air */
FP      Y_mix[4] ={0.,0.,0.,0.};  /* mixture */
const  char*  RMS_Name[11] = {"Rho","Rho*U","Rho*V","Rho*E","Rho*Yfu","Rho*Yox","Rho*Ycp","Rho*k","Rho*eps","Rho*omega","nu_t"};
int    useSwapFile=0;
char   RMSFileName[255];
char   MonitorsFileName[255];
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
ofstream*                                    pInputData;           // Output data stream
ofstream*                                    pRMS_OutFile;         // Output RMS stream
ofstream*                                    pMonitors_OutFile;    // Output Monitors stream
ofstream*                                    pHeatFlux_OutFile;    // Output HeatFlux stream
unsigned int                                 iter = 0;             // iteration number
unsigned int                                 last_iter=0;          // Global iteration number
int                                          isStop=0;             // Stop flag
InputData*                                   Data=NULL;            // Object data loader
UMatrix2D< FlowNode2D<FP,NUM_COMPONENTS> >*  J=NULL;               // Main computation area
UArray<Flow*>*                               FlowList=NULL;        // List of 'Flow' objects
UArray<Flow2D*>*                             Flow2DList=NULL;      // List of 'Flow2D' objects
UArray<Bound2D*>                             SingleBoundsList;     // Single Bounds List;

FP                                           dt;                   // time step
FP*                                          RoUx=NULL;
FP*                                          RoVy=NULL;
FP                                           GlobalTime=0.;
FP                                           CurrentTimePart=0;

ofstream*                                    pOutputData;          // output data stream (file)

int                                          I,NSaveStep;
unsigned int                                 MaxX=0;               // X dimension of computation area
unsigned int                                 MaxY=0;               // Y dimension of computation area

FP                                           dxdy,dx2,dy2;
FP                                           Gig,SigW,SigF,delta_bl;
FP                                           nrbc_beta0;

unsigned long                                FileSizeGas           = 0;
int                                          isVerboseOutput       = 0;
int                                          isTurbulenceReset     = 0;

void InitSharedData(InputData* _data,
                    void* CRM_data
#ifdef _MPI
                    ,int rank
#endif //_MPI
                    ) {
            
    ChemicalReactionsModelData2D* model_data = (ChemicalReactionsModelData2D*)CRM_data;

            fd_g=0;
            GasSwapData         = NULL;
            TurbMod             = 0;

            isVerboseOutput = _data->GetIntVal((char*)"isVerboseOutput");
            if ( _data->GetDataError()==-1 ) {
                Abort_OpenHyperFLOW2D();
            }

            bFF   = (BlendingFactorFunction)_data->GetIntVal((char*)"BFF");  // Blending factor function type
            if ( _data->GetDataError()==-1 ) {
                Abort_OpenHyperFLOW2D();
            }

            MaxX      = Data->GetIntVal((char*)"MaxX");          // X dimension of computation area
            if ( Data->GetDataError()==-1 ) {
                Abort_OpenHyperFLOW2D();
            }

            MaxY      = Data->GetIntVal((char*)"MaxY");          // Y dimension of computation area
            if ( Data->GetDataError()==-1 ) {
                Abort_OpenHyperFLOW2D();
            }

            dx        = Data->GetFloatVal((char*)"dx");          // x size of cell 
            if ( Data->GetDataError()==-1 ) {
                Abort_OpenHyperFLOW2D();
            }

            dy        = Data->GetFloatVal((char*)"dy");          // y size of sell
            if ( Data->GetDataError()==-1 ) {
                Abort_OpenHyperFLOW2D();
            }
#ifdef _MPI
            if(rank==0) {
#endif //_MPI
                *(_data->GetMessageStream()) << "X=" << MaxX << "  Y=" << MaxY << "  dx=" << dx << "  dy=" << dy <<"\n" << flush;
                _data->GetMessageStream()->flush();
#ifdef _MPI
            }
#endif //_MPI
            SigW     = _data->GetFloatVal((char*)"SigW");
            if ( _data->GetDataError()==-1 ) {
                Abort_OpenHyperFLOW2D();
            }

            SigF     = _data->GetFloatVal((char*)"SigF");
            if ( _data->GetDataError()==-1 ) {
                Abort_OpenHyperFLOW2D();
            }

            delta_bl  = _data->GetFloatVal((char*)"delta_bl");
            if ( _data->GetDataError()==-1 ) {
                Abort_OpenHyperFLOW2D();
            }

            TurbMod  = _data->GetIntVal((char*)"TurbulenceModel");
            if ( _data->GetDataError()==-1 ) {
                 Abort_OpenHyperFLOW2D();
            }

            TurbStartIter  = _data->GetIntVal((char*)"TurbStartIter");
            if ( _data->GetDataError()==-1 ) {
                 Abort_OpenHyperFLOW2D();
            }

            TurbExtModel = _data->GetIntVal((char*)"TurbExtModel");
            if ( _data->GetDataError()==-1 ) {
                 Abort_OpenHyperFLOW2D();
            }

            isTurbulenceReset = _data->GetIntVal((char*)"isTurbulenceReset");
            if ( _data->GetDataError()==-1 ) {
                 Abort_OpenHyperFLOW2D();
            }

            FlowNode2D<FP,NUM_COMPONENTS>::FT = (FlowType)(_data->GetIntVal((char*)"FlowType"));
            if ( _data->GetDataError()==-1 ) {
                Abort_OpenHyperFLOW2D();
            }
            
            ProblemType = (SolverMode)(_data->GetIntVal((char*)"ProblemType"));
            if ( _data->GetDataError()==-1 ) {
                Abort_OpenHyperFLOW2D();
            }

            CFL  = _data->GetFloatVal((char*)"CFL");               // Courant number
            if ( _data->GetDataError()==-1 ) {
                Abort_OpenHyperFLOW2D();
            }

            CFL_Scenario  = _data->GetTable((char*)"CFL_Scenario"); // Courant number scenario
            if ( _data->GetDataError()==-1 ) {
                Abort_OpenHyperFLOW2D();
            }

            NSaveStep = _data->GetIntVal((char*)"NSaveStep");
            if ( _data->GetDataError()==-1 ) Abort_OpenHyperFLOW2D();

            Nstep           = _data->GetIntVal((char*)"Nmax");     // Sync computation area each after NMax iterations  
            if ( _data->GetDataError()==-1 ) {
                Abort_OpenHyperFLOW2D();
            }
            NOutStep  = _data->GetIntVal((char*)"NOutStep");      // Output step
            if ( NOutStep<=0 )     NOutStep=1;

            if(NOutStep >= Nstep)
               Nstep = NOutStep + 1;

            isAlternateRMS = _data->GetIntVal((char*)"isAlternateRMS");
            if ( _data->GetDataError()==-1 ) {
                Abort_OpenHyperFLOW2D();
            }

            isIgnoreUnsetNodes = _data->GetIntVal((char*)"isIgnoreUnsetNodes");
            if ( _data->GetDataError()==-1 ) {
                Abort_OpenHyperFLOW2D();
            }

            MonitorIndex = _data->GetIntVal((char*)"MonitorIndex");
            if ( _data->GetDataError()==-1 ) {
                Abort_OpenHyperFLOW2D();
            }

            if(MonitorIndex > 5 ||
               MonitorIndex < 0) {
               MonitorIndex = 0;
            }

            ExitMonitorValue  = _data->GetFloatVal((char*)"ExitMonitorValue");   // Monitor value for exit
            if ( _data->GetDataError()==-1 ) {
                Abort_OpenHyperFLOW2D();
            }


            int NumMonitorPoints = _data->GetIntVal((char*)"NumMonitorPoints");

            if ( _data->GetDataError()==-1 ) {
                Abort_OpenHyperFLOW2D();
            }

            if (NumMonitorPoints > 0 ) {
#ifdef _MPI
            if(rank==0) {
#endif //_MPI
                *(_data->GetMessageStream()) << "Read monitor points...\n" << flush;
#ifdef _MPI
            }
#endif //_MPI
                MonitorPointsArray = new UArray< MonitorPoint >(); 
                MonitorPoint TmpPoint;
                for (int i=0;i<NumMonitorPoints;i++) {
                   char   point_str[256];
                   FP     point_xy;
                   snprintf(point_str,256,"Point-%i.X",i+1);
                   point_xy =  _data->GetFloatVal(point_str);
                   if ( _data->GetDataError()==-1 ) {
                       Abort_OpenHyperFLOW2D();
                   }
                   TmpPoint.MonitorXY.SetX(point_xy);
                   snprintf(point_str,256,"Point-%i.Y",i+1);
                   point_xy =  _data->GetFloatVal(point_str);
                   if ( _data->GetDataError()==-1 ) {
                       Abort_OpenHyperFLOW2D();
                   }
                   
                   TmpPoint.MonitorXY.SetY(point_xy);

                   if(TmpPoint.MonitorXY.GetX() < 0.0 ||
                      TmpPoint.MonitorXY.GetY() < 0.0 ||
                      TmpPoint.MonitorXY.GetX() > MaxX*dx ||
                      TmpPoint.MonitorXY.GetY() > MaxY*dy ) {
                      *(_data->GetMessageStream()) << "Point no " << i+1 << " X=" << TmpPoint.MonitorXY.GetX() << "m Y=" << TmpPoint.MonitorXY.GetY() << " m out of domain...monitor ignored." << endl;
                   } else {
                       MonitorPointsArray->AddElement(&TmpPoint);
#ifdef _MPI
                  if(rank==0) {
#endif //_MPI
                    *(_data->GetMessageStream()) << "Point no " << i+1 << " X=" << TmpPoint.MonitorXY.GetX() << "m Y=" << TmpPoint.MonitorXY.GetY() << " m" << endl;
#ifdef _MPI
                  }
#endif //_MPI
                }
               }
#ifdef _MPI

           if(rank==0) {
#endif //_MPI
                *(_data->GetMessageStream()) << "Load " << NumMonitorPoints << " monitor points...OK" << endl;
#ifdef _MPI
             }
#endif //_MPI
            
            }

            beta0 = _data->GetFloatVal((char*)"beta");       // Base blending factor.
            if ( _data->GetDataError()==-1 ) {
                Abort_OpenHyperFLOW2D();
            }
            
            nrbc_beta0 = _data->GetFloatVal((char*)"beta_NonReflectedBC");  // Base blending factor for non-reflected BC.
            if ( _data->GetDataError()==-1 ) {
                Abort_OpenHyperFLOW2D();
            }
            
            beta_Scenario  = _data->GetTable((char*)"beta_Scenario"); //  Base blending factor scenario
            if ( _data->GetDataError()==-1 ) {
                Abort_OpenHyperFLOW2D();
            }

            model_data->K0     = _data->GetFloatVal((char*)"K0");    // Stoichiometric ratio (OX/Fuel)
            if ( _data->GetDataError()==-1 ) {
                Abort_OpenHyperFLOW2D();
            }

            model_data->gamma  = _data->GetFloatVal((char*)"gamma"); // Factor completion of a chemical reaction
            if ( _data->GetDataError()==-1 ) {
                Abort_OpenHyperFLOW2D();
            }

            model_data->Tf     = _data->GetFloatVal((char*)"Tf");    // Ignition temperature
            if ( _data->GetDataError()==-1 ) {
                Abort_OpenHyperFLOW2D();
            }

            isAdiabaticWall     = _data->GetIntVal((char*)"isAdiabaticWall");  // is walls adiabatic ?
            if ( _data->GetDataError()==-1 ) {
                Abort_OpenHyperFLOW2D();
            }

/* Load combustion products properties */

            model_data->R_cp   = _data->GetFloatVal((char*)"R_cp");
            if ( _data->GetDataError()==-1 ) {
                Abort_OpenHyperFLOW2D();
            }
            model_data->H_cp   = _data->GetFloatVal((char*)"H_cp");
            if ( _data->GetDataError()==-1 ) {
                Abort_OpenHyperFLOW2D();
            }
            model_data->lam_cp  = _data->GetTable((char*)"lam_cp");
            if ( _data->GetDataError()==-1 ) {
                Abort_OpenHyperFLOW2D();
            }
            model_data->mu_cp   = _data->GetTable((char*)"mu_cp");
            if ( _data->GetDataError()==-1 ) {
                Abort_OpenHyperFLOW2D();
            }
            model_data->Cp_cp   = _data->GetTable((char*)"Cp_cp");
            if ( _data->GetDataError()==-1 ) {
                Abort_OpenHyperFLOW2D();
            }

/* Load fuel properties */

            model_data->R_Fuel = _data->GetFloatVal((char*)"R_Fuel");
            if ( _data->GetDataError()==-1 ) {
                Abort_OpenHyperFLOW2D();
            }
            model_data->H_Fuel = _data->GetFloatVal((char*)"H_Fuel");
            if ( _data->GetDataError()==-1 ) {
                Abort_OpenHyperFLOW2D();
            }
            model_data->lam_Fuel     = _data->GetTable((char*)"lam_Fuel");
            if ( _data->GetDataError()==-1 ) {
                Abort_OpenHyperFLOW2D();
            }
            model_data->mu_Fuel      = _data->GetTable((char*)"mu_Fuel");
            if ( _data->GetDataError()==-1 ) {
                Abort_OpenHyperFLOW2D();
            }
            model_data->Cp_Fuel      = _data->GetTable((char*)"Cp_Fuel");
            if ( _data->GetDataError()==-1 ) {
                Abort_OpenHyperFLOW2D();
            }

/* Load OX properties */

            model_data->R_OX   = _data->GetFloatVal((char*)"R_OX");
            if ( _data->GetDataError()==-1 ) {
                Abort_OpenHyperFLOW2D();
          }
            model_data->H_OX   = _data->GetFloatVal((char*)"H_OX");
            if ( _data->GetDataError()==-1 ) {
                Abort_OpenHyperFLOW2D();
            }
            model_data->lam_OX       = _data->GetTable((char*)"lam_OX");
            if ( _data->GetDataError()==-1 ) {
                Abort_OpenHyperFLOW2D();
            }
            model_data->mu_OX        = _data->GetTable((char*)"mu_OX");
            if ( _data->GetDataError()==-1 ) {
                Abort_OpenHyperFLOW2D();
            }
            model_data->Cp_OX        = _data->GetTable((char*)"Cp_OX");
            if ( _data->GetDataError()==-1 ) {
                Abort_OpenHyperFLOW2D();
            }

/* Load air properties */

            model_data->R_air   = _data->GetFloatVal((char*)"R_air");
            if ( _data->GetDataError()==-1 ) {
                Abort_OpenHyperFLOW2D();
            }
            model_data->H_air   = _data->GetFloatVal((char*)"H_air");
            if ( _data->GetDataError()==-1 ) {
                Abort_OpenHyperFLOW2D();
            }
            model_data->lam_air       = _data->GetTable((char*)"lam_air");
            if ( _data->GetDataError()==-1 ) {
                Abort_OpenHyperFLOW2D();
            }
            model_data->mu_air        = _data->GetTable((char*)"mu_air");
            if ( _data->GetDataError()==-1 ) {
                Abort_OpenHyperFLOW2D();
            }
            model_data->Cp_air        = _data->GetTable((char*)"Cp_air");
            if ( _data->GetDataError()==-1 ) {
                Abort_OpenHyperFLOW2D();
            }

#ifdef _UNIFORM_MESH_
            FlowNode2D<FP,NUM_COMPONENTS>::dx        = dx;
            FlowNode2D<FP,NUM_COMPONENTS>::dy        = dy;
#endif //_UNIFORM_MESH_
            FlowNode2D<FP,NUM_COMPONENTS>::Hu[h_fu]  = model_data->H_Fuel;
            FlowNode2D<FP,NUM_COMPONENTS>::Hu[h_ox]  = model_data->H_OX;
            FlowNode2D<FP,NUM_COMPONENTS>::Hu[h_cp]  = model_data->H_cp;
            FlowNode2D<FP,NUM_COMPONENTS>::Hu[h_air] = model_data->H_air;
};

#ifdef _MPI
struct Var_pack {
       FP   dt_min;
#ifdef __INTEL_COMPILER
       DD_pack  DD[4+NUM_COMPONENTS+2];
#else 
       DD_pack  DD[FlowNode2D<FP,NUM_COMPONENTS>::NumEq];
#endif //__INTEL_COMPILER
};
#endif // _MPI

void DEEPS2D_Run(ofstream* f_stream
#ifdef _MPI
                ,UMatrix2D< FlowNode2D<FP,NUM_COMPONENTS> >*     pJ,
                 UMatrix2D< FlowNodeCore2D<FP,NUM_COMPONENTS> >* pC,
                 int rank, int last_rank, FP x0
#endif //_MPI
                 ) {

   // local variables
#ifdef __ICC
    __declspec(align(_ALIGN)) FlowNode2D< FP,NUM_COMPONENTS >* CurrentNode=NULL;
    __declspec(align(_ALIGN)) FlowNodeCore2D< FP,NUM_COMPONENTS >* NextNode=NULL;
    __declspec(align(_ALIGN)) FlowNode2D< FP,NUM_COMPONENTS >* UpNode=NULL;          // neast
    __declspec(align(_ALIGN)) FlowNode2D< FP,NUM_COMPONENTS >* DownNode=NULL;        // nodes
    __declspec(align(_ALIGN)) FlowNode2D< FP,NUM_COMPONENTS >* LeftNode=NULL;        // references
    __declspec(align(_ALIGN)) FlowNode2D< FP,NUM_COMPONENTS >* RightNode=NULL;
    __declspec(align(_ALIGN)) FP   beta; 
    __declspec(align(_ALIGN)) FP   _beta;
    __declspec(align(_ALIGN)) FP   dXX;
    __declspec(align(_ALIGN)) FP   dYY;
    __declspec(align(_ALIGN)) FP   AAA;
    __declspec(align(_ALIGN)) FP   dtdx;
    __declspec(align(_ALIGN)) FP   dtdy;
    __declspec(align(_ALIGN)) FP   dyy;
    __declspec(align(_ALIGN)) FP   dxx;
    __declspec(align(_ALIGN)) FP   dx_1;       // 1/dx
    __declspec(align(_ALIGN)) FP   dy_1;       // 1/dy
    __declspec(align(_ALIGN)) FP   dx_1_n_n_1; // 1/(2*dx)
    __declspec(align(_ALIGN)) FP   dy_1_m_m_1; // 1/(2*dy)
    __declspec(align(_ALIGN)) FP   DD_local[FlowNode2D<FP,NUM_COMPONENTS>::NumEq];
    __declspec(align(_ALIGN)) FP   DD[FlowNode2D<FP,NUM_COMPONENTS>::NumEq];
    __declspec(align(_ALIGN)) FP   n_n;
    __declspec(align(_ALIGN)) FP   m_m;
    __declspec(align(_ALIGN)) FP   n_n_1;
    __declspec(align(_ALIGN)) FP   m_m_1;
    __declspec(align(_ALIGN)) FP   beta_Scenario_Val; 
    __declspec(align(_ALIGN)) FP   CFL_Scenario_Val;  
#else
    FlowNode2D< FP,NUM_COMPONENTS >* CurrentNode  __attribute__ ((aligned (_ALIGN)))=NULL;
    FlowNodeCore2D< FP,NUM_COMPONENTS >* NextNode __attribute__ ((aligned (_ALIGN)))=NULL;
    FlowNode2D< FP,NUM_COMPONENTS >* UpNode       __attribute__ ((aligned (_ALIGN)))=NULL;        // neast     
    FlowNode2D< FP,NUM_COMPONENTS >* DownNode     __attribute__ ((aligned (_ALIGN)))=NULL;        // nodes     
    FlowNode2D< FP,NUM_COMPONENTS >* LeftNode     __attribute__ ((aligned (_ALIGN)))=NULL;        // references
    FlowNode2D< FP,NUM_COMPONENTS >* RightNode    __attribute__ ((aligned (_ALIGN)))=NULL;                    
    FP   beta   __attribute__ ((aligned (_ALIGN)));  
    FP   _beta  __attribute__ ((aligned (_ALIGN)));
    FP   dXX    __attribute__ ((aligned (_ALIGN))); 
    FP   dYY    __attribute__ ((aligned (_ALIGN))); 
    FP   AAA    __attribute__ ((aligned (_ALIGN))); 
    FP   dtdx   __attribute__ ((aligned (_ALIGN)));
    FP   dtdy   __attribute__ ((aligned (_ALIGN)));
    FP   dyy    __attribute__ ((aligned (_ALIGN)));
    FP   dxx    __attribute__ ((aligned (_ALIGN)));
    FP   dx_1   __attribute__ ((aligned (_ALIGN)));
    FP   dy_1   __attribute__ ((aligned (_ALIGN)));
    FP   dx_1_n_n_1 __attribute__ ((aligned (_ALIGN))); // 1/(2*dx)
    FP   dy_1_m_m_1 __attribute__ ((aligned (_ALIGN))); // 1/(2*dy)
    FP   DD_local[FlowNode2D<FP,NUM_COMPONENTS>::NumEq] __attribute__ ((aligned (_ALIGN)));
    FP   DD[FlowNode2D<FP,NUM_COMPONENTS>::NumEq] __attribute__ ((aligned (_ALIGN)));
    FP   n_n   __attribute__ ((aligned (_ALIGN)));
    FP   m_m   __attribute__ ((aligned (_ALIGN)));
    FP   n_n_1 __attribute__ ((aligned (_ALIGN)));
    FP   m_m_1 __attribute__ ((aligned (_ALIGN)));
    FP   beta_Scenario_Val __attribute__ ((aligned (_ALIGN)));
    FP   CFL_Scenario_Val  __attribute__ ((aligned (_ALIGN)));
#endif //__ICC
    unsigned int Num_Eq,n1,n2,n3,n4;
    unsigned int j,k;
    unsigned int StartXLocal,MaxXLocal;
    FP   d_time;
    FP   t,VCOMP;
    timeval  start, stop, mark1, mark2;
    int  N1,N2,N3,N4;
#ifdef _OPEN_MP
    unsigned int i_max,j_max,k_max;
#endif //_OPEN_MP
#ifdef _MPI
    unsigned int i_max,j_max,k_max;
#endif //_MPI

#ifdef __ICC
    __declspec(align(_ALIGN)) FP   dt_min_local;
    __declspec(align(_ALIGN)) FP   max_RMS;
#else
    FP   dt_min_local __attribute__ ((aligned (_ALIGN)));
    FP   max_RMS __attribute__ ((aligned (_ALIGN)));
#endif // __ICC
    
    int  k_max_RMS;
#ifndef _MPI
    UMatrix2D< FlowNode2D<FP,NUM_COMPONENTS> >*     pJ=NULL;
    UMatrix2D< FlowNodeCore2D<FP,NUM_COMPONENTS> >* pC=NULL;
    FP*   dt_min;
    int*  i_c;
    int*  j_c;

#ifdef __ICC
    __declspec(align(_ALIGN)) FP    dtmin;
    __declspec(align(_ALIGN)) FP    sum_RMS[FlowNode2D<FP,NUM_COMPONENTS>::NumEq];
    __declspec(align(_ALIGN)) FP    sum_RMS_Div[FlowNode2D<FP,NUM_COMPONENTS>::NumEq];
#else
    FP    dtmin __attribute__ ((aligned (_ALIGN)));
    FP    sum_RMS[FlowNode2D<FP,NUM_COMPONENTS>::NumEq] __attribute__ ((aligned (_ALIGN)));
    FP    sum_RMS_Div[FlowNode2D<FP,NUM_COMPONENTS>::NumEq] __attribute__ ((aligned (_ALIGN)));
#endif // __ICC
    
    unsigned long  sum_iRMS[FlowNode2D<FP,NUM_COMPONENTS>::NumEq];
    UMatrix2D<FP>  RMS(FlowNode2D<FP,NUM_COMPONENTS>::NumEq,SubDomainArray->GetNumElements());
    UMatrix2D<FP>  sumDiv(FlowNode2D<FP,NUM_COMPONENTS>::NumEq,SubDomainArray->GetNumElements());
    UMatrix2D<int> iRMS(FlowNode2D<FP,NUM_COMPONENTS>::NumEq,SubDomainArray->GetNumElements());
    UMatrix2D<FP>  DD_max(FlowNode2D<FP,NUM_COMPONENTS>::NumEq,SubDomainArray->GetNumElements());
#else

#ifdef __ICC
    __declspec(align(_ALIGN)) FP   RMS[FlowNode2D<FP,NUM_COMPONENTS>::NumEq];
    __declspec(align(_ALIGN)) FP   sumDiv[FlowNode2D<FP,NUM_COMPONENTS>::NumEq];

#else
    FP   RMS[FlowNode2D<FP,NUM_COMPONENTS>::NumEq] __attribute__ ((aligned (_ALIGN)));
    FP   sumDiv[FlowNode2D<FP,NUM_COMPONENTS>::NumEq] __attribute__ ((aligned (_ALIGN)));
#endif // __ICC
    
    unsigned long iRMS[FlowNode2D<FP,NUM_COMPONENTS>::NumEq];
    Var_pack DD_max[last_rank+1];
#ifdef _MPI_NB
    MPI::Request  HaloExchange[4];
    MPI::Request  DD_Exchange[2*last_rank];
#endif //_MPI_NB
    unsigned int r_Overlap, l_Overlap;
#endif // _MPI
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
    dt_min = new FP[SubDomainArray->GetNumElements()];
    i_c    = new int[SubDomainArray->GetNumElements()];
    j_c    = new int[SubDomainArray->GetNumElements()];
#endif // _MPI

#ifndef _MPI
    
    for(int ii=0;ii<(int)SubDomainArray->GetNumElements();ii++) {
        dt_min[ii] = dtmin = dt;
        i_c[ii] = j_c[ii] = 0;
     }
    snprintf(RMSFileName,255,"RMS-%s",OutFileName);
    CutFile(RMSFileName);
    pRMS_OutFile = OpenData(RMSFileName);
    SaveRMSHeader(pRMS_OutFile);
    if(MonitorPointsArray && MonitorPointsArray->GetNumElements() > 0) {
        snprintf(MonitorsFileName,255,"Monitors-%s",OutFileName);
        CutFile(MonitorsFileName);
        pMonitors_OutFile = OpenData(MonitorsFileName);
        SaveMonitorsHeader(pMonitors_OutFile, MonitorPointsArray);
    }
#else
    MPI::COMM_WORLD.Barrier();
    MPI::COMM_WORLD.Bcast(&dt,1,MPI::DOUBLE,0);

    if( rank == 0 ) {
        snprintf(RMSFileName,255,"RMS-%s",OutFileName);
        CutFile(RMSFileName);
        pRMS_OutFile = OpenData(RMSFileName);
        SaveRMSHeader(pRMS_OutFile);
        if(MonitorPointsArray && MonitorPointsArray->GetNumElements() > 0) {
            snprintf(MonitorsFileName,255,"Monitors-%s",OutFileName);
            CutFile(MonitorsFileName);
            pMonitors_OutFile = OpenData(MonitorsFileName);
            SaveMonitorsHeader(pMonitors_OutFile, MonitorPointsArray);
        }
    }

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
                    
                    dx_1 = 1.0/dx;
                    dy_1 = 1.0/dy;
                    
             do {
                  gettimeofday(&mark2,NULL);
                  gettimeofday(&start,NULL);
                  
                  if( AddSrcStartIter < iter + last_iter){
                    FlowNode2D<FP,NUM_COMPONENTS>::isSrcAdd = 1;
                  } else {
                    FlowNode2D<FP,NUM_COMPONENTS>::isSrcAdd = 0;
                  }
                 
                  iter = 0;
                  
                  do {
                      
                      beta_Scenario_Val = beta_Scenario->GetVal(iter+last_iter); 
                      CFL_Scenario_Val  = CFL_Scenario->GetVal(iter+last_iter); 
#ifdef _MPI
// MPI version
                      if(rank == 0 ) {
                          max_RMS =  0;
                      } else {
                        if(MonitorIndex < 5)
                           max_RMS = 0;
                        else
                           max_RMS = 0.5*ExitMonitorValue;
                      }

                  k_max_RMS = -1;

                  for (int kk=0;kk<FlowNode2D<FP,NUM_COMPONENTS>::NumEq;kk++ ) {
                      RMS[kk]    = 0.;                                         // Clean sum residual
                      iRMS[kk]   = 0;                                          // num involved nodes
                      sumDiv[kk] = 0.;
                      
                      DD_max[rank].DD[kk].RMS    = 0.;                          // sum residual per rank
                      DD_max[rank].DD[kk].iRMS   = 0;                           // num involved nodes per rank
                      DD_max[rank].DD[kk].sumDiv = 0;                           // sum Div param
                      DD_max[rank].DD[kk].DD     = 0.;                          // max residual per rank
                      DD_max[rank].DD[kk].i      = 0;                           // max residual x-coord
                      DD_max[rank].DD[kk].j      = 0;                           // max residual y-coord
                      DD_local[kk]               = 0.;                          // local residual
                    }
#else
                   n_s = (int)SubDomainArray->GetNumElements();
#ifdef _OPENMP
#pragma omp parallel shared(f_stream,CoreSubDomainArray, SubDomainArray, chemical_reactions,Y_mix,sum_RMS, sum_iRMS, \
                            Cp,i_max,j_max,k_max,Tg,beta0,CurrentTimePart,DD,dx,dy,MaxX,MaxY,dt_min,RMS,iRMS,DD_max,i_c,j_c,n_s) \
                     private(iter,j,k,n1,n2,n3,n4,N1,N2,N3,N4,n_n,m_m,pC,pJ,err_i,err_j,\
                             dXX,dYY,DD_local,AAA,StartXLocal,MaxXLocal,\
                             dtdx,dtdy,dt) reduction(min: dtmin)
//#pragma omp single
#endif //_OPENMP
                   {

                    if(MonitorIndex < 5)
                      max_RMS =  0.5*ExitMonitorValue;
                    else
                      max_RMS = 0.;

                    k_max_RMS = -1;

                    for ( k=0;k<(int)FlowNode2D<FP,NUM_COMPONENTS>::NumEq;k++ ) {
                           for(int ii=0;ii<n_s;ii++) {
                               RMS(k,ii)    = 0.;    // Clean sum residual
                               iRMS(k,ii)   = 0;     // num involved nodes
                               sumDiv(k,ii) = 0;     // sum div params 
                           }
                           sum_iRMS[k] = 0;
                           sum_RMS_Div[k] = sum_RMS[k] = DD[k] = 0.;
                       }
                   }
#endif //_MPI
#ifdef _OPENMP
#pragma omp for private(CurrentNode,NextNode,UpNode,DownNode,LeftNode,RightNode) ordered nowait 
                for(int ii=0;ii<n_s;ii++) {  // OpenMP version
#endif //_OPENMP

#ifndef _MPI
#ifndef _OPENMP
                    int ii = 0;                                                    // Single thread version
                    dt = dt_min[0];
#endif //_OPENMP
#else
                     dt  = 1.;
                     for (int ii=0;(int)ii<last_rank+1;ii++) {                  // Choose
                         dt  = min(DD_max[ii].dt_min,dt);                       // minimal
                    }

#endif // _MPI

#ifdef _OPENMP
                   dt    = dtmin;
                   
                   i_c[ii] = j_c[ii] = 0;
#endif //_OPENMP

#ifndef _MPI
                    pJ = SubDomainArray->GetElement(ii);
                    pC = CoreSubDomainArray->GetElement(ii);
#endif // _MPI

#ifdef _MPI
                    MPI::COMM_WORLD.Bcast(&dt,1,MPI::DOUBLE,0);
#else
                    for ( k=0;k<(int)FlowNode2D<FP,NUM_COMPONENTS>::NumEq;k++ ) {
                         DD_max(k,ii) =  0.;
                      }

                    if( ii == 0)
                       StartXLocal=0;
                    else
                       StartXLocal=1;
                    
                    if( ii == (int)SubDomainArray->GetNumElements()-1) {
                        MaxXLocal=pJ->GetX();
                    } else {
                        MaxXLocal=pJ->GetX()-1;
                    }
#endif // _MPI
                    
                    dtdx = dt/dx;
                    dtdy = dt/dy;

#ifdef _OPENMP
                    dt_min[ii]   = 1.;
#else
#ifdef _MPI
                    DD_max[rank].dt_min = dt_min_local = 1.;
#endif // _MPI
#endif // _OPENMP
                    for (int i = StartXLocal;i<(int)MaxXLocal;i++ ) {
                       for (int j=0;j<(int)MaxY;j++ ) {
                          
                          err_i = i;
                          err_j = j;

                          CurrentNode = &(pJ->GetValue(i,j)); 

                          if ( CurrentNode->isCond2D(CT_NODE_IS_SET_2D) &&
                              !CurrentNode->isCond2D(CT_SOLID_2D)  &&
                              !CurrentNode->isCond2D(NT_FC_2D)) {

                                NextNode    = &(pC->GetValue(i,j)); 

                                CurrentNode->time=GlobalTime;

                                n1=CurrentNode->idXl;
                                n2=CurrentNode->idXr;
                                n3=CurrentNode->idYu;
                                n4=CurrentNode->idYd;

                                N1 = i - n1;
                                N2 = i + n2;
                                N3 = j + n3;
                                N4 = j - n4;

                                n_n = max(n1+n2,1);
                                m_m = max(n3+n4,1);

                                n_n_1 = 1./n_n;
                                m_m_1 = 1./m_m;

                                UpNode    = &(pJ->GetValue(i,N3));
                                DownNode  = &(pJ->GetValue(i,N4));
                                RightNode = &(pJ->GetValue(N2,j));
                                LeftNode  = &(pJ->GetValue(N1,j));

                                Num_Eq =  (FlowNode2D<FP,NUM_COMPONENTS>::NumEq-SetTurbulenceModel(CurrentNode));
                                
                                // Scan equation system ... k - number of equation
                                for (int k=0;k<Num_Eq;k++ ) {
                                        int      c_flag = 0;
                                        int      dx_flag, dx2_flag;
                                        int      dy_flag, dy2_flag;
                                        
                                        beta  = CurrentNode->beta[k];
                                        _beta = 1. - beta;

                                // Precomputed variables for current node ...
                                        c_flag  = dx_flag = dy_flag = dx2_flag = dy2_flag = 0;
                                    if ( k < 4 ) { // Make bit flags for future test for current equation
                                        c_flag   = CT_Rho_CONST_2D     << k; 
                                        dx_flag  = CT_dRhodx_NULL_2D   << k;
                                        dy_flag  = CT_dRhody_NULL_2D   << k;
                                        dx2_flag = CT_d2Rhodx2_NULL_2D << k;
                                        dy2_flag = CT_d2Rhody2_NULL_2D << k;
                                    } else if (k < (4+NUM_COMPONENTS)) {
                                        c_flag   = CT_Y_CONST_2D;
                                        dx_flag  = CT_dYdx_NULL_2D;
                                        dy_flag  = CT_dYdy_NULL_2D;
                                        dx2_flag = CT_d2Ydx2_NULL_2D;
                                        dy2_flag = CT_d2Ydy2_NULL_2D;
                                    } else if (ProblemType == SM_NS && 
                                               (CurrentNode->isTurbulenceCond2D(TCT_k_eps_Model_2D) ||
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
                                    } else if(ProblemType == SM_NS && 
                                              (CurrentNode->isTurbulenceCond2D(TCT_k_eps_Model_2D) ||
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
                                            NextNode->S[k] = CurrentNode->S[k]*beta+_beta*(dxx*(LeftNode->S[k]+RightNode->S[k])+dyy*(UpNode->S[k]+DownNode->S[k]))*0.5
                                                           - (dtdx*dXX+dtdy*(dYY+CurrentNode->F[k]/(j+1))) + (CurrentNode->Src[k])*dt+CurrentNode->SrcAdd[k];
                                        } else {
                                            NextNode->S[k] = CurrentNode->S[k]*beta+_beta*(dxx*(LeftNode->S[k]+RightNode->S[k])+dyy*(UpNode->S[k]+DownNode->S[k]))*0.5
                                                           - (dtdx*dXX+dtdy*dYY) + (CurrentNode->Src[k])*dt+CurrentNode->SrcAdd[k];
                                        }
                               }
                         }
                     }
                  }
               }
                   
               for (int i=StartXLocal;i<(int)MaxXLocal;i++ ) {
                 for ( j=0;j<MaxY;j++ ) {

                          CurrentNode = &(pJ->GetValue(i,j)); 

                          if (CurrentNode->isCond2D(CT_NODE_IS_SET_2D) &&
                              !CurrentNode->isCond2D(CT_SOLID_2D) &&
                              !CurrentNode->isCond2D(NT_FC_2D)) {

                              NextNode    = &(pC->GetValue(i,j)); 
                              
                              n1=CurrentNode->idXl;
                              n2=CurrentNode->idXr;
                              n3=CurrentNode->idYu;
                              n4=CurrentNode->idYd;

                              N1 = i - n1;
                              N2 = i + n2;
                              N3 = j + n3;
                              N4 = j - n4;

                              UpNode     = &(pJ->GetValue(i,N3));
                              DownNode   = &(pJ->GetValue(i,N4));
                              RightNode  = &(pJ->GetValue(N2,j));
                              LeftNode   = &(pJ->GetValue(N1,j));

                              n_n = max(n1+n2,1);
                              m_m = max(n3+n4,1);

                              dx_1_n_n_1 = dx_1/n_n; 
                              dy_1_m_m_1 = dy_1/m_m;
                              
                              Num_Eq = (FlowNode2D<FP,NUM_COMPONENTS>::NumEq-SetTurbulenceModel(CurrentNode));

                              for (int k=0;k<Num_Eq;k++ ) {

                                  int c_flag = 0;

                                  if ( k < 4 ) // Make bit flags for future test for current equation // FlowNode2D<FP,NUM_COMPONENTS>::NumEq-AddEq-NUM_COMPONENTS ?
                                      c_flag  = CT_Rho_CONST_2D   << k;
                                  else if (k<(4+NUM_COMPONENTS))  // 7 ?
                                      c_flag  = CT_Y_CONST_2D;
                                  else if( ProblemType == SM_NS  &&
                                          (CurrentNode->isTurbulenceCond2D(TCT_k_eps_Model_2D) ||
                                           CurrentNode->isTurbulenceCond2D(TCT_Spalart_Allmaras_Model_2D) )) 
                                      c_flag  = TCT_k_CONST_2D << (k-4-NUM_COMPONENTS); 

                                  if ( !CurrentNode->isCond2D((CondType2D)c_flag) && 
                                        CurrentNode->S[k] != 0. ) {

                                        FP Tmp = CurrentNode->S[k];
                                        FP beta_min;
                                        FP sqrt_RES = 0;
                                        FP absDD    = 0.;

                                        if(k == i2d_RhoU && k == i2d_RhoV ) {
                                            //Tmp = sqrt(CurrentNode->S[i2d_RhoU]*CurrentNode->S[i2d_RhoU]+
                                            //           CurrentNode->S[i2d_RhoV]*CurrentNode->S[i2d_RhoV]+1.e-30); // Flux
                                            Tmp = max(fabs(CurrentNode->S[i2d_RhoU]),fabs(CurrentNode->S[i2d_RhoV]));// max Flux
                                        }
                                        
                                        absDD       = NextNode->S[k]-CurrentNode->S[k];
                                        
                                        if(fabs(Tmp) > 1.e-15) {
                                            DD_local[k] = fabs(absDD/Tmp);
                                            sqrt_RES    = sqrt(DD_local[k]);
                                        } else {
                                            DD_local[k] = 1.0;
                                        }

                                        beta_min = min(beta0,beta_Scenario_Val);

                                        if(CurrentNode->isCond2D(CT_NONREFLECTED_2D)) {
                                           beta_min = nrbc_beta0;
                                        }

                                        if( bFF == BFF_L) {
                                        //LINEAR locally adopted blending factor function  (LLABFF)
                                          CurrentNode->beta[k] = min(beta_min,(beta_min*beta_min)/(beta_min+DD_local[k]));
                                        } else if( bFF == BFF_LR) {
                                        //LINEAR locally adopted blending factor function with relaxation (LLABFFR)
                                          CurrentNode->beta[k] = min((beta_min+CurrentNode->beta[k])*0.5,(beta_min*beta_min)/(beta_min+DD_local[k]));
                                        } else if( bFF == BFF_S) {
                                          //SQUARE locally adopted blending factor function (SLABF)
                                          CurrentNode->beta[k] = min(beta_min,(beta_min*beta_min)/(beta_min+DD_local[k]*DD_local[k]));
                                        } else if (bFF == BFF_SR) {
                                          //SQUARE locally adopted blending factor function with relaxation (SLABFFR)
                                          CurrentNode->beta[k] = min((beta_min+CurrentNode->beta[k])*0.5,(beta_min*beta_min)/(beta_min+DD_local[k]*DD_local[k]));
                                        } else if( bFF == BFF_SQR) {
                                        //SQRT() locally adopted blending factor function (SQRLABF) + most accurate & stable +
                                          CurrentNode->beta[k] = min(beta_min,(beta_min*beta_min)/(beta_min+sqrt_RES));
                                        } else if( bFF == BFF_SQRR) {
                                          CurrentNode->beta[k] = min((beta_min+CurrentNode->beta[k])*0.5,(beta_min*beta_min)/(beta_min+sqrt_RES)); 
                                        }
#ifdef _MPI
                                            DD_max[rank].DD[k].DD      = max(DD_max[rank].DD[k].DD,DD_local[k]);
                                            
                                            if (isAlternateRMS) {
                                                DD_max[rank].DD[k].RMS    += absDD*absDD;
                                            } else {
                                                DD_max[rank].DD[k].RMS    += DD_local[k]*DD_local[k];
                                            }
                                            
                                            DD_max[rank].DD[k].sumDiv += Tmp*Tmp;
                                            DD_max[rank].DD[k].iRMS++;

                                              if (DD_max[rank].DD[k].DD==DD_local[k] ) {
                                                  DD_max[rank].DD[k].i = i;
                                                  DD_max[rank].DD[k].j = j;
                                              }

#else
                                        
                                        if (isAlternateRMS) {
                                             sumDiv(k,ii) += Tmp*Tmp;
                                             RMS(k,ii)    += absDD;
                                         } else {
                                             RMS(k,ii)    += DD_local[k]*DD_local[k];
                                         }
                                        
                                         iRMS(k,ii)++;
                                        
                                         DD_max(k,ii) = max(DD_max(k,ii),DD_local[k]);

                                        if ( DD_max(k,ii) == DD_local[k] ) {
                                             i_c[ii] = i;
                                             j_c[ii] = j;
                                        }
#endif // _MPI

                                        }
                                        
                                        if (k<(4+NUM_COMPONENTS)) {
                                            if ( !CurrentNode->isCond2D((CondType2D)c_flag) )
                                                  CurrentNode->S[k]   = NextNode->S[k];
                                        } else if (ProblemType == SM_NS  &&
                                                   (CurrentNode->isTurbulenceCond2D(TCT_k_eps_Model_2D) ||
                                                    CurrentNode->isTurbulenceCond2D(TCT_Spalart_Allmaras_Model_2D)) ){
                                            if ( !CurrentNode->isTurbulenceCond2D((TurbulenceCondType2D)c_flag) )
                                                  CurrentNode->S[k]   =  NextNode->S[k];
                                        }
                                    }
                                    
                                    //CurrentNode->beta[i2d_RhoV] = CurrentNode->beta[i2d_RhoU] = max(CurrentNode->beta[i2d_RhoU],CurrentNode->beta[i2d_RhoV]);  // for symmetry keeping
                                    
                                    if(ProblemType == SM_NS) {

                                        FP  rhoY_air_Right = RightNode->S[i2d_Rho];
                                        FP  rhoY_air_Left  = LeftNode->S[i2d_Rho];
                                        FP  rhoY_air_Up    = UpNode->S[i2d_Rho];
                                        FP  rhoY_air_Down  = DownNode->S[i2d_Rho];
                                        
                                        CurrentNode->droYdx[NUM_COMPONENTS]=CurrentNode->droYdy[NUM_COMPONENTS]=0.;
                                        
                                        for (int k=4;k<FlowNode2D<FP,NUM_COMPONENTS>::NumEq-2;k++ ) {
                                            if ( !CurrentNode->isCond2D(CT_dYdx_NULL_2D) ) {
                                                CurrentNode->droYdx[k-4]=(RightNode->S[k]-LeftNode->S[k])*dx_1_n_n_1;
                                                rhoY_air_Right -= RightNode->S[k];
                                                rhoY_air_Left  -= LeftNode->S[k];
                                            }
                                            if ( !CurrentNode->isCond2D(CT_dYdy_NULL_2D) ) {
                                                  CurrentNode->droYdy[k-4]=(UpNode->S[k]-DownNode->S[k])*dy_1_m_m_1;
                                                  rhoY_air_Up    -= UpNode->S[k];
                                                  rhoY_air_Down  -= DownNode->S[k];
                                            }
                                        }
                                        
                                        if ( !CurrentNode->isCond2D(CT_dYdx_NULL_2D) ) {
                                            CurrentNode->droYdx[NUM_COMPONENTS]=(rhoY_air_Right - rhoY_air_Left)*dx_1_n_n_1;
                                        }
                                        
                                        if ( !CurrentNode->isCond2D(CT_dYdy_NULL_2D) ) {
                                            CurrentNode->droYdy[NUM_COMPONENTS]=(rhoY_air_Up - rhoY_air_Down)*dy_1_m_m_1;
                                        }
                                        
                                        if (CurrentNode->isCond2D(CT_WALL_NO_SLIP_2D) || CurrentNode->isCond2D(CT_WALL_LAW_2D) )  {
                                            CurrentNode->dUdx=(RightNode->U*n1-LeftNode->U*n2)*dx_1_n_n_1;
                                            CurrentNode->dVdx=(RightNode->V*n1-LeftNode->V*n2)*dx_1_n_n_1;

                                            CurrentNode->dUdy=(UpNode->U*n3-DownNode->U*n4)*dy_1_m_m_1;
                                            CurrentNode->dVdy=(UpNode->V*n3-DownNode->V*n4)*dy_1_m_m_1;

                                            if(CurrentNode->isTurbulenceCond2D(TCT_k_eps_Model_2D)){
                                              CurrentNode->dkdx   =(RightNode->S[i2d_k]*n1-LeftNode->S[i2d_k]*n2)*dx_1_n_n_1/CurrentNode->S[i2d_Rho];
                                              CurrentNode->depsdx =(RightNode->S[i2d_eps]*n1-LeftNode->S[i2d_eps]*n2)*dx_1_n_n_1/CurrentNode->S[i2d_Rho];

                                              CurrentNode->dkdy   =(UpNode->S[i2d_k]*n3-DownNode->S[i2d_k]*n4)*dy_1_m_m_1/CurrentNode->S[i2d_Rho];
                                              CurrentNode->depsdy =(UpNode->S[i2d_eps]*n3-DownNode->S[i2d_eps]*n4)*dy_1_m_m_1/CurrentNode->S[i2d_Rho];
                                            } else if (CurrentNode->isTurbulenceCond2D(TCT_Spalart_Allmaras_Model_2D)) {
                                                       turb_mod_name_index = 3;
                                                       CurrentNode->dkdx   =(RightNode->S[i2d_k]*n1-LeftNode->S[i2d_k]*n2)*dx_1_n_n_1/CurrentNode->S[i2d_Rho];
                                                       CurrentNode->dkdy   =(UpNode->S[i2d_k]*n3-DownNode->S[i2d_k]*n4)*dy_1_m_m_1/CurrentNode->S[i2d_Rho];
                                            }
                                        } else {
                                            CurrentNode->dUdx   =(RightNode->U-LeftNode->U)*dx_1_n_n_1;
                                            CurrentNode->dVdx   =(RightNode->V-LeftNode->V)*dx_1_n_n_1;

                                            CurrentNode->dUdy   =(UpNode->U-DownNode->U)*dy_1_m_m_1;
                                            CurrentNode->dVdy   =(UpNode->V-DownNode->V)*dy_1_m_m_1;
                                            if(CurrentNode->isTurbulenceCond2D(TCT_k_eps_Model_2D)){
                                              CurrentNode->dkdx   =(RightNode->S[i2d_k]-LeftNode->S[i2d_k])*dx_1_n_n_1/CurrentNode->S[i2d_Rho];
                                              CurrentNode->depsdx =(RightNode->S[i2d_eps]-LeftNode->S[i2d_eps])*dx_1_n_n_1/CurrentNode->S[i2d_Rho];

                                              CurrentNode->dkdy   =(UpNode->S[i2d_k]-DownNode->S[i2d_k])*dy_1_m_m_1/CurrentNode->S[i2d_Rho];
                                              CurrentNode->depsdy =(UpNode->S[i2d_eps]-DownNode->S[i2d_eps])*dy_1_m_m_1/CurrentNode->S[i2d_Rho];
                                            } else if (CurrentNode->isTurbulenceCond2D(TCT_Spalart_Allmaras_Model_2D)) {
                                                       turb_mod_name_index = 3;
                                                       CurrentNode->dkdx   =(RightNode->S[i2d_k]-LeftNode->S[i2d_k])*dx_1_n_n_1/CurrentNode->S[i2d_Rho];
                                                       CurrentNode->dkdy   =(UpNode->S[i2d_k]-DownNode->S[i2d_k])*dy_1_m_m_1/CurrentNode->S[i2d_Rho];
                                            }
                                        }

                                        CurrentNode->dTdx=(RightNode->Tg-LeftNode->Tg)*dx_1_n_n_1;
                                        CurrentNode->dTdy=(UpNode->Tg-DownNode->Tg)*dy_1_m_m_1;
                                    }
                                    
                                    if((int)(iter+last_iter) < TurbStartIter) {
                                       CurrentNode->FillNode2D(0,isTurbulenceReset,SigW,SigF,(TurbulenceExtendedModel)TurbExtModel,delta_bl,ProblemType);
                                    } else {
                                       CurrentNode->FillNode2D(1,0,SigW,SigF,(TurbulenceExtendedModel)TurbExtModel,delta_bl,ProblemType);
                                    }

                                    if( CurrentNode->Tg < 0. ) {
                                        *f_stream << "\nTg=" << CurrentNode->Tg << " K. p=" << CurrentNode->p <<" Pa dt=" << dt << " sec.\n" << flush;
#ifdef _MPI
                                        *f_stream << "\nERROR: Computational unstability in UMatrix2D< FlowNode2D<FP,NUM_COMPONENTS> >(" << (int)(x0/FlowNode2D<FP,NUM_COMPONENTS>::dx) + i <<","<< j << ") \nNode Conditions {\n";
                                         PrintCond(f_stream,CurrentNode);
                                        *f_stream  <<"} on iteration " << iter+last_iter<< "...\n";
#else
  #ifdef _OPENMP
                                        *f_stream << "\nERROR: Computational unstability in local (num_thread="<< omp_get_thread_num() << ") UMatrix2D< FlowNode2D<FP,NUM_COMPONENTS> >(" << i <<","<< j <<") \nNode Conditions {\n";
                                        PrintCond(f_stream,CurrentNode);
                                        *f_stream  <<"} on iteration " << iter+last_iter<< "...\n";
  #else
                                        *f_stream << "\nERROR: Computational unstability in UMatrix2D< FlowNode2D<FP,NUM_COMPONENTS> >(" << i <<","<< j <<") \nNode Conditions {\n";
                                        PrintCond(f_stream,CurrentNode);
                                        *f_stream  <<"} on iteration " << iter+last_iter<< "...\n";
  #endif // _OPENMP
#endif // _MPI

#ifdef _OPENMP
#pragma omp critical
                                          {
                                            char omp_ErrFileName[255];
                                            snprintf(omp_ErrFileName,255,"tid-%d-%s",omp_get_thread_num(),ErrFileName);  
                                            DataSnapshot(omp_ErrFileName);
#endif // _OPENMP
#ifdef _MPI
                                          if( rank == 0 ) {
                                            char mpi_ErrFileName[255];
                                            snprintf(mpi_ErrFileName,255,"rank-%d-%s",rank,ErrFileName);  
                                            DataSnapshot(mpi_ErrFileName);
#else  //  _MPI

                                            DataSnapshot(ErrFileName);
#endif // _MPI
#ifdef _OPENMP
                                           *f_stream << "Computation terminated. Error data saved in file "<< omp_ErrFileName <<" \n" ;
#else
#endif //  _OPENMP

#ifdef _MPI
#ifdef _OPENMP
                                            if(rank == 0)
                                              *f_stream << "Computation terminated. Error data saved in file "<< mpi_ErrFileName <<" \n";
#else  //  _MPI

                                            *f_stream << "Computation terminated. Error data saved in file "<< ErrFileName <<" \n" ;
#endif //  _OPENMP
#endif // _MPI

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
                                            }
#ifdef _OPENMP
                                         }
#endif // _OPENMP
#ifdef _MPI
                                         }
#endif // _MPI

                                      isRun = 0;
                                      Abort_OpenHyperFLOW2D();
                                    }  else {
                                            FP CFL_min      = min(CFL,CFL_Scenario_Val);
                                            AAA                 = sqrt(CurrentNode->k*CurrentNode->R*CurrentNode->Tg); 
                                            dt_min_local        = CFL_min*
                                                                  min(dx/(AAA+fabs(CurrentNode->U)),dy/(AAA+fabs(CurrentNode->V)));
#ifdef _MPI
                                            DD_max[rank].dt_min = min(DD_max[rank].dt_min, dt_min_local);

#else
                                            dt_min[ii]          = min(dt_min[ii],dt_min_local);
#endif // _MPI
                                            CalcChemicalReactions(CurrentNode,CRM_ZELDOVICH, (void*)(&chemical_reactions));
                                    }
                       } else if (CurrentNode->isCond2D(NT_FC_2D)) {
                         CurrentNode->FillNode2D(1,0,SigW,SigF,(TurbulenceExtendedModel)TurbExtModel,delta_bl,ProblemType);
                       }
                  }
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
                                      ,ii,(int)SubDomainArray->GetNumElements()-1 
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
 
 dtmin = min(dt_min[ii],dtmin);
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
          for(int ii=1;ii<last_rank+1;ii++) 
              DD_Exchange[last_rank+ii].Wait();  
        }
#endif //_MPI_NB
          
          if (MonitorPointsArray &&
              iter/NOutStep*NOutStep == iter ) {
              for(int ii_monitor=0;ii_monitor<(int)MonitorPointsArray->GetNumElements();ii_monitor++) {
                  if(rank == MonitorPointsArray->GetElement(ii_monitor).rank) {

                      int i_i = (MonitorPointsArray->GetElement(ii_monitor).MonitorXY.GetX() - x0 - FlowNode2D<FP,NUM_COMPONENTS>::dx*0.5)/FlowNode2D<FP,NUM_COMPONENTS>::dx;
                      int j_j = MonitorPointsArray->GetElement(ii_monitor).MonitorXY.GetY()/FlowNode2D<FP,NUM_COMPONENTS>::dy;
                      MonitorPointsArray->GetElement(ii_monitor).p  = pJ->GetValue(i_i,j_j).p;
                      MonitorPointsArray->GetElement(ii_monitor).T  = pJ->GetValue(i_i,j_j).Tg;
                  }

                  MPI::COMM_WORLD.Bcast(&MonitorPointsArray->GetElement(ii_monitor).p,
                                        sizeof(FP),
                                        MPI::DOUBLE,
                                        MonitorPointsArray->GetElement(ii_monitor).rank);

                  MPI::COMM_WORLD.Bcast(&MonitorPointsArray->GetElement(ii_monitor).T,
                                        sizeof(FP),
                                        MPI::DOUBLE,
                                        MonitorPointsArray->GetElement(ii_monitor).rank);
              }
        }
        
        if(rank == 0) {
            for(int ii=0,DD_max_var=0.;ii<last_rank+1;ii++) {
                for (k=0;k<(int)(FlowNode2D<FP,NUM_COMPONENTS>::NumEq);k++ ) {
                   
                   DD_max_var = max(DD_max[ii].DD[k].DD,DD_max_var);
                   
                   RMS[k]    += DD_max[ii].DD[k].RMS;
                   iRMS[k]   += DD_max[ii].DD[k].iRMS;
                   sumDiv[k] += DD_max[ii].DD[k].sumDiv;

                   if(DD_max[ii].DD[k].DD == DD_max_var) {
                    i_max = DD_max[ii].DD[k].i;
                    j_max = DD_max[ii].DD[k].j;
                    k_max = k;
                  }
              }
            }
            
            for (k=0;k<(int)(FlowNode2D<FP,NUM_COMPONENTS>::NumEq);k++ ) {
                
                
                if (isAlternateRMS) {
                    if(RMS[k] > 0.0 && sumDiv[k] > 0) { 
                       RMS[k] = sqrt(RMS[k]/sumDiv[k]);
                    }
                } else {
                    if(iRMS[k] > 0) {
                       RMS[k] = sqrt(RMS[k]/iRMS[k]);
                    }
                }
                  
                  if(MonitorIndex == 0 || MonitorIndex > 4) {
                     max_RMS = max(RMS[k],max_RMS);
                     if(max_RMS == RMS[k])
                        k_max_RMS = k;
                  } else {
                     max_RMS = max(RMS[MonitorIndex-1],max_RMS);
                     if(max_RMS == RMS[MonitorIndex-1])
                        k_max_RMS = k;
                  }

                 
                }
#else
#ifdef _OPENMP
}
                                       

#pragma omp single
#endif // _OPENMP
    {
        for(k=0;k<(int)(FlowNode2D<FP,NUM_COMPONENTS>::NumEq);k++ )     {
         
            for(int ii=0;ii<(int)SubDomainArray->GetNumElements();ii++) {
             if (isAlternateRMS) {
                 if(sumDiv(k,ii) > 0.0) {
                   sum_RMS[k]      += RMS(k,ii);
                   sum_RMS_Div[k]  +=sumDiv(k,ii);
                 }
                 if(sum_RMS_Div[k] > 0.0 && sum_RMS[k] > 0.0)
                    sum_RMS[k] = sqrt(sum_RMS[k]/sum_RMS_Div[k]);
                 else
                    sum_RMS[k] = 0.0;
             } else {
                 if(iRMS(k,ii) > 0) {
                     sum_RMS[k]  += RMS(k,ii);
                     sum_iRMS[k] += iRMS(k,ii);
                 }
                 if(sum_iRMS[k] != 0 && sum_RMS[k] > 0.0 )
                    sum_RMS[k] = sqrt(sum_RMS[k]/sum_iRMS[k]);
                 else
                    sum_RMS[k] = 0;
             }
             
           }

           if( MonitorIndex == 0 || MonitorIndex > 4) {
               max_RMS = max(sum_RMS[k],max_RMS);
               if(max_RMS == sum_RMS[k])
                  k_max_RMS = k;
           } else {
               max_RMS = max(sum_RMS[MonitorIndex-1],max_RMS);
               if(max_RMS == sum_RMS[MonitorIndex-1])
                  k_max_RMS = MonitorIndex-1;
           }

         }

        if (MonitorPointsArray &&
            iter/NOutStep*NOutStep == iter ) {
            for(int ii_monitor=0;ii_monitor<(int)MonitorPointsArray->GetNumElements();ii_monitor++) {
                    int i_i = (MonitorPointsArray->GetElement(ii_monitor).MonitorXY.GetX()
#ifdef _MPI
                               - x0
#endif //_MPI
                               )/FlowNode2D<FP,NUM_COMPONENTS>::dx;
                    int j_j = MonitorPointsArray->GetElement(ii_monitor).MonitorXY.GetY()/FlowNode2D<FP,NUM_COMPONENTS>::dy;
                    
                    MonitorPointsArray->GetElement(ii_monitor).p = 
#ifdef _MPI
                        pJ->GetValue(i_i,j_j).p;
#else
                        J->GetValue(i_i,j_j).p;
#endif //_MPI
                        MonitorPointsArray->GetElement(ii_monitor).T = 
#ifdef _MPI
                        pJ->GetValue(i_i,j_j).Tg;
#else
                        J->GetValue(i_i,j_j).Tg;
#endif //_MPI
            
            }
        }
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
             SaveRMS(pRMS_OutFile,last_iter+iter,
#ifdef  _MPI
             RMS);
#else
             sum_RMS);
#endif // _MPI
             
             if(MonitorPointsArray && MonitorPointsArray->GetNumElements() > 0) {
                SaveMonitors(pMonitors_OutFile,GlobalTime+CurrentTimePart,MonitorPointsArray);
             }


             if(k_max_RMS == i2d_nu_t)
                k_max_RMS +=turb_mod_name_index;

             if(k_max_RMS != -1 && (MonitorIndex == 0 || MonitorIndex == 5))
             *f_stream << "Step No " << iter+last_iter << " maxRMS["<< RMS_Name[k_max_RMS] << "]="<< (FP)(max_RMS*100.) \
                        <<  " % step_time=" << (FP)d_time << " sec (" << (FP)VCOMP <<" step/sec) dt="<< dt <<"\n" << flush;
             else if(MonitorIndex > 0 &&  MonitorIndex < 5 )
                 *f_stream << "Step No " << iter+last_iter << " maxRMS["<< RMS_Name[MonitorIndex-1] << "]="<< (FP)(max_RMS*100.) \
                  <<  " % step_time=" << (FP)d_time << " sec (" << (FP)VCOMP <<" step/sec) dt="<< dt <<"\n" << flush;
             else
             *f_stream << "Step No " << iter+last_iter << " maxRMS["<< k_max_RMS << "]="<< (FP)(max_RMS*100.) \
                        <<  " % step_time=" << (FP)d_time << " sec (" << (FP)VCOMP <<" step/sec) dt="<< dt <<"\n" << flush;
              f_stream->flush();
             }
//#ifdef _MPI
        }
//#endif // _MPI
     iter++;
   } while((int)iter < Nstep);
#ifdef _OPENMP
#pragma omp single 
          {
#endif //  _OPENMP

#ifdef _MPI
#ifdef _PARALLEL_RECALC_Y_PLUS_
   if(ProblemType == SM_NS) {
       
       if (rank == 0 ) {
           
           if(WallNodes && !WallNodesUw_2D && J) 
               WallNodesUw_2D = GetWallFrictionVelocityArray2D(J,WallNodes);

           if(WallNodes && 
              WallNodesUw_2D && J ) {
               *f_stream << "Recalc wall friction velocity in (" << WallNodesUw_2D->GetNumElements() <<") wall nodes..." ;
               RecalcWallFrictionVelocityArray2D(J,WallNodesUw_2D,WallNodes);
               *f_stream << "OK" << endl;
           }
           
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
   }
#endif // _PARALLEL_RECALC_Y_PLUS_
     // Collect all subdomain
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
           tmp_RecvPtr=(void*)((u_long)(ArraySubDomain->GetElement(ii)->GetMatrixPtr())+
                                        ArraySubDomain->GetElement(ii)->GetColSize()*(r_Overlap+
                                                                                      l_Overlap));
           tmp_RecvSize=ArraySubDomain->GetElement(ii)->GetMatrixSize()-
                        ArraySubDomain->GetElement(ii)->GetColSize()*(r_Overlap-l_Overlap);
#ifdef _IMPI_
           LongMatrixRecv(ii, tmp_RecvPtr, tmp_RecvSize);  // Low Mem Send subdomain
#else
           MPI::COMM_WORLD.Recv(tmp_RecvPtr,
                                tmp_RecvSize,
                                MPI::BYTE,ii,tag_Matrix); 
#endif // _IMPI_
           }
        
        }
        
            if( rank == 0 ) {
                
                if ( isGasSource && SrcList) {
                      *f_stream << "\nSet gas sources...";
                      SrcList->SetSources2D();
                      *f_stream << "OK" << endl;
                 }
#ifndef _PARALLEL_RECALC_Y_PLUS_
                *f_stream << "Recalc y+...";
                 Recalc_y_plus(J,WallNodes);
#endif // _IMPI_
            }
                
    
// Get Data back
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
    } else {
      void*  tmp_SendPtr;
      u_long tmp_SendSize;
      for(int ii=1;ii<last_rank+1;ii++) {
          tmp_SendPtr=(void*)((u_long)(ArraySubDomain->GetElement(ii)->GetMatrixPtr())+
          ArraySubDomain->GetElement(ii)->GetColSize()*(r_Overlap+l_Overlap));
          tmp_SendSize=ArraySubDomain->GetElement(ii)->GetMatrixSize()-
          ArraySubDomain->GetElement(ii)->GetColSize()*(r_Overlap-l_Overlap);
#ifdef _IMPI_
          LongMatrixSend(ii, tmp_SendPtr, tmp_SendSize);  // Low Mem Send subdomain
#else
          MPI::COMM_WORLD.Send(tmp_SendPtr,
                               tmp_SendSize,
                               MPI::BYTE,ii,tag_Matrix);
#endif //_IMPI_ 
      }
    }

    if( rank == 0) {
#endif //  _MPI
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

         if(is_Cx_calc) { // For Airfoils only
          *f_stream << "\nCx = " << Calc_Cx_2D(J,x0_body,y0_body,dx_body,dy_body,Flow2DList->GetElement(Cx_Flow_index-1)) << 
                       " Cy = "  << Calc_Cy_2D(J,x0_body,y0_body,dx_body,dy_body,Flow2DList->GetElement(Cx_Flow_index-1)) << 
                       " Fx = "  << CalcXForce2D(J,x0_body,y0_body,dx_body,dy_body) << " Fy = " << CalcYForce2D(J,x0_body,y0_body,dx_body,dy_body)  << endl;
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
     MPI::COMM_WORLD.Bcast(&GlobalTime,1,MPI::DOUBLE,0);
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
//#pragma omp barrier
#endif //  _OPENMP

          if(MonitorIndex < 5) {
              if( max_RMS > ExitMonitorValue )
                MonitorCondition = 1;
              else
                MonitorCondition = 0;
           } else {
              if( GlobalTime <  ExitMonitorValue )
                MonitorCondition = 1;
              else
                MonitorCondition = 0;
           }
#ifdef _MPI
 MPI::COMM_WORLD.Bcast(&MonitorCondition,1,MPI::INT,0);
#endif //  _MPI
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
                isRun = 0;
                //Abort_OpenHyperFLOW2D();
#ifdef _MPI
       }

        MPI::COMM_WORLD.Barrier();

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
                //Abort_OpenHyperFLOW2D();
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
                //Abort_OpenHyperFLOW2D();
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
                  //Abort_OpenHyperFLOW2D();
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
                //Abort_OpenHyperFLOW2D();
            }
            __end_except;
#endif // _DEBUG_0

#ifdef _MPI
           if(rank == 0)
#endif // _MPI
             *f_stream << "\nReady. Computation finished.\n"  << flush;

           Exit_OpenHyperFLOW2D();
};

u_long SetWallNodes(ofstream* f_str, ComputationalMatrix2D* pJ) {
u_long NumWallNodes=0;
    for (int j=0;j<(int)pJ->GetY();j++ ) {
       for (int i=0;i<(int)pJ->GetX();i++ ) {

                FlowNode2D< FP,NUM_COMPONENTS >* CurrentNode=NULL;
                FlowNode2D< FP,NUM_COMPONENTS >* UpNode=NULL;
                FlowNode2D< FP,NUM_COMPONENTS >* DownNode=NULL;
                FlowNode2D< FP,NUM_COMPONENTS >* LeftNode=NULL;
                FlowNode2D< FP,NUM_COMPONENTS >* RightNode=NULL;

                CurrentNode = &pJ->GetValue(i,j);

                if(!CurrentNode->isCond2D(CT_SOLID_2D) &&
                   !CurrentNode->isCond2D(NT_FC_2D)) {

                if(j < (int)pJ->GetY()-1)
                  UpNode    = &(pJ->GetValue(i,j+1));
                else
                  UpNode    = NULL;

                if(j>0)
                  DownNode  = &(pJ->GetValue(i,j-1));
                else
                  DownNode  = NULL;

                if(i>0)
                  LeftNode  = &(pJ->GetValue(i-1,j));
                else
                  LeftNode  = NULL;

                if(i < (int)pJ->GetX()-1)
                  RightNode = &(pJ->GetValue(i+1,j));
                else
                  RightNode = NULL;

               if(UpNode && UpNode->isCond2D(CT_SOLID_2D)) {
                   CurrentNode->SetCond2D(NT_WNS_2D);
                   NumWallNodes++;
               } else if(DownNode && DownNode->isCond2D(CT_SOLID_2D)) {
                   CurrentNode->SetCond2D(NT_WNS_2D);
                   NumWallNodes++;
               } else if(LeftNode && LeftNode->isCond2D(CT_SOLID_2D)) {
                   CurrentNode->SetCond2D(NT_WNS_2D);
                   NumWallNodes++;
               } else if(RightNode && RightNode->isCond2D(CT_SOLID_2D)) {
                   CurrentNode->SetCond2D(NT_WNS_2D);
                   NumWallNodes++;
               }
           }
       }
    }
 *f_str << "\nFound and restore " << NumWallNodes << " wall nodes" << endl;
 return  NumWallNodes;
}

UArray< XY<int> >* GetWallNodes(ofstream* f_str, ComputationalMatrix2D* pJ, int isPrint) {
    XY<int> ij;
    UArray< XY<int> >* WallNodes;
    WallNodes = new UArray< XY<int> >();
    if ( isPrint )
        *f_str << "Scan computation area for lookup wall nodes.\n" << flush;
    for (int j=0;j<(int)pJ->GetY();j++ ) {
       for (int i=0;i<(int)pJ->GetX();i++ ) {
           if (!pJ->GetValue(i,j).isCond2D(CT_SOLID_2D)){
            if ( pJ->GetValue(i,j).isCond2D(CT_WALL_LAW_2D) || 
                 pJ->GetValue(i,j).isCond2D(CT_WALL_NO_SLIP_2D)) {
                ij.SetXY(i,j);
                WallNodes->AddElement(&ij);
            }
           } 
        }
       if ( isPrint == 1 )  {
          *f_str << "Scan "<< 100*j/MaxY << "% nodes\r"<< flush; 
       } else if ( isPrint == 2 ) {
          *f_str << "Scan "<< 100*j/MaxY << "% nodes\n"<< flush; 
       }
    }
 return WallNodes;
}


UArray< XY<int> >*
ScanArea(ofstream* f_str,ComputationalMatrix2D* pJ ,int isPrint) {
    TurbulenceCondType2D TM = TCT_No_Turbulence_2D;
    timeval  time1, time2;
    unsigned int i,j,jj,num_threads=1; //ii,iw,jw,
    long num_active_nodes,active_nodes_per_SubDomain;
    UArray< XY<int> >* pWallNodes;
    UArray< XY<int> >* SubDomain;
    pWallNodes = new UArray< XY<int> >();
    SubDomain = new UArray< XY<int> >();
    XY<int> ij;
    XY<int> ijsm;

    num_active_nodes=0;

    if ( isPrint )
        *f_str << "Scan computation area for lookup wall nodes.\n" << flush;
    for ( j=0;j<MaxY;j++ ) {
       for ( i=0;i<MaxX;i++ ) {
           if (!pJ->GetValue(i,j).isCond2D(CT_SOLID_2D)){
               num_active_nodes++;
               pJ->GetValue(i,j).SetCond2D(CT_NODE_IS_SET_2D);
            if ( pJ->GetValue(i,j).isCond2D(CT_WALL_LAW_2D) || 
                 pJ->GetValue(i,j).isCond2D(CT_WALL_NO_SLIP_2D)) {
                ij.SetXY(i,j);
                pWallNodes->AddElement(&ij);
            }
           } 
        }
       if ( isPrint == 1 ) {
           *f_str << "Scan "<< 100*j/MaxY << "% nodes\r"<< flush;
       } else if ( isPrint == 2 ) {
           *f_str << "Scan "<< 100*j/MaxY << "% nodes\n"<< flush;
       }
    }

#ifdef _MPI
num_threads = MPI::COMM_WORLD.Get_size();
#endif // _MPI

#ifdef _OPENMP
#pragma omp parallel
{
 num_threads = omp_get_num_threads(); 
}
#endif // _OPENMP
active_nodes_per_SubDomain = num_active_nodes/num_threads; 
if ( isPrint )
        *f_str << "Found " << pWallNodes->GetNumElements() <<" wall nodes from " << num_active_nodes << " gas filled nodes (" << num_threads <<" threads, "<< active_nodes_per_SubDomain <<" active nodes per thread).\n" << flush;

    if ( isPrint )
        *f_str << "Scan computation area for lookup minimal distance from internal nodes to wall.\n" << flush;
#ifdef _DEBUG_0 // 2
gettimeofday(&time2,NULL);
#endif // _DEBUG_0 // 2
num_active_nodes=0;
ijsm.SetX(0);

if(isTurbulenceReset && ProblemType == SM_NS) {
   if ( TurbMod == 0 )
        TM = (TurbulenceCondType2D)(TCT_No_Turbulence_2D);
   else if ( TurbMod == 1 )
        TM = (TurbulenceCondType2D)(TCT_Integral_Model_2D);
   else if ( TurbMod == 2 )
        TM = (TurbulenceCondType2D)(TCT_Prandtl_Model_2D);
   else if ( TurbMod == 3 )
        TM = (TurbulenceCondType2D)(TCT_Spalart_Allmaras_Model_2D);
   else if ( TurbMod == 4 )
        TM = (TurbulenceCondType2D)(TCT_k_eps_Model_2D);
   else if ( TurbMod == 5 )
        TM = (TurbulenceCondType2D)(TCT_Smagorinsky_Model_2D);
    if ( isPrint )
        *f_str << "Reset turbulence model to " << PrintTurbCond(TurbMod) << "\n" << flush;
}
#ifdef _OPENMP
#pragma omp for
#endif // _OPENMP
for (int i=0;i<(int)MaxX;i++ ) {
        for (int j=0;j<(int)MaxY;j++ ) {

                    if(isTurbulenceReset && ProblemType == SM_NS) {
                       if(pJ->GetValue(i,j).isTurbulenceCond2D(TCT_Integral_Model_2D))
                          pJ->GetValue(i,j).CleanTurbulenceCond2D(TCT_Integral_Model_2D);
                       if(pJ->GetValue(i,j).isTurbulenceCond2D(TCT_Prandtl_Model_2D))
                          pJ->GetValue(i,j).CleanTurbulenceCond2D(TCT_Prandtl_Model_2D);
                       if(pJ->GetValue(i,j).isTurbulenceCond2D(TCT_Spalart_Allmaras_Model_2D))
                          pJ->GetValue(i,j).CleanTurbulenceCond2D(TCT_Spalart_Allmaras_Model_2D);
                       if(pJ->GetValue(i,j).isTurbulenceCond2D(TCT_k_eps_Model_2D))
                          pJ->GetValue(i,j).CleanTurbulenceCond2D(TCT_k_eps_Model_2D);
                       if(pJ->GetValue(i,j).isTurbulenceCond2D(TCT_Smagorinsky_Model_2D))
                          pJ->GetValue(i,j).CleanTurbulenceCond2D(TCT_Smagorinsky_Model_2D);

                       pJ->GetValue(i,j).SetTurbulenceCond2D(TM);
                       pJ->GetValue(i,j).dkdx = pJ->GetValue(i,j).dkdy = pJ->GetValue(i,j).depsdx = pJ->GetValue(i,j).depsdy =
                       pJ->GetValue(i,j).S[i2d_k] = pJ->GetValue(i,j).S[i2d_eps] =
                       pJ->GetValue(i,j).Src[i2d_k] = pJ->GetValue(i,j).Src[i2d_eps] =
                       pJ->GetValue(i,j).mu_t = pJ->GetValue(i,j).lam_t = 0.0;
                       pJ->GetValue(i,j).FillNode2D(0,1);
                    }

                if (pJ->GetValue(i,j).isCond2D(CT_NODE_IS_SET_2D) && !pJ->GetValue(i,j).isCond2D(CT_SOLID_2D)) {
                    num_active_nodes++;
                    
                    if(num_active_nodes >= active_nodes_per_SubDomain) {
                       ijsm.SetY(i+1);
                       SubDomain->AddElement(&ijsm);
                       ijsm.SetX(i);
                       num_active_nodes = 0;
                    }
                }
        }
#ifndef _OPENMP
   if ( isPrint == 1 && isVerboseOutput ) {
       *f_str << "Scan "<< 100.*i/MaxX << "% nodes   \r"<< flush;
   } else   if ( isPrint == 2 && isVerboseOutput ) {
       *f_str << "Scan "<< 100.*i/MaxX << "% nodes   \n"<< flush;
   }
#endif // _OPENMP

}
#ifdef _DEBUG_0 // 2
            gettimeofday(&time1,NULL);
            *f_str << " Time:" << (time1.tv_sec-time2.tv_sec)+(time1.tv_usec-time2.tv_usec)*1.e-6 << "sec.\n" << flush; 
#endif // _DEBUG_0 // 2
            delete WallNodes;
*f_str << "SubDomain decomposition was finished:\n";
for(jj=0;jj<SubDomain->GetNumElements();jj++) {
   *f_str << "SubDomain[" << jj << "]->["<< SubDomain->GetElementPtr(jj)->GetX() <<","<< SubDomain->GetElementPtr(jj)->GetY() <<"]\n";
}
if(isTurbulenceReset) {
   isTurbulenceReset = 0;
}
f_str->flush();
return SubDomain;
}

void SetInitBoundaryLayer(ComputationalMatrix2D* pJ, FP delta) {
    for (int i=0;i<(int)pJ->GetX();i++ ) {
           for (int j=0;j<(int)pJ->GetY();j++ ) {
                  if (pJ->GetValue(i,j).isCond2D(CT_NODE_IS_SET_2D) &&
                      !pJ->GetValue(i,j).isCond2D(CT_SOLID_2D) &&
                       pJ->GetValue(i,j).time == 0. &&
                       delta > 0) {
                       if(pJ->GetValue(i,j).l_min <= delta)
                          pJ->GetValue(i,j).S[i2d_RhoU] = pJ->GetValue(i,j).S[i2d_RhoU] * pJ->GetValue(i,j).l_min/delta;
                          pJ->GetValue(i,j).S[i2d_RhoV] = pJ->GetValue(i,j).S[i2d_RhoV] * pJ->GetValue(i,j).l_min/delta;
                          pJ->GetValue(i,j).FillNode2D(0,1);
                   }
           }
    }
}

#ifdef _PARALLEL_RECALC_Y_PLUS_
void RecalcWallFrictionVelocityArray2D(ComputationalMatrix2D* pJ,
                                       UArray<FP>* WallFrictionVelocityArray2D,
                                       UArray< XY<int> >* WallNodes2D) { 
    for(int ii=0;ii<(int)WallNodes2D->GetNumElements();ii++) {
        unsigned int iw,jw;
        FP tau_w;
        iw = WallNodes2D->GetElementPtr(ii)->GetX();
        jw = WallNodes2D->GetElementPtr(ii)->GetY();
        tau_w = (fabs(pJ->GetValue(iw,jw).dUdy)  +
                 fabs(pJ->GetValue(iw,jw).dVdx)) * pJ->GetValue(iw,jw).mu;
        WallFrictionVelocityArray2D->GetElement(ii) = sqrt(tau_w/pJ->GetValue(iw,jw).S[i2d_Rho]+1e-30);
    }
}

UArray<FP>* GetWallFrictionVelocityArray2D(ComputationalMatrix2D* pJ, 
                                           UArray< XY<int> >* WallNodes2D) { 
    UArray<FP>* WallFrictionVelocityArray2D;
    WallFrictionVelocityArray2D = new UArray<FP>();
    for(int ii=0;ii<(int)WallNodes2D->GetNumElements();ii++) {
        unsigned int iw,jw;
        FP tau_w, U_w;
        iw = WallNodes2D->GetElementPtr(ii)->GetX();
        jw = WallNodes2D->GetElementPtr(ii)->GetY();
        tau_w = (fabs(pJ->GetValue(iw,jw).dUdy)  +
                 fabs(pJ->GetValue(iw,jw).dVdx)) * pJ->GetValue(iw,jw).mu;
        U_w   = sqrt(tau_w/pJ->GetValue(iw,jw).S[i2d_Rho]+1e-30);
        WallFrictionVelocityArray2D->AddElement(&U_w);
    }
  return WallFrictionVelocityArray2D;
}

void ParallelRecalc_y_plus(ComputationalMatrix2D* pJ, 
                           UArray< XY<int> >* WallNodes,
                           UArray<FP>* WallFrictionVelocity2D,
                           FP x0) {
#ifndef _OLD_Y_PLUS_

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
                             
                             if(iw == pJ->GetValue(i,j).i_wall && 
                                jw == pJ->GetValue(i,j).j_wall ) {
                                 pJ->GetValue(i,j).y_plus = fabs(U_w*pJ->GetValue(i,j).l_min*pJ->GetValue(i,j).S[i2d_Rho]/pJ->GetValue(i,j).mu);
                             }
                        }
                 }
            }
        }

#else 
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
                                                                    (FlowNode2D<FP,NUM_COMPONENTS>::dy))*pJ->GetValue(i,j).S[i2d_Rho]/pJ->GetValue(i,j).mu;
                             } else {
                                 if(pJ->GetValue(i,j).l_min == min(pJ->GetValue(i,j).l_min,
                                                                   sqrt((x-wx)*(x-wx) + (y-wy)*(y-wy))))
                                 pJ->GetValue(i,j).y_plus = U_w*pJ->GetValue(i,j).l_min*pJ->GetValue(i,j).S[i2d_Rho]/pJ->GetValue(i,j).mu;
                             }
                        }
                 }
            }
        }

#endif // _OLD_Y_PLUS_
}
#else
void Recalc_y_plus(ComputationalMatrix2D* pJ, UArray< XY<int> >* WallNodes) {
    unsigned int iw,jw;
    FP tau_w, U_w;
    
    for (int i=0;i<(int)pJ->GetX();i++ ) {
           for (int j=0;j<(int)pJ->GetY();j++ ) {
               if (pJ->GetValue(i,j).isCond2D(CT_NODE_IS_SET_2D) &&
                   !pJ->GetValue(i,j).isCond2D(CT_SOLID_2D)) {
                   iw = pJ->GetValue(i,j).i_wall;
                   jw = pJ->GetValue(i,j).j_wall;
                   tau_w = (fabs(pJ->GetValue(iw,jw).dUdy) + fabs(pJ->GetValue(iw,jw).dVdx)) * pJ->GetValue(iw,jw).mu;
                   if (pJ->GetValue(iw,jw).S[i2d_Rho] > 0.0 &&
                       tau_w > 0.0) {
                       U_w   = sqrt(tau_w/pJ->GetValue(iw,jw).S[i2d_Rho]+1e-30);
                       pJ->GetValue(i,j).y_plus = fabs(U_w*min(dx,dy)*pJ->GetValue(i,j).S[i2d_Rho]/pJ->GetValue(i,j).mu);
                   } else {
                       pJ->GetValue(i,j).y_plus = 0.0;
                   }
               
               }
           }
    }
               
}
#endif // _PARALLEL_RECALC_Y_PLUS_

const char* PrintTurbCond(int TM) {
   if ( TM == 0 )
          return "TCT_No_Turbulence_2D";
   else if ( TM == 1 )
          return "TCT_Integral_Model_2D";
   else if ( TM == 2 )
          return "TCT_Prandtl_Model_2D";
   else if ( TM == 3 )
          return "TCT_Spalart_Allmaras_Model_2D";
   else if ( TM == 4 )
          return "TCT_k_eps_Model_2D";
   else if ( TM == 5 )
          return "TCT_Smagorinsky_Model_2D";
   return "Unknown turbulence model";
}

void PrintCond(ofstream* OutputData, FlowNode2D<FP,NUM_COMPONENTS>* fn) {
    //Const conditions

    if ( fn->isCond2D(CT_Rho_CONST_2D) )
        *OutputData << " CT_Rho_CONST_2D" << endl << flush;
    if ( fn->isCond2D(CT_U_CONST_2D) )
        *OutputData << " CT_U_CONST_2D" << endl << flush;
    if ( fn->isCond2D(CT_V_CONST_2D) )
        *OutputData << " CT_V_CONST_2D" << endl << flush;
    if ( fn->isCond2D(CT_T_CONST_2D) )
        *OutputData << " CT_T_CONST_2D" << endl << flush;
    if ( fn->isCond2D(CT_Y_CONST_2D) )
        *OutputData << " CT_Y_CONST_2D" << endl << flush;

    // dF/dx = 0
    if ( fn->isCond2D(CT_dRhodx_NULL_2D) )
        *OutputData << " CT_dRhodx_NULL_2D" << endl << flush;
    if ( fn->isCond2D(CT_dUdx_NULL_2D) )
        *OutputData << " CT_dUdx_NULL_2D" << endl << flush;
    if ( fn->isCond2D(CT_dVdx_NULL_2D) )
        *OutputData << " CT_dVdx_NULL_2D" << endl << flush;
    if ( fn->isCond2D(CT_dTdx_NULL_2D) )
        *OutputData << " CT_dTdx_NULL_2D" << endl << flush;
    if ( fn->isCond2D(CT_dYdx_NULL_2D) )
        *OutputData << " CT_dYdx_NULL_2D" << endl << flush;

    // dF/dy = 0

    if ( fn->isCond2D(CT_dRhody_NULL_2D) )
        *OutputData << " CT_dRhody_NULL_2D"<< endl << flush;
    if ( fn->isCond2D(CT_dUdy_NULL_2D) )
        *OutputData << " CT_dUdy_NULL_2D" << endl << flush;
    if ( fn->isCond2D(CT_dVdy_NULL_2D) )
        *OutputData << " CT_dVdy_NULL_2D" << endl << flush;
    if ( fn->isCond2D(CT_dTdy_NULL_2D) )
        *OutputData << " CT_dTdy_NULL_2D" << endl << flush;
    if ( fn->isCond2D(CT_dYdy_NULL_2D) )
        *OutputData << " CT_dYdy_NULL_2D" << endl << flush;

    // d2F/dx2 =
    if ( fn->isCond2D(CT_d2Rhodx2_NULL_2D) )
        *OutputData << " CT_d2Rhodx2_NULL_2D" << endl << flush;
    if ( fn->isCond2D(CT_d2Udx2_NULL_2D) )
        *OutputData << " CT_d2Udx2_NULL_2D" << endl << flush;
    if ( fn->isCond2D(CT_d2Vdx2_NULL_2D) )
        *OutputData << " CT_d2Vdx2_NULL_2D" << endl << flush;
    if ( fn->isCond2D(CT_d2Tdx2_NULL_2D) )
        *OutputData << " CT_d2Tdx2_NULL_2D" << endl << flush;
    if ( fn->isCond2D(CT_d2Ydx2_NULL_2D) )
        *OutputData << " CT_d2Ydx2_NULL_2D" << endl << flush;

    // d2F/dy2 = 0

    if ( fn->isCond2D(CT_d2Rhody2_NULL_2D) )
        *OutputData << " CT_d2Rhody2_NULL_2D"<< endl << flush;
    if ( fn->isCond2D(CT_d2Udy2_NULL_2D) )
        *OutputData << " CT_d2Udy2_NULL_2D" << endl << flush;
    if ( fn->isCond2D(CT_d2Vdy2_NULL_2D) )
        *OutputData << " CT_d2Vdy2_NULL_2D" << endl << flush;
    if ( fn->isCond2D(CT_d2Tdy2_NULL_2D) )
        *OutputData << " CT_d2Tdy2_NULL_2D" << endl << flush;
    if ( fn->isCond2D(CT_d2Ydy2_NULL_2D) )
        *OutputData << " CT_d2Ydy2_NULL_2D" << endl << flush;

    // Wall conditions:
    // - Slip
    if ( fn->isCond2D(CT_WALL_LAW_2D) )
        *OutputData << " CT_WALL_LAW_2D" << endl << flush;
    // - Glide
    if ( fn->isCond2D(CT_WALL_NO_SLIP_2D) )
        *OutputData << " CT_WALL_NO_SLIP_2D"  << endl << flush;
    // is solid body ?
    if ( fn->isCond2D(CT_SOLID_2D) )
        *OutputData << " CT_SOLID_2D" << endl << flush;
    // is node activate ?
    if ( fn->isCond2D(CT_NODE_IS_SET_2D) )
        *OutputData << " CT_NODE_IS_SET_2D" << endl << flush;
    // boundary layer refinement
    if ( fn->isCond2D(CT_BL_REFINEMENT_2D) )
        *OutputData << " CT_BL_REFINEMENT_2D" << endl << flush;

    // Special BC
    // non-reflected BC
    if ( fn->isCond2D(CT_NONREFLECTED_2D) )
        *OutputData << " CT_NONREFLECTED_2D" << endl << flush;
}


void DataSnapshot(char* filename, WRITE_MODE ioMode) {
#ifdef _DEBUG_0
    static SysException E=0;
    ___try {
#endif  // _DEBUG_0
        if ( ioMode ) // 0 - append(TecPlot), 1- rewrite(GNUPlot)
            CutFile(filename);
        pOutputData = OpenData(filename);
        SaveData2D(pOutputData,ioMode);
#ifdef _DEBUG_0
    } __except(SysException e) {
        E=e;
    }
    __end_except;
#endif // _DEBUG_0
    pOutputData ->close();
#ifdef _DEBUG_0
    if ( E )
        throw(E);
    else
#endif // _DEBUG_0
        return;
}

ofstream* OpenData(char* outputDataFile) {
    static ofstream* pOutputData;
    pOutputData = new ofstream();
    pOutputData->open(outputDataFile, ios::app);
    return(pOutputData);
}

void CutFile(char* cutFile) {
    static ofstream* pCutFile;
    pCutFile = new ofstream();
    pCutFile->open(cutFile, ios::trunc);
    pCutFile->close();
}

void SaveMonitorsHeader(ofstream* MonitorsFile,UArray< MonitorPoint >* MonitorPtArray) {
    char  TecPlotTitle[1024]={'\0'};
    char  MonitorStr[256];
    snprintf(TecPlotTitle,1024,"#VARIABLES = Time");
    
    for (int i=0;i<(int)MonitorPtArray->GetNumElements();i++) {
         // Point-%i.U, Point-%i.V, 
         snprintf(MonitorStr,256,", Point-%i.p, Point-%i.T",i+1,i+1);
         strcat(TecPlotTitle,MonitorStr);
    }
   *MonitorsFile << TecPlotTitle << endl;
}

void SaveRMSHeader(ofstream* OutputData) {
        char  TecPlotTitle[1024];
        char  TmpData[256]={'\0'};

        //if(is_Cx_calc)
        //   snprintf(TmpData,256,", Cx(N), Cy(N), Fx(N), Fy(N)");

        if(is_Cd_calc)
           snprintf(TmpData,256,", Cd(N), Cv(N)");

        snprintf(TecPlotTitle,1024,"#VARIABLES = N, RMS_Ro(N), RMS_RoU(N), RMS_RoV(N), RMS_RoE(N), RMS_RoY_fu(N), RMS_RoY_ox(N), RMS_RoY_cp(N), RMS_k(N), RMS_eps(N)%s",TmpData);

        *OutputData <<  TecPlotTitle << endl;
    }

void SaveMonitors(ofstream* MonitorsFile, 
                  FP t, 
                  UArray< MonitorPoint >* MonitorPtArray) {
    *MonitorsFile  << t << " ";
    for (int i=0;i<(int)MonitorPtArray->GetNumElements();i++) {
        // MonitorPtArray->GetElement(i).MonitorNode.U << " " << MonitorPtArray->GetElement(i).MonitorNode.V << " " <<
        *MonitorsFile << MonitorPtArray->GetElement(i).p << " " << MonitorPtArray->GetElement(i).T << " ";
    }
    *MonitorsFile << endl;
}

void SaveRMS(ofstream* OutputData,unsigned int n, FP* outRMS) {
         *OutputData <<  n  << " ";
         for(int i=0;i<FlowNode2D<FP,NUM_COMPONENTS>::NumEq;i++) {
             *OutputData <<  outRMS[i] << " " << flush;
         }

        //if(is_Cx_calc) {
        //  *OutputData << Calc_Cx_2D(J,x0_body,y0_body,dx_body,dy_body,Flow2DList->GetElement(Cx_Flow_index-1)) << " " <<  Calc_Cy_2D(J,x0_body,y0_body,dx_body,dy_body,Flow2DList->GetElement(Cx_Flow_index-1)) << " ";
        //  *OutputData << CalcXForce2D(J,x0_body,y0_body,dx_body,dy_body) << " " <<  CalcYForce2D(J,x0_body,y0_body,dx_body,dy_body);
        //}
        
        if(is_Cd_calc) {
          *OutputData << " " << Calc_Cd(J,x0_nozzle,y0_nozzle,dy_nozzle,Flow2DList->GetElement(Cd_Flow_index-1)) << " " <<  Calc_Cv(J,x0_nozzle,y0_nozzle,dy_nozzle,p_ambient,Flow2DList->GetElement(Cd_Flow_index-1)) << " ";
        }

        *OutputData << endl;
    }

    void SaveData2D(ofstream* OutputData, int type) { // type = 1 - GNUPLOT
        int    i,j;
        char   TechPlotTitle1[1024]={0};
        char   TechPlotTitle2[256]={0};
        char   YR[2];
        FP Mach,A,W,Re,Re_t,dx_out,dy_out;
        char   RT[10];
        if(is_p_asterisk_out)
          snprintf(RT,10,"p*");
        else
          snprintf(RT,10,"mu_t/mu");

        if(FlowNode2D<FP,NUM_COMPONENTS>::FT == 1) // FT_FLAT
          snprintf(YR,2,"R");
        else
          snprintf(YR,2,"Y");

        snprintf(TechPlotTitle1,1024,"VARIABLES = X, %s, U, V, T, p, Rho, Y_fuel, Y_ox, Y_cp, Y_i, %s, Mach, l_min, y+, Cp"
                                     "\n",YR, RT); 
        snprintf(TechPlotTitle2,256,"ZONE T=\"Time: %g sec.\" I= %i J= %i F=POINT\n",GlobalTime, MaxX, MaxY);

        if ( type ) {
            *OutputData <<  TechPlotTitle1;
            *OutputData <<  TechPlotTitle2;
        } else {
            *OutputData <<  TechPlotTitle1;
            *OutputData <<  TechPlotTitle2;
        }
        
        dx_out = (dx*MaxX)/(MaxX-1); // dx
        dy_out = (dy*MaxY)/(MaxY-1); // dy
        
        for ( j=0;j<(int)MaxY;j++ ) {
            for ( i=0;i<(int)MaxX;i++ ) {

                *OutputData << i*dx_out*1.e3                    << "  "; // 1
                *OutputData << dy_out*j*1.e3                    << "  "; // 2
                Mach = Re = Re_t = 0;
                if ( !J->GetValue(i,j).isCond2D(CT_SOLID_2D) ) {
                    *OutputData << J->GetValue(i,j).U           << "  "; // 3
                    *OutputData << J->GetValue(i,j).V           << "  "; // 4
                    *OutputData << J->GetValue(i,j).Tg          << "  "; // 5
                    *OutputData << J->GetValue(i,j).p           << "  "; // 6
                    *OutputData << J->GetValue(i,j).S[0]        << "  "; // 7
                    
                    A = sqrt(J->GetValue(i,j).k*J->GetValue(i,j).R*J->GetValue(i,j).Tg+1.e-30);
                    W = sqrt(J->GetValue(i,j).U*J->GetValue(i,j).U+J->GetValue(i,j).V*J->GetValue(i,j).V+1.e-30);
                    Mach = W/A;
                    
                    if ( J->GetValue(i,j).S[0] != 0. ) {
                        *OutputData << J->GetValue(i,j).S[4]/J->GetValue(i,j).S[0] << "  ";  // 8
                        *OutputData << J->GetValue(i,j).S[5]/J->GetValue(i,j).S[0] << "  ";  // 9
                        *OutputData << J->GetValue(i,j).S[6]/J->GetValue(i,j).S[0] << "  ";  // 10
                        *OutputData << fabs(1-J->GetValue(i,j).S[4]/J->GetValue(i,j).S[0]-J->GetValue(i,j).S[5]/J->GetValue(i,j).S[0]-J->GetValue(i,j).S[6]/J->GetValue(i,j).S[0]) << "  "; //11

                        if(is_p_asterisk_out)
                          *OutputData << p_asterisk(&(J->GetValue(i,j))) << "  ";            // 12
                        else
                          *OutputData << J->GetValue(i,j).mu_t/J->GetValue(i,j).mu << "  ";  // 12
                        
                    } else {
                        *OutputData << " +0. +0  +0  +0  +0  "; /* 8 9 10 11 12 */
                    }
                } else {
                    *OutputData << "  0  0  ";                     /* 3 4 */
                    *OutputData << J->GetValue(i,j).Tg;            /* 5 */
                    *OutputData << "  0  0  0  0  0  0  0";        /* 6 7 8 9 10 11 12 */
                }
                if(!J->GetValue(i,j).isCond2D(CT_SOLID_2D)) {
                    if( Mach > 1.e-30) 
                      *OutputData << Mach  << "  " << J->GetValue(i,j).l_min << " " << J->GetValue(i,j).y_plus;  
                    else
                      *OutputData << "  0  0  0  ";
                } else {
                    *OutputData << "  0  0  0  ";
                }
                if(is_Cx_calc)
                   *OutputData << " " << Calc_Cp(&J->GetValue(i,j),Flow2DList->GetElement(Cx_Flow_index-1)) << "\n" ;
                else
                  *OutputData << " 0\n"; 
            }
            if ( type )
                *OutputData <<  "\n" ;
        }
    }

    inline FP kg(FP Cp, FP R) {
        return(Cp/(Cp-R));
    }

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

                if(!CurrentNode->isCond2D(CT_SOLID_2D)) { //  && CurrentNode->isCleanSources == 1
                    
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
                    
                    //if( CurrentNode->isCleanSources ) 
                    //    CurrentNode->SrcAdd[i2d_RhoE] = 0.; 
                    
                    
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
                        
                        CurrentNode->SrcAdd[i2d_RhoE] = -dt*DownNode->Q_conv/dy;
                        //CurrentNode->SrcAdd[i2d_RhoE] += dt*lam_eff*(DownNode->Tg - CurrentNode->Tg)/dy2;
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
                        
                        CurrentNode->SrcAdd[i2d_RhoE] = -dt*UpNode->Q_conv/dy;
                        //CurrentNode->SrcAdd[i2d_RhoE] += dt*lam_eff*(UpNode->Tg - CurrentNode->Tg)/dy2;
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
                        
                        CurrentNode->SrcAdd[i2d_RhoE] = -dt*LeftNode->Q_conv/dx;
                        //CurrentNode->SrcAdd[i2d_RhoE] += dt*lam_eff*(LeftNode->Tg - CurrentNode->Tg)/dx2;
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
                        
                        CurrentNode->SrcAdd[i2d_RhoE] = -dt*RightNode->Q_conv/dx;
                        //CurrentNode->SrcAdd[i2d_RhoE] += dt*lam_eff*(RightNode->Tg - CurrentNode->Tg)/dx2;
                    }
                }
            }
         }
    }
/*  Init -DEEPS2D- solver */
void* InitDEEPS2D(void* lpvParam)
    {
        unsigned int    i_last=0,i_err=0,j_err=0;
        unsigned int    k,j=0;
        ofstream*       f_stream=(ofstream*)lpvParam;
        char            FlowStr[256];
        Flow*           TmpFlow;
        Flow2D*         TmpFlow2D;
        Flow*           pTestFlow=NULL;
        Flow2D*         pTestFlow2D=NULL;
        char            NameBound[128];
        char            NameContour[128];
        char*           BoundStr;
        int             FlowIndex;
        Table*          ContourTable;
        CondType2D      TmpCT = CT_NO_COND_2D;
        TurbulenceCondType2D TmpTurbulenceCT = TCT_No_Turbulence_2D;
        BoundContour2D*   BC; 
        char            ErrorMessage[255];
        u_long          FileSizeGas=0;

        SubDomainArray     = new UArray<UMatrix2D< FlowNode2D<FP,NUM_COMPONENTS> >* >();
        CoreSubDomainArray = new UArray<UMatrix2D< FlowNodeCore2D<FP,NUM_COMPONENTS> >* >();
#ifdef _DEBUG_0
        ___try {
#endif  // _DEBUG_0

            unsigned int NumFlow = Data->GetIntVal((char*)"NumFlow");       // Number of "Flow" objects (for additional information see libFlow library source)
            if ( Data->GetDataError()==-1 ) {
                Abort_OpenHyperFLOW2D();
            }
            unsigned int NumFlow2D = Data->GetIntVal((char*)"NumFlow2D");   // Number of "Flow2D" objects (for additional information see libFlow library source)
            if ( Data->GetDataError()==-1 ) {
                Abort_OpenHyperFLOW2D();
            }
            unsigned int NumArea = Data->GetIntVal((char*)"NumArea");       // Number of "Area" objects 
            if ( Data->GetDataError()==-1 ) {
                Abort_OpenHyperFLOW2D();
            }

            NumContour = Data->GetIntVal((char*)"NumContour");              // Number of "Contour" objects 
            if ( Data->GetDataError()==-1 ) {
                Abort_OpenHyperFLOW2D();
            }
//---------------------    File Names   ------------------------------
            ProjectName = Data->GetStringVal((char*)"ProjectName");         // Project Name
            if ( Data->GetDataError()==-1 ) Abort_OpenHyperFLOW2D();
            sprintf(GasSwapFileName,"%s%s",ProjectName,Data->GetStringVal((char*)"GasSwapFile")) ; // Swap File name for gas area...
            if ( Data->GetDataError()==-1 ) Abort_OpenHyperFLOW2D();
            sprintf(OutFileName,"%s%s",ProjectName,Data->GetStringVal((char*)"OutputFile"));
            if ( Data->GetDataError()==-1 ) Abort_OpenHyperFLOW2D();
            sprintf(TecPlotFileName,"tp-%s",OutFileName);
            sprintf(ErrFileName,"%s%s",ProjectName,Data->GetStringVal((char*)"ErrorFile"));
            if ( Data->GetDataError()==-1 ) Abort_OpenHyperFLOW2D();
//-----------------------------------------------------------------------
            Ts0 = Data->GetFloatVal((char*)"Ts0");
            if ( Data->GetDataError()==-1 ) Abort_OpenHyperFLOW2D();

            isOutHeatFluxX = Data->GetIntVal((char*)"isOutHeatFluxX");
            if ( Data->GetDataError()==-1 ) Abort_OpenHyperFLOW2D();
            if(isOutHeatFluxX) { 
                Cp_Flow_index = Data->GetIntVal((char*)"Cp_Flow_Index");
                if(Data->GetDataError()==-1) Abort_OpenHyperFLOW2D();
                y_max = Data->GetIntVal((char*)"y_max");
                if(Data->GetDataError()==-1) Abort_OpenHyperFLOW2D();
                y_min = Data->GetIntVal((char*)"y_min");
                if(Data->GetDataError()==-1) Abort_OpenHyperFLOW2D();
            }

            isOutHeatFluxY = Data->GetIntVal((char*)"isOutHeatFluxY");
            if ( Data->GetDataError()==-1 ) Abort_OpenHyperFLOW2D();

            is_p_asterisk_out  = Data->GetIntVal((char*)"is_p_asterisk_out");
            if ( Data->GetDataError()==-1 ) Abort_OpenHyperFLOW2D();

            // Clear Flow list
            if ( FlowList ) {
                for (int i=0;i<(int)FlowList->GetNumElements();i++ ) {
                    delete  FlowList->GetElement(i);
                }
                delete FlowList;
                FlowList=NULL;
            }

            // Load Flow list
            FlowList=new UArray<Flow*>();
            for (int i=0;i<(int)NumFlow;i++ ) {
                snprintf(FlowStr,256,"Flow%i.p",i+1);
                Pg = Data->GetFloatVal(FlowStr);
                if ( Data->GetDataError()==-1 ) {
                    Abort_OpenHyperFLOW2D();
                }
                snprintf(FlowStr,256,"Flow%i.T",i+1);
                Tg = Data->GetFloatVal(FlowStr);
                if ( Data->GetDataError()==-1 ) {
                    Abort_OpenHyperFLOW2D();
                }
                snprintf(FlowStr,256,"Flow%i.CompIndex",i+1);
                CompIndex = Data->GetIntVal(FlowStr);
                if ( Data->GetDataError()==-1 ) {
                     Abort_OpenHyperFLOW2D();
                }
                    if ( CompIndex==0 ) {   // Fuel
                        Cp =chemical_reactions.Cp_Fuel->GetVal(Tg);
                        lam=chemical_reactions.lam_Fuel->GetVal(Tg);
                        mu =chemical_reactions.mu_Fuel->GetVal(Tg);
                        Rg =chemical_reactions.R_Fuel;
                    } else if ( CompIndex==1 ) {  // ox
                        Cp =chemical_reactions.Cp_OX->GetVal(Tg);
                        lam=chemical_reactions.lam_OX->GetVal(Tg);
                        mu =chemical_reactions.mu_OX->GetVal(Tg);
                        Rg =chemical_reactions.R_OX;
                    } else if ( CompIndex==2 ) {  // Combustion products
                        Cp =chemical_reactions.Cp_cp->GetVal(Tg);
                        lam=chemical_reactions.lam_cp->GetVal(Tg);
                        mu =chemical_reactions.mu_cp->GetVal(Tg);
                        Rg =chemical_reactions.R_cp;
                    } else if ( CompIndex==3 ) {  // Air
                        Cp =chemical_reactions.Cp_air->GetVal(Tg);
                        lam=chemical_reactions.lam_air->GetVal(Tg);
                        mu =chemical_reactions.mu_air->GetVal(Tg);
                        Rg =chemical_reactions.R_air;
                    }  else if ( CompIndex==4 ) {  // Mixture
                       snprintf(FlowStr,256,"Flow%i.Y_fuel",i+1);
                       if ( Data->GetDataError()==-1 ) {
                            Abort_OpenHyperFLOW2D();
                       } 
                       Y_mix[0] =  Data->GetFloatVal(FlowStr);
                       snprintf(FlowStr,256,"Flow%i.Y_ox",i+1);
                       if ( Data->GetDataError()==-1 ) {
                            Abort_OpenHyperFLOW2D();
                       } 
                       Y_mix[1] =  Data->GetFloatVal(FlowStr);
                       snprintf(FlowStr,256,"Flow%i.Y_cp",i+1);
                       if ( Data->GetDataError()==-1 ) {
                            Abort_OpenHyperFLOW2D();
                       } 
                       Y_mix[2] =  Data->GetFloatVal(FlowStr);
                       snprintf(FlowStr,256,"Flow%i.Y_air",i+1);
                       if ( Data->GetDataError()==-1 ) {
                            Abort_OpenHyperFLOW2D();
                       } 
                       Y_mix[3] = 1 - Y_mix[0] + Y_mix[1] + Y_mix[2];
                       Cp = Y_mix[0]*chemical_reactions.Cp_Fuel->GetVal(Tg) + Y_mix[1]*chemical_reactions.Cp_OX->GetVal(Tg) + Y_mix[2]*chemical_reactions.Cp_cp->GetVal(Tg)+Y_mix[3]*chemical_reactions.Cp_air->GetVal(Tg);
                       lam= Y_mix[0]*chemical_reactions.lam_Fuel->GetVal(Tg)+ Y_mix[1]*chemical_reactions.lam_OX->GetVal(Tg)+ Y_mix[2]*chemical_reactions.lam_cp->GetVal(Tg)+Y_mix[3]*chemical_reactions.lam_air->GetVal(Tg);
                       mu = Y_mix[0]*chemical_reactions.mu_Fuel->GetVal(Tg) + Y_mix[1]*chemical_reactions.mu_OX->GetVal(Tg) + Y_mix[2]*chemical_reactions.mu_cp->GetVal(Tg)+Y_mix[3]*chemical_reactions.mu_air->GetVal(Tg);
                       Rg = Y_mix[0]*chemical_reactions.R_Fuel + Y_mix[1]*chemical_reactions.R_OX + Y_mix[2]*chemical_reactions.R_cp+Y_mix[3]*chemical_reactions.R_air;
                    } else {
                    *f_stream << "\n";
                    *f_stream << "Bad component index \""<< CompIndex <<"\" use in Flow"<< i+1 <<"\n" << flush;
                    isRun=0;
                    f_stream->flush();
                    Abort_OpenHyperFLOW2D();
                }

                TmpFlow = new Flow(Cp,Tg,Pg,Rg,lam,mu);

                snprintf(FlowStr,256,"Flow%i.Type",i+1);
                FP  LamF;
                int F_Type  = Data->GetIntVal(FlowStr);
                if ( (int)Data->GetDataError()==-1 ) {
                    Abort_OpenHyperFLOW2D();
                }
                if ( F_Type == 0 ) {
                    snprintf(FlowStr,256,"Flow%i.Lam",i+1);
                    LamF = Data->GetFloatVal(FlowStr);
                    if ( Data->GetDataError()==-1 ) {
                        Abort_OpenHyperFLOW2D();
                    }
                    TmpFlow->LAM(LamF);
                } else {
                    snprintf(FlowStr,256,"Flow%i.W",i+1);
                    Wg = Data->GetFloatVal(FlowStr);
                    if ( Data->GetDataError()==-1 ) {
                        Abort_OpenHyperFLOW2D();
                    }
                    TmpFlow->Wg(Wg);
                }
                FlowList->AddElement(&TmpFlow);
                *f_stream << "Add object \"Flow" << i+1 <<  "\"...OK\n" << flush;
                f_stream->flush();
            }

            // Clear Flow2D list
            if ( Flow2DList ) {
                for (int i=0;i<(int)Flow2DList->GetNumElements();i++ ) {
                    delete  Flow2DList->GetElement(i);
                }
                delete Flow2DList;
                Flow2DList=NULL;
            }

            // Load Flow2D list
            Flow2DList=new UArray<Flow2D*>();
            for (int i=0;i<(int)NumFlow2D;i++ ) {
                snprintf(FlowStr,256,"Flow2D-%i.CompIndex",i+1);
                CompIndex = Data->GetIntVal(FlowStr);
                if ( Data->GetDataError()==-1 ) {
                    Abort_OpenHyperFLOW2D();
                }

                snprintf(FlowStr,256,"Flow2D-%i.p",i+1);
                Pg = Data->GetFloatVal(FlowStr);
                if ( Data->GetDataError()==-1 ) {
                    Abort_OpenHyperFLOW2D();
                }
                snprintf(FlowStr,256,"Flow2D-%i.T",i+1);
                Tg = Data->GetFloatVal(FlowStr);
                if ( Data->GetDataError()==-1 ) {
                    Abort_OpenHyperFLOW2D();
                }

                if ( CompIndex==0 ) {   // Fuel
                    Cp =chemical_reactions.Cp_Fuel->GetVal(Tg);
                    lam=chemical_reactions.lam_Fuel->GetVal(Tg);
                    mu =chemical_reactions.mu_Fuel->GetVal(Tg);
                    Rg =chemical_reactions.R_Fuel;
                } else if ( CompIndex==1 ) {  // OX
                    Cp =chemical_reactions.Cp_OX->GetVal(Tg);
                    lam=chemical_reactions.lam_OX->GetVal(Tg);
                    mu =chemical_reactions.mu_OX->GetVal(Tg);
                    Rg =chemical_reactions.R_OX;
                } else if ( CompIndex==2 ) {  // Combustion products
                    Cp =chemical_reactions.Cp_cp->GetVal(Tg);
                    lam=chemical_reactions.lam_cp->GetVal(Tg);
                    mu =chemical_reactions.mu_cp->GetVal(Tg);
                    Rg =chemical_reactions.R_cp;
                } else if ( CompIndex==3 ) {  // Air
                    Cp =chemical_reactions.Cp_air->GetVal(Tg);
                    lam=chemical_reactions.lam_air->GetVal(Tg);
                    mu =chemical_reactions.mu_air->GetVal(Tg);
                    Rg =chemical_reactions.R_air;
                }  else if ( CompIndex==4 ) {  // Mixture
                    snprintf(FlowStr,256,"Flow2D-%i.Y_fuel",i+1);
                    if ( Data->GetDataError()==-1 ) {
                         Abort_OpenHyperFLOW2D();
                    } 
                    Y_mix[0] =  Data->GetFloatVal(FlowStr);
                    snprintf(FlowStr,256,"Flow2D-%i.Y_ox",i+1);
                    if ( Data->GetDataError()==-1 ) {
                         Abort_OpenHyperFLOW2D();
                    } 
                    Y_mix[1] =  Data->GetFloatVal(FlowStr);
                    snprintf(FlowStr,256,"Flow2D-%i.Y_cp",i+1);
                    if ( Data->GetDataError()==-1 ) {
                         Abort_OpenHyperFLOW2D();
                    } 
                    Y_mix[2] =  Data->GetFloatVal(FlowStr);
                    snprintf(FlowStr,256,"Flow2D-%i.Y_air",i+1);
                    if ( Data->GetDataError()==-1 ) {
                         Abort_OpenHyperFLOW2D();
                    }
                    Y_mix[3] = 1 - Y_mix[0] + Y_mix[1] + Y_mix[2];
                    Cp = Y_mix[0]*chemical_reactions.Cp_Fuel->GetVal(Tg) + Y_mix[1]*chemical_reactions.Cp_OX->GetVal(Tg) + Y_mix[2]*chemical_reactions.Cp_cp->GetVal(Tg)+Y_mix[3]*chemical_reactions.Cp_air->GetVal(Tg);
                    lam= Y_mix[0]*chemical_reactions.lam_Fuel->GetVal(Tg)+ Y_mix[1]*chemical_reactions.lam_OX->GetVal(Tg)+ Y_mix[2]*chemical_reactions.lam_cp->GetVal(Tg)+Y_mix[3]*chemical_reactions.lam_air->GetVal(Tg);
                    mu = Y_mix[0]*chemical_reactions.mu_Fuel->GetVal(Tg) + Y_mix[1]*chemical_reactions.mu_OX->GetVal(Tg) + Y_mix[2]*chemical_reactions.mu_cp->GetVal(Tg)+Y_mix[3]*chemical_reactions.mu_air->GetVal(Tg);
                    Rg = Y_mix[0]*chemical_reactions.R_Fuel + Y_mix[1]*chemical_reactions.R_OX + Y_mix[2]*chemical_reactions.R_cp + Y_mix[3]*chemical_reactions.R_air;
                }else {
                    *f_stream << "\n";
                    *f_stream << "Bad component index \""<< CompIndex <<"\" use in Flow2D-"<< i+1 <<"\n" << flush;
                    f_stream->flush();
                    Abort_OpenHyperFLOW2D();
                }

                snprintf(FlowStr,256,"Flow2D-%i.U",i+1);
                Ug = Data->GetFloatVal(FlowStr);
                if ( Data->GetDataError()==-1 ) {
                    Abort_OpenHyperFLOW2D();
                }
                snprintf(FlowStr,256,"Flow2D-%i.V",i+1);
                Vg = Data->GetFloatVal(FlowStr);
                if ( Data->GetDataError()==-1 ) {
                    Abort_OpenHyperFLOW2D();
                }

                snprintf(FlowStr,256,"Flow2D-%i.Mode",i+1); 
                int FlowMode = Data->GetIntVal(FlowStr);
                if ( Data->GetDataError()==-1 ) {
                    Abort_OpenHyperFLOW2D();
                }
                // 0 - U + V + static p, T
                // 1 - U + V + total p*, T*
                // 2 - Mach + Angle + static p*, T*
                // 3 - Mach + Angle + total p, T

                if(FlowMode == 2) {
                  Ug = Vg = 0;
                }

                TmpFlow2D = new Flow2D(mu,lam,Cp,Tg,Pg,Rg,Ug,Vg);

                if(FlowMode == 0) {
                    TmpFlow2D->CorrectFlow(Tg,Pg,sqrt(Ug*Ug+Vg*Vg+1.e-30),FV_VELOCITY);

                } if(FlowMode == 2 || FlowMode == 3) {
                    snprintf(FlowStr,256,"Flow2D-%i.Mach",i+1);
                    FP Mach = Data->GetFloatVal(FlowStr);
                    if ( Data->GetDataError()==-1 ) {
                        Abort_OpenHyperFLOW2D();
                    }

                    snprintf(FlowStr,256,"Flow2D-%i.Angle",i+1); // Angle (deg.)
                    FP Angle = Data->GetFloatVal(FlowStr);
                    if ( Data->GetDataError()==-1 ) {
                        Abort_OpenHyperFLOW2D();
                    }

                    if(FlowMode == 2 ) {
                       TmpFlow2D->CorrectFlow(Tg,Pg,Mach);
                    }

                    TmpFlow2D->MACH(Mach);

                    Wg = TmpFlow2D->Wg();
                    Ug = cos(Angle*M_PI/180)*Wg; 
                    Vg = sin(Angle*M_PI/180)*Wg;
                    TmpFlow2D->Wg(Ug,Vg);
                }

                Flow2DList->AddElement(&TmpFlow2D);
                *f_stream << "Add object \"Flow2D-" << i+1 << " Mach=" << TmpFlow2D->MACH()
                                                           << " U="    << TmpFlow2D->U()   << " m/sec"
                                                           << " V="    << TmpFlow2D->V()   << " m/sec"
                                                           << " Wg="   << TmpFlow2D->Wg()  << " m/sec"
                                                           << " T="    << TmpFlow2D->Tg()  << " K"
                                                           << " p="    << TmpFlow2D->Pg()  << " Pa"
                                                           << " p*="   << TmpFlow2D->P0()  << " Pa"
                                                           << " T*="   << TmpFlow2D->T0()  << " K\"...OK\n" << flush;
                f_stream->flush(); 
            }

            // Load X cuts list
             XCutArray   = new UArray<XCut>();
             XCut        TmpXCut;
             unsigned int NumXCut = Data->GetIntVal((char*)"NumXCut");
             if ( Data->GetDataError()==-1 ) {
                    Abort_OpenHyperFLOW2D();
             }

            for (int i=0;i<(int)NumXCut;i++ ) {
                snprintf(FlowStr,256,"CutX-%i.x0",i+1);
                TmpXCut.x0 = Data->GetFloatVal(FlowStr);
                if ( Data->GetDataError()==-1 ) {
                    Abort_OpenHyperFLOW2D();
                }

                snprintf(FlowStr,256,"CutX-%i.y0",i+1);
                TmpXCut.y0 = Data->GetFloatVal(FlowStr);
                if ( Data->GetDataError()==-1 ) {
                    Abort_OpenHyperFLOW2D();
                }

                snprintf(FlowStr,256,"CutX-%i.dy",i+1);
                TmpXCut.dy = Data->GetFloatVal(FlowStr);
                if ( Data->GetDataError()==-1 ) {
                    Abort_OpenHyperFLOW2D();
                }
                *f_stream << "Add test XCut No" << i+1 << " X=" << TmpXCut.x0 << " Y=" << TmpXCut.y0 << " dY=" << TmpXCut.dy;
                XCutArray->AddElement(&TmpXCut);
                *f_stream << "...OK\n" << flush;
            }

            FileSizeGas =  MaxX*MaxY*sizeof(FlowNode2D<FP,NUM_COMPONENTS>);


            GasSwapData   = LoadSwapFile2D(GasSwapFileName,
                                           (int)MaxX,
                                           (int)MaxY,
                                           sizeof(FlowNode2D<FP,NUM_COMPONENTS>),
                                           &p_g,
                                           &fd_g,
                                           f_stream);
            PreloadFlag = p_g;

            if ( GasSwapFileName!=NULL )
                if ( GasSwapData == 0 ) {
                    PreloadFlag = 0;
                    *f_stream << "\n";
                    *f_stream << "Error mapping swap file...Start without swap file.\n" << flush;
                    f_stream->flush();
                    if ( GasSwapData != 0 ) {
                        munmap(GasSwapData,FileSizeGas);
                        close(fd_g);
                        unlink(GasSwapFileName);
                    }
                }
#ifdef _DEBUG_0
            ___try {
#endif  // _DEBUG_0
                if ( GasSwapData!=0 ) {
                    *f_stream << "Mapping computation area..." << flush;
                    f_stream->flush();
                    J = new UMatrix2D< FlowNode2D<FP,NUM_COMPONENTS> >((FlowNode2D<FP,NUM_COMPONENTS>*)GasSwapData,MaxX,MaxY);
                    useSwapFile=1;
                    sprintf(OldSwapFileName,"%s",GasSwapFileName);
                    OldSwapData = GasSwapData;
                    OldFileSizeGas = FileSizeGas;
                    Old_fd = fd_g;
                } else {
                    *f_stream << "Allocate computation area..." << flush;
                    f_stream->flush();
                    if ( J ) {
                        delete J;J=NULL;
                    }
                    J = new UMatrix2D< FlowNode2D<FP,NUM_COMPONENTS> >(MaxX,MaxY);
                }
#ifdef _DEBUG_0
            } __except( ComputationalMatrix2D*  m) {  // ExceptLib know bug...
#endif  // _DEBUG_0
                if ( J->GetMatrixState()!=MXS_OK ) {
                    *f_stream << "\n";
                    *f_stream << " Memory allocation error !\n" << flush;
                    f_stream->flush();
                    Abort_OpenHyperFLOW2D();
                }
#ifdef _DEBUG_0
            } __end_except;
#endif  // _DEBUG_0

            *f_stream << "OK\n" << flush;
            f_stream->flush();


            sprintf(ErrorMessage, "\nComputation terminated.\nInternal error in OpenHyperFLOW2D/DEEPS.\n");
            static int   isFlow2D=0;
            static int   is_reset;
            unsigned int ix,iy;

#ifdef _DEBUG_0
            ___try {
#endif // _DEBUG_0

 /* Start SingleBound loading */
 //if(!PreloadFlag) {
   Table*       SingleBoundTable;
   Bound2D*     SingleBound;
   static int   BoundMaterialID;
   unsigned int s_x,s_y,e_x,e_y;
   unsigned int numSingleBounds=Data->GetIntVal((char*)"NumSingleBounds");
                if(Data->GetDataError()==-1) Abort_OpenHyperFLOW2D();
                    for (int i=1;i<(int)numSingleBounds+1;i++ ) {
                            sprintf(NameContour,"SingleBound%i",i);
                            sprintf(NameBound,"%s.Points",NameContour);
                            SingleBoundTable = Data->GetTable(NameBound);
                            if ( Data->GetDataError()==-1 ) Abort_OpenHyperFLOW2D();

                            s_x = max((unsigned int)(SingleBoundTable->GetX(0)/dx),0);
                            s_y = max((unsigned int)(SingleBoundTable->GetY(0)/dy),0);
                            e_x = max((unsigned int)(SingleBoundTable->GetX(1)/dx),0);
                            e_y = max((unsigned int)(SingleBoundTable->GetY(1)/dy),0);

                            sprintf(NameBound,"%s.Cond",NameContour);
                            BoundStr = Data->GetStringVal(NameBound);
                            if ( Data->GetDataError()==-1 ) Abort_OpenHyperFLOW2D();

                            sprintf(NameBound,"%s.TurbulenceModel",NameContour);
                            int TurbMod=Data->GetIntVal(NameBound);
                            if ( Data->GetDataError()==-1 ) Abort_OpenHyperFLOW2D();

                            TmpCT           = CT_NO_COND_2D;
                            TmpTurbulenceCT = TCT_No_Turbulence_2D;

                            if ( TurbMod == 0 )
                                TmpTurbulenceCT = (TurbulenceCondType2D)(TmpTurbulenceCT | TCT_No_Turbulence_2D);
                            else if ( TurbMod == 1 )
                                TmpTurbulenceCT = (TurbulenceCondType2D)(TmpTurbulenceCT | TCT_Integral_Model_2D);
                            else if ( TurbMod == 2 )
                                TmpTurbulenceCT = (TurbulenceCondType2D)(TmpTurbulenceCT | TCT_Prandtl_Model_2D);
                            else if ( TurbMod == 3 )
                                TmpTurbulenceCT = (TurbulenceCondType2D)(TmpTurbulenceCT | TCT_Spalart_Allmaras_Model_2D);
                            else if ( TurbMod == 4 )
                                TmpTurbulenceCT = (TurbulenceCondType2D)(TmpTurbulenceCT | TCT_k_eps_Model_2D);
                            else if ( TurbMod == 5 )
                                TmpTurbulenceCT = (TurbulenceCondType2D)(TmpTurbulenceCT | TCT_Smagorinsky_Model_2D);

                            // Atomic conditions
                            if ( strstr(BoundStr,"CT_Rho_CONST_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_Rho_CONST_2D);
                            if ( strstr(BoundStr,"CT_U_CONST_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_U_CONST_2D);
                            if ( strstr(BoundStr,"CT_V_CONST_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_V_CONST_2D);
                            if ( strstr(BoundStr,"CT_T_CONST_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_T_CONST_2D);
                            if ( strstr(BoundStr,"CT_Y_CONST_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_Y_CONST_2D);
                            if ( strstr(BoundStr,"CT_WALL_LAW_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_WALL_LAW_2D);
                            if ( strstr(BoundStr,"CT_WALL_NO_SLIP_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_WALL_NO_SLIP_2D);
                            if ( strstr(BoundStr,"CT_dRhodx_NULL_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_dRhodx_NULL_2D);
                            if ( strstr(BoundStr,"CT_dUdx_NULL_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_dUdx_NULL_2D);
                            if ( strstr(BoundStr,"CT_dVdx_NULL_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_dVdx_NULL_2D);
                            if ( strstr(BoundStr,"CT_dTdx_NULL_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_dTdx_NULL_2D);
                            if ( strstr(BoundStr,"CT_dYdx_NULL_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_dYdx_NULL_2D);
                            if ( strstr(BoundStr,"CT_dRhody_NULL_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_dRhody_NULL_2D);
                            if ( strstr(BoundStr,"CT_dUdy_NULL_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_dUdy_NULL_2D);
                            if ( strstr(BoundStr,"CT_dVdy_NULL_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_dVdy_NULL_2D);
                            if ( strstr(BoundStr,"CT_dTdy_NULL_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_dTdy_NULL_2D);
                            if ( strstr(BoundStr,"CT_dYdy_NULL_2D_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_dYdy_NULL_2D);
                            if ( strstr(BoundStr,"CT_d2Rhodx2_NULL_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_d2Rhodx2_NULL_2D);
                            if ( strstr(BoundStr,"CT_d2Udx2_NULL_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_d2Udx2_NULL_2D);
                            if ( strstr(BoundStr,"CT_d2Vdx2_NULL_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_d2Vdx2_NULL_2D);
                            if ( strstr(BoundStr,"CT_d2Tdx2_NULL_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_d2Tdx2_NULL_2D);
                            if ( strstr(BoundStr,"CT_d2Ydx2_NULL_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_d2Ydx2_NULL_2D);
                            if ( strstr(BoundStr,"CT_d2Rhody2_NULL_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_d2Rhody2_NULL_2D);
                            if ( strstr(BoundStr,"CT_d2Udy2_NULL_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_d2Udy2_NULL_2D);
                            if ( strstr(BoundStr,"CT_d2Vdy2_NULL_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_d2Vdy2_NULL_2D);
                            if ( strstr(BoundStr,"CT_d2Tdy2_NULL_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_d2Tdy2_NULL_2D);
                            if ( strstr(BoundStr,"CT_d2Ydy2_NULL_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_d2Ydy2_NULL_2D);
                            if ( strstr(BoundStr,"CT_SOLID_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_SOLID_2D);
                            if ( strstr(BoundStr,"CT_BL_REFINEMENT_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_BL_REFINEMENT_2D);
                            if ( strstr(BoundStr,"CT_NONREFLECTED_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_NONREFLECTED_2D);

                            if ( strstr(BoundStr,"TCT_k_eps_Model_2D") )
                                TmpTurbulenceCT = (TurbulenceCondType2D)(TmpTurbulenceCT | TCT_k_eps_Model_2D);
                            else if (strstr(BoundStr,"TCT_Smagorinsky_Model_2D") )
                                TmpTurbulenceCT = (TurbulenceCondType2D)(TmpTurbulenceCT | TCT_Smagorinsky_Model_2D);
                            else if (strstr(BoundStr,"TCT_Spalart_Allmaras_Model_2D") )
                                TmpTurbulenceCT = (TurbulenceCondType2D)(TmpTurbulenceCT | TCT_Spalart_Allmaras_Model_2D);
                            else if (strstr(BoundStr,"TCT_Prandtl_Model_2D") )
                                TmpTurbulenceCT = (TurbulenceCondType2D)(TmpTurbulenceCT | TCT_Prandtl_Model_2D);
                            else if (strstr(BoundStr,"TCT_Integral_Model_2D") )
                                TmpTurbulenceCT = (TurbulenceCondType2D)(TmpTurbulenceCT | TCT_Integral_Model_2D);

                            if (isTurbulenceCond2D(TmpTurbulenceCT,TCT_k_eps_Model_2D) || 
                                isTurbulenceCond2D(TmpTurbulenceCT,TCT_Spalart_Allmaras_Model_2D)) {
                                if ( strstr(BoundStr,"TCT_k_CONST_2D") )
                                    TmpTurbulenceCT = (TurbulenceCondType2D)(TmpTurbulenceCT | TCT_k_CONST_2D);
                                if ( strstr(BoundStr,"TCT_eps_CONST_2D") )
                                    TmpTurbulenceCT = (TurbulenceCondType2D)(TmpTurbulenceCT | TCT_eps_CONST_2D);
                                if ( strstr(BoundStr,"TCT_dkdx_NULL_2D") )
                                    TmpTurbulenceCT = (TurbulenceCondType2D)(TmpTurbulenceCT | TCT_dkdx_NULL_2D);
                                if ( strstr(BoundStr,"TCT_depsdx_NULL_2D") )
                                    TmpTurbulenceCT = (TurbulenceCondType2D)(TmpTurbulenceCT | TCT_depsdx_NULL_2D);
                                if ( strstr(BoundStr,"TCT_dkdy_NULL_2D") )
                                    TmpTurbulenceCT = (TurbulenceCondType2D)(TmpTurbulenceCT | TCT_dkdy_NULL_2D);
                                if ( strstr(BoundStr,"TCT_depsdy_NULL_2D") )
                                    TmpTurbulenceCT = (TurbulenceCondType2D)(TmpTurbulenceCT | TCT_depsdy_NULL_2D);
                                if ( strstr(BoundStr,"TCT_d2kdx2_NULL_2D") )
                                    TmpTurbulenceCT = (TurbulenceCondType2D)(TmpTurbulenceCT | TCT_d2kdx2_NULL_2D);
                                if ( strstr(BoundStr,"TCT_d2epsdx2_NULL_2D") )
                                    TmpTurbulenceCT = (TurbulenceCondType2D)(TmpTurbulenceCT | TCT_d2epsdx2_NULL_2D);
                                if ( strstr(BoundStr,"TCT_d2kdy2_NULL_2D") )
                                    TmpTurbulenceCT = (TurbulenceCondType2D)(TmpTurbulenceCT | TCT_d2kdy2_NULL_2D);
                                if ( strstr(BoundStr,"TCT_d2epsdy2_NULL_2D") )
                                    TmpTurbulenceCT = (TurbulenceCondType2D)(TmpTurbulenceCT | TCT_d2epsdy2_NULL_2D);

                                if ( strstr(BoundStr,"TCT_eps_mud2kdx2_WALL_2D") )
                                    TmpTurbulenceCT = (TurbulenceCondType2D)(TmpTurbulenceCT | TCT_eps_mud2kdx2_WALL_2D);
                                if ( strstr(BoundStr,"TCT_eps_mud2kdy2_WALL_2D") )
                                    TmpTurbulenceCT = (TurbulenceCondType2D)(TmpTurbulenceCT | TCT_eps_mud2kdy2_WALL_2D);
                                if ( strstr(BoundStr,"TCT_eps_Cmk2kXn_WALL_2D") )
                                    TmpTurbulenceCT = (TurbulenceCondType2D)(TmpTurbulenceCT | TCT_eps_Cmk2kXn_WALL_2D);
                            }

                            // Macro conditions
                            if ( strstr(BoundStr,"NT_AX_2D") )
                                TmpCT = (CondType2D)(TmpCT | NT_AX_2D);
                            else if ( strstr(BoundStr,"NT_AY_2D") )
                                TmpCT = (CondType2D)(TmpCT | NT_AY_2D);

                            if ( strstr(BoundStr,"NT_D0X_2D") )
                                TmpCT = (CondType2D)(TmpCT | NT_D0X_2D);
                            if ( strstr(BoundStr,"NT_D0Y_2D") )
                                TmpCT = (CondType2D)(TmpCT | NT_D0Y_2D);
                            if ( strstr(BoundStr,"NT_D2X_2D") )
                                TmpCT = (CondType2D)(TmpCT | NT_D2X_2D);
                            if ( strstr(BoundStr,"NT_D2Y_2D") )
                                TmpCT = (CondType2D)(TmpCT | NT_D2Y_2D);
                            if ( strstr(BoundStr,"NT_WALL_LAW_2D") )
                                TmpCT = (CondType2D)(TmpCT | NT_WALL_LAW_2D);
                            else if ( strstr(BoundStr,"NT_WNS_2D") )
                                TmpCT = (CondType2D)(TmpCT | NT_WNS_2D);
                            if ( strstr(BoundStr,"NT_FC_2D") )
                                TmpCT = (CondType2D)(TmpCT | NT_FC_2D);
                            if ( strstr(BoundStr,"NT_FARFIELD_2D") )
                                TmpCT = (CondType2D)(TmpCT | NT_FARFIELD_2D);
                            if ( strstr(BoundStr,"NT_S_2D") )
                                TmpCT = (CondType2D)(TmpCT | NT_S_2D);
                            if ( strstr(BoundStr,"NT_FALSE_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_NODE_IS_SET_2D);

                            if ( TmpCT==CT_NO_COND_2D ) {
                                *f_stream << "\n";
                                *f_stream << "Unknown condition type "<< BoundStr << " in "<< NameContour <<" \n" << flush;
                                f_stream->flush();
                                if ( GasSwapData!=0 ) {
                                    CloseSwapFile(GasSwapFileName,
                                                  GasSwapData,
                                                  FileSizeGas,
                                                  fd_g,
                                                  1);
                                    useSwapFile=0;
                                } else {
                                    delete J;
                                }
                                J=NULL;
                                Abort_OpenHyperFLOW2D();
                            }
                            // Check Flow2D at first ...
                            sprintf(NameBound,"%s.Flow2D",NameContour);
                            FlowIndex = Data->GetIntVal(NameBound);

                            if ( FlowIndex < 1 ) {
                                *f_stream << "\n";
                                sprintf(NameBound,"%s.Flow",NameContour);
                                FlowIndex = Data->GetIntVal(NameBound);
                                if ( FlowIndex < 1 ) {
                                    *f_stream << "\n";
                                    *f_stream << "Bad Flow index [" << FlowIndex << "]\n"<< flush;
                                    f_stream->flush();
                                    Abort_OpenHyperFLOW2D();
                                }
                                if ( Data->GetDataError()==-1 ) Abort_OpenHyperFLOW2D();
                                pTestFlow = FlowList->GetElement(FlowIndex-1);
                                sprintf(FlowStr,"Flow%i.CompIndex",FlowIndex);
                                isFlow2D=0;
                            } else if (FlowIndex <= Flow2DList->GetNumElements()){
                                pTestFlow2D = Flow2DList->GetElement(FlowIndex-1);
                                sprintf(FlowStr,"Flow2D-%i.CompIndex",FlowIndex);
                                isFlow2D=1;
                            } else {
                                *f_stream << "\n";
                                *f_stream << "Bad Flow index [" << FlowIndex << "]\n"<< flush;
                                Abort_OpenHyperFLOW2D();
                            }

                            CompIndex = Data->GetIntVal(FlowStr);
                            if ( CompIndex==0 )      Y=Y_fuel;
                            else if ( CompIndex==1 ) Y=Y_ox;
                            else if ( CompIndex==2 ) Y=Y_cp;
                            else if ( CompIndex==3 ) Y=Y_air;
                            else if ( CompIndex==4 ) Y=Y_mix;

                            sprintf(NameBound,"%s.isReset",NameContour);
                            is_reset = Data->GetIntVal(NameBound);
                            if(!p_g) { 
                                is_reset      = 1; // If swapfile not exist, bounds reset always
                                isFirstStart  = 1;
                            }

                            sprintf(NameBound,"%s.MaterialID",NameContour);
                            BoundMaterialID = Data->GetIntVal(NameBound);
                            //}
                            *f_stream << "\nAdd object \"SingleBound" << i << "\"  ["<< s_x << ";"<< s_y <<"]" << flush;

                                if(!is_reset)
                                    pTestFlow = pTestFlow2D=NULL;

                                if ( !isFlow2D )
                                    SingleBound = new Bound2D(NameContour,J,s_x,s_y,e_x,e_y,TmpCT, pTestFlow,Y,TmpTurbulenceCT);
                                else
                                    SingleBound = new Bound2D(NameContour,J,s_x,s_y,e_x,e_y,TmpCT, pTestFlow2D,Y,TmpTurbulenceCT);

                                *f_stream << "-["<< e_x << ";"<< e_y <<"]" << flush;
                                if(is_reset)
                                   *f_stream << "...Reset parameters\n" << flush;
                                else
                                   *f_stream << "\n" << flush;
                                SingleBound->SetBound(BoundMaterialID);
                                delete SingleBound;
                                f_stream->flush();
                            }
                /* end SingleBound loading */
                /* Set bounds */
                //if ( !PreloadFlag )
                   *f_stream << "\nNum contours: " <<  NumContour << endl;
                   for (int jc=0;jc<NumContour;jc++ ) {
                       sprintf(NameContour,"Contour%i",jc+1);
                        ContourTable = Data->GetTable(NameContour);
                        if ( Data->GetDataError()==-1 ) Abort_OpenHyperFLOW2D();

                        ix=max((int)(ContourTable->GetX(0)/dx),0);
                        iy=max((int)(ContourTable->GetY(0)/dy-1),0);
                        BC = new BoundContour2D(NameContour,J,ix,iy);

                        *f_stream << "Add object \""<< NameContour << "\"...\n" << flush;
                        f_stream->flush();

                        BoundMaterialID = 0;
                        sprintf(NameBound,"%s.MaterialID",NameContour);
                        BoundMaterialID = Data->GetIntVal(NameBound);

                        for (int i=1;i<(int)ContourTable->GetNumNodes()+1;i++ ) {
                            i_last = ContourTable->GetNumNodes()+1;
                            sprintf(NameBound,"%s.Bound%i.Cond",NameContour,i);
                            BoundStr = Data->GetStringVal(NameBound);
                            if ( Data->GetDataError()==-1 ) Abort_OpenHyperFLOW2D();

                            sprintf(NameBound,"%s.Bound%i.TurbulenceModel",NameContour,i);
                            int TurbMod=Data->GetIntVal(NameBound);
                            if ( Data->GetDataError()==-1 ) Abort_OpenHyperFLOW2D();

                            TmpCT           = CT_NO_COND_2D;
                            TmpTurbulenceCT = TCT_No_Turbulence_2D;

                            if ( TurbMod == 0 )
                                TmpTurbulenceCT = (TurbulenceCondType2D)(TmpTurbulenceCT | TCT_No_Turbulence_2D);
                            else if ( TurbMod == 1 )
                                TmpTurbulenceCT = (TurbulenceCondType2D)(TmpTurbulenceCT | TCT_Integral_Model_2D);
                            else if ( TurbMod == 2 )
                                TmpTurbulenceCT = (TurbulenceCondType2D)(TmpTurbulenceCT | TCT_Prandtl_Model_2D);
                            else if ( TurbMod == 3 )
                                TmpTurbulenceCT = (TurbulenceCondType2D)(TmpTurbulenceCT | TCT_Spalart_Allmaras_Model_2D);
                            else if ( TurbMod == 4 )
                                TmpTurbulenceCT = (TurbulenceCondType2D)(TmpTurbulenceCT | TCT_k_eps_Model_2D);
                            else if ( TurbMod == 5 )
                                TmpTurbulenceCT = (TurbulenceCondType2D)(TmpTurbulenceCT | TCT_Smagorinsky_Model_2D);

                            // Simple conditions
                            if ( strstr(BoundStr,"CT_Rho_CONST_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_Rho_CONST_2D);
                            if ( strstr(BoundStr,"CT_U_CONST_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_U_CONST_2D);
                            if ( strstr(BoundStr,"CT_V_CONST_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_V_CONST_2D);
                            if ( strstr(BoundStr,"CT_T_CONST_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_T_CONST_2D);
                            if ( strstr(BoundStr,"CT_Y_CONST_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_Y_CONST_2D);
                            if ( strstr(BoundStr,"CT_dRhodx_NULL_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_dRhodx_NULL_2D);
                            if ( strstr(BoundStr,"CT_dUdx_NULL_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_dUdx_NULL_2D);
                            if ( strstr(BoundStr,"CT_dVdx_NULL_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_dVdx_NULL_2D);
                            if ( strstr(BoundStr,"CT_dTdx_NULL_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_dTdx_NULL_2D);
                            if ( strstr(BoundStr,"CT_dYdx_NULL_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_dYdx_NULL_2D);
                            if ( strstr(BoundStr,"CT_dRhody_NULL_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_dRhody_NULL_2D);
                            if ( strstr(BoundStr,"CT_dUdy_NULL_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_dUdy_NULL_2D);
                            if ( strstr(BoundStr,"CT_dVdy_NULL_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_dVdy_NULL_2D);
                            if ( strstr(BoundStr,"CT_dTdy_NULL_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_dTdy_NULL_2D);
                            if ( strstr(BoundStr,"CT_dYdy_NULL_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_dYdy_NULL_2D);
                            if ( strstr(BoundStr,"CT_d2Rhodx2_NULL_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_d2Rhodx2_NULL_2D);
                            if ( strstr(BoundStr,"CT_d2Udx2_NULL_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_d2Udx2_NULL_2D);
                            if ( strstr(BoundStr,"CT_d2Vdx2_NULL_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_d2Vdx2_NULL_2D);
                            if ( strstr(BoundStr,"CT_d2Tdx2_NULL_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_d2Tdx2_NULL_2D);
                            if ( strstr(BoundStr,"CT_d2Ydx2_NULL_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_d2Ydx2_NULL_2D);
                            if ( strstr(BoundStr,"CT_d2Rhody2_NULL_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_d2Rhody2_NULL_2D);
                            if ( strstr(BoundStr,"CT_d2Udy2_NULL_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_d2Udy2_NULL_2D);
                            if ( strstr(BoundStr,"CT_d2Vdy2_NULL_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_d2Vdy2_NULL_2D);
                            if ( strstr(BoundStr,"CT_d2Tdy2_NULL_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_d2Tdy2_NULL_2D);
                            if ( strstr(BoundStr,"CT_d2Ydy2_NULL_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_d2Ydy2_NULL_2D);
                            if ( strstr(BoundStr,"CT_SOLID_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_SOLID_2D);
                            if ( strstr(BoundStr,"CT_BL_REFINEMENT_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_BL_REFINEMENT_2D);
                            if ( strstr(BoundStr,"CT_NONREFLECTED_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_NONREFLECTED_2D);


                            if ( strstr(BoundStr,"TCT_k_eps_Model_2D") )
                                TmpTurbulenceCT = (TurbulenceCondType2D)(TmpTurbulenceCT | TCT_k_eps_Model_2D);
                            else if (strstr(BoundStr,"TCT_Smagorinsky_Model_2D") )
                                TmpTurbulenceCT = (TurbulenceCondType2D)(TmpTurbulenceCT | TCT_Smagorinsky_Model_2D);
                            else if ( strstr(BoundStr,"TCT_Spalart_Allmaras_Model_2D") )
                                TmpTurbulenceCT = (TurbulenceCondType2D)(TmpTurbulenceCT | TCT_Spalart_Allmaras_Model_2D);
                            else if (strstr(BoundStr,"TCT_Prandtl_Model_2D") )
                                TmpTurbulenceCT = (TurbulenceCondType2D)(TmpTurbulenceCT | TCT_Prandtl_Model_2D);
                            else if (strstr(BoundStr,"TCT_Integral_Model_2D") )
                                TmpTurbulenceCT = (TurbulenceCondType2D)(TmpTurbulenceCT | TCT_Integral_Model_2D);

                            if (isTurbulenceCond2D(TmpTurbulenceCT,TCT_k_eps_Model_2D) ||
                                isTurbulenceCond2D(TmpTurbulenceCT,TCT_Spalart_Allmaras_Model_2D)) {
                                if ( strstr(BoundStr,"TCT_k_CONST_2D") )
                                    TmpTurbulenceCT = (TurbulenceCondType2D)(TmpTurbulenceCT | TCT_k_CONST_2D);
                                if ( strstr(BoundStr,"TCT_eps_CONST_2D") )
                                    TmpTurbulenceCT = (TurbulenceCondType2D)(TmpTurbulenceCT | TCT_eps_CONST_2D);
                                if ( strstr(BoundStr,"TCT_dkdx_NULL_2D") )
                                    TmpTurbulenceCT = (TurbulenceCondType2D)(TmpTurbulenceCT | TCT_dkdx_NULL_2D);
                                if ( strstr(BoundStr,"TCT_depsdx_NULL_2D") )
                                    TmpTurbulenceCT = (TurbulenceCondType2D)(TmpTurbulenceCT | TCT_depsdx_NULL_2D);
                                if ( strstr(BoundStr,"TCT_dkdy_NULL_2D") )
                                    TmpTurbulenceCT = (TurbulenceCondType2D)(TmpTurbulenceCT | TCT_dkdy_NULL_2D);
                                if ( strstr(BoundStr,"TCT_depsdy_NULL_2D") )
                                    TmpTurbulenceCT = (TurbulenceCondType2D)(TmpTurbulenceCT | TCT_depsdy_NULL_2D);
                                if ( strstr(BoundStr,"TCT_d2kdx2_NULL_2D") )
                                    TmpTurbulenceCT = (TurbulenceCondType2D)(TmpTurbulenceCT | TCT_d2kdx2_NULL_2D);
                                if ( strstr(BoundStr,"TCT_d2epsdx2_NULL_2D") )
                                    TmpTurbulenceCT = (TurbulenceCondType2D)(TmpTurbulenceCT | TCT_d2epsdx2_NULL_2D);
                                if ( strstr(BoundStr,"TCT_d2kdy2_NULL_2D") )
                                    TmpTurbulenceCT = (TurbulenceCondType2D)(TmpTurbulenceCT | TCT_d2kdy2_NULL_2D);
                                if ( strstr(BoundStr,"TCT_d2epsdy2_NULL_2D") )
                                    TmpTurbulenceCT = (TurbulenceCondType2D)(TmpTurbulenceCT | TCT_d2epsdy2_NULL_2D);

                                if ( strstr(BoundStr,"TCT_eps_mud2kdx2_WALL_2D") )
                                    TmpTurbulenceCT = (TurbulenceCondType2D)(TmpTurbulenceCT | TCT_eps_mud2kdx2_WALL_2D);
                                if ( strstr(BoundStr,"TCT_eps_mud2kdy2_WALL_2D") )
                                    TmpTurbulenceCT = (TurbulenceCondType2D)(TmpTurbulenceCT | TCT_eps_mud2kdy2_WALL_2D);
                                if ( strstr(BoundStr,"TCT_eps_Cmk2kXn_WALL_2D") )
                                    TmpTurbulenceCT = (TurbulenceCondType2D)(TmpTurbulenceCT | TCT_eps_Cmk2kXn_WALL_2D);
                            }

                            // Macro conditions
                            if ( strstr(BoundStr,"NT_AX_2D") )
                                TmpCT = (CondType2D)(TmpCT | NT_AX_2D);
                            else if ( strstr(BoundStr,"NT_AY_2D") )
                                TmpCT = (CondType2D)(TmpCT | NT_AY_2D);

                            if ( strstr(BoundStr,"NT_D0X_2D") )
                                TmpCT = (CondType2D)(TmpCT | NT_D0X_2D);
                            if ( strstr(BoundStr,"NT_D0Y_2D") )
                                TmpCT = (CondType2D)(TmpCT | NT_D0Y_2D);
                            if ( strstr(BoundStr,"NT_D2X_2D") )
                                TmpCT = (CondType2D)(TmpCT | NT_D2X_2D);
                            if ( strstr(BoundStr,"NT_D2Y_2D") )
                                TmpCT = (CondType2D)(TmpCT | NT_D2Y_2D);
                            if ( strstr(BoundStr,"NT_WALL_LAW_2D") )
                                TmpCT = (CondType2D)(TmpCT | NT_WALL_LAW_2D);
                            else if ( strstr(BoundStr,"NT_WNS_2D") )
                                TmpCT = (CondType2D)(TmpCT | NT_WNS_2D);
                            if ( strstr(BoundStr,"NT_FC_2D") )
                                TmpCT = (CondType2D)(TmpCT | NT_FC_2D);
                            if ( strstr(BoundStr,"NT_FARFIELD_2D") )
                                TmpCT = (CondType2D)(TmpCT | NT_FARFIELD_2D);
                            if ( strstr(BoundStr,"NT_S_2D") )
                                TmpCT = (CondType2D)(TmpCT | NT_S_2D);
                            if ( strstr(BoundStr,"NT_FALSE_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_NODE_IS_SET_2D);


                            if ( TmpCT==CT_NO_COND_2D  &&  TmpTurbulenceCT == 0) {
                                *f_stream << "\n";
                                *f_stream << "Unknown condition type "<< BoundStr << " in "<< NameBound <<" \n" << flush;
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
                                J=NULL;
                                Abort_OpenHyperFLOW2D();
                            }
                            // Check Flow2D at first ...
                            sprintf(NameBound,"%s.Bound%i.Flow2D",NameContour,i);
                            FlowIndex = Data->GetIntVal(NameBound);

                            if ( FlowIndex < 1 ) {
                                *f_stream << "\n";
                                sprintf(NameBound,"%s.Flow",NameContour);
                                FlowIndex = Data->GetIntVal(NameBound);
                                if ( FlowIndex < 1 ) {
                                    *f_stream << "\n";
                                    *f_stream << "Bad Flow index [" << FlowIndex << "]\n"<< flush;
                                    f_stream->flush();
                                    Abort_OpenHyperFLOW2D();
                                }
                                if ( Data->GetDataError()==-1 ) Abort_OpenHyperFLOW2D();
                                pTestFlow = FlowList->GetElement(FlowIndex-1);
                                sprintf(FlowStr,"Flow%i.CompIndex",FlowIndex);
                                isFlow2D=0;
                            } else if (FlowIndex <= Flow2DList->GetNumElements()){
                                pTestFlow2D = Flow2DList->GetElement(FlowIndex-1);
                                sprintf(FlowStr,"Flow2D-%i.CompIndex",FlowIndex);
                                isFlow2D=1;
                            } else {
                                *f_stream << "\n";
                                *f_stream << "Bad Flow index [" << FlowIndex << "]\n"<< flush;
                                Abort_OpenHyperFLOW2D();
                            }

                            CompIndex = Data->GetIntVal(FlowStr);
                            if ( CompIndex==0 )      Y=Y_fuel;
                            else if ( CompIndex==1 ) Y=Y_ox;
                            else if ( CompIndex==2 ) Y=Y_cp;
                            else if ( CompIndex==3 ) Y=Y_air;
                            else if ( CompIndex==4 ) Y=Y_mix;

                            sprintf(NameBound,"%s.Bound%i.isReset",NameContour,i);
                            is_reset = Data->GetIntVal(NameBound);
                            if(!p_g) is_reset = 1; // If swapfile not exist, bounds reset always
                            *f_stream << "Add object \"Bound" << i << "\"  ["<< BC->GetCurrentX() << ";"<< BC->GetCurrentY() <<"]" << flush;
                            if ( i < (int)ContourTable->GetNumNodes() ) {
                                ix=max((int)(ContourTable->GetX(i)/dx),0);
                                iy=max((int)(ContourTable->GetY(i)/dy-1),0);

                                if(!is_reset)
                                    pTestFlow = pTestFlow2D=NULL;
                                sprintf(NameBound,"%s.Bound%i",NameContour,i);
                                if ( !isFlow2D )
                                    BC->AddBound2D(NameBound,ix,iy,TmpCT,pTestFlow,NULL,Y,TmpTurbulenceCT);
                                else
                                    BC->AddBound2D(NameBound,ix,iy,TmpCT,NULL,pTestFlow2D,Y,TmpTurbulenceCT);

                                *f_stream << "-["<< BC->GetCurrentX() << ";"<< BC->GetCurrentY() <<"]" << flush;

                                if(is_reset)
                                   *f_stream << "...Reset parameters" << flush;
                                *f_stream << "\n";
                                f_stream->flush();
                            }
                        }

                        if(!is_reset)
                            pTestFlow = pTestFlow2D=NULL;

                        sprintf(NameBound,"%s.Bound%i",NameContour,i_last);

                        if ( !isFlow2D )
                            BC->CloseContour2D(NameBound,TmpCT,pTestFlow,NULL,Y,TmpTurbulenceCT);
                        else
                            BC->CloseContour2D(NameBound,TmpCT,NULL,pTestFlow2D,Y,TmpTurbulenceCT);
                        *f_stream << "-["<< BC->GetCurrentX() << ";"<< BC->GetCurrentY() <<"]" << flush;
                        if(is_reset)
                           *f_stream << "...Reset parameters" << flush;
                        *f_stream << "\nEnd bounds..." << flush;
                        f_stream->flush();

                        if ( !BC->IsContourClosed() ) {
                            *f_stream << "\n";
                            *f_stream << "Contour is not looped.\n" << flush;
                            f_stream->flush();
                            Abort_OpenHyperFLOW2D();
                        }
                        *f_stream << "OK\nSet bounds for \""<< NameContour << "\"..." << flush;
                        f_stream->flush();
                        BC->SetBounds(BoundMaterialID);
                        *f_stream << "OK\n" << flush;
                        f_stream->flush();
                        delete BC;
                   }
#ifdef _DEBUG_0
            }__except(int num_bound) {
                *f_stream << "\nERROR:";
                *f_stream << "Set Bound error (bound No. "<< num_bound+1 << ") in \""<< NameContour <<"\"\n" << flush;
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
                J=NULL;
                Abort_OpenHyperFLOW2D();
            }__except(void*) {
                *f_stream << " \nERROR: ";
                *f_stream << "Unknown Bound error in \""<< NameContour <<"\"\n" << flush;
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
                J=NULL;
                Abort_OpenHyperFLOW2D();
            }
            __end_except;
#endif // _DEBUG_0

            dt = 1;

            for (int i = 0;i<(int)FlowList->GetNumElements();i++ ) {
                FP CFL_min  = min(CFL,CFL_Scenario->GetVal(iter+last_iter));
                dt = min(dt,CFL_min*min(dx/(FlowList->GetElement(i)->Asound()+FlowList->GetElement(i)->Wg()),
                                        dy/(FlowList->GetElement(i)->Asound()+FlowList->GetElement(i)->Wg())));
            }

            for (int i = 0;i<(int)Flow2DList->GetNumElements();i++ ) {
                FP CFL_min  = min(CFL,CFL_Scenario->GetVal(iter+last_iter));
                dt = min(dt,CFL_min*min(dx/(Flow2DList->GetElement(i)->Asound()+Flow2DList->GetElement(i)->Wg()),
                                        dy/(Flow2DList->GetElement(i)->Asound()+Flow2DList->GetElement(i)->Wg())));
            }

            if ( !PreloadFlag ) {
                CutFile(OutFileName);  
                *f_stream << "Init computation area..." << flush;
                f_stream->flush();
#ifdef _DEBUG_0
           ___try {
#endif  // _DEBUG_0
                    for (int j=0;j<(int)MaxY;j++ )
                        for (int i=0;i<(int)MaxX;i++ ) {
                          i_err = i;
                          j_err = j;
#ifndef _UNIFORM_MESH_
                          J->GetValue(i,j).dx  = dx;
                          J->GetValue(i,j).dy  = dy;
#endif //_UNIFORM_MESH_

                          //if ( FlowNode2D<FP,NUM_COMPONENTS>::FT == FT_AXISYMMETRIC )
                          //     J->GetValue(i,j).r     = (j+1)*dy;
                          J->GetValue(i,j).x     = (i+0.5)*dx;
                          J->GetValue(i,j).y     = (j+0.5)*dy;
                          J->GetValue(i,j).Tf    = chemical_reactions.Tf;
                          J->GetValue(i,j).BGX   = 1.;
                          J->GetValue(i,j).BGY   = 1.;
                          J->GetValue(i,j).NGX   = 0;
                          J->GetValue(i,j).NGY   = 0;

                          for ( k=0;k<(int)FlowNode2D<FP,NUM_COMPONENTS>::NumEq;k++ )
                              J->GetValue(i,j).Src[k]= J->GetValue(i,j).SrcAdd[k] = 0;
                        }
#ifdef _DEBUG_0
                }__except(SysException e) {
                    //int isDel=0;
                    if ( e == SIGINT ) {
                        *f_stream << "\n";
                        *f_stream << "Interrupted by user\n" <<  flush;
                    } else {
                        *f_stream << "\n";
                        *f_stream << "Error init computation area ("<<i_err<<","<<j_err<<")..." << flush;
                        *f_stream << "\n";
                        *f_stream << "System error " << e << " \n" << flush;
                    }
                    f_stream->flush();
                    if ( GasSwapData!=0 ) {
#ifdef _REMOVE_SWAPFILE_
                        if ( e != SIGINT ) isDel = 1;
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
                      J=NULL;
                    Abort_OpenHyperFLOW2D();
                } __except( ComputationalMatrix2D*  m) {
                    if ( m->GetMatrixState()==MXS_ERR_OUT_OF_INDEX ) {
                        *f_stream << "\n";
                        *f_stream << "Error init computation area..("<<i_err<<","<<j_err<<")..." << flush;
                        *f_stream << "\n";
                        *f_stream << "Matrix out of index error\n" << flush;
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
                        J=NULL;
                        Abort_OpenHyperFLOW2D();
                    }
                } __except(void*) {
                    *f_stream << "\n";
                    *f_stream << "Error init computation area..("<<i_err<<","<<j_err<<")..." << flush;
                    *f_stream << "Unknown error\n" << flush;
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
                    J=NULL;
                    Abort_OpenHyperFLOW2D();
                }
                __end_except;
#endif  // _DEBUG_0

                *f_stream << "OK\n" << flush;
                f_stream->flush();
                /* End set bounds */
            }
            /* -- if use swapfile */

            //Cx,Cy,Cd,Cv calc
            is_Cx_calc = Data->GetIntVal((char*)"is_Cx_calc");
            if(Data->GetDataError()==-1) Abort_OpenHyperFLOW2D();
            
            if (is_Cx_calc) {
             x0_body=Data->GetFloatVal((char*)"x_body");
             if(Data->GetDataError()==-1) Abort_OpenHyperFLOW2D();
             y0_body=Data->GetFloatVal((char*)"y_body");
             if(Data->GetDataError()==-1) Abort_OpenHyperFLOW2D();
             dx_body=Data->GetFloatVal((char*)"dx_body");
             if(Data->GetDataError()==-1) Abort_OpenHyperFLOW2D();
             dy_body=Data->GetFloatVal((char*)"dy_body");
             if(Data->GetDataError()==-1) Abort_OpenHyperFLOW2D();
             Cx_Flow_index = Data->GetIntVal((char*)"Cx_Flow_Index");
            }
            
            is_Cd_calc = Data->GetIntVal((char*)"is_Cd_calc");
            if(Data->GetDataError()==-1) Abort_OpenHyperFLOW2D();
            
            if (is_Cd_calc) {
             x0_nozzle=Data->GetFloatVal((char*)"x_nozzle");
             if(Data->GetDataError()==-1) Abort_OpenHyperFLOW2D();
             y0_nozzle=Data->GetFloatVal((char*)"y_nozzle");
             if(Data->GetDataError()==-1) Abort_OpenHyperFLOW2D();
             dy_nozzle=Data->GetFloatVal((char*)"dy_nozzle");
             if(Data->GetDataError()==-1) Abort_OpenHyperFLOW2D();
             Cd_Flow_index = Data->GetIntVal((char*)"Cd_Flow_Index");
             if(Data->GetDataError()==-1) Abort_OpenHyperFLOW2D();
             p_ambient = Data->GetFloatVal((char*)"p_ambient");
             if(Data->GetDataError()==-1) Abort_OpenHyperFLOW2D();
            }

            //Set bound Primitives
            // Rects
                FP        Xstart,Ystart,X_0,Y_0;
                unsigned int  numRects=Data->GetIntVal((char*)"NumRects");
                if(Data->GetDataError()==-1) Abort_OpenHyperFLOW2D();
                //unsigned int ix0,iy0;
                SolidBoundRect2D* SBR;
                GlobalTime=Data->GetFloatVal((char*)"InitTime");
                if(Data->GetDataError()==-1) Abort_OpenHyperFLOW2D();

                if(numRects)
                  {
                   for(j=0;j<numRects;j++) {

                       TurbulenceCondType2D TM;

                       sprintf(NameContour,"Rect%i",j+1);
                       *f_stream << "Add object \""<< NameContour << "\"..." << flush;

                       sprintf(NameContour,"Rect%i.Xstart",j+1);
                       Xstart=Data->GetFloatVal(NameContour);
                       if(Data->GetDataError()==-1) Abort_OpenHyperFLOW2D();

                       sprintf(NameContour,"Rect%i.Ystart",j+1);
                       Ystart=Data->GetFloatVal(NameContour);
                       if(Data->GetDataError()==-1) Abort_OpenHyperFLOW2D();

                       sprintf(NameContour,"Rect%i.DX",j+1);
                       X_0=Data->GetFloatVal(NameContour);
                       if(Data->GetDataError()==-1) Abort_OpenHyperFLOW2D();

                       sprintf(NameContour,"Rect%i.DY",j+1);
                       Y_0=Data->GetFloatVal(NameContour);
                       if(Data->GetDataError()==-1) Abort_OpenHyperFLOW2D();

                       chord = max(chord,X_0);

                       sprintf(NameContour,"Rect%i.Flow2D",j+1);
                       FlowIndex=Data->GetIntVal(NameContour);
                       if(Data->GetDataError()==-1) Abort_OpenHyperFLOW2D();

                       sprintf(NameContour,"Rect%i.TurbulenceModel",j+1);
                       int TurbMod=Data->GetIntVal(NameContour);
                       if ( Data->GetDataError()==-1 ) Abort_OpenHyperFLOW2D();

                       TM = TCT_No_Turbulence_2D;
                       if ( TurbMod == 0 )
                           TM = (TurbulenceCondType2D)(TM | TCT_No_Turbulence_2D);
                       else if ( TurbMod == 1 )
                           TM = (TurbulenceCondType2D)(TM | TCT_Integral_Model_2D);
                       else if ( TurbMod == 2 )
                           TM = (TurbulenceCondType2D)(TM | TCT_Prandtl_Model_2D);
                       else if ( TurbMod == 3 )
                           TM = (TurbulenceCondType2D)(TM | TCT_Spalart_Allmaras_Model_2D);
                       else if ( TurbMod == 4 )
                           TM = (TurbulenceCondType2D)(TM | TCT_k_eps_Model_2D);
                       else if ( TurbMod == 5 )
                           TM = (TurbulenceCondType2D)(TM | TCT_Smagorinsky_Model_2D);
                       if(FlowIndex < 1 || FlowIndex > Flow2DList->GetNumElements() )
                         {
                           *f_stream << "\nBad Flow index [" << FlowIndex << "]\n"<< flush;
                           f_stream->flush();
                           Abort_OpenHyperFLOW2D();
                         }

                       pTestFlow2D = Flow2DList->GetElement(FlowIndex-1);
                       snprintf(FlowStr,256,"Flow2D-%i.CompIndex",FlowIndex);

                       CompIndex = Data->GetIntVal(FlowStr);
                       if(CompIndex==0)      Y=Y_fuel;
                       else if(CompIndex==1) Y=Y_ox;
                       else if(CompIndex==2) Y=Y_cp;
                       else if(CompIndex==3) Y=Y_air;
                       else if(CompIndex==4) Y=Y_mix;
                       else {
                           *f_stream << "\nBad component index [" << CompIndex << "]\n"<< flush;
                           f_stream->flush();
                           Abort_OpenHyperFLOW2D();
                       }
                       sprintf(NameContour,"Rect%i",j+1);

                       if(!PreloadFlag) {
                         SBR = new SolidBoundRect2D(NameContour,J,Xstart,Ystart,X_0,Y_0,dx,dy,(CondType2D)NT_WNS_2D,pTestFlow2D,Y,TM);
                          delete SBR;
                       }

                       *f_stream << "OK\n" << flush;
                     }
                  }
            // Bound Circles
            unsigned int numCircles=Data->GetIntVal((char*)"NumCircles");
            TurbulenceCondType2D TM;
            BoundCircle2D* SBC;
                if ( numCircles ) {
                    for ( j=0;j<numCircles;j++ ) {

                        sprintf(NameContour,"Circle%i",j+1);
                        *f_stream << "Add object \""<< NameContour << "\"..." << flush;

                        sprintf(NameContour,"Circle%i.Xstart",j+1);
                        Xstart=Data->GetFloatVal(NameContour);
                        if ( Data->GetDataError()==-1 ) Abort_OpenHyperFLOW2D();

                        sprintf(NameContour,"Circle%i.Ystart",j+1);
                        Ystart=Data->GetFloatVal(NameContour);
                        if ( Data->GetDataError()==-1 ) Abort_OpenHyperFLOW2D();

                        sprintf(NameContour,"Circle%i.X0",j+1);
                        X_0=Data->GetFloatVal(NameContour);
                        if ( Data->GetDataError()==-1 ) Abort_OpenHyperFLOW2D();

                        sprintf(NameContour,"Circle%i.Y0",j+1);
                        Y_0=Data->GetFloatVal(NameContour);
                        if ( Data->GetDataError()==-1 ) Abort_OpenHyperFLOW2D();

                        chord = max(chord,sqrt((X_0-Xstart)*(X_0-Xstart) + (Y_0-Ystart*(Y_0-Ystart) )+1.e-60));

                        sprintf(NameContour,"Circle%i.MaterialID",j+1);
                        int MaterialID=Data->GetIntVal(NameContour);
                        if ( Data->GetDataError()==-1 ) Abort_OpenHyperFLOW2D();

                        sprintf(NameContour,"Circle%i.TurbulenceModel",j+1);
                        int TurbMod=Data->GetIntVal(NameContour);
                        if ( Data->GetDataError()==-1 ) Abort_OpenHyperFLOW2D();

                        TM = TCT_No_Turbulence_2D;
                        if ( TurbMod == 0 )
                            TM = (TurbulenceCondType2D)(TM | TCT_No_Turbulence_2D);
                        else if ( TurbMod == 1 )
                            TM = (TurbulenceCondType2D)(TM | TCT_Integral_Model_2D);
                        else if ( TurbMod == 2 )
                            TM = (TurbulenceCondType2D)(TM | TCT_Prandtl_Model_2D);
                        else if ( TurbMod == 3 )
                            TM = (TurbulenceCondType2D)(TM | TCT_Spalart_Allmaras_Model_2D);
                        else if ( TurbMod == 4 )
                            TM = (TurbulenceCondType2D)(TM | TCT_k_eps_Model_2D);
                        else if ( TurbMod == 5 )
                            TM = (TurbulenceCondType2D)(TM | TCT_Smagorinsky_Model_2D);

                        sprintf(NameContour,"Circle%i.Flow2D",j+1);
                        FlowIndex=Data->GetIntVal(NameContour);
                        if ( Data->GetDataError()==-1 ) Abort_OpenHyperFLOW2D();

                        if ( FlowIndex < 1 || FlowIndex > Flow2DList->GetNumElements() ) {
                            *f_stream << "\n";
                            *f_stream << "Bad Flow index [" << FlowIndex << "]\n"<< flush;
                            f_stream->flush();
                            Abort_OpenHyperFLOW2D();
                        }

                        pTestFlow2D = Flow2DList->GetElement(FlowIndex-1);
                        snprintf(FlowStr,256,"Flow2D-%i.CompIndex",FlowIndex);

                        CompIndex = Data->GetIntVal(FlowStr);
                        if ( CompIndex==0 )      Y=Y_fuel;
                        else if ( CompIndex==1 ) Y=Y_ox;
                        else if ( CompIndex==2 ) Y=Y_cp;
                        else if ( CompIndex==3 ) Y=Y_air;
                        else if(CompIndex==4) Y=Y_mix;
                        else {
                            *f_stream << "\n";
                            *f_stream << "Bad component index [" << CompIndex << "]\n"<< flush;
                            f_stream->flush();
                            Abort_OpenHyperFLOW2D();
                        }
                        sprintf(NameContour,"Circle%i",j+1);
                        if( !PreloadFlag ) {
                        if(MaterialID) { // Solid
                            SBC = new BoundCircle2D(NameContour,J,Xstart,Ystart,X_0,Y_0,dx,dy,(CondType2D)NT_WNS_2D,MaterialID,pTestFlow2D,Y,TM
#ifdef _DEBUG_1
                                                    ,f_stream
#endif //_DEBUG_1 
                                                    );
                        } else {         // GAS
                            SBC = new BoundCircle2D(NameContour,J,Xstart,Ystart,X_0,Y_0,dx,dy,(CondType2D)CT_NODE_IS_SET_2D,MaterialID,pTestFlow2D,Y,TM
#ifdef _DEBUG_1
                                                    ,f_stream
#endif //_DEBUG_1 
                                                    );
                        }
                        delete SBC;
                        }
                        *f_stream << "OK\n" << flush;
                    }
                }
// Solid Bound Airfoils
            unsigned int  numAirfoils=Data->GetIntVal((char*)"NumAirfoils");
                if ( numAirfoils ) {
                    FP            mm,pp,thick,scale,attack_angle;
                    int           Airfoil_Type;
                    InputData*    AirfoilInputData;
                    SolidBoundAirfoil2D* SBA;
                    for ( j=0;j<numAirfoils;j++ ) {
                        sprintf(NameContour,"Airfoil%i",j+1);    
                        *f_stream << "Add object \""<< NameContour << "\"..." << flush;

                        sprintf(NameContour,"Airfoil%i.Xstart",j+1);
                        Xstart=Data->GetFloatVal(NameContour);
                        if ( Data->GetDataError()==-1 ) Abort_OpenHyperFLOW2D();

                        sprintf(NameContour,"Airfoil%i.Ystart",j+1);
                        Ystart=Data->GetFloatVal(NameContour);
                        if ( Data->GetDataError()==-1 ) Abort_OpenHyperFLOW2D();


                        sprintf(NameContour,"Airfoil%i.Type",j+1);
                        Airfoil_Type=Data->GetFloatVal(NameContour);
                        if ( Data->GetDataError()==-1 ) Abort_OpenHyperFLOW2D();

                        if (Airfoil_Type == 0) { // NACA Airfoil
                            sprintf(NameContour,"Airfoil%i.pp",j+1);
                            pp=Data->GetFloatVal(NameContour);
                            if ( Data->GetDataError()==-1 ) Abort_OpenHyperFLOW2D();

                            sprintf(NameContour,"Airfoil%i.mm",j+1);
                            mm=Data->GetFloatVal(NameContour);
                            if ( Data->GetDataError()==-1 ) Abort_OpenHyperFLOW2D();

                            sprintf(NameContour,"Airfoil%i.thick",j+1);
                            thick=Data->GetFloatVal(NameContour);
                            if ( Data->GetDataError()==-1 ) Abort_OpenHyperFLOW2D();
                        } else { // TsAGI Airfoil
                           char* AirfoilInputDataFileName;
                           sprintf(NameContour,"Airfoil%i.InputData",j+1);
                           AirfoilInputDataFileName=Data->GetStringVal(NameContour);
                           if ( Data->GetDataError()==-1 ) Abort_OpenHyperFLOW2D();
                           AirfoilInputData = new InputData(AirfoilInputDataFileName,DS_FILE,f_stream);
                        }

                        sprintf(NameContour,"Airfoil%i.scale",j+1);
                        scale=Data->GetFloatVal(NameContour);
                        if ( Data->GetDataError()==-1 ) Abort_OpenHyperFLOW2D();

                        chord = max(chord,scale);

                        sprintf(NameContour,"Airfoil%i.attack_angle",j+1);    
                        attack_angle=Data->GetFloatVal(NameContour);
                        if ( Data->GetDataError()==-1 ) Abort_OpenHyperFLOW2D();

                        sprintf(NameContour,"Airfoil%i.Flow2D",j+1);    
                        FlowIndex=Data->GetIntVal(NameContour);
                        if ( Data->GetDataError()==-1 ) Abort_OpenHyperFLOW2D();

                        sprintf(NameContour,"Airfoil%i.TurbulenceModel",j+1);    
                        int TurbMod=Data->GetIntVal(NameContour);
                        if ( Data->GetDataError()==-1 ) Abort_OpenHyperFLOW2D();

                        TM = TCT_No_Turbulence_2D;
                        if ( TurbMod == 0 )
                            TM = (TurbulenceCondType2D)(TM | TCT_No_Turbulence_2D);
                        else if ( TurbMod == 1 )
                            TM = (TurbulenceCondType2D)(TM | TCT_Integral_Model_2D);
                        else if ( TurbMod == 2 )
                            TM = (TurbulenceCondType2D)(TM | TCT_Prandtl_Model_2D);
                        else if ( TurbMod == 3 )
                            TM = (TurbulenceCondType2D)(TM | TCT_Spalart_Allmaras_Model_2D);
                        else if ( TurbMod == 4 )
                            TM = (TurbulenceCondType2D)(TM | TCT_k_eps_Model_2D);
                        else if ( TurbMod == 5 )
                            TM = (TurbulenceCondType2D)(TM | TCT_Smagorinsky_Model_2D);

                        if ( FlowIndex < 1 || FlowIndex > Flow2DList->GetNumElements()) {
                            *f_stream << "\n";
                            *f_stream << "Bad Flow index [" << FlowIndex << "]\n"<< flush;
                            f_stream->flush();
                            Abort_OpenHyperFLOW2D();
                        }

                        pTestFlow2D = Flow2DList->GetElement(FlowIndex-1);
                        snprintf(FlowStr,256,"Flow2D-%i.CompIndex",FlowIndex);

                        CompIndex = Data->GetIntVal(FlowStr);
                        if ( CompIndex==0 )      Y=Y_fuel;
                        else if ( CompIndex==1 ) Y=Y_ox;
                        else if ( CompIndex==2 ) Y=Y_cp;
                        else if ( CompIndex==3 ) Y=Y_air;
                        else if ( CompIndex==4 ) Y=Y_mix;
                        else {
                            *f_stream << "\n";
                            *f_stream << "Bad component index [" << CompIndex << "]\n"<< flush;
                            f_stream->flush();
                            Abort_OpenHyperFLOW2D();
                        }

                        *f_stream  << "Airfoil Re = "  << Re_Airfoil(chord,Flow2DList->GetElement(Cx_Flow_index-1)) << endl;

                        sprintf(NameContour,"Airfoil%i",j+1);
                     if (!PreloadFlag) {
                        if (Airfoil_Type == 0) { // NACA Airfoil
                           SBA = new SolidBoundAirfoil2D(NameContour,J,Xstart,Ystart,mm,pp,thick,dx,dy,(CondType2D)NT_WNS_2D,pTestFlow2D,Y,TM,scale,attack_angle,f_stream);
                        } else {                // TsAGI Airfoil
                           SBA = new SolidBoundAirfoil2D(NameContour,J,Xstart,Ystart,AirfoilInputData,dx,dy,(CondType2D)NT_WNS_2D,pTestFlow2D,Y,TM,scale,attack_angle,f_stream);
                        }
                        delete SBA;
                     }
                        *f_stream << "OK\n" << flush;
                  }
                }
               //  Areas
            if ( !PreloadFlag ) {
#ifdef _DEBUG_0
           ___try {
#endif  // _DEBUG_0
                    static Area2D*  TmpArea=NULL; 
                    static char     AreaName[256];
                    static int      AreaType;
                    static Table*   AreaPoint=NULL;
                    static int      AreaMaterialID;

                    for (int i=0;i<(int)NumArea;i++ ) {
                        snprintf(AreaName,256,"Area%i",i+1);

                        TmpArea = new Area2D(AreaName,J);
                        pTestFlow = pTestFlow2D = NULL;

                        *f_stream << "Add object \""<< AreaName << "\":" << flush;
                        f_stream->flush();
                        AreaPoint = Data->GetTable(AreaName);
                        if ( Data->GetDataError()==-1 ) {
                            Abort_OpenHyperFLOW2D();
                        }
                        snprintf(AreaName,256,"Area%i.Type",i+1);
                        AreaType = Data->GetIntVal(AreaName);
                        if ( Data->GetDataError()==-1 ) {
                            Abort_OpenHyperFLOW2D();
                        }
                        snprintf(AreaName,256,"initial Area point (%d,%d)...",(unsigned int)AreaPoint->GetX(0),
                                                                              (unsigned int)AreaPoint->GetY(0));
                        *f_stream << AreaName; 

                        snprintf(AreaName,256,"Area%i.MaterialID",i+1);
                        AreaMaterialID = Data->GetIntVal(AreaName);

                        if ( AreaType==0 ) {
                            if ( Data->GetDataError()==-1 ) {
                                Abort_OpenHyperFLOW2D();
                            }
                            TmpArea->FillArea2D((unsigned int)AreaPoint->GetX(0),
                                                (unsigned int)AreaPoint->GetY(0),
                                                CT_SOLID_2D,
                                                TCT_No_Turbulence_2D,
                                                AreaMaterialID);
                        } else if ( AreaType==1 ) {
                            //AreaMaterialID = GAS_ID;
                            snprintf(AreaName,256,"Area%i.Flow2D",i+1);
                            FlowIndex = Data->GetIntVal(AreaName);

                            if ( FlowIndex < 1 || FlowIndex > Flow2DList->GetNumElements() ) {
                                *f_stream << "\n";
                                *f_stream << "Bad Flow index [" << FlowIndex << "] \n"<< flush;
                                f_stream->flush();
                                Abort_OpenHyperFLOW2D();
                            }

                            if ( Data->GetDataError()==-1 ) {
                                snprintf(AreaName,256,"Area%i.Flow",i+1);
                                FlowIndex = Data->GetIntVal(AreaName);
                                if ( FlowIndex < 1 || FlowIndex > Flow2DList->GetNumElements() ) {
                                    *f_stream << "\n";
                                    *f_stream << "Bad Flow index [" << FlowIndex << "]\n"<< flush;
                                    f_stream->flush();
                                    Abort_OpenHyperFLOW2D();
                                }
                                if ( Data->GetDataError()==-1 ) {
                                    Abort_OpenHyperFLOW2D();
                                }
                                pTestFlow = FlowList->GetElement(FlowIndex-1);
                            } else {
                                pTestFlow2D = Flow2DList->GetElement(FlowIndex-1);
                            }

                            if ( pTestFlow )
                                snprintf(FlowStr,256,"Flow%i.CompIndex",FlowIndex);
                            else if ( pTestFlow2D )
                                snprintf(FlowStr,256,"Flow2D-%i.CompIndex",FlowIndex);

                            CompIndex = Data->GetIntVal(FlowStr);
                            if ( Data->GetDataError()==-1 ) {
                                Abort_OpenHyperFLOW2D();
                            }

                            if ( CompIndex==0 )      Y=Y_fuel;
                            else if ( CompIndex==1 ) Y=Y_ox;
                            else if ( CompIndex==2 ) Y=Y_cp;
                            else if ( CompIndex==3 ) Y=Y_air;
                            else if ( CompIndex==4 ) Y=Y_mix;
                            else {
                                *f_stream << "Bad component index [" << CompIndex << "]\n" << flush;
                                Abort_OpenHyperFLOW2D();
                            }

                            sprintf(AreaName,"Area%i.TurbulenceModel",i+1);    
                            int TurbMod=Data->GetIntVal(AreaName);
                            if ( Data->GetDataError()==-1 ) Abort_OpenHyperFLOW2D();

                            TM = TCT_No_Turbulence_2D;
                            if ( TurbMod == 0 )
                                TM = (TurbulenceCondType2D)(TM | TCT_No_Turbulence_2D);
                            else if ( TurbMod == 1 )
                                TM = (TurbulenceCondType2D)(TM | TCT_Integral_Model_2D);
                            else if ( TurbMod == 2 )
                                TM = (TurbulenceCondType2D)(TM | TCT_Prandtl_Model_2D);
                            else if ( TurbMod == 3 )
                                TM = (TurbulenceCondType2D)(TM | TCT_Spalart_Allmaras_Model_2D);
                            else if ( TurbMod == 4 )
                                TM = (TurbulenceCondType2D)(TM | TCT_k_eps_Model_2D);
                            else if ( TurbMod == 5 )
                                TM = (TurbulenceCondType2D)(TM | TCT_Smagorinsky_Model_2D);

                            if ( pTestFlow )
                                TmpArea->FillArea2D((unsigned int)AreaPoint->GetX(0),
                                                    (unsigned int)AreaPoint->GetY(0),
                                                    CT_NO_COND_2D,pTestFlow,Y,TM,AreaMaterialID);
                            else if ( pTestFlow2D )
                                TmpArea->FillArea2D((unsigned int)AreaPoint->GetX(0),
                                                    (unsigned int)AreaPoint->GetY(0),
                                                    CT_NO_COND_2D,pTestFlow2D,Y,TM,AreaMaterialID);
                        } else {
                            *f_stream << "\n";
                            *f_stream << "Bad Area type index \""<< AreaType <<"\" use in \"Area"<< i+1 <<"\"\n" << flush;
                            isRun=0;
                            f_stream->flush();
                            Abort_OpenHyperFLOW2D();
                        }
                        delete TmpArea;
                        *f_stream << "OK\n";
                        f_stream->flush();
                        TmpArea = NULL;
                        AreaPoint=NULL;
                    }
#ifdef _DEBUG_0
                }__except(Area2D* errAreaPtr) {
                *f_stream << "\n";
                *f_stream << "Error init computation area...\n";
                if ( errAreaPtr->GetAreaState() == AS_ERR_OUT_OF_RANGE ) {
                    *f_stream << "\n";
                    *f_stream << "Init Area point ["<< errAreaPtr->GetStartX()  <<"," << errAreaPtr->GetStartY() <<"] out of range.\n" << flush;
                } else if ( errAreaPtr->GetAreaState() == AS_ERR_INIT_POINT ) {
                    *f_stream << "\n";
                    *f_stream << "Init Area point ["<< errAreaPtr->GetStartX()  <<"," << errAreaPtr->GetStartY() <<"] already in initialized node.\n" << flush;
//-- debug ----
//                     DataSnapshot(ErrFileName);
//-- debug ----
                } else {
                    *f_stream << "\n";
                    *f_stream << "Unknown Init Area error.\n" << flush;
                }
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
                J = NULL;
                f_stream->flush();
                Abort_OpenHyperFLOW2D();
            }__except( ComputationalMatrix2D*  m) {
                if ( m->GetMatrixState()==MXS_ERR_OUT_OF_INDEX ) {
                    *f_stream << "\n";
                    *f_stream << "Error init computation area..." << flush;
                    *f_stream << "\n";
                    *f_stream << "Matrix out of index error\n" << flush;
                    f_stream->flush();
                    isRun=0;
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
                    J=NULL;
                    Abort_OpenHyperFLOW2D();
                }
            }__except(void*) {
                *f_stream << "\n";
                *f_stream << "Error init computation area..." << flush;
                *f_stream << "Unknown error\n" << flush;
                f_stream->flush();
                isRun=0;
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
                J=NULL;
                Abort_OpenHyperFLOW2D();
            }
            __end_except;
#endif  // _DEBUG_0
            }

            if ( !PreloadFlag )
#ifdef _DEBUG_0
              ___try {
#endif  // _DEBUG_0

                    *f_stream << "\nFirst initialization computation area...";
                    for (int i=0;i<(int)MaxX;i++ )
                        for (int j=0;j<(int)MaxY;j++ ) {

                                i_err = i;
                                j_err = j;

                                J->GetValue(i,j).idXl  = 1;
                                J->GetValue(i,j).idXr  = 1;
                                J->GetValue(i,j).idYu  = 1;
                                J->GetValue(i,j).idYd  = 1;
                                J->GetValue(i,j).l_min = min(dx*MaxX,dy*MaxY);

                                for (int k=0;k < FlowNode2D<FP,NUM_COMPONENTS>::NumEq;k++) {
                                   J->GetValue(i,j).beta[k]  = beta0;
                                }

                                if ( j==0 || J->GetValue(i,j-1).isCond2D(CT_SOLID_2D) ) {           // is down node present ? (0 or 1)
                                    J->GetValue(i,j).idYd=0;/*J->GetValue(i,j).NGY =2;*/
                                }
                                if ( j==(int)MaxY-1 || J->GetValue(i,j+1).isCond2D(CT_SOLID_2D) ) { // is up node present ? (0 or 1)
                                    J->GetValue(i,j).idYu=0;/*J->GetValue(i,j).NGY =2;*/
                                }
                                if ( i==0 || J->GetValue(i-1,j).isCond2D(CT_SOLID_2D) ) {           // is left node present ? (0 or 1)
                                    J->GetValue(i,j).idXl=0;/*J->GetValue(i,j).NGX =2;*/
                                }
                                if ( i==(int)MaxX-1 ||J->GetValue(i+1,j).isCond2D(CT_SOLID_2D) ) {  // is right node present ? (0 or 1)
                                    J->GetValue(i,j).idXr=0;/*J->GetValue(i,j).NGX =2;*/
                                }

                                if (J->GetValue(i,j).isCond2D(CT_WALL_NO_SLIP_2D) || J->GetValue(i,j).isCond2D(CT_WALL_LAW_2D)) {
                                    J->GetValue(i,j).NGX   = J->GetValue(i,j).idXl - J->GetValue(i,j).idXr + J->GetValue(i,j).idXl*J->GetValue(i,j).idXr;
                                    J->GetValue(i,j).NGY   = J->GetValue(i,j).idYd - J->GetValue(i,j).idYu + J->GetValue(i,j).idYd*J->GetValue(i,j).idYu;
                                }

                                if (!isIgnoreUnsetNodes && !J->GetValue(i,j).isCond2D(CT_NODE_IS_SET_2D) ) {
                                    *f_stream << "\n";
                                    *f_stream << "Node ("<< i_err << "," << j_err <<") has not CT_NODE_IS_SET flag." << flush;
                                    *f_stream << "\n";
                                    *f_stream << "Possible some \"Area\" objects not defined.\n" << flush;
                                    f_stream->flush();
//-- debug ----
                                    DataSnapshot(ErrFileName);
//-- debug ----
                                    Abort_OpenHyperFLOW2D();
                                }

                                if (J->GetValue(i,j).isCond2D(CT_SOLID_2D) ) {
                                    J->GetValue(i,j).Tg = Ts0;
                                }   else {
                                    J->GetValue(i,j).FillNode2D(0,1);
                                }

                                if(J->GetValue(i,j).p == 0.) {
                                   J->GetValue(i,j).Tg = Ts0;
                                }
                            }
                 *f_stream << "OK" << endl;
#ifdef _DEBUG_0
                }__except( ComputationalMatrix2D*  m) {
                *f_stream << "\n";
                *f_stream << "Error set computation area state ("<< i_err << "," << j_err <<")." << "\n" << flush;
                *f_stream << "\n";
                *f_stream << "UMatrix2D< FlowNode2D<FP,NUM_COMPONENTS> >" << "\n" << flush;
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
                J=NULL;
                Abort_OpenHyperFLOW2D();
#ifdef _DEBUG_0
            }__except(Area2D*) {
#endif  // _DEBUG_0
                *f_stream << "\n";
                *f_stream << "Error set computation area state ("<< i_err << "," << j_err <<")." << "\n" << flush;
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
                J=NULL;
                Abort_OpenHyperFLOW2D();
            }
            __end_except;
#endif  // _DEBUG_0

            if(GlobalTime > 0.)
               J->GetValue(0,0).time = GlobalTime;
            else
               GlobalTime=J->GetValue(0,0).time;

            *f_stream << "\nInitial dt=" << dt << "sec." << endl;
//------> place here <------*
           if(ProblemType == SM_NS) { 
               SetWallNodes(f_stream, J);
           }
            
            GlobalSubDomain = ScanArea(f_stream,J, isVerboseOutput);
/* Load additional sources */
             isGasSource  = Data->GetIntVal((char*)"NumSrc");
             if ( Data->GetDataError()==-1 ) Abort_OpenHyperFLOW2D();

             if ( isGasSource ) {
                  SrcList = new SourceList2D(J,Data);
                  SrcList->SetSources2D();
             }

             if(!PreloadFlag) {
                *f_stream << "Seek nodes with non-reflected BC..." << endl;
                *f_stream << "Found " << SetNonReflectedBC(J,nrbc_beta0,f_stream) << " nodes with non-reflected BC...OK" << endl;
             }

             if(ProblemType == SM_NS) {
                 if(!PreloadFlag) {
                    *f_stream << "Set initial boundary layer...";
                    SetInitBoundaryLayer(J,delta_bl);
                    *f_stream << "OK" << endl; 
                   }
             }

#ifdef _DEBUG_0
        }__except(SysException e) {
            {
                const char ConstErrorMessage[]=" handled by <LibExcept> module in <InitSolver>.\n";
                if ( e == SIGSEGV )
                    *f_stream << "\nSIGSEGV" << ConstErrorMessage << flush;
                else if ( e == SIGBUS )
                    *f_stream << "\nSIGBUS" << ConstErrorMessage << flush;
                else if ( e == SIGFPE )
                    *f_stream << "\nSIGFPE" << ConstErrorMessage << flush;
                else if ( e == SIGINT )
                    *f_stream << "\nSIGINT" << ConstErrorMessage << flush;
                else if ( e == SIGABRT )
                    *f_stream << "\nSIGABRT" << ConstErrorMessage << flush;
                else if ( e == SIGIO )
                    *f_stream << "\nSIGIO" << ConstErrorMessage << flush;
                else if ( e == SIGTRAP )
                    *f_stream << "\nSIGTRAP" << ConstErrorMessage << flush;
                else if ( e == SIGKILL )
                    *f_stream << "\nSIGKILL" << ConstErrorMessage << flush;
                f_stream->flush();
            }
           Abort_OpenHyperFLOW2D();
        } __end_except;
//-- debug ----
//           DataSnapshot(ErrFileName);
//-- debug ----

#endif  // _DEBUG_0
        return(NULL);
};
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

inline int CalcChemicalReactions(FlowNode2D<FP,NUM_COMPONENTS>* CalcNode,
                                 ChemicalReactionsModel cr_model, void* CRM_data) {
    ChemicalReactionsModelData2D* model_data = (ChemicalReactionsModelData2D*)CRM_data;
    FP   Y0,Yfu,Yox,Ycp,Yair;

    Yfu  = CalcNode->S[i2d_Yfu]/CalcNode->S[0]; // Fuel
    Yox  = CalcNode->S[i2d_Yox]/CalcNode->S[0]; // OX
    Ycp  = CalcNode->S[i2d_Ycp]/CalcNode->S[0]; // cp
    Yair = 1. - (Yfu+Yox+Ycp);                  // air

    if(cr_model==CRM_ZELDOVICH) {
//--- chemical reactions (Zeldovich model) -------------------------------------------------->
      if ( !CalcNode->isCond2D(CT_Y_CONST_2D) ) {
          Y0   = 1./(Yfu+Yox+Ycp+Yair);
          Yfu  = Yfu*Y0;
          Yox  = Yox*Y0;
          Ycp  = Ycp*Y0;

          if ( CalcNode->Tg > CalcNode->Tf ) {
              if ( Yox > Yfu*model_data->K0 ) { // Yo2 > Yfuel
                   Yox = Yox - Yfu*model_data->K0;
                   Yfu = 0.;
                   Ycp = 1.-Yox-Yair;
              } else {                          // Yo2 < Yfuel
                   Yfu = Yfu - Yox/model_data->K0;
                   Yox = 0.;
                   Ycp = 1. -Yfu-Yair;
              }
           }
        }
//--- chemical reactions (Zeldovich model) -------------------------------------------------->
    }

    CalcNode->R   = model_data->R_Fuel*Yfu+
                            model_data->R_OX*Yox+
                            model_data->R_cp*Ycp+
                            model_data->R_air*Yair;
    CalcNode->CP  = model_data->Cp_Fuel->GetVal(CalcNode->Tg)*Yfu+
                            model_data->Cp_OX->GetVal(CalcNode->Tg)*Yox+
                            model_data->Cp_cp->GetVal(CalcNode->Tg)*Ycp+
                            model_data->Cp_air->GetVal(CalcNode->Tg)*Yair;
    if( ProblemType == SM_NS ) {
        CalcNode->lam = model_data->lam_Fuel->GetVal(CalcNode->Tg)*Yfu+
                                model_data->lam_OX->GetVal(CalcNode->Tg)*Yox+
                                model_data->lam_cp->GetVal(CalcNode->Tg)*Ycp+
                                model_data->lam_air->GetVal(CalcNode->Tg)*Yair;
        CalcNode->mu  = model_data->mu_Fuel->GetVal(CalcNode->Tg)*Yfu+
                                model_data->mu_OX->GetVal(CalcNode->Tg)*Yox+
                                model_data->mu_cp->GetVal(CalcNode->Tg)*Ycp+
                                model_data->mu_air->GetVal(CalcNode->Tg)*Yair;
    }

    if ( Yair<1.e-5 ) {
         Yair =0.;
      }
    if ( Ycp<1.e-8 ) {
         Ycp =0.;
      }
    if ( Yox<1.e-8 ) {
         Yox =0.;
      }
    if ( Yfu<1.e-8 ) {
         Yfu =0.;
      }

     Y0   = 1./(Yfu+Yox+Ycp+Yair);
     Yfu  = Yfu*Y0;
     Yox  = Yox*Y0;
     Ycp  = Ycp*Y0;
     Yair = Yair*Y0;


    CalcNode->Y[0] = Yfu;
    CalcNode->Y[1] = Yox;
    CalcNode->Y[2] = Ycp;
    CalcNode->Y[3] = Yair;

    if ( !CalcNode->isCond2D(CT_Y_CONST_2D) ) {
          CalcNode->S[i2d_Yfu] = fabs(Yfu*CalcNode->S[0]);
          CalcNode->S[i2d_Yox] = fabs(Yox*CalcNode->S[0]);
          CalcNode->S[i2d_Ycp] = fabs(Ycp*CalcNode->S[0]);
     }
 return 1;
}


void SetMinDistanceToWall2D(ComputationalMatrix2D* pJ2D,
                            UArray< XY<int> >* WallNodes2D, 
                            FP x0 ) {

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
                             
                             FP _l_min = sqrt((x-wx)*(x-wx) + (y-wy)*(y-wy));

                             pJ2D->GetValue(i,j).l_min = min(pJ2D->GetValue(i,j).l_min,_l_min);

                             if (pJ2D->GetValue(i,j).l_min == _l_min) {
                                 pJ2D->GetValue(i,j).i_wall = iw;
                                 pJ2D->GetValue(i,j).j_wall = jw;
                             }
                             pJ2D->GetValue(i,j).l_min =  max(min_l_min,pJ2D->GetValue(i,j).l_min);
                        }
             }
         }
      }
   }
}


int SetNonReflectedBC(UMatrix2D< FlowNode2D<FP,NUM_COMPONENTS> >* OutputMatrix2D,
                      FP beta_nrbc,
                      ofstream* f_stream) {
    int nr_nodes = 0;
// Clean walls nodes
#ifdef _OPEN_MP
#pragma omp parallel for \
        reduction(+:nr_nodes)
#endif //_OPEN_MP

    for ( int ii=0;ii<OutputMatrix2D->GetX();ii++ ) {
        for ( int jj=0;jj<OutputMatrix2D->GetY();jj++ ) {
                if ( OutputMatrix2D->GetValue(ii,jj).isCond2D(NT_FARFIELD_2D) ) {
                    nr_nodes++;
                    if ( ii > 0 &&
                         OutputMatrix2D->GetValue(ii-1,jj).isCond2D(CT_NODE_IS_SET_2D) &&
                         !OutputMatrix2D->GetValue(ii-1,jj).isCond2D(CT_WALL_NO_SLIP_2D) &&
                         !OutputMatrix2D->GetValue(ii-1,jj).isCond2D(CT_SOLID_2D) &&
                         !OutputMatrix2D->GetValue(ii-1,jj).isCond2D(NT_FC_2D) ) {

                        OutputMatrix2D->GetValue(ii-1,jj).SetCond2D(CT_NONREFLECTED_2D);
                        nr_nodes++;
                    }
                    if ( ii < OutputMatrix2D->GetX() - 1 &&
                         OutputMatrix2D->GetValue(ii+1,jj).isCond2D(CT_NODE_IS_SET_2D) &&
                         !OutputMatrix2D->GetValue(ii+1,jj).isCond2D(CT_WALL_NO_SLIP_2D) &&
                         !OutputMatrix2D->GetValue(ii+1,jj).isCond2D(CT_SOLID_2D) &&
                         !OutputMatrix2D->GetValue(ii+1,jj).isCond2D(NT_FC_2D) ) {

                        OutputMatrix2D->GetValue(ii+1,jj).SetCond2D(CT_NONREFLECTED_2D);
                        nr_nodes++;
                    }
                    if ( jj > 0 &&
                         OutputMatrix2D->GetValue(ii,jj-1).isCond2D(CT_NODE_IS_SET_2D) &&
                         !OutputMatrix2D->GetValue(ii,jj-1).isCond2D(CT_WALL_NO_SLIP_2D) &&
                         !OutputMatrix2D->GetValue(ii,jj-1).isCond2D(CT_SOLID_2D) &&
                         !OutputMatrix2D->GetValue(ii,jj-1).isCond2D(NT_FC_2D) ) {

                        OutputMatrix2D->GetValue(ii,jj-1).SetCond2D(CT_NONREFLECTED_2D);
                        nr_nodes++;
                    }
                    if ( jj < OutputMatrix2D->GetY() - 1 &&
                         OutputMatrix2D->GetValue(ii,jj+1).isCond2D(CT_NODE_IS_SET_2D) &&
                         !OutputMatrix2D->GetValue(ii,jj+1).isCond2D(CT_WALL_NO_SLIP_2D) &&
                         !OutputMatrix2D->GetValue(ii,jj+1).isCond2D(CT_SOLID_2D) &&
                         !OutputMatrix2D->GetValue(ii,jj+1).isCond2D(NT_FC_2D) ) {

                        OutputMatrix2D->GetValue(ii,jj+1).SetCond2D(CT_NONREFLECTED_2D);
                        nr_nodes++;
                    }
                }
            }
        }
    return nr_nodes;
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



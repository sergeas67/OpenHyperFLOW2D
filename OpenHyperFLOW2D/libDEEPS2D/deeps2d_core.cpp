/*******************************************************************************
*   OpenHyperFLOW2D                                                            *
*                                                                              *
*   Transient, Density based Effective Explicit Parallel Solver (TDEEPS2D)     *
*                                                                              *
*   Version  1.0.0                                                             *
*   Copyright (C)  1995-2013 by Serge A. Suchkov                               *
*   Copyright policy: LGPL V3                                                  *
*                                                                              *
*  deeps2d_core.cpp: OpenHyperFLOW2D solver core code....                      *
*                                                                              *
*  last update: 15/12/2013                                                     *
********************************************************************************/
#include "deeps2d_core.hpp"

#include <sys/time.h>
#include <sys/timeb.h>
#include <sys/file.h>

int NumContour;
int start_iter = 5;
int n_s;

SourceList2D*  SrcList = NULL;
int            isGasSource=0;
int            TurbStartIter;
int            TurbExtModel;
int            err_i, err_j;
int            turb_mod_name_index = 0;
double         Ts0,A,W,Mach;

UArray< XY<int> >* GlobalSubmatrix;
UArray< XY<int> >* WallNodes;

UArray<UMatrix2D< FlowNode2D<double,NUM_COMPONENTS> >*>*     SubmatrixArray;
UArray<UMatrix2D< FlowNodeCore2D<double,NUM_COMPONENTS> >*>* CoreSubmatrixArray;
UArray<XCut>*                                                XCutArray;

//double  makeZero;

int     Nstep;
double  ExitMonitorValue;
int     MonitorNumber;
int     MonitorCondition; // 0 - equal
                          // 1 - less than
                          // 2 - great than

unsigned int     AddSrcStartIter = 0;
double  beta;
double  beta0;
double  CFL;
Table*  CFL_Scenario;
Table*  beta_Scenario;
#ifdef _OPEN_MP
double  DD_max_var;
#endif // OPEN_MP

//------------------------------------------
// Cx,Cy
//------------------------------------------
int     is_Cx_calc;
double  x0_body,y0_body,dx_body,dy_body;
int     Cx_Flow_index;
int     SigmaFi_Flow_index;
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

double* Y=NULL;
double  Cp=0.;
double  lam=0.;
double  mu=0.;
double  Tg=0.;
double  Rg=0.;
double  Pg=0.;
double  Wg=0.;
double  Ug,Vg;
int     CompIndex;

double Y_fuel[4]={1.,0.,0.,0.};  /* fuel */
double Y_ox[4]  ={0.,1.,0.,0.};  /* OX */
double Y_cp[4]  ={0.,0.,1.,0.};  /* cp */
double Y_air[4] ={0.,0.,0.,1.};  /* air */
double Y_mix[4] ={0.,0.,0.,0.};  /* mixture */
const  char*  RMS_Name[11] = {"Rho","Rho*U","Rho*V","Rho*E","Rho*Yfu","Rho*Yox","Rho*Ycp","Rho*k","Rho*eps","Rho*omega","nu_t"};
int    useSwapFile=0;
char   RMSFileName[255];
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
ofstream*                                    pInputData;      // Output data stream
ofstream*                                    pRMS_OutFile;    // Output RMS stream
ofstream*                                    pHeatFlux_OutFile; // Output HeatFlux stream
unsigned int                                 iter = 0;        // iteration number
unsigned int                                 last_iter=0;     // Global iteration number
int                                          isStop=0;        // Stop flag
InputData*                                   Data=NULL;       // Object data loader
UMatrix2D< FlowNode2D<double,NUM_COMPONENTS> >*  J=NULL;      // Main computation area
UArray<Flow*>*                               FlowList=NULL;   // List of 'Flow' objects
UArray<Flow2D*>*                             Flow2DList=NULL; // List of 'Flow2D' objects
UArray<Bound2D*>                             SingleBoundsList;// Single Bounds List;

double                                       dt;              // time step
double*                                      RoUx=NULL;
double*                                      RoVy=NULL;
double                                       GlobalTime=0.;
double                                       CurrentTimePart=0;

ofstream*                                    pOutputData;     // output data stream (file)

int                                          I,NSaveStep;
unsigned int                                 MaxX=0;          // X dimension of computation area
unsigned int                                 MaxY=0;          // Y dimension of computation area

double                                       dxdy,dx2,dy2;
double                                       Gig,SigW,SigF,delta_bl;

unsigned long FileSizeGas      = 0;
int      isVerboseOutput       = 0;
int      isTurbulenceReset     = 0;

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

            FlowNode2D<double,NUM_COMPONENTS>::FT = (FlowType)(_data->GetIntVal((char*)"FlowType"));
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

            MonitorNumber = _data->GetIntVal((char*)"MonitorNumber");
            if ( _data->GetDataError()==-1 ) {
                Abort_OpenHyperFLOW2D();
            }

            if(MonitorNumber > 5 ||
               MonitorNumber < 0) {
               MonitorNumber = 0;
            }

            ExitMonitorValue  = _data->GetFloatVal((char*)"ExitMonitorValue");   // Monitor value for exit
            if ( _data->GetDataError()==-1 ) {
                Abort_OpenHyperFLOW2D();
            }

            beta0 = beta  = _data->GetFloatVal((char*)"beta");     // Base blending factor.
            if ( _data->GetDataError()==-1 ) {
                Abort_OpenHyperFLOW2D();
            }

            beta_Scenario  = _data->GetTable((char*)"beta_Scenario"); //  Base blending factor scenario
            if ( _data->GetDataError()==-1 ) {
                Abort_OpenHyperFLOW2D();
            }

            model_data->K0     = _data->GetFloatVal((char*)"K0"); // Stoichiometric ratio (OX/Fuel)
            if ( _data->GetDataError()==-1 ) {
                Abort_OpenHyperFLOW2D();
            }

            model_data->gamma  = _data->GetFloatVal((char*)"gamma"); // Factor completion of a chemical reaction
            if ( _data->GetDataError()==-1 ) {
                Abort_OpenHyperFLOW2D();
            }

            model_data->Tf     = _data->GetFloatVal((char*)"Tf");   // Ignition temperature
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
            FlowNode2D<double,NUM_COMPONENTS>::dx    = dx;
            FlowNode2D<double,NUM_COMPONENTS>::dy    = dy;
#endif //_UNIFORM_MESH_
            FlowNode2D<double,NUM_COMPONENTS>::Hu[h_fu] = model_data->H_Fuel;
            FlowNode2D<double,NUM_COMPONENTS>::Hu[h_ox] = model_data->H_OX;
            FlowNode2D<double,NUM_COMPONENTS>::Hu[h_cp] = model_data->H_cp;
            FlowNode2D<double,NUM_COMPONENTS>::Hu[h_air] = model_data->H_air;
};

#ifdef _MPI
struct Var_pack {
       double   dt_min;
#ifdef __INTEL_COMPILER
       DD_pack  DD[4+NUM_COMPONENTS+2];
#else 
       DD_pack  DD[FlowNode2D<double,NUM_COMPONENTS>::NumEq];
#endif //__INTEL_COMPILER
};
#endif // _MPI

void DEEPS2D_Run(ofstream* f_stream
#ifdef _MPI
                ,UMatrix2D< FlowNode2D<double,NUM_COMPONENTS> >*     pJ,
                 UMatrix2D< FlowNodeCore2D<double,NUM_COMPONENTS> >* pC,
                 int rank, int last_rank
#endif //_MPI
                 ) {

   // local variables
    double   n_n,m_m,n_n_1,m_m_1;
    double   dXX,dYY,AAA;
    int      n1,n2,n3,n4;
    unsigned int j,k;
    unsigned int StartXLocal,MaxXLocal;
    double   dtdx;
    double   dtdy;
    double   dyy;
    double   dxx;
    double   dx_1,dy_1; // 1/dx, 1/dy
    double   d_time;
    double   t,VCOMP;
    timeval  start, stop, mark1, mark2;
    int      N1,N2,N3,N4;
#ifdef _OPEN_MP
    unsigned int i_max,j_max,k_max;
#endif //_OPEN_MP
#ifdef _MPI
    unsigned int i_max,j_max,k_max;
#endif //_MPI
    double   dt_min_local;
    double   _beta;
    double   DD_local[FlowNode2D<double,NUM_COMPONENTS>::NumEq];
    double   max_RMS;
    int      k_max_RMS;

#ifndef _MPI
    UMatrix2D< FlowNode2D<double,NUM_COMPONENTS> >*     pJ=NULL;
    UMatrix2D< FlowNodeCore2D<double,NUM_COMPONENTS> >* pC=NULL;
    double*   dt_min;
    int*      i_c;
    int*      j_c;
    double    dtmin=1.0;
    double    DD[FlowNode2D<double,NUM_COMPONENTS>::NumEq];
    double    sum_RMS[FlowNode2D<double,NUM_COMPONENTS>::NumEq];
    unsigned long sum_iRMS[FlowNode2D<double,NUM_COMPONENTS>::NumEq];
    UMatrix2D<double> RMS(FlowNode2D<double,NUM_COMPONENTS>::NumEq,SubmatrixArray->GetNumElements());
    UMatrix2D<int>    iRMS(FlowNode2D<double,NUM_COMPONENTS>::NumEq,SubmatrixArray->GetNumElements());
    UMatrix2D<double> DD_max(FlowNode2D<double,NUM_COMPONENTS>::NumEq,SubmatrixArray->GetNumElements());
#else
    double   RMS[FlowNode2D<double,NUM_COMPONENTS>::NumEq];
    unsigned long iRMS[FlowNode2D<double,NUM_COMPONENTS>::NumEq];
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
    dt_min = new double[SubmatrixArray->GetNumElements()];
    i_c    = new int[SubmatrixArray->GetNumElements()];
    j_c    = new int[SubmatrixArray->GetNumElements()];
#endif // _MPI

#ifndef _MPI
    for(int ii=0;ii<(int)SubmatrixArray->GetNumElements();ii++) {
        dt_min[ii] = dtmin = dt;
        i_c[ii] = j_c[ii] = 0;
     }
    snprintf(RMSFileName,255,"RMS-%s",OutFileName);
    CutFile(RMSFileName);
    pRMS_OutFile = OpenData(RMSFileName);
    SaveRMSHeader(pRMS_OutFile);
#else
    MPI::COMM_WORLD.Barrier();
    MPI::COMM_WORLD.Bcast(&dt,1,MPI::DOUBLE,0);

    if( rank == 0 ) {
        char    RMSFileName[255];
        snprintf(RMSFileName,255,"RMS-%s",OutFileName);
        CutFile(RMSFileName);
        pRMS_OutFile = OpenData(RMSFileName);
        SaveRMSHeader(pRMS_OutFile);
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
                    
                 
             do {

                  gettimeofday(&mark2,NULL);
                  gettimeofday(&start,NULL);
                  if( AddSrcStartIter < iter + last_iter){
                   FlowNode2D<double,NUM_COMPONENTS>::isSrcAdd = 1;
                  } else {
                   FlowNode2D<double,NUM_COMPONENTS>::isSrcAdd = 0;
                  }
#ifdef _OPENMP
                  n_s = (int)SubmatrixArray->GetNumElements();
#pragma omp parallel shared(f_stream,CoreSubmatrixArray, SubmatrixArray, chemical_reactions,Y_mix,sum_RMS, sum_iRMS, \
                            Cp,i_max,j_max,k_max,Tg,beta0,CurrentTimePart,DD,dx,dy,MaxX,MaxY,dt_min, dtmin,RMS,iRMS,DD_max,i_c,j_c,n_s) \
                     private(iter,j,k,n1,n2,n3,n4,N1,N2,N3,N4,n_n,m_m,pC,pJ,err_i,err_j,\
                             beta,_beta,AddEq,dXX,dYY,DD_local,AAA,StartXLocal,MaxXLocal,\
                             dtdx,dtdy,dt)
#endif //_OPENMP
                  iter = 0;
                  do {
                      FlowNode2D< double,NUM_COMPONENTS >* CurrentNode=NULL;
                      FlowNodeCore2D< double,NUM_COMPONENTS >* NextNode=NULL;
                      FlowNode2D< double,NUM_COMPONENTS >* UpNode=NULL;          // neast
                      FlowNode2D< double,NUM_COMPONENTS >* DownNode=NULL;        // nodes
                      FlowNode2D< double,NUM_COMPONENTS >* LeftNode=NULL;        // references
                      FlowNode2D< double,NUM_COMPONENTS >* RightNode=NULL;
#ifdef _MPI
// MPI version
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

                  for (int kk=0;kk<FlowNode2D<double,NUM_COMPONENTS>::NumEq;kk++ ) {
                      RMS[kk] = 0.;                                           // Clean sum residual
                      iRMS[kk]= 0;                                            // num involved nodes
                      DD_max[rank].DD[kk].RMS  = 0.;                          // sum residual per rank
                      DD_max[rank].DD[kk].iRMS = 0;                           // num involved nodes per rank
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
                       max_RMS = 0.5 * ExitMonitorValue;
                       k_max_RMS = -1;
                       for ( k=0;k<(int)FlowNode2D<double,NUM_COMPONENTS>::NumEq;k++ ) {
                           for(int ii=0;ii<(int)SubmatrixArray->GetNumElements();ii++) {
                               sum_RMS[k]  = RMS(k,ii)  = 0.;    // Clean sum residual
                               sum_iRMS[k] = iRMS(k,ii) = 0;     // num involved nodes
                           }
                           DD[k]      = 0.;
                       }
                   }
#endif //_MPI

#ifdef _OPENMP
#pragma omp barrier
#pragma omp for private(CurrentNode,NextNode,UpNode,DownNode,LeftNode,RightNode)
// schedule(dynamic) ordered nowait
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
#pragma omp critical
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
                    for ( k=0;k<(int)FlowNode2D<double,NUM_COMPONENTS>::NumEq;k++ ) {
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
                    for (int i = StartXLocal;i<(int)MaxXLocal;i++ ) {
                       for (int j=0;j<(int)MaxY;j++ ) {

                          err_i = i;
                          err_j = j;

                          CurrentNode = &(pJ->GetValue(i,j)); 

                          if ( CurrentNode->isCond2D(CT_NODE_IS_SET_2D) &&
                                 !CurrentNode->isCond2D(CT_SOLID_2D)
                              && !CurrentNode->isCond2D(NT_FC_2D)
                              ) {

                              NextNode    = &(pC->GetValue(i,j)); 

                              CurrentNode->time=GlobalTime;

                                if(CurrentNode->isTurbulenceCond2D(TCT_Spalart_Allmaras_Model_2D))
                                   turb_mod_name_index = 3;

#ifdef _MPI
                                CurrentNode->NodeID = rank;
#else
                                CurrentNode->NodeID = ii;
#endif // _MPI

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


                                AddEq = SetTurbulenceModel(pJ,i,j);

                                beta  = CurrentNode->beta;
                                _beta = 1. - beta;

                                // Scan equation system ... k - number of equation
                                for ( k=0;k<(FlowNode2D<double,NUM_COMPONENTS>::NumEq-AddEq);k++ ) {
                                        int      c_flag = 0;
                                        int      dx_flag, dx2_flag;
                                        int      dy_flag, dy2_flag;

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
                                        
                                        /*
                                        if(((iter + last_iter) < (int)start_iter)
                                            && (k == 1 || k == 2 ))
                                           makeZero = 0.;
                                        else
                                           makeZero = 1.;
                                        */
                                        if ( CurrentNode->FT ) {
                                            NextNode->S[k] = CurrentNode->S[k]*beta+_beta*(dxx*(LeftNode->S[k]+RightNode->S[k])+dyy*(UpNode->S[k]+DownNode->S[k]))*0.5
                                                           - /*makeZero*/(dtdx*dXX+dtdy*(dYY+CurrentNode->F[k]/(j+1))) + (CurrentNode->Src[k])*dt+CurrentNode->SrcAdd[k];
                                        } else {
                                            NextNode->S[k] = CurrentNode->S[k]*beta+_beta*(dxx*(LeftNode->S[k]+RightNode->S[k])+dyy*(UpNode->S[k]+DownNode->S[k]))*0.5
                                                           - /*makeZero*/(dtdx*dXX+dtdy*dYY) + (CurrentNode->Src[k])*dt+CurrentNode->SrcAdd[k];
                                        }
                                
                               }
                         }
                     }
                  }
               }
#ifdef _OPENMP
           }
#endif // _OPENMP

#ifdef _OPENMP
#pragma omp barrier
#pragma omp for private(CurrentNode,NextNode,UpNode,DownNode,LeftNode,RightNode) 
//schedule(dynamic) ordered nowait
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

                              n_n_1 = 1./n_n;
                              m_m_1 = 1./m_m;

                              AddEq = SetTurbulenceModel(pJ,i,j);

                              for ( k=0;k<(FlowNode2D<double,NUM_COMPONENTS>::NumEq-AddEq);k++ ) {

                                  int c_flag = 0;

                                  if ( k < 4 ) // Make bit flags for future test for current equation // FlowNode2D<double,NUM_COMPONENTS>::NumEq-AddEq-NUM_COMPONENTS ?
                                      c_flag  = CT_Ro_CONST_2D   << k;
                                  else if (k<(4+NUM_COMPONENTS))  // 7 ?
                                      c_flag  = CT_Y_CONST_2D;
                                  else if((CurrentNode->isTurbulenceCond2D(TCT_k_eps_Model_2D) ||
                                           CurrentNode->isTurbulenceCond2D(TCT_Spalart_Allmaras_Model_2D) )) 
                                      c_flag  = TCT_k_CONST_2D << (k-4-NUM_COMPONENTS); 

                                  DD_local[k] = 0;

                                  if ( !CurrentNode->isCond2D((CondType2D)c_flag) && 
                                        CurrentNode->S[k] != 0. ) {

                                        double Tmp = CurrentNode->S[k];
                                        double beta_min;

                                        if(k == i2d_RoU && k == i2d_RoV ) {
                                            Tmp = sqrt(CurrentNode->S[i2d_RoU]*CurrentNode->S[i2d_RoU]+
                                                       CurrentNode->S[i2d_RoV]*CurrentNode->S[i2d_RoV]+1.e-30); // Flux
                                        }

                                        if(fabs(Tmp) > 1.e-15)
                                           DD_local[k] = fabs((NextNode->S[k]-CurrentNode->S[k])/Tmp);
                                        else
                                           DD_local[k] = 0.0;

                                        beta_min = min(beta0,beta_Scenario->GetVal(iter+last_iter));

                                        if( bFF == BFF_L) {
                                        //LINEAR locally adopted blending factor function  (LLABFF)
                                          CurrentNode->beta = min(beta_min,(beta_min*beta_min)/(beta_min+DD_local[k]));
                                        } else if( bFF == BFF_LR) {
                                        //LINEAR locally adopted blending factor function with relaxation (LLABFFR)
                                          CurrentNode->beta = min((beta_min+CurrentNode->beta)*0.5,(beta_min*beta_min)/(beta_min+DD_local[k]));
                                        } else if( bFF == BFF_S) {
                                          //SQUARE locally adopted blending factor function (SLABF)
                                          CurrentNode->beta = min(beta_min,(beta_min*beta_min)/(beta_min+DD_local[k]*DD_local[k]));
                                        } else if (bFF == BFF_SR) {
                                          //SQUARE locally adopted blending factor function with relaxation (SLABFFR)
                                          CurrentNode->beta = min((beta_min+CurrentNode->beta)*0.5,(beta_min*beta_min)/(beta_min+DD_local[k]*DD_local[k]));
                                        } else if( bFF == BFF_SQR) {
                                        //SQRT() locally adopted blending factor function (SQRLABF) + most accurate & stable +
                                          CurrentNode->beta = min(beta_min,(beta_min*beta_min)/(beta_min+sqrt(DD_local[k])));
                                        } else if( bFF == BFF_SQRR) {
                                          CurrentNode->beta = min((beta_min+CurrentNode->beta)*0.5,(beta_min*beta_min)/(beta_min+sqrt(DD_local[k]))); 
                                        }
#ifdef _MPI
                                            DD_max[rank].DD[k].DD   = max(DD_max[rank].DD[k].DD,DD_local[k]);
                                            DD_max[rank].DD[k].RMS += DD_local[k]*DD_local[k];
                                            DD_max[rank].DD[k].iRMS++;


                                              if (DD_max[rank].DD[k].DD==DD_local[k] ) {
                                                  DD_max[rank].DD[k].i = i;
                                                  DD_max[rank].DD[k].j = j;
                                              }

#else
                                        RMS(k,ii) += DD_local[k];
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
                                        } else if ((CurrentNode->isTurbulenceCond2D(TCT_k_eps_Model_2D) ||
                                                    CurrentNode->isTurbulenceCond2D(TCT_Spalart_Allmaras_Model_2D)) ){
                                            if ( !CurrentNode->isTurbulenceCond2D((TurbulenceCondType2D)c_flag) )
                                                  CurrentNode->S[k]   =  NextNode->S[k];
                                        }
                                    }

                                    CurrentNode->droYdx[NUM_COMPONENTS]=CurrentNode->droYdy[NUM_COMPONENTS]=0.;

                                    for ( k=4;k<FlowNode2D<double,NUM_COMPONENTS>::NumEq-2;k++ ) {
                                        if ( !CurrentNode->isCond2D(CT_dYdx_NULL_2D) ) {
                                            CurrentNode->droYdx[k-4]=(RightNode->S[k]-LeftNode->S[k])*dx_1*0.5;
                                            CurrentNode->droYdx[NUM_COMPONENTS]+=(RightNode->S[k]-LeftNode->S[k])*dx_1*0.5;
                                        }
                                        if ( !CurrentNode->isCond2D(CT_dYdy_NULL_2D) ) {
                                              CurrentNode->droYdy[k-4]=(UpNode->S[k]-DownNode->S[k])*dy_1*0.5;
                                              CurrentNode->droYdy[NUM_COMPONENTS]+=(DownNode->S[k]-UpNode->S[k])*dy_1*0.5;
                                        }
                                    }

                                    if (CurrentNode->isCond2D(CT_WALL_NO_SLIP_2D) || CurrentNode->isCond2D(CT_WALL_SLIP_2D) )  {
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
                                    CalcChemicalReactions(CurrentNode,CRM_ZELDOVICH, (void*)(&chemical_reactions));

                                    if((int)(iter+last_iter) < TurbStartIter) {
                                       CurrentNode->FillNode2D(0,1,SigW,SigF,(TurbulenceExtendedModel)TurbExtModel,delta_bl);
                                    } else {
                                       CurrentNode->FillNode2D(1,0,SigW,SigF,(TurbulenceExtendedModel)TurbExtModel,delta_bl);
                                    }

                                    if( CurrentNode->Tg < 0. ) {
                                        *f_stream << "\nTg=" << CurrentNode->Tg << " K. p=" << CurrentNode->p <<" Pa dt=" << dt << " sec.\n" << flush;
#ifdef _MPI
                                        *f_stream << "\nERROR: Computational unstability in UMatrix2D< FlowNode2D<double,NUM_COMPONENTS> >(" << (int)(x0/FlowNode2D<double,NUM_COMPONENTS>::dx) + i <<","<< j << ") on iteration " << iter+last_iter<< "...\n";
#else
  #ifdef _OPENMP
                                        *f_stream << "\nERROR: Computational unstability in local (num_thread="<< omp_get_thread_num() << ") UMatrix2D< FlowNode2D<double,NUM_COMPONENTS> >(" << i <<","<< j <<") on iteration " << iter+last_iter << "...\n";
  #else
                                        *f_stream << "\nERROR: Computational unstability in UMatrix2D< FlowNode2D<double,NUM_COMPONENTS> >(" << i <<","<< j <<") on iteration " << iter+last_iter<< "...\n";
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
                                            double CFL_min      = min(CFL,CFL_Scenario->GetVal(iter+last_iter));
                                            AAA                 = sqrt(CurrentNode->k*CurrentNode->R*CurrentNode->Tg); 
                                            dt_min_local        = CFL_min*
                                                                  min(dx/(AAA+fabs(CurrentNode->U)),dy/(AAA+fabs(CurrentNode->V)));
#ifdef _MPI
                                            DD_max[rank].dt_min = min(DD_max[rank].dt_min, dt_min_local); 
#else
                                            dt_min[ii]          = min(dt_min[ii],dt_min_local);
#endif // _MPI

                                    }
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
                                      ,ii,(int)SubmatrixArray->GetNumElements()-1 
#endif // _MPI
                                   );
#ifdef _OPENMP
 for (DD_max_var = 1.,k=0;k<(int)FlowNode2D<double,NUM_COMPONENTS>::NumEq;k++ ) {
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
               for (k=0;k<(int)(FlowNode2D<double,NUM_COMPONENTS>::NumEq);k++ ) {
                  DD_max_var = max(DD_max[ii].DD[k].DD,DD_max_var);
                  RMS[k] += DD_max[ii].DD[k].RMS;
                  iRMS[k]+= DD_max[ii].DD[k].iRMS;
                  if(DD_max[ii].DD[k].DD == DD_max_var) {
                   i_max = DD_max[ii].DD[k].i;
                   j_max = DD_max[ii].DD[k].j;
                   k_max = k;
                 }
             }
           }
               for (k=0;k<(int)(FlowNode2D<double,NUM_COMPONENTS>::NumEq);k++ ) {
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
#else
#ifdef _OPENMP
 }

#pragma omp single
#endif // _OPENMP
    {
        for(k=0;k<(int)(FlowNode2D<double,NUM_COMPONENTS>::NumEq);k++ ) {
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
#endif // _MPI
        CurrentTimePart += dt;
         if ( isVerboseOutput && iter/NOutStep*NOutStep == iter ) {
             gettimeofday(&mark1,NULL);
             d_time = (double)(mark1.tv_sec-mark2.tv_sec)+(double)(mark1.tv_usec-mark2.tv_usec)*1.e-6; 

             if(d_time > 0.)
                 VCOMP = (double)(NOutStep)/d_time;
             else
                 VCOMP = 0.;
             memcpy(&mark2,&mark1,sizeof(mark1));
             SaveRMS(pRMS_OutFile,last_iter+iter,
#ifdef  _MPI
             RMS);
#else
             sum_RMS);
#endif // _MPI

             if(k_max_RMS == i2d_nu_t)
                k_max_RMS +=turb_mod_name_index;

             if(k_max_RMS != -1 && (MonitorNumber == 0 || MonitorNumber == 5))
             *f_stream << "Step No " << iter+last_iter << " maxRMS["<< RMS_Name[k_max_RMS] << "]="<< (double)(max_RMS*100.) \
                        <<  " % step_time=" << (double)d_time << " sec (" << (double)VCOMP <<" step/sec) dt="<< dt <<"\n" << flush;
             else if(MonitorNumber > 0 &&  MonitorNumber < 5 )
                 *f_stream << "Step No " << iter+last_iter << " maxRMS["<< RMS_Name[MonitorNumber-1] << "]="<< (double)(max_RMS*100.) \
                  <<  " % step_time=" << (double)d_time << " sec (" << (double)VCOMP <<" step/sec) dt="<< dt <<"\n" << flush;
             else
             *f_stream << "Step No " << iter+last_iter << " maxRMS["<< k_max_RMS << "]="<< (double)(max_RMS*100.) \
                        <<  " % step_time=" << (double)d_time << " sec (" << (double)VCOMP <<" step/sec) dt="<< dt <<"\n" << flush;
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
#endif // _PARALLEL_RECALC_Y_PLUS_

#ifdef _PARALLEL_RECALC_Y_PLUS_
    ParallelRecalc_y_plus(pJ,WallNodes,WallNodesUw_2D,x0);
#endif // _PARALLEL_RECALC_Y_PLUS_

    if( rank == 0) {
#ifndef _PARALLEL_RECALC_Y_PLUS_
         Recalc_y_plus(J,WallNodes);
#endif // _PARALLEL_RECALC_Y_PLUS_
#else
#ifdef _PARALLEL_RECALC_Y_PLUS_
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
         d_time = (double)(stop.tv_sec-start.tv_sec)+(double)(stop.tv_usec-start.tv_usec)*1.e-6; 
         *f_stream << "HyperFLOW/DEEPS computation cycle time=" << (double)d_time << " sec ( average  speed " << (double)(Nstep/d_time) <<" step/sec).       \n" << flush;
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
            WallNodesUw_2D = new UArray<double>(NumWallNodes,-1);

        MPI::COMM_WORLD.Recv(WallNodesUw_2D->GetArrayPtr(),
                             NumWallNodes*sizeof(double),
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
                               NumWallNodes*sizeof(double),
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
            }__except( UMatrix2D<double>*  m) {
                *f_stream << "\n";
                *f_stream << "Error in UMatrix2D<double> ("<< err_i << "," << err_j <<")." << "\n" << flush;
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
                *f_stream << "Error in UMatrix2D< FlowNode2D<double,NUM_COMPONENTS> ("<< err_i << "," << err_j <<")->";
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

u_long SetWallNodes(ofstream* f_str, ComputationalMatrix2D* pJ) {
u_long NumWallNodes=0;
    for (int j=0;j<(int)pJ->GetY();j++ ) {
       for (int i=0;i<(int)pJ->GetX();i++ ) {

                FlowNode2D< double,NUM_COMPONENTS >* CurrentNode=NULL;
                FlowNode2D< double,NUM_COMPONENTS >* UpNode=NULL;
                FlowNode2D< double,NUM_COMPONENTS >* DownNode=NULL;
                FlowNode2D< double,NUM_COMPONENTS >* LeftNode=NULL;
                FlowNode2D< double,NUM_COMPONENTS >* RightNode=NULL;

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
            if ( pJ->GetValue(i,j).isCond2D(CT_WALL_SLIP_2D) || 
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
    long num_active_nodes,active_nodes_per_submatrix;
    UArray< XY<int> >* pWallNodes;
    UArray< XY<int> >* Submatrix;
    pWallNodes = new UArray< XY<int> >();
    Submatrix = new UArray< XY<int> >();
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
            if ( pJ->GetValue(i,j).isCond2D(CT_WALL_SLIP_2D) || 
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
active_nodes_per_submatrix = num_active_nodes/num_threads; 
if ( isPrint )
        *f_str << "Found " << pWallNodes->GetNumElements() <<" wall nodes from " << num_active_nodes << " gas filled nodes (" << num_threads <<" threads, "<< active_nodes_per_submatrix <<" active nodes per thread).\n" << flush;

    if ( isPrint )
        *f_str << "Scan computation area for lookup minimal distance from internal nodes to wall.\n" << flush;
#ifdef _DEBUG_0 // 2
gettimeofday(&time2,NULL);
#endif // _DEBUG_0 // 2
num_active_nodes=0;
ijsm.SetX(0);

if(isTurbulenceReset) {
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

                    if(isTurbulenceReset) {
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
                    
                    if(num_active_nodes >= active_nodes_per_submatrix) {
                       ijsm.SetY(i+1);
                       Submatrix->AddElement(&ijsm);
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
*f_str << "SubMatrix decomposition was finished:\n";
for(jj=0;jj<Submatrix->GetNumElements();jj++) {
   *f_str << "SubMatrix[" << jj << "]->["<< Submatrix->GetElementPtr(jj)->GetX() <<","<< Submatrix->GetElementPtr(jj)->GetY() <<"]\n";
}
if(isTurbulenceReset) {
   isTurbulenceReset = 0;
}
f_str->flush();
return Submatrix;
}

void SetInitBoundaryLayer(ComputationalMatrix2D* pJ, double delta) {
    for (int i=0;i<(int)pJ->GetX();i++ ) {
           for (int j=0;j<(int)pJ->GetY();j++ ) {
                  if (pJ->GetValue(i,j).isCond2D(CT_NODE_IS_SET_2D) &&
                      !pJ->GetValue(i,j).isCond2D(CT_SOLID_2D) &&
                       pJ->GetValue(i,j).time == 0. &&
                       delta > 0) {
                       if(pJ->GetValue(i,j).l_min <= delta)
                          pJ->GetValue(i,j).S[i2d_RoU] = pJ->GetValue(i,j).S[i2d_RoU] * pJ->GetValue(i,j).l_min/delta;
                          pJ->GetValue(i,j).S[i2d_RoV] = pJ->GetValue(i,j).S[i2d_RoV] * pJ->GetValue(i,j).l_min/delta;
                          pJ->GetValue(i,j).FillNode2D(0,1);
                   }
           }
    }
}

#ifdef _PARALLEL_RECALC_Y_PLUS_
void RecalcWallFrictionVelocityArray2D(ComputationalMatrix2D* pJ,
                                       UArray<double>* WallFrictionVelocityArray2D,
                                       UArray< XY<int> >* WallNodes2D) { 
    for(int ii=0;ii<(int)WallNodes2D->GetNumElements();ii++) {
        unsigned int iw,jw;
        double tau_w;
        iw = WallNodes2D->GetElementPtr(ii)->GetX();
        jw = WallNodes2D->GetElementPtr(ii)->GetY();
        tau_w = (fabs(pJ->GetValue(iw,jw).dUdy)  +
                 fabs(pJ->GetValue(iw,jw).dVdx)) * pJ->GetValue(iw,jw).mu;
        WallFrictionVelocityArray2D->GetElement(ii) = sqrt(tau_w/pJ->GetValue(iw,jw).S[i2d_Ro]+1e-30);
    }
}

UArray<double>* GetWallFrictionVelocityArray2D(ComputationalMatrix2D* pJ, 
                                               UArray< XY<int> >* WallNodes2D) { 
    UArray<double>* WallFrictionVelocityArray2D;
    WallFrictionVelocityArray2D = new UArray<double>();
    for(int ii=0;ii<(int)WallNodes2D->GetNumElements();ii++) {
        unsigned int iw,jw;
        double tau_w, U_w;
        iw = WallNodes2D->GetElementPtr(ii)->GetX();
        jw = WallNodes2D->GetElementPtr(ii)->GetY();
        tau_w = (fabs(pJ->GetValue(iw,jw).dUdy)  +
                 fabs(pJ->GetValue(iw,jw).dVdx)) * pJ->GetValue(iw,jw).mu;
        U_w   = sqrt(tau_w/pJ->GetValue(iw,jw).S[i2d_Ro]+1e-30);
        WallFrictionVelocityArray2D->AddElement(&U_w);
    }
  return WallFrictionVelocityArray2D;
}

void ParallelRecalc_y_plus(ComputationalMatrix2D* pJ, 
                           UArray< XY<int> >* WallNodes,
                           UArray<double>* WallFrictionVelocity2D,
                           double x0) {
#ifdef _OPEN_MP
#pragma omp parallel for
#endif //_OPEN_MP
        for (int i=0;i<(int)pJ->GetX();i++ ) {
            for (int j=0;j<(int)pJ->GetY();j++ ) {
                 if ( pJ->GetValue(i,j).isCond2D(CT_NODE_IS_SET_2D) &&
                      !pJ->GetValue(i,j).isCond2D(CT_SOLID_2D)) {
                        
                        for (int ii=0;ii<(int)WallNodes->GetNumElements();ii++) {

                             unsigned int iw,jw;
                             double U_w,x,y,wx,wy;

                             iw = WallNodes->GetElementPtr(ii)->GetX();
                             jw = WallNodes->GetElementPtr(ii)->GetY();

                             U_w   =  WallFrictionVelocity2D->GetElement(ii);

                             x = x0 + i * FlowNode2D<double,NUM_COMPONENTS>::dx;
                             y = j * FlowNode2D<double,NUM_COMPONENTS>::dy;

                             wx = iw * FlowNode2D<double,NUM_COMPONENTS>::dx;
                             wy = jw * FlowNode2D<double,NUM_COMPONENTS>::dy;

                             if(x == wx && y == wy ) {
                                 pJ->GetValue(i,j).y_plus = U_w*min((FlowNode2D<double,NUM_COMPONENTS>::dx),
                                                                    (FlowNode2D<double,NUM_COMPONENTS>::dy))/pJ->GetValue(i,j).mu;
                             } else {
                                 if(pJ->GetValue(i,j).l_min == min(pJ->GetValue(i,j).l_min,
                                                                   sqrt((x-wx)*(x-wx) + (y-wy)*(y-wy))))
                                 pJ->GetValue(i,j).y_plus = U_w*pJ->GetValue(i,j).l_min/pJ->GetValue(i,j).mu;
                             }
                        }
                 }
            }
        }
}
#else
void Recalc_y_plus(ComputationalMatrix2D* pJ, UArray< XY<int> >* WallNodes) {
    unsigned int iw,jw;
    double tau_w, U_w;

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
                                      pJ->GetValue(i,j).y_plus = U_w*pJ->GetValue(i,j).l_min/pJ->GetValue(i,j).mu;
                               }
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

void PrintCond(ofstream* OutputData, FlowNode2D<double,NUM_COMPONENTS>* fn) {
    //Const conditions

    if ( fn->isCond2D(CT_Ro_CONST_2D) )
        *OutputData << "-CT_Ro_CONST_2D" << flush;
    if ( fn->isCond2D(CT_U_CONST_2D) )
        *OutputData << "-CT_U_CONST_2D" << flush;
    if ( fn->isCond2D(CT_V_CONST_2D) )
        *OutputData << "-CT_V_CONST_2D" << flush;
    if ( fn->isCond2D(CT_T_CONST_2D) )
        *OutputData << "-CT_T_CONST_2D" << flush;
    if ( fn->isCond2D(CT_Y_CONST_2D) )
        *OutputData << "-CT_Y_CONST_2D" << flush;

    // dF/dx = 0
    if ( fn->isCond2D(CT_dRodx_NULL_2D) )
        *OutputData << "-CT_dRodx_NULL_2D" << flush;
    if ( fn->isCond2D(CT_dUdx_NULL_2D) )
        *OutputData << "-CT_dUdx_NULL_2D" << flush;
    if ( fn->isCond2D(CT_dVdx_NULL_2D) )
        *OutputData << "-CT_dVdx_NULL_2D" << flush;
    if ( fn->isCond2D(CT_dTdx_NULL_2D) )
        *OutputData << "-CT_dTdx_NULL_2D" << flush;
    if ( fn->isCond2D(CT_dYdx_NULL_2D) )
        *OutputData << "-CT_dYdx_NULL_2D" << flush;

    // dF/dy = 0

    if ( fn->isCond2D(CT_dRody_NULL_2D) )
        *OutputData << "-CT_dRody_NULL_2D"<< flush;
    if ( fn->isCond2D(CT_dUdy_NULL_2D) )
        *OutputData << "-CT_dUdy_NULL_2D" << flush;
    if ( fn->isCond2D(CT_dVdy_NULL_2D) )
        *OutputData << "-CT_dVdy_NULL_2D" << flush;
    if ( fn->isCond2D(CT_dTdy_NULL_2D) )
        *OutputData << "-CT_dTdy_NULL_2D" << flush;
    if ( fn->isCond2D(CT_dYdy_NULL_2D) )
        *OutputData << "-CT_dYdy_NULL_2D" << flush;

    // d2F/dx2 =
    if ( fn->isCond2D(CT_d2Rodx2_NULL_2D) )
        *OutputData << "-CT_d2Rodx2_NULL_2D" << flush;
    if ( fn->isCond2D(CT_d2Udx2_NULL_2D) )
        *OutputData << "-CT_d2Udx2_NULL_2D" << flush;
    if ( fn->isCond2D(CT_d2Vdx2_NULL_2D) )
        *OutputData << "-CT_d2Vdx2_NULL_2D" << flush;
    if ( fn->isCond2D(CT_d2Tdx2_NULL_2D) )
        *OutputData << "-CT_d2Tdx2_NULL_2D" << flush;
    if ( fn->isCond2D(CT_d2Ydx2_NULL_2D) )
        *OutputData << "-CT_d2Ydx2_NULL_2D" << flush;

    // d2F/dy2 = 0

    if ( fn->isCond2D(CT_d2Rody2_NULL_2D) )
        *OutputData << "-CT_d2Rody2_NULL_2D"<< flush;
    if ( fn->isCond2D(CT_d2Udy2_NULL_2D) )
        *OutputData << "-CT_d2Udy2_NULL_2D" << flush;
    if ( fn->isCond2D(CT_d2Vdy2_NULL_2D) )
        *OutputData << "-CT_d2Vdy2_NULL_2D" << flush;
    if ( fn->isCond2D(CT_d2Tdy2_NULL_2D) )
        *OutputData << "-CT_d2Tdy2_NULL_2D" << flush;
    if ( fn->isCond2D(CT_d2Ydy2_NULL_2D) )
        *OutputData << "-CT_d2Ydy2_NULL_2D" << flush;

    // Wall conditions:
    // - Slip
    if ( fn->isCond2D(CT_WALL_SLIP_2D) )
        *OutputData << "-CT_WALL_SLIP_2D" << flush;
    // - Glide
    if ( fn->isCond2D(CT_WALL_NO_SLIP_2D) )
        *OutputData << "-CT_WALL_NO_SLIP_2D"  << flush;
    // is solid body ?
    if ( fn->isCond2D(CT_SOLID_2D) )
        *OutputData << "-CT_SOLID_2D" << flush;
    // is node activate ?
    if ( fn->isCond2D(CT_NODE_IS_SET_2D) )
        *OutputData << "-CT_NODE_IS_SET_2D" << flush;
    // boundary layer refinement
    if ( fn->isCond2D(CT_BL_REFINEMENT_2D) )
        *OutputData << "-CT_BL_REFINEMENT_2D" << flush;
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

void SaveRMSHeader(ofstream* OutputData) {
        char  TechPlotTitle[1024];
        char  TmpData[256]={'\0'};

        if(is_Cx_calc)
           snprintf(TmpData,256,", Cx(N), Cy(N), Fx(N), Fy(N)\n");

        snprintf(TechPlotTitle,1024,"#VARIABLES = N, RMS_Ro(N), RMS_RoU(N), RMS_RoV(N), RMS_RoE(N), RMS_RoY_fu(N), RMS_RoY_ox(N), RMS_RoY_cp(N), RMS_k(N), RMS_eps(N)%s",TmpData);

        *OutputData <<  TechPlotTitle << endl;
    }

    void SaveRMS(ofstream* OutputData,unsigned int n, double* outRMS) {
         *OutputData <<  n  << " ";
         for(int i=0;i<FlowNode2D<double,NUM_COMPONENTS>::NumEq;i++) {
             *OutputData <<  outRMS[i]  << " ";
         }

        if(is_Cx_calc) {
          *OutputData << Calc_Cx_2D(J,x0_body,y0_body,dx_body,dy_body,Flow2DList->GetElement(Cx_Flow_index-1)) << " " <<  Calc_Cy_2D(J,x0_body,y0_body,dx_body,dy_body,Flow2DList->GetElement(Cx_Flow_index-1)) << " ";
          *OutputData << CalcXForce2D(J,x0_body,y0_body,dx_body,dy_body) << " " <<  CalcYForce2D(J,x0_body,y0_body,dx_body,dy_body);
        }

        *OutputData << endl;
    }

    void SaveData2D(ofstream* OutputData, int type) { // type = 1 - GNUPLOT
        int    i,j;
        char   TechPlotTitle1[1024]={0};
        char   TechPlotTitle2[256]={0};
        char   YR[2];
        double Mach,A,W,Re,Re_t,dx_out,dy_out;
        char   RT[10];
        if(is_p_asterisk_out)
          snprintf(RT,10,"p*");
        else
          snprintf(RT,10,"mu_t/mu");

        if(FlowNode2D<double,NUM_COMPONENTS>::FT == 1) // FT_FLAT
          snprintf(YR,2,"R");
        else
          snprintf(YR,2,"Y");

        snprintf(TechPlotTitle1,1024,"VARIABLES = X, %s, U, V, T, p, Rho, Y_fuel, Y_ox, Y_cp, Y_i, %s, Mach"
                                     "\n",YR, RT); 
        snprintf(TechPlotTitle2,256,"ZONE T=\"Time: %g sec.\" I= %i J= %i F=POINT\n",GlobalTime, MaxX, MaxY);

        if ( type ) {
            *OutputData <<  TechPlotTitle1;
            *OutputData <<  TechPlotTitle2;
        } else {
            *OutputData <<  TechPlotTitle1;
            *OutputData <<  TechPlotTitle2;
        }
        
        dx_out =(dx*MaxX)/(MaxX-1); 
        dy_out =(dy*MaxY)/(MaxY-1); 
        
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
                      *OutputData << Mach  << "  ";  
                    else
                      *OutputData << "  0  ";
                } else {
                    *OutputData << "  0  ";
                }
                *OutputData <<  "\n" ;
            }
            if ( type )
                *OutputData <<  "\n" ;
        }
    }

    inline double kg(double Cp, double R) {
        return(Cp/(Cp-R));
    }

inline  void CalcHeatOnWallSources(UMatrix2D< FlowNode2D<double,NUM_COMPONENTS> >* F, double dx, double dy, double dt, int rank, int last_rank) {

        unsigned int StartXLocal,MaxXLocal;
        double dx_local, dy_local;

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

                FlowNode2D< double,NUM_COMPONENTS >* CurrentNode=NULL;
                FlowNode2D< double,NUM_COMPONENTS >* UpNode=NULL;
                FlowNode2D< double,NUM_COMPONENTS >* DownNode=NULL;
                FlowNode2D< double,NUM_COMPONENTS >* LeftNode=NULL;
                FlowNode2D< double,NUM_COMPONENTS >* RightNode=NULL;


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

                if ( CurrentNode->isCond2D(CT_WALL_SLIP_2D) || 
                     CurrentNode->isCond2D(CT_WALL_NO_SLIP_2D)) { 
                    double lam_eff = 0;
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
        static int      PreloadFlag = 0,p_g=0;

        SubmatrixArray     = new UArray<UMatrix2D< FlowNode2D<double,NUM_COMPONENTS> >* >();
        CoreSubmatrixArray = new UArray<UMatrix2D< FlowNodeCore2D<double,NUM_COMPONENTS> >* >();
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

            NumContour = Data->GetIntVal((char*)"NumContour"); // Number of "Contour" objects 
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
                double  LamF;
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

                TmpFlow2D = new Flow2D(lam,mu,Cp,Tg,Pg,Rg,Ug,Vg);
                TmpFlow2D->CorrectFlow(Tg,Pg,sqrt(Ug*Ug+Vg*Vg+1.e-30),FV_VELOCITY);

                Flow2DList->AddElement(&TmpFlow2D);
                *f_stream << "Add object \"Flow2D-" << i+1 <<  "\"...OK\n" << flush;
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

            FileSizeGas =  MaxX*MaxY*sizeof(FlowNode2D<double,NUM_COMPONENTS>);
            

            GasSwapData   = LoadSwapFile2D(GasSwapFileName,
                                           (int)MaxX,
                                           (int)MaxY,
                                           sizeof(FlowNode2D<double,NUM_COMPONENTS>),
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
                    J = new UMatrix2D< FlowNode2D<double,NUM_COMPONENTS> >((FlowNode2D<double,NUM_COMPONENTS>*)GasSwapData,MaxX,MaxY);
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
                    J = new UMatrix2D< FlowNode2D<double,NUM_COMPONENTS> >(MaxX,MaxY);
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

            dt=1.;
            for (int i = 0;i<(int)FlowList->GetNumElements();i++ ) {
                 double CFL_min      = min(CFL,CFL_Scenario->GetVal(iter+last_iter));
                 dt = min(dt,CFL_min*dx*dy/(dy*(FlowList->GetElement(i)->Asound()*2.)+
                                            dx*FlowList->GetElement(i)->Asound()));
            }

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
                            if ( strstr(BoundStr,"CT_Ro_CONST_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_Ro_CONST_2D);
                            if ( strstr(BoundStr,"CT_U_CONST_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_U_CONST_2D);
                            if ( strstr(BoundStr,"CT_V_CONST_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_V_CONST_2D);
                            if ( strstr(BoundStr,"CT_T_CONST_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_T_CONST_2D);
                            if ( strstr(BoundStr,"CT_Y_CONST_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_Y_CONST_2D);
                            if ( strstr(BoundStr,"CT_WALL_SLIP_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_WALL_SLIP_2D);
                            if ( strstr(BoundStr,"CT_WALL_NO_SLIP_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_WALL_NO_SLIP_2D);
                            if ( strstr(BoundStr,"CT_dRodx_NULL_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_dRodx_NULL_2D);
                            if ( strstr(BoundStr,"CT_dUdx_NULL_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_dUdx_NULL_2D);
                            if ( strstr(BoundStr,"CT_dVdx_NULL_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_dVdx_NULL_2D);
                            if ( strstr(BoundStr,"CT_dTdx_NULL_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_dTdx_NULL_2D);
                            if ( strstr(BoundStr,"CT_dYdx_NULL_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_dYdx_NULL_2D);
                            if ( strstr(BoundStr,"CT_dRody_NULL_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_dRody_NULL_2D);
                            if ( strstr(BoundStr,"CT_dUdy_NULL_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_dUdy_NULL_2D);
                            if ( strstr(BoundStr,"CT_dVdy_NULL_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_dVdy_NULL_2D);
                            if ( strstr(BoundStr,"CT_dTdy_NULL_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_dTdy_NULL_2D);
                            if ( strstr(BoundStr,"CT_dYdy_NULL_2D_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_dYdy_NULL_2D);
                            if ( strstr(BoundStr,"CT_d2Rodx2_NULL_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_d2Rodx2_NULL_2D);
                            if ( strstr(BoundStr,"CT_d2Udx2_NULL_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_d2Udx2_NULL_2D);
                            if ( strstr(BoundStr,"CT_d2Vdx2_NULL_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_d2Vdx2_NULL_2D);
                            if ( strstr(BoundStr,"CT_d2Tdx2_NULL_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_d2Tdx2_NULL_2D);
                            if ( strstr(BoundStr,"CT_d2Ydx2_NULL_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_d2Ydx2_NULL_2D);
                            if ( strstr(BoundStr,"CT_d2Rody2_NULL_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_d2Rody2_NULL_2D);
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
                            if ( strstr(BoundStr,"NT_WS_2D") )
                                TmpCT = (CondType2D)(TmpCT | NT_WS_2D);
                            else if ( strstr(BoundStr,"NT_WNS_2D") )
                                TmpCT = (CondType2D)(TmpCT | NT_WNS_2D);
                            if ( strstr(BoundStr,"NT_FC_2D") )
                                TmpCT = (CondType2D)(TmpCT | NT_FC_2D);
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
                            } else {
                                pTestFlow2D = Flow2DList->GetElement(FlowIndex-1);
                                sprintf(FlowStr,"Flow2D-%i.CompIndex",FlowIndex);
                                isFlow2D=1;
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
                            if ( strstr(BoundStr,"CT_Ro_CONST_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_Ro_CONST_2D);
                            if ( strstr(BoundStr,"CT_U_CONST_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_U_CONST_2D);
                            if ( strstr(BoundStr,"CT_V_CONST_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_V_CONST_2D);
                            if ( strstr(BoundStr,"CT_T_CONST_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_T_CONST_2D);
                            if ( strstr(BoundStr,"CT_Y_CONST_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_Y_CONST_2D);
                            if ( strstr(BoundStr,"CT_dRodx_NULL_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_dRodx_NULL_2D);
                            if ( strstr(BoundStr,"CT_dUdx_NULL_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_dUdx_NULL_2D);
                            if ( strstr(BoundStr,"CT_dVdx_NULL_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_dVdx_NULL_2D);
                            if ( strstr(BoundStr,"CT_dTdx_NULL_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_dTdx_NULL_2D);
                            if ( strstr(BoundStr,"CT_dYdx_NULL_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_dYdx_NULL_2D);
                            if ( strstr(BoundStr,"CT_dRody_NULL_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_dRody_NULL_2D);
                            if ( strstr(BoundStr,"CT_dUdy_NULL_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_dUdy_NULL_2D);
                            if ( strstr(BoundStr,"CT_dVdy_NULL_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_dVdy_NULL_2D);
                            if ( strstr(BoundStr,"CT_dTdy_NULL_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_dTdy_NULL_2D);
                            if ( strstr(BoundStr,"CT_dYdy_NULL_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_dYdy_NULL_2D);
                            if ( strstr(BoundStr,"CT_d2Rodx2_NULL_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_d2Rodx2_NULL_2D);
                            if ( strstr(BoundStr,"CT_d2Udx2_NULL_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_d2Udx2_NULL_2D);
                            if ( strstr(BoundStr,"CT_d2Vdx2_NULL_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_d2Vdx2_NULL_2D);
                            if ( strstr(BoundStr,"CT_d2Tdx2_NULL_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_d2Tdx2_NULL_2D);
                            if ( strstr(BoundStr,"CT_d2Ydx2_NULL_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_d2Ydx2_NULL_2D);
                            if ( strstr(BoundStr,"CT_d2Rody2_NULL_2D") )
                                TmpCT = (CondType2D)(TmpCT | CT_d2Rody2_NULL_2D);
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
                            if ( strstr(BoundStr,"NT_WS_2D") )
                                TmpCT = (CondType2D)(TmpCT | NT_WS_2D);
                            else if ( strstr(BoundStr,"NT_WNS_2D") )
                                TmpCT = (CondType2D)(TmpCT | NT_WNS_2D);
                            if ( strstr(BoundStr,"NT_FC_2D") )
                                TmpCT = (CondType2D)(TmpCT | NT_FC_2D);
                            if ( strstr(BoundStr,"NT_S_2D") )
                                TmpCT = (CondType2D)(TmpCT | NT_S_2D);

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
                                sprintf(NameBound,"%s.Bound%i.Flow",NameContour,i);
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
                            } else {
                                pTestFlow2D = Flow2DList->GetElement(FlowIndex-1);
                                sprintf(FlowStr,"Flow2D-%i.CompIndex",FlowIndex);
                                isFlow2D=1;
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

            for (int i = 0;i<(int)Flow2DList->GetNumElements();i++ ) {
                double CFL_min  = min(CFL,CFL_Scenario->GetVal(iter+last_iter));
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
                            if ( meshType ) {
                                J->GetValue(i,j).dx  = dx;
                                J->GetValue(i,j).dy  = dy;
                            }
#endif //_UNIFORM_MESH_

                            if ( FlowNode2D<double,NUM_COMPONENTS>::FT == FT_AXISYMMETRIC )
                                 J->GetValue(i,j).r     = (j+1)*dy;
                            J->GetValue(i,j).Tf    = chemical_reactions.Tf;
                            J->GetValue(i,j).BGX   = 1.;
                            J->GetValue(i,j).BGY   = 1.;
                            J->GetValue(i,j).isCleanSources = 1;
                            J->GetValue(i,j).NGX   = 0;
                            J->GetValue(i,j).NGY   = 0;

                            for ( k=0;k<(int)FlowNode2D<double,NUM_COMPONENTS>::NumEq;k++ )
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

            //Cx,Cy calc
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
            //Set bound Primitives
            // Rects
                double        Xstart,Ystart,X_0,Y_0;
                unsigned int  numRects=Data->GetIntVal((char*)"NumRects");
                if(Data->GetDataError()==-1) Abort_OpenHyperFLOW2D();
                //unsigned int ix0,iy0;
                SolidBoundRect2D* SBR;
                GlobalTime=Data->GetFloatVal((char*)"InitTime");
                if(Data->GetDataError()==-1) Abort_OpenHyperFLOW2D();

                if(p_g==0)
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
                       if(FlowIndex < 1)
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
                       SBR = new SolidBoundRect2D(NameContour,J,Xstart,Ystart,X_0,Y_0,dx,dy,(CondType2D)NT_WNS_2D,pTestFlow2D,Y,TM);
                       *f_stream << "OK\n" << flush;
                       delete SBR;
                     }
                  }
            // Solid Bound Circles
            unsigned int numCircles=Data->GetIntVal((char*)"NumCircles");
            TurbulenceCondType2D TM;
            SolidBoundCircle2D* SBC;
            if ( p_g==0 )
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

                        if ( FlowIndex < 1 ) {
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
                        SBC = new SolidBoundCircle2D(NameContour,J,Xstart,Ystart,X_0,Y_0,dx,dy,(CondType2D)NT_WNS_2D,pTestFlow2D,Y,TM);
                        *f_stream << "OK\n" << flush;
                        delete SBC;
                    }
                }
// Solid Bound Airfoils
/*
            unsigned int numAirfoils=Data->GetIntVal((char*)"NumAirfoils");
            double       mm,pp,thick,scale,attack_angle;
            SolidBoundAirfoil2D* SBA;
            if ( p_g==0 )
                if ( numAirfoils ) {
                    for ( j=0;j<numAirfoils;j++ ) {
                        sprintf(NameContour,"Airfoil%i",j+1);    
                        *f_stream << "Add object \""<< NameContour << "\"..." << flush;

                        sprintf(NameContour,"Airfoil%i.Xstart",j+1);
                        Xstart=Data->GetFloatVal(NameContour);
                        if ( Data->GetDataError()==-1 ) Abort_OpenHyperFLOW2D();

                        sprintf(NameContour,"Airfoil%i.Ystart",j+1);
                        Ystart=Data->GetFloatVal(NameContour);
                        if ( Data->GetDataError()==-1 ) Abort_OpenHyperFLOW2D();

                        sprintf(NameContour,"Airfoil%i.pp",j+1);
                        pp=Data->GetFloatVal(NameContour);
                        if ( Data->GetDataError()==-1 ) Abort_OpenHyperFLOW2D();

                        sprintf(NameContour,"Airfoil%i.mm",j+1);
                        mm=Data->GetFloatVal(NameContour);
                        if ( Data->GetDataError()==-1 ) Abort_OpenHyperFLOW2D();

                        sprintf(NameContour,"Airfoil%i.thick",j+1);
                        thick=Data->GetFloatVal(NameContour);
                        if ( Data->GetDataError()==-1 ) Abort_OpenHyperFLOW2D();

                        sprintf(NameContour,"Airfoil%i.scale",j+1);
                        scale=Data->GetFloatVal(NameContour);
                        if ( Data->GetDataError()==-1 ) Abort_OpenHyperFLOW2D();

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

                        if ( FlowIndex < 1 ) {
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
                        sprintf(NameContour,"Airfoil%i",j+1);
                        SBA = new SolidBoundAirfoil2D(NameContour,J,Xstart,Ystart,mm,pp,thick,dx,dy,(CondType2D)NT_WNS_2D,pTestFlow2D,Y,TM,scale,attack_angle,f_stream);
                        *f_stream << "OK\n" << flush;
                        delete SBA;
                    }
                }
                */
                //  Areas
            if ( !PreloadFlag )
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

                            if ( FlowIndex < 1 ) {
                                *f_stream << "\n";
                                *f_stream << "Bad Flow index [" << FlowIndex << "] \n"<< flush;
                                f_stream->flush();
                                Abort_OpenHyperFLOW2D();
                            }

                            if ( Data->GetDataError()==-1 ) {
                                snprintf(AreaName,256,"Area%i.Flow",i+1);
                                FlowIndex = Data->GetIntVal(AreaName);
                                if ( FlowIndex < 1 ) {
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
                                J->GetValue(i,j).beta  = beta0;
                                
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

                                if (J->GetValue(i,j).isCond2D(CT_WALL_NO_SLIP_2D) || J->GetValue(i,j).isCond2D(CT_WALL_SLIP_2D)) {
                                    J->GetValue(i,j).NGX   = J->GetValue(i,j).idXl - J->GetValue(i,j).idXr + J->GetValue(i,j).idXl*J->GetValue(i,j).idXr;
                                    J->GetValue(i,j).NGY   = J->GetValue(i,j).idYd - J->GetValue(i,j).idYu + J->GetValue(i,j).idYd*J->GetValue(i,j).idYu;
                                }

                                if ( !J->GetValue(i,j).isCond2D(CT_NODE_IS_SET_2D) ) {
                                    *f_stream << "\n";
                                    *f_stream << "Node ("<< i_err << "," << j_err <<") has not CT_NODE_IS_SET flag." << flush;
                                    *f_stream << "\n";
                                    *f_stream << "Possible some \"Area\" objects not defined.\n" << flush;
                                    f_stream->flush();
//-- debug ----
//                                    DataSnapshot(ErrFileName);
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
                *f_stream << "UMatrix2D< FlowNode2D<double,NUM_COMPONENTS> >" << "\n" << flush;
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
            SetWallNodes(f_stream, J);
            GlobalSubmatrix = ScanArea(f_stream,J, isVerboseOutput);
/* Load additional sources */
             isGasSource  = Data->GetIntVal((char*)"NumSrc");
             if ( Data->GetDataError()==-1 ) Abort_OpenHyperFLOW2D();

             if ( isGasSource ) {
                  SrcList = new SourceList2D(J,Data);
                  SrcList->SetSources2D();
             }
             if(!PreloadFlag) {
                *f_stream << "Set initial boundary layer...";
                SetInitBoundaryLayer(J,delta_bl);
                *f_stream << "OK" << endl; 
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

inline int SetTurbulenceModel(UMatrix2D< FlowNode2D<double,NUM_COMPONENTS> >* pJ, int i, int j) {
int      AddEq = 2;
       if (pJ->GetValue(i,j).isTurbulenceCond2D(TCT_Prandtl_Model_2D)) {
           AddEq = 2;
       } else if ( pJ->GetValue(i,j).isTurbulenceCond2D(TCT_k_eps_Model_2D)) {
           AddEq = 0;
       } else if ( pJ->GetValue(i,j).isTurbulenceCond2D(TCT_Spalart_Allmaras_Model_2D)) {
           AddEq = 1;
       } else {
           AddEq = 2;
       }
 return AddEq;
}

inline int CalcChemicalReactions(FlowNode2D<double,NUM_COMPONENTS>* CalcNode,
                                 ChemicalReactionsModel cr_model, void* CRM_data) {
    ChemicalReactionsModelData2D* model_data = (ChemicalReactionsModelData2D*)CRM_data;
    double   Y0,Yfu,Yox,Ycp,Yair;

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
    CalcNode->mu  = model_data->mu_Fuel->GetVal(CalcNode->Tg)*Yfu+
                            model_data->mu_OX->GetVal(CalcNode->Tg)*Yox+
                            model_data->mu_cp->GetVal(CalcNode->Tg)*Ycp+
                            model_data->mu_air->GetVal(CalcNode->Tg)*Yair;
    CalcNode->CP  = model_data->Cp_Fuel->GetVal(CalcNode->Tg)*Yfu+
                            model_data->Cp_OX->GetVal(CalcNode->Tg)*Yox+
                            model_data->Cp_cp->GetVal(CalcNode->Tg)*Ycp+
                            model_data->Cp_air->GetVal(CalcNode->Tg)*Yair;
    CalcNode->lam = model_data->lam_Fuel->GetVal(CalcNode->Tg)*Yfu+
                            model_data->lam_OX->GetVal(CalcNode->Tg)*Yox+
                            model_data->lam_cp->GetVal(CalcNode->Tg)*Ycp+
                            model_data->lam_air->GetVal(CalcNode->Tg)*Yair;

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
                            UArray< XY<int> >* WallNodes2D
                            ,double x0
                            ) {

double  min_l_min = min((FlowNode2D<double,NUM_COMPONENTS>::dx),
                        (FlowNode2D<double,NUM_COMPONENTS>::dy));
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
                        double wx, wy;
                        double x, y;

                        pJ2D->GetValue(i,j).l_min = max((x0+FlowNode2D<double,NUM_COMPONENTS>::dx*pJ2D->GetX()),
                                                           (FlowNode2D<double,NUM_COMPONENTS>::dy*pJ2D->GetY()));
                        
                        x = x0 + i * FlowNode2D<double,NUM_COMPONENTS>::dx;
                        y = j * FlowNode2D<double,NUM_COMPONENTS>::dy;

                        for (int ii=0;ii<(int)WallNodes2D->GetNumElements();ii++) {

                             iw = WallNodes2D->GetElementPtr(ii)->GetX();
                             jw = WallNodes2D->GetElementPtr(ii)->GetY();

                             wx = iw * FlowNode2D<double,NUM_COMPONENTS>::dx;
                             wy = jw * FlowNode2D<double,NUM_COMPONENTS>::dy;

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



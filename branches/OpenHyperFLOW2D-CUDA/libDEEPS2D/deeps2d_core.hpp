/*******************************************************************************
*   OpenHyperFLOW2D-CUDA                                                       *
*                                                                              *
*   Transient, Density based Effective Explicit Parallel Hybrid Solver         *
*   TDEEPHS (CUDA+MPI)                                                         *
*   Version  1.0.1                                                             *
*   Copyright (C)  1995-2014 by Serge A. Suchkov                               *
*   Copyright policy: LGPL V3                                                  *
*                                                                              *
* Common function declarations file.                                           *
*                                                                              *
*  last update: 16/01/2014                                                     *
********************************************************************************/
#ifdef _MPI
#include <mpi.h>
#define _PARALLEL_ONLY
#endif // _MPI

#ifdef _OPENMP
#include <omp.h>
#define _PARALLEL_ONLY
#endif //OPENMP

#ifdef _CUDA_
#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#endif //_CUDA_

#ifdef _NO_MMAP_SHARED_
#ifndef _WRITE_LARGE_FILE_
#define _WRITE_LARGE_FILE_
#endif //_WRITE_LARGE_FILE_
#else
#ifdef _NO_MMAP_
#ifndef _WRITE_LARGE_FILE_
#define _WRITE_LARGE_FILE_
#endif //_WRITE_LARGE_FILE_
#else
#ifndef _MSYNC_LARGE_FILE_
#define _MSYNC_LARGE_FILE_
#endif //_MSYNC_LARGE_FILE_
#endif //_NO_MMAP_
#endif // _NO_MMAP_SHARED_

#include "utl/umatrix2d.hpp"
#include "obj_data/obj_data.hpp"
#include "libOpenHyperFLOW2D/hyper_flow2d.hpp"
#include "libOutCFD/out_cfd_param.hpp"

#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/mman.h>

#ifdef _DEBUG_0
#include "libExcept/except.hpp"
#endif // _DEBUG_0

#ifndef ComputationalMatrix2D
#define ComputationalMatrix2D  UMatrix2D< FlowNode2D<double,NUM_COMPONENTS> >
#endif //ComputationalMatrix2D

enum BlendingFactorFunction {
     BFF_L,   // LINEAR locally adopted blending factor function (LLABF)
     BFF_LR,  // LINEAR locally adopted blending factor function with relaxation (LLABFR)
     BFF_S,   // SQUARE locally adopted blending factor function (SLABFR)
     BFF_SR,  // SQUARE locally adopted blending factor function with relaxation (SLABFR)
     BFF_SQR, // SQRT() locally adopted blending factor function (SRLABF)
     BFF_SQRR,// SQRT() locally adopted blending factor function with relaxation (SRLABFR)
     BFF_MACH,// Mach number depended locally adapted blending factor function (MLABF)
     BFF_LG,  // Local gradient  adopted blending factor function (LGABF)
     BFF_MIXED, // BFF_SQR + BFF_LG
     BFF_HYBRID, // BFF_SQR + BFF_MACH + BFF_LG
     BFF_SQR_PRESSURE, // BFF_SQR + pressure check
     BFF_SR_LIMITED,
};

enum WRITE_MODE {
     WM_APPEND,
     WM_REWRITE
};

enum data_tag {
    tag_MaxX=4000,
    tag_MaxY,
    tag_X0,
    tag_Matrix,
    tag_MatrixRow,
    tag_MatrixTail,
    tag_MatrixHead,
    tag_NumWallNodes,
    tag_WallNodesArray,
    tag_WallFrictionVelocity,
    tag_DD
};

#ifdef _MPI_
struct DD_pack {
          double         DD;   // max residual
          double         RMS;  // sum residual
          unsigned long  iRMS; // num involved nodes
          unsigned int   i,j;  // x,y -coordinates
};
#endif // _MPI_

extern int    fd_g;
extern int    fd_s;
extern int    fd_l;
extern double delta_bl;
extern void*  SolidSwapData;
extern void*  GasSwapData;
extern UArray< UMatrix2D< FlowNode2D<double,NUM_COMPONENTS> >* >* ArraySubmatrix;
extern ChemicalReactionsModelData2D chemical_reactions;
extern UArray< XY<int> >* WallNodes;
extern UArray< XY<int> >* GetWallNodes(ofstream* f_str, ComputationalMatrix2D* pJ, int isPrint);

extern UArray< double >*     WallNodesUw_2D;
extern int                   NumWallNodes;
extern double                x0;
// External variables
extern SolverMode ProblemType;
extern UArray< XY<int> >*                                           GlobalSubmatrix;
extern UArray<UMatrix2D< FlowNode2D<double,NUM_COMPONENTS> >*>*     SubmatrixArray;
extern UArray<UMatrix2D< FlowNodeCore2D<double,NUM_COMPONENTS> >*>* CoreSubmatrixArray;
extern InputData*                            Data;                     // Object data loader
extern UMatrix2D< FlowNode2D<double,NUM_COMPONENTS> >*     J;          // Main computation area
extern UMatrix2D< FlowNodeCore2D<double,NUM_COMPONENTS> >* C;          // Main core area
extern UArray<Flow*>*                        FlowList;                 // Flow list
extern unsigned int                          MaxX;
extern unsigned int                          MaxY;
extern double                                dx;
extern double                                dy;
extern int                                   EndFLAG;
extern int                                   PrintFLAG;
extern int                                   isVerboseOutput;
extern unsigned int                          iter;
extern SourceList2D*                         SrcList;
extern int                                   TurbMod;
extern BlendingFactorFunction                bFF;
extern SourceList2D*                         SrcList;
extern int                                   isGasSource;
extern int                                   TurbStartIter;
extern int                                   TurbExtModel;
extern int                                   err_i, err_j;
extern int                                   turb_mod_name_index;
extern double                                Ts0,A,W,Mach;
extern double                                SigW,SigF,delta_bl;
extern int                                   isVerboseOutput;
extern int                                   isTurbulenceReset;
extern double                                CFL;
extern Table*                                CFL_Scenario;
extern Table*                                beta_Scenario;
extern int                                   NSaveStep;
extern int                                   NOutStep;
extern int                                   isFirstStart;
extern int                                   ScaleFLAG;
extern int                                   isAdiabaticWall;
extern int                                   isOutHeatFluxX;
extern int                                   isOutHeatFluxY;
extern int                                   is_p_asterisk_out;
extern int                                   Nstep;
extern double                                ExitMonitorValue;
extern int                                   MonitorNumber;
extern int                                   MonitorCondition;
extern unsigned int                          AddSrcStartIter;
extern double                                beta[6+NUM_COMPONENTS];
extern double                                beta0;
extern int                                   NumContour;
extern char*                                 ProjectName;
extern char                                  GasSwapFileName[255];
extern char                                  ErrFileName[255];
extern char                                  OutFileName[255];
extern char                                  TecPlotFileName[255];
extern int                                   fd_g;
extern void*                                 GasSwapData;
extern int                                   TurbMod;
extern int                                   EndFLAG;
extern int                                   PrintFLAG;
extern ofstream*                             pInputData;        // Output data stream
extern ofstream*                             pHeatFlux_OutFile; // Output HeatFlux stream
extern unsigned int                          iter;              // iteration number
extern unsigned int                          last_ite;          // Global iteration number

extern int                                   isStop;            // Stop flag
extern InputData*                            Data;              // Object data loader
extern UArray<Flow*>*                        FlowList;          // List of 'Flow' objects
extern UArray<Flow2D*>*                      Flow2DList;        // List of 'Flow2D' objects
extern UArray<Bound2D*>                      SingleBoundsList;  // Single Bounds List;

extern double                                dt;                // time step
extern double*                               RoUx;
extern double*                               RoVy;
extern double                                GlobalTime;
extern double                                CurrentTimePart;

extern double*                               Y;
extern double                                Cp;
extern double                                lam;
extern double                                mu;
extern double                                Tg;
extern double                                Rg;
extern double                                Pg;
extern double                                Wg;
extern double                                Ug,Vg;
extern int                                   CompIndex;

extern double                                Y_fuel[4];  /* fuel */
extern double                                Y_ox[4]  ;  /* OX */
extern double                                Y_cp[4]  ;  /* cp */
extern double                                Y_air[4] ;  /* air */
extern double                                Y_mix[4] ;  /* mixture */
                                              
#ifdef _RMS_
extern char                                  RMSFileName[255];
extern ofstream*                             pRMS_OutFile;    // Output RMS stream
#endif //_RMS_
extern int                                   useSwapFile;
extern char                                  OldSwapFileName[255];
extern void*                                 OldSwapData;
extern u_long                                OldFileSizeGas;
extern int                                   Old_fd;
extern int                                   isScan;
extern int                                   MemLen;
extern int                                   isRun;
extern int                                   isDraw;
extern int                                   isSnapshot;
extern UArray<XCut>*                         XCutArray;
extern unsigned int                          last_iter;     // Global iteration number
extern int                                   is_Cx_calc;
extern double                                x0_body,y0_body,dx_body,dy_body;
extern int                                   Cx_Flow_index;
extern int                                   SigmaFi_Flow_index;
extern int                                   num_threads;
extern int                                   num_blocks;
extern ofstream*                             pOutputData;     // output data stream (file)
//extern double                                MemFullness;
extern unsigned int*                         dt_min_host;
extern unsigned int*                         dt_min_device;
extern UArray< unsigned int* >*              dt_min_host_Array;
extern UArray< unsigned int* >*              dt_min_device_Array;
extern double*                               cudaHu;  
extern UArray< double* >*                    cudaHuArray;  
extern UArray< FlowNode2D<double,NUM_COMPONENTS>* >*     cudaArraySubmatrix;
extern UArray< FlowNodeCore2D<double,NUM_COMPONENTS>* >* cudaArrayCoreSubmatrix ;
extern UArray< XY<int> >*                                cudaDimArray;
extern XY<int>*                                          cudaWallNodes;
extern FlowNode2D<double,NUM_COMPONENTS>*                cudaSubmatrix;
extern FlowNodeCore2D<double,NUM_COMPONENTS>*            cudaCoreSubmatrix;
extern int                                               warp_size;
extern int                                               max_num_threads;
extern int                                               num_gpus;  // number of CUDA GPUs

// External functions

#ifdef _CUDA_
extern void                                  CopyDeviceToHost(void* src, void* dst, size_t length, cudaStream_t stream=0);
extern void                                  CopyHostToDevice(void* src, void* dst, size_t length, cudaStream_t stream=0);
extern void                                  CopyDeviceToDevice(void* src, void* dst, size_t length, cudaStream_t stream=0);
extern void                                  CopyDeviceToDeviceP2P(void* src, int src_dev,void* dst, int dst_dev,size_t length, cudaStream_t stream=0);
extern void                                  SetP2PAccess(int dev1, int dev2);
extern void                                  DisableP2PAccess(int dev1, int dev2);
 
#endif //_CUDA_ 
extern const char*                           PrintTurbCond(int TM);
extern void*                                 InitDEEPS2D(void*);
extern void                                  SaveData2D(ofstream* OutputData, int);
extern ofstream*                             OpenData(char* outputDataFile);
extern void                                  SaveRMSHeader(ofstream* OutputData);
extern void                                  SaveRMS(ofstream* OutputData,unsigned int n, double* outRMS);
extern void                                  CutFile(char* cutFile);
extern u_long                                SetWallNodes(ofstream* f_str, ComputationalMatrix2D* pJ);
#ifdef _IMPI_
extern void                                  LongMatrixSend(int rank, void* src,  size_t len);
extern void                                  LongMatrixRecv(int rank, void* dst,  size_t len);
#endif // _IMPI_
extern void                                  SetInitialSources(UMatrix2D< FlowNode2D<double,NUM_COMPONENTS> >* pJ);
extern void                                  SetInitBoundaryLayer(ComputationalMatrix2D* pJ, double delta);

#ifdef _CUDA_
 __host__ __device__
#endif //_CUDA_ 
extern int SetTurbulenceModel(FlowNode2D<double,NUM_COMPONENTS>* pJ);
extern void DataSnapshot(char* filename, WRITE_MODE ioMode=WM_REWRITE);
extern void CalcHeatOnWallSources(UMatrix2D< FlowNode2D<double,NUM_COMPONENTS> >* F, double dx, double dr, double dt, int rank, int last_rank);
extern UArray< XY<int> >* ScanArea(ofstream* f_str,ComputationalMatrix2D* pJ ,int isPrint);


int CalcChemicalReactions(FlowNode2D<double,NUM_COMPONENTS>* CalcNode,
                          ChemicalReactionsModel cr_model, void* CRM_data);

void RecalcWallFrictionVelocityArray2D(ComputationalMatrix2D* pJ,
                                       UArray<double>* WallFrictionVelocityArray2D,
                                       UArray< XY<int> >* WallNodes2D);

UArray<double>* GetWallFrictionVelocityArray2D(ComputationalMatrix2D* pJ, 
                                               UArray< XY<int> >* WallNodes2D);

#ifdef _PARALLEL_RECALC_Y_PLUS_
void ParallelRecalc_y_plus(ComputationalMatrix2D* pJ, 
                           UArray< XY<int> >* WallNodes,
                           UArray<double>* WallFrictionVelocity2D,
                           double x0);
#else
extern void Recalc_y_plus(ComputationalMatrix2D* pJ, UArray< XY<int> >* WallNodes);
#endif //_PARALLEL_RECALC_Y_PLUS_
 
extern void SetMinDistanceToWall2D(ComputationalMatrix2D* pJ2D,UArray< XY<int> >* WallNodes2D, double x0=0.0);
extern int BuildMesh(int mode);
extern void InitSharedData(InputData*, void*
#ifdef _MPI
                           ,int
#endif // _MPI
                           );

extern void DEEPS2D_Run(ofstream* o_stream 
#ifdef _MPI
                        ,UMatrix2D< FlowNode2D<double,NUM_COMPONENTS> >*    pJ,
                        UMatrix2D< FlowNodeCore2D<double,NUM_COMPONENTS> >* pC,
                        int rank , int last_rank
#endif // _MPI
                        );

#ifdef _CUDA_
extern __host__ int LoadTable2GPU(Table* Src, Table*& Dst, int i_dev);


extern void DEEPS2D_Run(ofstream* f_stream, 
                        UMatrix2D<FlowNode2D<double,NUM_COMPONENTS> >*     pJ,
                        UMatrix2D<FlowNodeCore2D<double,NUM_COMPONENTS> >* pC,
                        UArray< FlowNode2D<double,NUM_COMPONENTS>* >*      cudaSubmatrix,
                        UArray< FlowNodeCore2D<double,NUM_COMPONENTS>* >*  cudaCoreSubmatrix,
                        UArray< XY<int> >*                                 cudaDimArray,
                        UArray< ChemicalReactionsModelData2D* >*           cudaCRM2D,
                        int                                                num_mp,
                        cudaStream_t*                                      cuda_streams,
                        cudaEvent_t*                                       cuda_events,
                        int                                                max_thread_block);

extern __global__  void  
cuda_DEEPS2D_Stage1(FlowNode2D<double,NUM_COMPONENTS>*     pLJ,
                    FlowNodeCore2D<double,NUM_COMPONENTS>* pLC,
                    unsigned long int index_limit,
                    int MAX_X, int MAX_Y, 
                    unsigned long r_limit,
                    unsigned long l_limit,
                    double dxx, double dyy,
                    double dtdx, double dtdy,
                    double _dt,
                    int _FT, int Num_Eq, SolverMode pt);

extern __global__  void 
cuda_DEEPS2D_Stage2(FlowNode2D<double,NUM_COMPONENTS>*     pLJ,
                    FlowNodeCore2D<double,NUM_COMPONENTS>* pLC,
                    unsigned long int index_limit,
                    int MAX_X, int MAX_Y,
                    unsigned long r_limit,
                    unsigned long l_limit,
                    double beta_init, double  beta0, 
                    int b_FF, double CFL0,
                    ChemicalReactionsModelData2D* pCRMD,
                    int noTurbCond,
                    double SigW, double SigF, double dx_1, double dy_1, double delta_bl,
                    FlowType _FT, int Num_Eq,
#ifdef _RMS_             
                    double*  RMS, 
                    int*     iRMS,
                    double   DD_max,
                    int*     i_c,
                    int*     j_c,
#endif // _RMS_
                    double* _Hu,
                    int _isSrcAdd,
                    unsigned int* dt_global, double int2float_scale,
                    TurbulenceExtendedModel TurbExtModel, 
                    SolverMode pt);

extern __global__ void 
cuda_SetInitBoundaryLayer(FlowNode2D<double,NUM_COMPONENTS>* pJ2D,
                          unsigned long int index_limit,
                          int X0, int MAX_Y,
                          double delta,
                          double sig_w, 
                          double sig_f,
                          TurbulenceExtendedModel etm,
                          double _dx, double _dy,
                          double* _Hu,
                          int _isSrcAdd,
                          FlowType _FT,
                          SolverMode sm);

extern __global__ void
cuda_SetMinDistanceToWall2D(FlowNode2D<double,NUM_COMPONENTS>* pJ2D,
                            unsigned long int index_limit,
                            XY<int>* WallNodes2D, 
                            int NumWallNodes2D,
                            double min_l_min,
                            double max_l_min,
                            double _dx, double _dy);

extern __global__ void 
cuda_Recalc_y_plus(FlowNode2D<double,NUM_COMPONENTS>* pJ2D,
                   unsigned long int index_limit,
                   XY<int>* WallNodes2D, 
                   int NumWallNodes2D,
                   double min_l_min,
                   double max_l_min,
                   double _dx, 
                   double _dy,
                   int max_y);

extern void CUDA_BARRIER(char* KernelName);
#endif // _CUDA_

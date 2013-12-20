/*******************************************************************************
*   OpenHyperFLOW2D                                                            *
*                                                                              *
*   Transient, Density based Effective Explicit Parallel Solver (T-DEEPS2D)    *
*                                                                              *
*   Version  1.0.0                                                             *
*   Copyright (C)  1995-2013 by Serge A. Suchkov                               *
*   Copyright policy: LGPL V3                                                  *
*                                                                              *
* Common function declarations file.                                           *
*                                                                              *
*  last update: 15/12/2013                                                     *
********************************************************************************/
#ifdef _MPI
#include <mpi.h>
#define _PARALLEL_ONLY
#endif // _MPI

#ifdef _OPENMP
#include <omp.h>
#define _PARALLEL_ONLY
#endif //OPENMP

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

struct DD_pack {
          double         DD;   // max residual
          double         RMS;  // sum residual
          unsigned long  iRMS; // num involved nodes
          unsigned int   i,j;  // x,y -coordinates
};
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
extern UArray< XY<int> >*                                           GlobalSubmatrix;
extern UArray<UMatrix2D< FlowNode2D<double,NUM_COMPONENTS> >*>*     SubmatrixArray;
extern UArray<UMatrix2D< FlowNodeCore2D<double,NUM_COMPONENTS> >*>* CoreSubmatrixArray;
extern InputData*                            Data;                     // Object data loader
extern UMatrix2D< FlowNode2D<double,NUM_COMPONENTS> >*     J;          // Main computation area
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
// External functions
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
extern void SetInitialSources(UMatrix2D< FlowNode2D<double,NUM_COMPONENTS> >* pJ);
extern void SetInitBoundaryLayer(ComputationalMatrix2D* pJ, double delta);
extern int  SetTurbulenceModel(UMatrix2D< FlowNode2D<double,NUM_COMPONENTS> >* pJ, int i, int j);
extern void DataSnapshot(char* filename, WRITE_MODE ioMode=WM_REWRITE);
extern void CalcHeatOnWallSources(UMatrix2D< FlowNode2D<double,NUM_COMPONENTS> >* F, double dx, double dr, double dt, int rank, int last_rank);
extern UArray< XY<int> >* ScanArea(ofstream* f_str,ComputationalMatrix2D* pJ ,int isPrint);
extern int CalcChemicalReactions(FlowNode2D<double,NUM_COMPONENTS>* CalcNode,
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
                        ,UMatrix2D< FlowNode2D<double,NUM_COMPONENTS> >*     pJ,
                        UMatrix2D< FlowNodeCore2D<double,NUM_COMPONENTS> >* pC,
                        int rank , int last_rank
#endif // _MPI
                        );


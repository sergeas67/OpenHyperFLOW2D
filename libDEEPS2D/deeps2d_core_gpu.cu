/*******************************************************************************
*   OpenHyperFLOW2D-CUDA                                                       *
*                                                                              *
*   Transient, Density based Effective Explicit Parallel Hybrid Solver         *
*   TDEEPHS (CUDA+MPI)                                                         *
*   Version  2.0.1                                                             *
*   Copyright (C)  1995-2016 by Serge A. Suchkov                               *
*   Copyright policy: LGPL V3                                                  *
*   http://github.com/sergeas67/openhyperflow2d                                *
*                                                                              *
*   deeps2d_core.cpp: OpenHyperFLOW2D solver core code....                     *
*                                                                              *
*  last update: 14/04/2016                                                     *
********************************************************************************/
#include "deeps2d_core.hpp"

#include <sys/time.h>
#include <sys/timeb.h>
#include <sys/file.h>

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
                             FP&  opt_round_trip,
                             FP   round_trip) {
   
   if(ThreadBlockSize == 0) {
       if(opt_round_trip == 0.0) {
            opt_block_size[0] = opt_block_size[1] =  cur_block_size;
            opt_round_trip    =  round_trip;
       } else if (opt_round_trip >= round_trip) {
            opt_block_size[0] = cur_block_size;
            opt_block_size[1] = cur_block_size/2;
            opt_round_trip    = round_trip;
       } else {
            opt_block_size[1] = cur_block_size/2; 
       }
         return opt_block_size[1];
   } else {
         return ThreadBlockSize;
   }
};


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
// CFL_Scenario->GetVal(iter+last_iter)

#ifdef _CUDA_
void DEEPS2D_Run(ofstream* f_stream, 
                 UMatrix2D<FlowNode2D<FP,NUM_COMPONENTS> >*     pJ,
                 UMatrix2D<FlowNodeCore2D<FP,NUM_COMPONENTS> >* pC,
                 UArray< FlowNode2D<FP,NUM_COMPONENTS>* >*      cudaSubDomainArray,
                 UArray< FlowNodeCore2D<FP,NUM_COMPONENTS>* >*  cudaCoreSubDomainArray,
                 UArray< XY<int> >*                             cudaDimArray,
                 UArray< XY<int>* >*                            cudaWallNodesArray,
                 UArray< ChemicalReactionsModelData2D* >*       cudaCRM2D,
                 int num_mp,
                 cudaStream_t *cuda_streams,
                 cudaEvent_t  *cuda_events) {

   // local variables
    char Eff_warp_size[64];
    char Time_step[64];
    FlowNode2D<FP,NUM_COMPONENTS>* TmpMatrixPtr;
    //FlowNodeCore2D<FP,NUM_COMPONENTS>* TmpCoreMatrixPtr;
    //int    SubMaxX, SubStartIndex; 
    int      TmpMaxX;
    int      current_div_stage1=1;
    int      current_div_stage2=1;
    int      current_div_heat=1;
    int      I;
    int      n_s;
    unsigned int   int2float_var;

    int      opt_thread_block_size_stage1[2];
    int      opt_thread_block_size_stage2[2];
    int      opt_thread_block_size_heat[2];
    
    FP       opt_round_trip_stage1 = 0.0;
    FP       opt_round_trip_stage2 = 0.0;
    FP       opt_round_trip_heat   = 0.0;

    int       r_Overlap;
    int       l_Overlap;
    
    int      num_cuda_threads;
    int      num_cuda_blocks; 

    FP       dtdx;
    FP       dtdy;
    FP       dyy;
    FP       dxx;
    FP       dx_1,dy_1; // 1/dx, 1/dy
    FP       d_time_stage1;
    FP       d_time_stage2;
    FP       d_time_heat;
    FP       VCOMP;
    FP       d_time;
    timeval  start, stop, mark1, mark2, start1, end;
    FP       int2float_scale;
#ifdef _RMS_
    FP       int2float_RMS_scale;
#endif // _RMS_
    FP*      dt_min;
    int*     i_c;
    int*     j_c;
    FP       dtmin=1.0;
    //FP     DD[FlowNode2D<FP,NUM_COMPONENTS>::NumEq];

#ifndef _P2P_ACCESS_
    void* host_HalloBuff[8];
#endif //_P2P_ACCESS_ 
    isScan = 0;
    dyy    = dx/(dx+dy);
    dxx    = dy/(dx+dy);

    d_time = d_time_stage1  = d_time_stage2 = d_time_heat = 0.;
    
    dt_min = new FP[cudaArraySubDomain->GetNumElements()];
    i_c    = new int[cudaArraySubDomain->GetNumElements()];
    j_c    = new int[cudaArraySubDomain->GetNumElements()];

    for(int ii=0;ii<(int)cudaArraySubDomain->GetNumElements();ii++) {
        dt_min[ii] = dtmin = dt;
        i_c[ii] = j_c[ii] = 0;
#ifndef _P2P_ACCESS_
        *f_stream << "\nAlloc " << pJ->GetColSize() << " bytes in host memory for Hallo Exchange Buffer...";
        host_HalloBuff[ii]   =  malloc(pJ->GetColSize());
        if(host_HalloBuff[ii]) {
          *f_stream << "OK" << endl;
        } else {
          *f_stream << "ERROR!" << endl;
          Exit_OpenHyperFLOW2D(n_s);
        }
#endif // _P2P_ACCESS_
     }
#ifdef _RMS_    
    FP                    max_RMS;   
    int                   k_max_RMS; 
    cudaError_t           cudaState;
    
    UArray<RMS_pack*>     device_RMS_Array;
    UArray<RMS_pack*>     host_RMS_Array;
    
    RMS_pack*             cuda_RMS_pack;
    RMS_pack*             host_RMS_pack;
    RMS_pack              zero_RMS_pack;
    RMS_pack              tmp_RMS_pack;
    
    memset(&zero_RMS_pack,0,sizeof(RMS_pack));
    n_s = cudaDimArray->GetNumElements();

    for(int ii=0;ii<n_s;ii++) {    
        
        int igpu;
        if(isSingleGPU) {
            igpu = ActiveSingleGPU;
        } else {
            igpu = ii;
        }

        if(cudaSetDevice(igpu) != cudaSuccess ) {
           *f_stream << "\nError set CUDA device no: "<< igpu << endl;
           Exit_OpenHyperFLOW2D(n_s);
        }
    
        cudaState = cudaMalloc((void**)&cuda_RMS_pack, sizeof(RMS_pack));
        if (cudaState == cudaErrorMemoryAllocation) {
            *f_stream << "\nError allocate GPU memory for sum_RMS[] "<<  endl;
            Exit_OpenHyperFLOW2D(n_s);
        }
        
        host_RMS_pack   =  new RMS_pack();  
        
        device_RMS_Array.AddElement(&cuda_RMS_pack);
        host_RMS_Array.AddElement(&host_RMS_pack);
    }

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
                    current_div_heat = current_div_stage1 = current_div_stage2  =  dprop.multiProcessorCount*warp_size;
                    if(ThreadBlockSize == 0)
                        isCalibrate = 1;
            do {

                  gettimeofday(&start,NULL);

                  if( AddSrcStartIter < iter + last_iter) {
                   FlowNode2D<FP,NUM_COMPONENTS>::isSrcAdd = 1;
                  } else {
                   FlowNode2D<FP,NUM_COMPONENTS>::isSrcAdd = 0;
                  }
                  
                  n_s = cudaDimArray->GetNumElements();

                  iter = 0;

                  do {
                      if ( isVerboseOutput && iter/NOutStep*NOutStep == iter ) {
                           gettimeofday(&start1,NULL);
                      }
#ifdef _RMS_
                       max_RMS   = 0;
                       k_max_RMS = -1;
                       
                       int2float_RMS_scale = (FP)((INT_MAX) >> sizeof(unsigned int)*4)/dt;
#endif // _RMS_
                       if (isLocalTimeStep) {
                           int2float_scale  = 1;
                       } else {
                           int2float_scale  = (FP)((INT_MAX) >> sizeof(unsigned int)*4)/dt;
                       }
                       
                       unsigned int dtest_int = (unsigned int)(int2float_scale*10);

                       int iX0=0;

                       dtmin = 10*dt;
                       
                       if (isCalibrate) {
                           gettimeofday(&mark2,NULL);
                       }
#pragma unroll
                       for(int ii=0;ii<n_s;ii++) {  // CUDA version

                       int max_X = cudaDimArray->GetElement(ii).GetX();
                       int max_Y = cudaDimArray->GetElement(ii).GetY();
                       
                       int i_gpu;
                       if(isSingleGPU) {
                           i_gpu = ActiveSingleGPU;
                       } else {
                           i_gpu = ii;
                       }

                       if(cudaSetDevice(i_gpu) != cudaSuccess ) {
                          *f_stream << "\nError set CUDA device no: "<< i_gpu << endl;
                          Exit_OpenHyperFLOW2D(n_s);
                       }
#ifdef _DEVICE_MMAP_
                       dt_min_host = dt_min_host_Array->GetElement(ii);
                       *dt_min_host = dtest_int;
#else
                       dt_min_device = dt_min_device_Array->GetElement(ii);
                       CopyHostToDevice(&dtest_int,dt_min_device,sizeof(unsigned int));
#endif //_DEVICE_MMAP_
                       
                       cudaJ = cudaSubDomainArray->GetElement(ii);
                       cudaC = cudaCoreSubDomainArray->GetElement(ii);

                       if (ThreadBlockSize == 0) {
                           num_cuda_threads =  dprop.multiProcessorCount*warp_size/max(1,current_div_stage1);
                       } else {
                           num_cuda_threads = ThreadBlockSize;
                       }
                       
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
                       
                       dx_1 = 1.0/dx;
                       dy_1 = 1.0/dy;

                       dtdx = dt/dx;
                       dtdy = dt/dy;

#ifdef _RMS_
                       CopyHostToDevice(&zero_RMS_pack,device_RMS_Array.GetElement(ii),sizeof(RMS_pack),cuda_streams[ii]);
                       memset(&tmp_RMS_pack,0,sizeof(RMS_pack));
#endif // _RMS_
                       cuda_DEEPS2D_Stage1<<<num_cuda_blocks,num_cuda_threads, 0, cuda_streams[ii]>>>(cudaJ,cudaC,
                                                                                                      max_X*max_Y,
                                                                                                      iX0 + max_X - r_Overlap,
                                                                                                      iX0 + l_Overlap,
                                                                                                      dxx,dyy,dtdx,dtdy,dt,
                                                                                                      FlowNode2D<FP,NUM_COMPONENTS>::FT,
                                                                                                      FlowNode2D<FP,NUM_COMPONENTS>::NumEq-ConvertTurbMod(TurbMod),
                                                                                                      ProblemType,
                                                                                                      isLocalTimeStep);
                       iX0 += max_X;
                   }
#ifdef _BARRIER_
                   CUDA_BARRIER((char*)"cuda_DEEPS2D_Stage1");
#endif // _BARRIER_
                   if (isCalibrate) {
                       gettimeofday(&mark1,NULL);
                       d_time_stage1 = (FP)(mark1.tv_sec-mark2.tv_sec)+(FP)(mark1.tv_usec-mark2.tv_usec)*1.e-6; 
                       gettimeofday(&mark2,NULL);
                   }
                   iX0 = 0;
#pragma unroll
                   for(int ii=0;ii<n_s;ii++) {  // CUDA version

                       int i_gpu;
                       if(isSingleGPU) {
                           i_gpu = ActiveSingleGPU;
                       } else {
                           i_gpu = ii;
                       }

                       if(cudaSetDevice(i_gpu) != cudaSuccess ) {
                          *f_stream << "\nError set CUDA device no: "<< i_gpu << endl;
                          Exit_OpenHyperFLOW2D(n_s);
                       }
                       int max_X = cudaDimArray->GetElement(ii).GetX();
                       int max_Y = cudaDimArray->GetElement(ii).GetY();

                       cudaJ = cudaSubDomainArray->GetElement(ii);
                       cudaC = cudaCoreSubDomainArray->GetElement(ii);
                       
                       if (ThreadBlockSize == 0) {
                           num_cuda_threads = min(max_num_threads,dprop.multiProcessorCount*warp_size/max(1,current_div_stage2));
                       } else {
                           num_cuda_threads = min(max_num_threads,ThreadBlockSize);
                       }
                       
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
                                                                                                       iX0 + max_X - r_Overlap,
                                                                                                       iX0 + l_Overlap,
                                                                                                       beta_Scenario->GetVal(iter+last_iter),beta0,
                                                                                                       bFF, min(CFL0,CFL_Scenario->GetVal(iter+last_iter)),
                                                                                                       nrbc_beta0,
                                                                                                       cudaCRM2D->GetElement(ii),
                                                                                                       noTurbCond,
                                                                                                       SigW,SigF,dx_1,dy_1,delta_bl,dx,dy,
                                                                                                       FlowNode2D<FP,NUM_COMPONENTS>::FT,
                                                                                                       FlowNode2D<FP,NUM_COMPONENTS>::NumEq-ConvertTurbMod(TurbMod),
#ifdef _RMS_
                                                                                                       isAlternateRMS,
                                                                                                       device_RMS_Array.GetElement(ii),
                                                                                                       int2float_RMS_scale,    
#endif //_RMS_
                                                                                                       cudaHu,
                                                                                                       FlowNode2D<FP,NUM_COMPONENTS>::isSrcAdd,
                                                                                                       dt_min_device, int2float_scale,
                                                                                                       (TurbulenceExtendedModel)TurbExtModel,
                                                                                                       ProblemType,
                                                                                                       isLocalTimeStep);
                        
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
#ifdef _BARRIER_
                   CUDA_BARRIER((char*)"cuda_DEEPS2D_Stage2");
#endif // _BARRIER_
                    
                   if(isCalibrate) {
                        gettimeofday(&mark1,NULL);
                        d_time_stage2 = (FP)(mark1.tv_sec-mark2.tv_sec)+(FP)(mark1.tv_usec-mark2.tv_usec)*1.e-6; 
                    }
#pragma unroll
                   for(int ii=0;ii<n_s;ii++) {  // CUDA version

                       int i_gpu;
                       
                       if(isSingleGPU) {
                           i_gpu = ActiveSingleGPU;
                       } else {
                           i_gpu = ii;
                       }

                       if(cudaSetDevice(i_gpu) != cudaSuccess ) {
                          *f_stream << "\nError set CUDA device no: "<< i_gpu << endl;
                          Exit_OpenHyperFLOW2D(n_s);
                       
                       }

#ifdef _RMS_           
                       CopyDeviceToHost(device_RMS_Array.GetElement(ii),host_RMS_Array.GetElement(ii),sizeof(RMS_pack),cuda_streams[ii]);
#endif //_RMS_
// --- Halo exchange ---
                       
                       if(ii == n_s-1)
                         r_Overlap = 0;
                       else
                         r_Overlap = 1;
                       if(ii == 0)
                         l_Overlap = 0;
                       else
                         l_Overlap = 1;

                       if(r_Overlap > 0) {
                       
                           
                       int max_X = cudaDimArray->GetElement(ii).GetX();
                       int max_Y = cudaDimArray->GetElement(ii).GetY();

                       cudaJ = cudaSubDomainArray->GetElement(ii);
                       size_t cuda_HalloSize = pJ->GetColSize();

                       void*  cuda_Src  = (void*)((ulong)cudaJ+(max_X*max_Y*sizeof(FlowNode2D<FP,NUM_COMPONENTS>) - 2*cuda_HalloSize));
                       void*  cuda_Dst  = (void*)(cudaArraySubDomain->GetElement(ii+1));
#ifdef _P2P_ACCESS_

                       CopyDeviceToDeviceP2P(cuda_Src,
                                             ii,
                                             cuda_Dst,
                                             ii+1,
                                             cuda_HalloSize,
                                             NULL,
                                             cuda_streams[ii]);

#else
                       CopyDeviceToDeviceP2P(cuda_Src,
                                             ii,
                                             cuda_Dst,
                                             ii+1,
                                             cuda_HalloSize,
                                             host_HalloBuff[ii],
                                             cuda_streams[ii]);

#endif // _P2P_ACCESS_
                       cuda_Dst  = (void*)((ulong)cuda_Src + cuda_HalloSize);
                       cuda_Src  = (void*)((ulong)cudaArraySubDomain->GetElement(ii+1)+cuda_HalloSize);

                       if(cudaSetDevice(ii+1) != cudaSuccess ) {
                          *f_stream << "\nError set CUDA device no: "<< ii << endl;
                          Exit_OpenHyperFLOW2D(n_s);
                       }
#ifdef _P2P_ACCESS_

                       CopyDeviceToDeviceP2P(cuda_Src,
                                             ii+1,
                                             cuda_Dst,
                                             ii,
                                             cuda_HalloSize,
                                             NULL,
                                             cuda_streams[ii+1]);
#else
                       CopyDeviceToDeviceP2P(cuda_Src,
                                             ii+1,
                                             cuda_Dst,
                                             ii,
                                             cuda_HalloSize,
                                             host_HalloBuff[ii],
                                             cuda_streams[ii+1]);

#endif // _P2P_ACCESS_
// --- Halo exchange ---
                       }
#ifdef _DEVICE_MMAP_
                      dt_min_host = dt_min_host_Array->GetElement(ii);
                      unsigned int  int2float_var = *dt_min_host;
#else
                      
                      dt_min_device = dt_min_device_Array->GetElement(ii);
                      CopyDeviceToHost(dt_min_device,&int2float_var,sizeof(unsigned int),cuda_streams[ii]);
#endif //_DEVICE_MMAP_
                      if (isLocalTimeStep) {
                          dtmin  = int2float_var;
                      } else {
                          dtmin  =  min(dtmin,(FP)(int2float_var)/int2float_scale);
                      }
                     
                      if(dtmin == 0 ) {
                         *f_stream << "\nERROR: Computational unstability on iteration " << iter+last_iter << " in domain " << ii <<  endl;
                         Abort_OpenHyperFLOW2D(n_s);
                       }
#ifdef _RMS_           
#pragma unroll
                       for (int n=0; n < FlowNode2D<FP,NUM_COMPONENTS>::NumEq; n++) {
                           //tmp_RMS_pack.DD_max[n]    = max(tmp_RMS_pack.DD_max[n],host_RMS_Array.GetElement(ii)->iDD_max[n]/int2float_RMS_scale);
                           //tmp_RMS_pack.sum_RMS[n]  += host_RMS_Array.GetElement(ii)->sum_RMS[n];
                           tmp_RMS_pack.RMS[n]      += host_RMS_Array.GetElement(ii)->RMS[n];
                           tmp_RMS_pack.sum_iRMS[n] += host_RMS_Array.GetElement(ii)->sum_iRMS[n];
                           tmp_RMS_pack.sumDiv[n]   += host_RMS_Array.GetElement(ii)->sumDiv[n];
                       }
#endif //_RMS_
                   }
#ifdef _BARRIER_                  
                   if(n_s > 1) {
                      CUDA_BARRIER((char*)"Halo/dt exchange");
                  }
#endif // _BARRIER_
             if(!isAdiabaticWall) {

                 iX0 = 0;
                 
                 if(isCalibrate) {
                     gettimeofday(&mark2,NULL);
                 }
 #pragma unroll
                   for(int ii=0;ii<n_s;ii++) {  // CUDA version

                       int i_gpu;
                       if(isSingleGPU) {
                           i_gpu = ActiveSingleGPU;
                       } else {
                           i_gpu = ii;
                       }

                       if(cudaSetDevice(i_gpu) != cudaSuccess ) {
                          *f_stream << "\nError set CUDA device no: "<< i_gpu << endl;
                          Exit_OpenHyperFLOW2D(n_s);
                       }
                       int max_X = cudaDimArray->GetElement(ii).GetX();
                       int max_Y = cudaDimArray->GetElement(ii).GetY();

                       cudaJ = cudaSubDomainArray->GetElement(ii);
                       
                       if (ThreadBlockSize == 0) {
                           num_cuda_threads = min(max_num_threads,dprop.multiProcessorCount*warp_size/max(1,current_div_heat));
                       } else {
                           num_cuda_threads = min(max_num_threads,ThreadBlockSize);
                       }
                       
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

                       cuda_CalcHeatOnWallSources<<<num_cuda_blocks,num_cuda_threads, 0, cuda_streams[ii]>>>(cudaJ,
                                                                                                             max_X*max_Y,
                                                                                                             max_X,max_Y,
                                                                                                             iX0 + max_X - r_Overlap,
                                                                                                             iX0 + l_Overlap,
                                                                                                             dx,dy,dt);
                       iX0 += max_X;
                   }
#ifdef _BARRIER_               
                   CUDA_BARRIER((char*)"Calc heat sources on wall");
#endif // _BARRIER_
                   if(isCalibrate) {
                       gettimeofday(&mark1,NULL);
                       d_time_heat = (FP)(mark1.tv_sec-mark2.tv_sec)+(FP)(mark1.tv_usec-mark2.tv_usec)*1.e-6; 
                   }
             }
        
        CurrentTimePart += dt;
        
        dt = dtmin;
         

       if ( isVerboseOutput && iter/NOutStep*NOutStep == iter ) {

             gettimeofday(&end,NULL);
             
             d_time = (FP)(end.tv_sec-start1.tv_sec)+(FP)(end.tv_usec-start1.tv_usec)*1.e-6; 
             
             if (ThreadBlockSize == 0 && isCalibrate  ) {

                    current_div_stage1 = max(1,CalibrateThreadBlockSize(current_div_stage1,
                                                                        opt_thread_block_size_stage1,
                                                                        opt_round_trip_stage1,
                                                                        d_time_stage1));
                    *f_stream << "Calibrate DEEPS_Stage1: thread_block_size=" << min(max_num_threads,dprop.multiProcessorCount*warp_size/max(1,current_div_stage1))  
                              << " time=" << d_time_stage1 <<   endl;
                    current_div_stage2 = max(1,CalibrateThreadBlockSize(current_div_stage2,
                                                                        opt_thread_block_size_stage2,
                                                                        opt_round_trip_stage2,
                                                                        d_time_stage2));
                    *f_stream << "Calibrate DEEPS_Stage2: thread_block_size=" << min(max_num_threads,dprop.multiProcessorCount*warp_size/max(1,current_div_stage2))  
                              << " time=" << d_time_stage2 <<   endl;
                    if(!isAdiabaticWall) {
                        current_div_heat = max(1,CalibrateThreadBlockSize(current_div_heat,
                                                                          opt_thread_block_size_heat,
                                                                          opt_round_trip_heat,
                                                                          d_time_heat));
                        *f_stream << "Calibrate HeatOnWall: thread_block_size=" << min(max_num_threads,dprop.multiProcessorCount*warp_size/max(1,current_div_heat))  
                                  << " time=" << d_time_heat <<   endl;
                    }                 

             }

             if(isCalibrate && ThreadBlockSize == 0 && iter+last_iter > (NOutStep+1)*8 ) {
                current_div_stage1 = opt_thread_block_size_stage1[0];
                current_div_stage2 = opt_thread_block_size_stage2[0];
                current_div_heat   = opt_thread_block_size_heat[0];
                isCalibrate = 0;
            }

             if(d_time > 0.)
                 VCOMP = (FP)(1.0/d_time);
             else
                 VCOMP = 0.;

             if(!isAdiabaticWall) {
                 snprintf(Eff_warp_size,64,"(%u/%u/%u)",
                          min(max_num_threads,dprop.multiProcessorCount*warp_size/max(1,current_div_stage1)),
                          min(max_num_threads,dprop.multiProcessorCount*warp_size/max(1,current_div_stage2)),
                          min(max_num_threads,dprop.multiProcessorCount*warp_size/max(1,current_div_heat)));
             } else {
                 snprintf(Eff_warp_size,64,"(%u/%u)",
                          min(max_num_threads,dprop.multiProcessorCount*warp_size/max(1,current_div_stage1)),
                          min(max_num_threads,dprop.multiProcessorCount*warp_size/max(1,current_div_stage2)));
             }
             if (!isLocalTimeStep) {
                 snprintf(Time_step,64,"dt=%g ",dt);
             } else {
                 snprintf(Time_step,64," ");
             }


#ifdef _RMS_
             for (int k=0;k<(int)(FlowNode2D<FP,NUM_COMPONENTS>::NumEq);k++ ) {
                 
                 if (isAlternateRMS) {
                     if(tmp_RMS_pack.RMS[k] > 0.0 && tmp_RMS_pack.sumDiv[k] > 0) { 
                        tmp_RMS_pack.RMS[k] = sqrt(tmp_RMS_pack.RMS[k]/tmp_RMS_pack.sumDiv[k]);
                     }
                 } else {
                     if(tmp_RMS_pack.sum_iRMS[k] > 0) {
                        tmp_RMS_pack.RMS[k] = sqrt(tmp_RMS_pack.RMS[k]/tmp_RMS_pack.sum_iRMS[k]);
                     }
                 }

                 if(MonitorIndex == 0 || MonitorIndex > 4) {
                    max_RMS = max(tmp_RMS_pack.RMS[k],max_RMS);
                    if(max_RMS == tmp_RMS_pack.RMS[k])
                       k_max_RMS = k;
                 } else {
                    max_RMS = max(tmp_RMS_pack.RMS[MonitorIndex-1],max_RMS);
                    if(max_RMS == tmp_RMS_pack.RMS[MonitorIndex-1])
                       k_max_RMS = k;
                 }
             }
             
             SaveRMS(pRMS_OutFile,last_iter+iter, tmp_RMS_pack.RMS);

             if(k_max_RMS == i2d_nu_t)
                k_max_RMS +=turb_mod_name_index;

             if(k_max_RMS != -1 && (MonitorIndex == 0 || MonitorIndex == 5))
             *f_stream << "Step No " << iter+last_iter <<  "/" << Nstep << " maxRMS["<< RMS_Name[k_max_RMS] << "]="<< (FP)(max_RMS*100.) \
                        <<  " % step_time=" << (FP)(NOutStep*d_time) << " sec (" << (FP)VCOMP <<" step/sec) "<< Time_step << Eff_warp_size <<endl;
             else if(MonitorIndex > 0 &&  MonitorIndex < 5 )
                 *f_stream << "Step No " << iter+last_iter <<  "/" << Nstep <<  " maxRMS["<< RMS_Name[MonitorIndex-1] << "]="<< (FP)(max_RMS*100.) \
                  <<  " % step_time=" << (FP)(NOutStep*d_time) << " sec (" << (FP)VCOMP <<" step/sec) " << Time_step << Eff_warp_size <<endl;
             else
             *f_stream << "Step No " << iter+last_iter <<  "/" << Nstep << " maxRMS["<< k_max_RMS << "]="<< (FP)(max_RMS*100.) \
                        <<  " % step_time=" << (FP)(NOutStep*d_time)<< " sec (" << (FP)VCOMP <<" step/sec) "<< Time_step << Eff_warp_size <<endl;
#else                     
             *f_stream << "Step No " << iter+last_iter <<  "/" << Nstep <<" step_time=" << (FP)(NOutStep*d_time) << " sec (" << (FP)VCOMP <<" step/sec) "<< Time_step << Eff_warp_size <<endl;
#endif // _RMS_
             if(MonitorPointsArray && MonitorPointsArray->GetNumElements() > 0) {
                SaveMonitors(pMonitors_OutFile,GlobalTime+CurrentTimePart,MonitorPointsArray);
             }
        }
   iter++;

  } while((int)iter < Nstep);

for (unsigned int i=0;i<GlobalSubDomain->GetNumElements();i++) {

      int SubStartIndex = GlobalSubDomain->GetElementPtr(i)->GetX();
      int SubMaxX = GlobalSubDomain->GetElementPtr(i)->GetY();

      int i_gpu;
      if(isSingleGPU) {
          i_gpu = ActiveSingleGPU;
      } else {
          i_gpu = i;
      }

      if(cudaSetDevice(i_gpu) != cudaSuccess ) {
         *f_stream << "\nError set CUDA device no: "<< i_gpu << endl;
         Exit_OpenHyperFLOW2D(n_s);
      }

      if(i == GlobalSubDomain->GetNumElements()-1)
        r_Overlap = 0;
      else
        r_Overlap = 1;
      if(i == 0)
        l_Overlap = 0;
      else
        l_Overlap = 1;

       TmpMaxX = (SubMaxX-SubStartIndex) - r_Overlap;
       TmpMatrixPtr = (FlowNode2D<FP,NUM_COMPONENTS>*)((ulong)J->GetMatrixPtr()+(ulong)(sizeof(FlowNode2D<FP,NUM_COMPONENTS>)*(SubStartIndex + r_Overlap)*MaxY)); // !!!

#ifdef _PARALLEL_RECALC_Y_PLUS_
       if(ProblemType == SM_NS &&
          (TurbExtModel == TEM_Spalart_Allmaras ||
           TurbExtModel == TEM_vanDriest ||
           TurbExtModel == TEM_k_eps_Chien)) {

           FP  x0 = SubStartIndex*FlowNode2D<FP,NUM_COMPONENTS>::dx;
           
           if(isSingleGPU) {
              *f_stream << "Parallel recalc y+ on CUDA device No  " << ActiveSingleGPU  << endl;
           } else {
              *f_stream << "Parallel recalc y+ on CUDA device No  " << i  << endl;
           }
           
           int max_X = cudaDimArray->GetElement(i).GetX();
           int max_Y = cudaDimArray->GetElement(i).GetY();
           
           if (ThreadBlockSize == 0) {
               num_cuda_threads = min(max_num_threads,dprop.multiProcessorCount*warp_size/max(1,current_div_heat));
           } else {
               num_cuda_threads = min(max_num_threads,ThreadBlockSize);
           }

           num_cuda_blocks  = (max_X*max_Y)/num_cuda_threads;

           if(num_cuda_blocks*num_cuda_threads != max_X*max_Y)
              num_cuda_blocks++;

           cuda_Recalc_y_plus<<<num_cuda_blocks,num_cuda_threads, 0, cuda_streams[i]>>>(cudaSubDomainArray->GetElement(i),
                                                                                        max_X*max_Y,
                                                                                        cudaWallNodesArray->GetElement(i),
                                                                                        NumWallNodes,
                                                                                        min(dx,dy),
                                                                                        max((x0+FlowNode2D<FP,NUM_COMPONENTS>::dx*max_X),
                                                                                               (FlowNode2D<FP,NUM_COMPONENTS>::dy*max_Y)),
                                                                                        FlowNode2D<FP,NUM_COMPONENTS>::dx,
                                                                                        FlowNode2D<FP,NUM_COMPONENTS>::dy,
                                                                                        max_Y);
           CUDA_BARRIER((char*)"y+ recalc");
       }

#endif // _PARALLEL_RECALC_Y_PLUS_

       CopyDeviceToHost(cudaArraySubDomain->GetElement(i),TmpMatrixPtr,(sizeof(FlowNode2D<FP,NUM_COMPONENTS>))*(TmpMaxX*MaxY),cuda_streams[i]);
  }
       CUDA_BARRIER((char*)"Data collection");
       
      if ( isGasSource && SrcList) {
           *f_stream << "\nSet gas sources...";
            SrcList->SetSources2D();
            *f_stream << "OK" << endl;
            for (unsigned int i=0;i<GlobalSubDomain->GetNumElements();i++) {

                 int SubStartIndex = GlobalSubDomain->GetElementPtr(i)->GetX();
                 int SubMaxX = GlobalSubDomain->GetElementPtr(i)->GetY();

                 int i_gpu;
                 if(isSingleGPU) {
                     i_gpu = ActiveSingleGPU;
                 } else {
                     i_gpu = i;
                 }

                 if(cudaSetDevice(i_gpu) != cudaSuccess ) {
                    *f_stream << "\nError set CUDA device no: "<< i_gpu << endl;
                    Exit_OpenHyperFLOW2D(n_s);
                 }

                 if(i == GlobalSubDomain->GetNumElements()-1)
                    r_Overlap = 0;
                 else
                    r_Overlap = 1;
                 if(i == 0)
                    l_Overlap = 0;
                 else
                    l_Overlap = 1;

                  TmpMaxX = (SubMaxX-SubStartIndex) - r_Overlap;
                  TmpMatrixPtr = (FlowNode2D<FP,NUM_COMPONENTS>*)((ulong)J->GetMatrixPtr()+(ulong)(sizeof(FlowNode2D<FP,NUM_COMPONENTS>)*(SubStartIndex)*MaxY));
                  CopyHostToDevice(TmpMatrixPtr,cudaArraySubDomain->GetElement(i),(sizeof(FlowNode2D<FP,NUM_COMPONENTS>))*(TmpMaxX*MaxY),cuda_streams[i]);
            }
        CUDA_BARRIER((char*)"Update sources in to device");
      }

      
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

         if(is_Cx_calc) { // For Airfoils only
          *f_stream << "\nCx = " << Calc_Cx_2D(J,x0_body,y0_body,dx_body,dy_body,Flow2DList->GetElement(Cx_Flow_index-1)) << 
                       " Cy = "  << Calc_Cy_2D(J,x0_body,y0_body,dx_body,dy_body,Flow2DList->GetElement(Cx_Flow_index-1)) << 
                       " Fx = "  << CalcXForce2D(J,x0_body,y0_body,dx_body,dy_body) << " Fy = " << CalcYForce2D(J,x0_body,y0_body,dx_body,dy_body) << endl;
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

if(MonitorIndex < 5) {
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

#ifndef _P2P_ACCESS_
for(int ii=0;ii<(int)cudaArraySubDomain->GetNumElements();ii++)
    free(host_HalloBuff[ii]);
#endif // _P2P_ACCESS_

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




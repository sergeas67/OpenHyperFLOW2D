/*******************************************************************************
*   OpenHyperFLOW2D-CUDA                                                       *
*                                                                              *
*   Version  1.0.2                                                             *
*   Copyright (C)  1995-2015 by Serge A. Suchkov                               *
*   Copyright policy: LGPL V3                                                  *
*   http://openhyperflow2d.googlecode.com                                      *
*                                                                              *
*   hf2d_start.cpp: OpenHyperFLOW2D solver init code....                       *
*                                                                              *
*  last update: 01/02/2015                                                     *
********************************************************************************/
#ifdef _CUDA_
#define _PARALLEL_ONLY

#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>

#include <sys/time.h>
#include <sys/timeb.h>
#include "libDEEPS2D/deeps2d_core.hpp"
 
SolverMode ProblemType;
int rank;
int last_rank;
int warp_size = 0;
FP _dt_test=1.0;

timeval mark1, mark2;

// Arrays for multiGPU
UArray< FlowNode2D<FP,NUM_COMPONENTS>* >*     cudaArraySubDomain      = NULL;
UArray< FlowNodeCore2D<FP,NUM_COMPONENTS>* >* cudaArrayCoreSubDomain  = NULL;
UArray< XY<int>  >*                           cudaDimArray            = NULL;
UArray< XY<int>* >*                           cudaWallNodesArray      = NULL;
UArray< ChemicalReactionsModelData2D* >*      cudaCRM2DArray          = NULL;
XY<int>*                                      cudaWallNodes           = NULL;
FlowNode2D<FP,NUM_COMPONENTS>*                cudaSubDomain           = NULL;
FlowNodeCore2D<FP,NUM_COMPONENTS>*            cudaCoreSubDomain       = NULL;
ChemicalReactionsModelData2D*                 cudaCRM2D               = NULL;

UArray< MonitorPoint >*                       MonitorPointsArray      = NULL;

FP*                                           cudaHu;
UArray<FP*>*                                  cudaHuArray;
FP                                            x0;
FP                                            dx;
FP                                            dy;

int num_gpus = 0;   // number of CUDA GPUs
int max_num_threads = 0;
size_t max_gpu_memsize = 0;

size_t task_size;
cudaError_t cudaState;


UArray< FP >*                                 WallNodesUw_2D = NULL;
int                                           NumWallNodes;
int                                           isSingleGPU = 1;
                                                 
/*                                               
FP GPUConf {                                 

}

struct NodeConfig {
 char*  HostName;
 int    numCPUcores;
 UArray<FP> GPUrate;
};

UArray<NodeConfig*>* clusterArray = NULL;

int CalibrateNodes(UArray<NodeConfig*>* CA, int rank) {

NodeConfig TmpNConf;  
    
 if (rank == 0) {
  
 } else {
 
 }

}
*/

int main( int argc, char **argv )
{
    const  FP                      ver=_VER;
    char                           inFile[256];
    ChemicalReactionsModelData2D   TmpCRM2D;

    FlowNode2D<FP,NUM_COMPONENTS>* TmpMatrixPtr;
    //FlowNodeCore2D<FP,NUM_COMPONENTS>* TmpCoreMatrixPtr;

    int TmpMaxX;
    
    ostream*     o_stream = &cout;
#ifdef    _DEBUG_0
    __ExceptLib.SetSystemException(SIGINT, AT_HANDLED);
    __ExceptLib.SetSystemException(SIGFPE, AT_HANDLED);
    __ExceptLib.SetSystemException(SIGSEGV,AT_HANDLED);
    __ExceptLib.SetSystemException(SIGBUS, AT_HANDLED);
#endif //_DEBUG_0

#ifdef _DEBUG_0
    ___try {
#endif  // _DEBUG_0
        if (argc < 2) {
            printf("OpenHyperFLOW2D/DEEPS/FP%u solver v %'.2f ",(unsigned int)(sizeof(FP)*8),ver);
            printf(" (parallel CUDA version)");
            printf("\nCopyright (C) 1995-2015 by Serge A. Suchkov\nCopyright policy: LGPL V3\nUsage: %s [{input_data_file}]\n",argv[0]);
            printf("\n\t* Density-based 2D-Navier-Stokes solver for ");
#ifdef _UNIFORM_MESH_
            printf("uniform cartesian mesh");
#else
            printf("non-uniform mesh");
#endif //_UNIFORM_MESH_

            printf("\n\n\tFlowNode2D size = %d bytes\n\n",sizeof(FlowNode2D<FP,NUM_COMPONENTS>));
            Exit_OpenHyperFLOW2D();
        } else {
            sprintf(inFile,"%s",argv[1]);
            Data = new InputData(inFile,DS_FILE,o_stream,0,10);
            if (Data->GetDataError()!=0) {
                *o_stream << "\nInput data error.\n" ;
                 o_stream->flush();
                 Exit_OpenHyperFLOW2D(0);
            }
        }
//------------------------- CUDA version ------------------------------------

        printf("\n\n\t");
        cudaGetDeviceCount(&num_gpus);

        if (num_gpus < 1) {
         printf("no CUDA capable devices were detected\n");
         Exit_OpenHyperFLOW2D(0);
        }

        printf("Number of CUDA devices:\t%d\n", num_gpus);
        //int num_dev;
        for (int i = 0; i < num_gpus; i++)
        {
         cudaDeviceProp dprop;
         cudaGetDeviceProperties(&dprop, i);
         printf("\t   %d: %s\n", i, dprop.name);
         max_num_threads  = dprop.maxThreadsPerBlock;
         printf("\t   Max threads per block: %d:\n", max_num_threads);
         max_gpu_memsize = dprop.totalGlobalMem;
         printf("\t   Max GPU memory size: %lu\n", max_gpu_memsize);
         printf("\t   Number of  multiprocessors: %d\n",dprop.multiProcessorCount);
         printf("\t   Is kernels concurrent: %d\n",dprop.concurrentKernels);
         warp_size = dprop.warpSize;
         printf("\t   Warp size: %d\n",warp_size);
         printf("\t   Enable timeout: %i\n\n",dprop.kernelExecTimeoutEnabled);
        }
        
        // Init shared data
        InitSharedData(Data,&chemical_reactions);        
        
        if(isSingleGPU)
           num_threads = num_blocks = num_gpus = 1;
        else
           num_threads = num_blocks = num_gpus;

        if(num_gpus > 1) {
           printf("Using multi GPU mode.\n\n");
        } else {
           printf("Using single GPU mode.\n\n");
        }
       
       //Create arrays  
        cudaArraySubDomain      =  new UArray< FlowNode2D<FP,NUM_COMPONENTS>* >();
        cudaArrayCoreSubDomain  =  new UArray< FlowNodeCore2D<FP,NUM_COMPONENTS>* >();
        cudaDimArray            =  new UArray< XY<int> >();
        cudaHuArray             =  new UArray< FP* >();
        dt_min_host_Array       =  new UArray<unsigned int*>();   
        dt_min_device_Array     =  new UArray<unsigned int*>();
        cudaWallNodesArray      =  new UArray< XY<int>* >();
        cudaCRM2DArray          =  new UArray< ChemicalReactionsModelData2D* >();          

        // Init solver (run on host)
        InitDEEPS2D((void*)o_stream);

        if(ProblemType == SM_NS) {                                          
            *o_stream << "\nSolver Mode: Navier-Stokes/FP" << 8*sizeof(FP) <<"\n" << endl;
        } else {
            *o_stream << "\nSolver Mode: Euler/FP" << 8*sizeof(FP) <<"\n" << endl;
        }

        cudaStream_t *cuda_streams = (cudaStream_t *) malloc(num_gpus * sizeof(cudaStream_t));
        cudaEvent_t  *cuda_events = (cudaEvent_t *) malloc(num_gpus * sizeof(cudaEvent_t));

        for (int i = 0; i < num_gpus; i++) {

           if(cudaSetDevice(i) != cudaSuccess ) {
              *o_stream << "\nError set CUDA device no: "<< i << endl;
              Exit_OpenHyperFLOW2D(num_gpus);
           }

           cudaState = cudaStreamCreate(&(cuda_streams[i]));
           if(cudaState != cudaSuccess ) {
              *o_stream << "\nError create stream no: "<< i << endl;
              printf("%s\n", cudaGetErrorString( cudaGetLastError() ) );
              Exit_OpenHyperFLOW2D(num_gpus);
           }

           cudaState = cudaEventCreate(&(cuda_events[i]));
           if(cudaState != cudaSuccess ) {
              *o_stream << "\nError create event no: "<< i << endl;
              printf("%s\n", cudaGetErrorString( cudaGetLastError() ) );
              Exit_OpenHyperFLOW2D(num_gpus);
           }
           
           // Load components properties
           LoadTable2GPU(chemical_reactions.Cp_OX,   TmpCRM2D.Cp_OX,   i);
           LoadTable2GPU(chemical_reactions.Cp_Fuel, TmpCRM2D.Cp_Fuel, i);
           LoadTable2GPU(chemical_reactions.Cp_cp,   TmpCRM2D.Cp_cp,   i);
           LoadTable2GPU(chemical_reactions.Cp_air,  TmpCRM2D.Cp_air,  i);
                                                     
           LoadTable2GPU(chemical_reactions.mu_air,  TmpCRM2D.mu_air,  i);
           LoadTable2GPU(chemical_reactions.mu_cp,   TmpCRM2D.mu_cp,   i);
           LoadTable2GPU(chemical_reactions.mu_Fuel, TmpCRM2D.mu_Fuel, i);
           LoadTable2GPU(chemical_reactions.mu_OX,   TmpCRM2D.mu_OX,   i);
           
           if(ProblemType == SM_NS) {                                          
               LoadTable2GPU(chemical_reactions.lam_air, TmpCRM2D.lam_air, i);
               LoadTable2GPU(chemical_reactions.lam_cp,  TmpCRM2D.lam_cp,  i);
               LoadTable2GPU(chemical_reactions.lam_Fuel,TmpCRM2D.lam_Fuel,i);
               LoadTable2GPU(chemical_reactions.lam_OX,  TmpCRM2D.lam_OX,  i);
           }

           TmpCRM2D.H_air  = chemical_reactions.H_air;
           TmpCRM2D.H_cp   = chemical_reactions.H_cp;
           TmpCRM2D.H_Fuel = chemical_reactions.H_Fuel;
           TmpCRM2D.H_OX   = chemical_reactions.H_OX;
           
           TmpCRM2D.R_air  = chemical_reactions.R_air;
           TmpCRM2D.R_cp   = chemical_reactions.R_cp;
           TmpCRM2D.R_Fuel = chemical_reactions.R_Fuel;
           TmpCRM2D.R_OX   = chemical_reactions.R_OX;

           TmpCRM2D.K0     = chemical_reactions.K0;
           TmpCRM2D.Tf     = chemical_reactions.Tf;
           TmpCRM2D.gamma  = chemical_reactions.gamma;
           
           if(cudaMalloc( (void**)&cudaCRM2D, sizeof(ChemicalReactionsModelData2D) ) == cudaErrorMemoryAllocation) {
              *o_stream << "\nError allocate GPU memory for CRM2D on device no:" << i << endl;
              Exit_OpenHyperFLOW2D(num_gpus);
           }
           
           CopyHostToDevice(&TmpCRM2D,cudaCRM2D,sizeof(ChemicalReactionsModelData2D));
           
           cudaCRM2DArray->AddElement(&cudaCRM2D);

       }

        for (int i_dev=0; i_dev < num_gpus;i_dev++) {

#ifdef _P2P_ACCESS_
            if((num_gpus>1)&&(i_dev+1 < num_gpus)) {
                SetP2PAccess(i_dev,i_dev+1);
            }
#endif // _P2P_ACCESS_

            if(cudaSetDevice(i_dev) != cudaSuccess ) {
               *o_stream << "\nError set CUDA device no: "<< i_dev << endl;
               Exit_OpenHyperFLOW2D(num_gpus);
            }

            cudaState = cudaMalloc( (void**)&cudaHu, sizeof(FP)*(NUM_COMPONENTS+1) );
            if(cudaState == cudaErrorMemoryAllocation) {
               *o_stream << "\nError allocate GPU memory for Hu[]" << endl;
               Exit_OpenHyperFLOW2D(num_gpus);
            }

            CopyHostToDevice(FlowNode2D<FP,NUM_COMPONENTS>::Hu,cudaHu,sizeof(FP)*(NUM_COMPONENTS+1));

            cudaHuArray->AddElement(&cudaHu);

#ifdef _DEVICE_MMAP_
            cudaState = cudaHostAlloc( (void**)&dt_min_host, sizeof(unsigned int), cudaHostAllocMapped ); 

            dt_min_host_Array->AddElement(&dt_min_host);  

            if(cudaState == cudaErrorMemoryAllocation) {
               *o_stream << "\nError allocate GPU memory for dt_min" << endl;
               Exit_OpenHyperFLOW2D(num_gpus);
            }

            cudaState = cudaHostGetDevicePointer( &dt_min_device, dt_min_host, 0 );

            dt_min_device_Array->AddElement(&dt_min_device);

            if(cudaState == cudaErrorMemoryAllocation) {
               *o_stream << "\nError mapped GPU memory for dt_min" << endl;
               Exit_OpenHyperFLOW2D(num_gpus);
            }
#else
            cudaState = cudaMalloc( (void**)&dt_min_device, sizeof(unsigned int));

            dt_min_device_Array->AddElement(&dt_min_device);

            if(cudaState == cudaErrorMemoryAllocation) {
               *o_stream << "\nError allocate GPU memory for dt_min "<<  endl;
               Exit_OpenHyperFLOW2D(num_gpus);
            }
#endif // _DEVICE_MMAP_
        }
        
        if(ProblemType == SM_NS) {                                          

            // Scan area for seek wall nodes      
            WallNodes = GetWallNodes((ofstream*)o_stream,J,Data->GetIntVal((char*)"isVerboseOutput")); 

            NumWallNodes = WallNodes->GetNumElements();
            *o_stream << NumWallNodes << " wall nodes found" << endl; 
        }
        
        gettimeofday(&mark2,NULL);

        for (int i_dev=0; i_dev < num_gpus;i_dev++) {

            if(cudaSetDevice(i_dev) != cudaSuccess ) {
               *o_stream << "\nError set CUDA device no: "<< i_dev << endl;
               Exit_OpenHyperFLOW2D(num_gpus);
            }

            if(ProblemType == SM_NS && NumWallNodes > 0) {
                
                cudaState = cudaMalloc( (void**)&cudaWallNodes, sizeof(XY<int>)*NumWallNodes );
                if(cudaState == cudaErrorMemoryAllocation) {
                   *o_stream << "\nError allocate GPU memory for WallNodes array." << endl;
                   Exit_OpenHyperFLOW2D(num_gpus);
                }

                if(cudaState != cudaSuccess ) {
                   *o_stream << "\nError set CUDA device no: "<< i_dev << endl;
                   Exit_OpenHyperFLOW2D(num_gpus);
                }

                CopyHostToDevice(WallNodes->GetArrayPtr(),cudaWallNodes,sizeof(XY<int>)*NumWallNodes);

                cudaWallNodesArray->AddElement(&cudaWallNodes);
            }
        }

        TmpMatrixPtr=J->GetMatrixPtr();
        int SubStartIndex, SubMaxX,r_Overlap=0,l_Overlap=0;
        SubStartIndex = 0;
        int iX0=0;

        *o_stream << "Allocate SubDomain:\n";

        //Allocate GPU buffers

        for (unsigned int i=0;i<GlobalSubDomain->GetNumElements();i++) {

            if(cudaSetDevice(i) != cudaSuccess ) {
               *o_stream << "\nError set CUDA device no: "<< i << endl;
               Exit_OpenHyperFLOW2D(num_gpus);
            }

            SubStartIndex = GlobalSubDomain->GetElementPtr(i)->GetX();  
            SubMaxX = GlobalSubDomain->GetElementPtr(i)->GetY();

            if(i == GlobalSubDomain->GetNumElements()-1)
              r_Overlap = 0;
            else
              r_Overlap = 1;
            if(i == 0)
              l_Overlap = 0;
            else
              l_Overlap = 1;

            TmpMaxX = (SubMaxX-SubStartIndex) + r_Overlap;

            // Allocate FlowNode2D<FP,NUM_COMPONENTS> subdomain
            cudaState = cudaMalloc( (void**)&cudaSubDomain, (sizeof(FlowNode2D<FP,NUM_COMPONENTS>))*(TmpMaxX*MaxY) );

            if(cudaState == cudaErrorMemoryAllocation) {
               *o_stream << "\nError allocate GPU memory for SubDomain"<< endl;
               Exit_OpenHyperFLOW2D(num_gpus);
            }

           cudaArraySubDomain->AddElement(&cudaSubDomain);

           // Allocate FlowNodeCore2D<FP,NUM_COMPONENTS> subdomain
           cudaState = cudaMalloc( (void**)&cudaCoreSubDomain, (sizeof(FlowNodeCore2D<FP,NUM_COMPONENTS>))*(TmpMaxX*MaxY) );

           if(cudaState == cudaErrorMemoryAllocation) {
              *o_stream << "\nError allocate GPU memory for CoreSubDomain"<< endl;
              Exit_OpenHyperFLOW2D(num_gpus);
           }
           cudaArrayCoreSubDomain->AddElement(&cudaCoreSubDomain);
        }

        for (unsigned int i=0;i<GlobalSubDomain->GetNumElements();i++) {

            XY<int> TmpDim;

            if(cudaSetDevice(i) != cudaSuccess ) {
               *o_stream << "\nError set CUDA device no: "<< i << endl;
               Exit_OpenHyperFLOW2D(num_gpus);
            }


            if(i == GlobalSubDomain->GetNumElements()-1)
              r_Overlap = 0;
            else
              r_Overlap = 1;
            if(i == 0)
              l_Overlap = 0;
            else
              l_Overlap = 1;

            SubStartIndex = GlobalSubDomain->GetElementPtr(i)->GetX();  
            SubMaxX = GlobalSubDomain->GetElementPtr(i)->GetY();

            TmpMaxX = (SubMaxX-SubStartIndex) - r_Overlap;
            TmpMatrixPtr = (FlowNode2D<FP,NUM_COMPONENTS>*)((ulong)J->GetMatrixPtr()+(ulong)(sizeof(FlowNode2D<FP,NUM_COMPONENTS>)*(SubStartIndex)*MaxY));

            int num_cuda_threads = warp_size;
            int num_cuda_blocks  = (TmpMaxX*MaxY)/num_cuda_threads;

            if (num_cuda_blocks*num_cuda_threads != TmpMaxX*MaxY)
                num_cuda_blocks++;

            x0 = SubStartIndex*FlowNode2D<FP,NUM_COMPONENTS>::dx;

            *o_stream << "SubDomain("<<i<<")[" << TmpMaxX << "x" << MaxY << "]  Size=" << (ulong)(sizeof(FlowNode2D<FP,NUM_COMPONENTS>)*TmpMaxX*MaxY)/(1024*1024) << " Mb\n"; 

            cudaSubDomain = cudaArraySubDomain->GetElement(i);

            CopyHostToDevice(TmpMatrixPtr,cudaSubDomain,(sizeof(FlowNode2D<FP,NUM_COMPONENTS>))*(TmpMaxX*MaxY));

            cudaHu =  cudaHuArray->GetElement(i);

            if(NumWallNodes > 0 &&
               ProblemType == SM_NS ) {

                cudaWallNodes = cudaWallNodesArray->GetElement(i);

                *o_stream << "GPU no: " << i << endl; 
                *o_stream << "CUDA threads: " << num_cuda_threads << endl;
                *o_stream << "CUDA thread blocks : " << num_cuda_blocks << endl;

                *o_stream << "\nParallel calc min distance to wall..." << endl;
                *o_stream << "Run cuda_SetMinDistanceToWall2D kernel..." << flush;

                cuda_SetMinDistanceToWall2D<<<num_cuda_blocks,num_cuda_threads, 0, cuda_streams[i]>>>(cudaSubDomain,
                                                                                                      TmpMaxX*MaxY,
                                                                                                      cudaWallNodes,
                                                                                                      NumWallNodes,
                                                                                                      min(dx,dy),
                                                                                                      max((x0+FlowNode2D<FP,NUM_COMPONENTS>::dx*TmpMaxX), 
                                                                                                      (FlowNode2D<FP,NUM_COMPONENTS>::dy*MaxY)),
                                                                                                      FlowNode2D<FP,NUM_COMPONENTS>::dx,
                                                                                                      FlowNode2D<FP,NUM_COMPONENTS>::dy,
                                                                                                      x0);

                 CUDA_BARRIER((char*)"cuda_SetMinDistanceToWall2D");
                 *o_stream << "OK" << endl;


                 if( TurbExtModel == TEM_Spalart_Allmaras ||
                     TurbExtModel == TEM_vanDriest ||
                     TurbExtModel == TEM_k_eps_Chien ) {

                      *o_stream << "Run cuda_Recalc_y_plus kernel on CUDA device No " << i << flush;


                      cuda_Recalc_y_plus<<<num_cuda_blocks,num_cuda_threads, 0, cuda_streams[i]>>>(cudaSubDomain,
                                                                                                   TmpMaxX*MaxY,
                                                                                                   cudaWallNodes,
                                                                                                   NumWallNodes,
                                                                                                   min(dx,dy),
                                                                                                   max((x0+FlowNode2D<FP,NUM_COMPONENTS>::dx*TmpMaxX),
                                                                                                       (FlowNode2D<FP,NUM_COMPONENTS>::dy*MaxY)),
                                                                                                   FlowNode2D<FP,NUM_COMPONENTS>::dx,
                                                                                                   FlowNode2D<FP,NUM_COMPONENTS>::dy,
                                                                                                   MaxY);

                      CUDA_BARRIER((char*)"cuda_Recalc_y_plus");
                      *o_stream << "...OK" << endl;
                     }
            }

            *o_stream << "Run cuda_SetInitBoundaryLayer kernel..." << flush;


            cuda_SetInitBoundaryLayer<<<num_cuda_blocks,num_cuda_threads, 0, cuda_streams[i]>>>(cudaSubDomain,
                                                                                                TmpMaxX*MaxY, iX0, MaxY,
                                                                                                delta_bl,
                                                                                                SigW,SigF,(TurbulenceExtendedModel)TurbExtModel, 
                                                                                                FlowNode2D<FP,NUM_COMPONENTS>::dx,
                                                                                                FlowNode2D<FP,NUM_COMPONENTS>::dy,
                                                                                                cudaHu,
                                                                                                FlowNode2D<FP,NUM_COMPONENTS>::isSrcAdd,
                                                                                                FlowNode2D<FP,NUM_COMPONENTS>::FT,
                                                                                                ProblemType);

            CUDA_BARRIER((char*)"cuda_SetInitBoundaryLayer");
            *o_stream << "OK" << endl;


            CopyDeviceToHost(cudaSubDomain,TmpMatrixPtr,(sizeof(FlowNode2D<FP,NUM_COMPONENTS>))*(TmpMaxX*MaxY));

            iX0 += TmpMaxX;
            TmpDim.SetXY(TmpMaxX,MaxY);
            cudaDimArray->AddElement(&TmpDim);
            o_stream->flush();
            if(MonitorPointsArray) {
                for(int ii_monitor=0;ii_monitor<(int)MonitorPointsArray->GetNumElements();ii_monitor++) {
                        if(MonitorPointsArray->GetElement(ii_monitor).MonitorXY.GetX() >= x0 &&
                           MonitorPointsArray->GetElement(ii_monitor).MonitorXY.GetX() < x0 + FlowNode2D<FP,NUM_COMPONENTS>::dx*TmpMaxX) {
                           MonitorPointsArray->GetElement(ii_monitor).rank = i; 
                        }
                }
            }
        }

        gettimeofday(&mark1,NULL);
        *o_stream << "OK\n" << "Time: " << (FP)(mark1.tv_sec-mark2.tv_sec)+(FP)(mark1.tv_usec-mark2.tv_usec)*1.e-6 << " sec." << endl; 



       DEEPS2D_Run((ofstream*)o_stream,        // ofstream* o_stream,
                   J,                          // UMatrix2D< FlowNode2D< FP,NUM_COMPONENTS > >*
                   C,                          // UMatrix2D< FlowNodeCore2D< FP,NUM_COMPONENTS > >*
                   cudaArraySubDomain,         // UArray< FlowNode2D< FP,NUM_COMPONENTS >* >*
                   cudaArrayCoreSubDomain,     // UArray< FlowNodeCore2D< FP,NUM_COMPONENTS >* >*
                   cudaDimArray,               // UArray< XY<int> >*
                   cudaWallNodesArray,         // UArray< XY<int>* >*
                   cudaCRM2DArray,
                   num_gpus,
                   cuda_streams,
                   cuda_events);          
       
       for (int i = 0; i < num_gpus; i++)  {

           if(cudaSetDevice(i) != cudaSuccess ) {
              *o_stream << "\nError set CUDA device no: "<< i << endl;
              Exit_OpenHyperFLOW2D(num_gpus);
           }

           cudaState = cudaStreamDestroy(cuda_streams[i]);
           if(cudaState != cudaSuccess ) {
              *o_stream << "\nError destroy stream no: "<< i << endl;
              printf("%s\n", cudaGetErrorString( cudaGetLastError() ) );
              Exit_OpenHyperFLOW2D(num_gpus);
           }

           cudaState = cudaEventDestroy(cuda_events[i]);
           if(cudaState != cudaSuccess ) {
              *o_stream << "\nError destroy event no: "<< i << endl;
              printf("%s\n", cudaGetErrorString( cudaGetLastError() ) );
              Exit_OpenHyperFLOW2D(num_gpus);
           }
       }

       free(cuda_streams);
       free(cuda_events);

       for (int i_dev=0; i_dev < num_gpus;i_dev++) {

           if(cudaSetDevice(i_dev) != cudaSuccess ) {
              *o_stream << "\nError set CUDA device no: "<< i_dev << endl;
              Exit_OpenHyperFLOW2D(num_gpus);
           }
#ifdef _P2P_ACCESS_
           if((num_gpus>1)&&(i_dev+1 < num_gpus)){
               DisableP2PAccess(i_dev,i_dev+1);
           }
#endif // _P2P_ACCESS_

           cudaSubDomain = cudaArraySubDomain->GetElement(i_dev);

           cudaState = cudaFree(cudaSubDomain);

           if(cudaState != cudaSuccess ) {
              *o_stream << "\nError free SubDomain  from GPU memory." << endl;
              Exit_OpenHyperFLOW2D(num_gpus);
           }

           cudaCoreSubDomain = cudaArrayCoreSubDomain->GetElement(i_dev);

           cudaState = cudaFree(cudaCoreSubDomain);

           if(cudaState != cudaSuccess ) {
              *o_stream << "\nError free CoreSubDomain from GPU memory." << endl;
              Exit_OpenHyperFLOW2D(num_gpus);
           }

           if(cudaWallNodesArray->GetNumElements()) {
              cudaWallNodes = cudaWallNodesArray->GetElement(i_dev);

              cudaState = cudaFree(cudaWallNodes);

              if(cudaState != cudaSuccess ) {
                 *o_stream << "\nError free WallNodes array from GPU memory." << endl;
                 Exit_OpenHyperFLOW2D(num_gpus);
              }
           }
#ifdef _DEVICE_MMAP_
           dt_min_host = dt_min_host_Array->GetElement(i_dev);

           cudaState = cudaFreeHost(dt_min_host);
           if(cudaState != cudaSuccess ) {
              *o_stream << "\nError free dt_min from GPU memory." << endl;
              Exit_OpenHyperFLOW2D(num_gpus);
           }
#else
           dt_min_device = dt_min_device_Array->GetElement(i_dev);

           cudaState = cudaFree(dt_min_device);

           if(cudaState != cudaSuccess ) {
              *o_stream << "\nError free dt_min from GPU memory." << endl;
              Exit_OpenHyperFLOW2D(num_gpus);
           }
#endif //_DEVICE_MMAP_

           cudaHu = cudaHuArray->GetElement(i_dev);

           cudaState = cudaFree(cudaHu);

           if(cudaState != cudaSuccess ) {
              *o_stream << "\nError free cudaHu from GPU memory." << endl;
              Exit_OpenHyperFLOW2D(num_gpus);
           }
       }
       //>>>>>>>>>>>>>>>>>>>>>>
       //DataSnapshot(OutFileName,WM_REWRITE);
       //>>>>>>>>>>>>>>>>>>>>>>

       delete cudaArraySubDomain;
       delete cudaArrayCoreSubDomain;
       delete cudaDimArray;
       delete cudaHuArray;
       delete dt_min_host_Array;
       delete dt_min_device_Array;
       delete cudaWallNodesArray;

#ifdef _DEBUG_0
    }__except(SysException e) {{
            const char ConstErrorMessage[]=" handled by <LibExcept> module in <InputData>.\n";
            *o_stream << "SIG";
            if (e == SIGSEGV)
                *o_stream << "SEGV" << ConstErrorMessage ;
            else if (e == SIGBUS)
                *o_stream << "BUS" << ConstErrorMessage ;
            else if (e == SIGFPE)
                *o_stream << "FPE" << ConstErrorMessage ;
            else if (e == SIGINT)
                *o_stream << "INT" << ConstErrorMessage ;
            else if (e == SIGABRT)
                *o_stream << "ABRT" << ConstErrorMessage ;
            else if (e == SIGIO)
                *o_stream << "IO" << ConstErrorMessage ;
            else if (e == SIGTRAP)
                *o_stream << "TRAP" << ConstErrorMessage ;
            else
                *o_stream << " No: "<< e << ConstErrorMessage ;
        }
    } __end_except;
#endif  // _DEBUG_0
        *o_stream << "Computation stopped.\n";
        o_stream->flush();
        Exit_OpenHyperFLOW2D(num_gpus);
}
#endif //_CUDA_


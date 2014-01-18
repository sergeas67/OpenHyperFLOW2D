/*******************************************************************************
*   OpenHyperFLOW2D-CUDA                                                       *
*                                                                              *
*   Version  1.0.1                                                             *
*   Copyright (C)  1995-2014 by Serge A. Suchkov                               *
*   Copyright policy: LGPL V3                                                  *
*                                                                              *
*  hf2d_start.cpp: OpenHyperFLOW2D solver init code....                        *
*                                                                              *
*  last update: 16/01/2014                                                     *
********************************************************************************/
#ifdef _CUDA_
#define _PARALLEL_ONLY

#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>


#include <sys/time.h>
#include <sys/timeb.h>
#include "libDEEPS2D/deeps2d_core.hpp"
 
int rank;
int last_rank;
int warp_size = 0;
float _dt_test=1.0;


timeval mark1, mark2;
UArray< FlowNode2D<double,NUM_COMPONENTS>* >*     cudaArraySubmatrix      = NULL;
UArray< FlowNodeCore2D<double,NUM_COMPONENTS>* >* cudaArrayCoreSubmatrix  = NULL;
UArray< XY<int> >*                                cudaDimArray            = NULL;
XY<int>*                                          cudaWallNodes           = NULL;
FlowNode2D<double,NUM_COMPONENTS>*                cudaSubmatrix           = NULL;
FlowNodeCore2D<double,NUM_COMPONENTS>*            cudaCoreSubmatrix       = NULL;

double*  cudaHu;
double   x0;
double   dx;
double   dy;

UArray< double >*                                WallNodesUw_2D = NULL;
int                                              NumWallNodes;

int main( int argc, char **argv )
{
    const  float   ver=_VER;
    static char    inFile[256];

    FlowNode2D<double,NUM_COMPONENTS>* TmpMatrixPtr;
   // FlowNodeCore2D<double,NUM_COMPONENTS>* TmpCoreMatrixPtr;

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
            printf("OpenHyperFLOW2D/DEEPS solver v %'.2f ",ver);
            printf(" (parallel CUDA version)\n");
            printf("Copyright (C) 1995-2014 by Serge A. Suchkov\nCopyright policy: LGPL V3\nUsage: %s [{input_data_file}]\n",argv[0]);
            printf("\n\t* Density-based 2D-Navier-Stokes solver for ");
#ifdef _UNIFORM_MESH_
            printf("uniform cartesian mesh");
#else
            printf("non-uniform mesh");
#endif //_UNIFORM_MESH_

            printf("\n");
            Exit_OpenHyperFLOW2D();
        } else {
            sprintf(inFile,"%s",argv[1]);
            Data = new InputData(inFile,DS_FILE,o_stream,0,10);
            if (Data->GetDataError()!=0) {
                *o_stream << "\nInput data error.\n" ;
                 o_stream->flush();
                 Exit_OpenHyperFLOW2D();
            }
        }
//------------------------- CUDA version ------------------------------------
        int num_gpus = 0;   // number of CUDA GPUs
        int max_num_threads = 0;
        int max_thread_block = 0;
        size_t max_gpu_memsize = 0;
        size_t task_size;
        cudaError_t cudaState;

        printf("\n\n\t");
        cudaGetDeviceCount(&num_gpus);

        if (num_gpus < 1) {
         printf("no CUDA capable devices were detected\n");
         Exit_OpenHyperFLOW2D();
        }

        printf("Number of CUDA devices:\t%d\n", num_gpus);

        for (int i = 0; i < num_gpus; i++)
        {
         cudaDeviceProp dprop;
         cudaGetDeviceProperties(&dprop, i);
         printf("\t   %d: %s\n", i, dprop.name);
         // Use only first device
         if( i == 0) {
            max_num_threads  = dprop.maxThreadsPerBlock;
            printf("\t   Max threads per block: %d:\n", max_num_threads);
            max_thread_block = dprop.maxThreadsDim[0];
            printf("\t   Max block size: %d:\n", max_thread_block);
            max_gpu_memsize = dprop.memPitch;
            printf("\t   Max GPU memory size: %d\n", max_gpu_memsize);
            warp_size = dprop.warpSize;
            printf("\t   Warp size: %d\n",warp_size);
            printf("\t   Enable timeout: %i\n\n",dprop.kernelExecTimeoutEnabled);

            break;
         }
        }
        
        // Init shared data
        InitSharedData(Data,&chemical_reactions);        
                                                                                             
        task_size = MaxX*MaxY*(sizeof(FlowNode2D<double,NUM_COMPONENTS>)+sizeof(FlowNodeCore2D<double,NUM_COMPONENTS>))*1.1; // +10% task size

        num_blocks = task_size/(max_gpu_memsize);    // Num mem blocks
        
        if(MaxX*MaxY*sizeof(FlowNode2D<double,NUM_COMPONENTS>) != num_blocks*max_gpu_memsize)
           num_blocks++;

        num_threads = num_blocks;
        
        cudaArraySubmatrix      =  new UArray<  FlowNode2D<double,NUM_COMPONENTS>* >();
        cudaArrayCoreSubmatrix  =  new UArray<  FlowNodeCore2D<double,NUM_COMPONENTS>* >();
        cudaDimArray            =  new UArray<  XY<int> >();
        // Load components properties
        // TODO
        
        
        
        // Init solver (run on host)
        InitDEEPS2D((void*)o_stream);  
        
         cudaState = cudaSetDevice(0);
        
        if(cudaState != cudaSuccess ) {
           *o_stream << "\nError set CUDA device no: "<< 0 << endl;
           Exit_OpenHyperFLOW2D();
        }
        
        cudaState = cudaMalloc( (void**)&cudaHu, sizeof(double)*(NUM_COMPONENTS+1) );
        if(cudaState == cudaErrorMemoryAllocation) {
           *o_stream << "\nError allocate GPU memory for Hu[]" << endl;
           Exit_OpenHyperFLOW2D();
        }
                
        CopyHostToDevice(FlowNode2D<double,NUM_COMPONENTS>::Hu,cudaHu,sizeof(double)*(NUM_COMPONENTS+1));
        
        cudaState = cudaMalloc( (void**)&dt_global, sizeof(float));

        if(cudaState == cudaErrorMemoryAllocation) {
           *o_stream << "\nError allocate GPU memory for dt_global. "<<  endl;
           Exit_OpenHyperFLOW2D();
        }

        // Scan area for seek wall nodes      
        WallNodes = GetWallNodes((ofstream*)o_stream,J,Data->GetIntVal((char*)"isVerboseOutput")); 
        
            
        NumWallNodes = WallNodes->GetNumElements();
        *o_stream << NumWallNodes << " wall nodes found" << endl; 
        
        
        gettimeofday(&mark2,NULL);

            if(NumWallNodes > 0) {
                cudaState = cudaMalloc( (void**)&cudaWallNodes, sizeof(XY<int>)*NumWallNodes );
                if(cudaState == cudaErrorMemoryAllocation) {
                   *o_stream << "\nError allocate GPU memory for WallNodes array." << endl;
                   Exit_OpenHyperFLOW2D();
                }

                CopyHostToDevice(WallNodes->GetArrayPtr(),cudaWallNodes,sizeof(XY<int>)*NumWallNodes);
            }

            TmpMatrixPtr=J->GetMatrixPtr();
            int SubStartIndex, SubMaxX,r_Overlap=0,l_Overlap=0;
            *o_stream << "Allocate SubMatrix:\n";
            SubStartIndex = 0;
            
            int iX0=0;
            
            size_t max_CudaSubmatrixSize     = 0L;
            size_t max_CudaCoreSubmatrixSize = 0L;

            //Allocate max GPU buffers
            for (unsigned int i=0;i<GlobalSubmatrix->GetNumElements();i++) {
                SubStartIndex = GlobalSubmatrix->GetElementPtr(i)->GetX();  
                SubMaxX = GlobalSubmatrix->GetElementPtr(i)->GetY();
                
                if(i == GlobalSubmatrix->GetNumElements()-1)
                  r_Overlap = 0;
                else
                  r_Overlap = 1;
                if(i == 0)
                  l_Overlap = 0;
                else
                  l_Overlap = 1;

                TmpMaxX = (SubMaxX-SubStartIndex)+r_Overlap;
                
                max_CudaSubmatrixSize = max(max_CudaSubmatrixSize,(sizeof(FlowNode2D<double,NUM_COMPONENTS>))*(TmpMaxX*MaxY));
                max_CudaCoreSubmatrixSize = max(max_CudaCoreSubmatrixSize,(sizeof(FlowNodeCore2D<double,NUM_COMPONENTS>))*(TmpMaxX*MaxY));

            }
            
            // Allocate FlowNode2D<double,NUM_COMPONENTS> subdomain
            cudaState = cudaMalloc( (void**)&cudaSubmatrix, max_CudaSubmatrixSize );

            if(cudaState == cudaErrorMemoryAllocation) {
               *o_stream << "\nError allocate GPU memory for Submatrix"<< endl;
               Exit_OpenHyperFLOW2D();
            }
            
            // Allocate FlowNodeCore2D<double,NUM_COMPONENTS> subdomain
            cudaState = cudaMalloc( (void**)&cudaCoreSubmatrix, max_CudaCoreSubmatrixSize );

            if(cudaState == cudaErrorMemoryAllocation) {
               *o_stream << "\nError allocate GPU memory for CoreSubmatrix"<< endl;
               Exit_OpenHyperFLOW2D();
            }

            for (unsigned int i=0;i<GlobalSubmatrix->GetNumElements();i++) {
                
                XY<int> TmpDim;

                SubStartIndex = GlobalSubmatrix->GetElementPtr(i)->GetX();  
                SubMaxX = GlobalSubmatrix->GetElementPtr(i)->GetY();
                
                if(i == GlobalSubmatrix->GetNumElements()-1)
                  r_Overlap = 0;
                else
                  r_Overlap = 1;
                if(i == 0)
                  l_Overlap = 0;
                else
                  l_Overlap = 1;

                TmpMaxX = (SubMaxX-SubStartIndex)+r_Overlap;
                TmpMatrixPtr = (FlowNode2D<double,NUM_COMPONENTS>*)((ulong)J->GetMatrixPtr()+(ulong)(sizeof(FlowNode2D<double,NUM_COMPONENTS>)*(SubStartIndex)*MaxY));                      
                
                int num_cuda_threads = warp_size/4;
                int num_cuda_blocks  = (TmpMaxX*MaxY)/num_cuda_threads;
                if (num_cuda_blocks*num_cuda_threads != TmpMaxX*MaxY)
                    num_cuda_blocks++;

                x0 = SubStartIndex*FlowNode2D<double,NUM_COMPONENTS>::dx;
                
                *o_stream << "SubMatrix("<<i<<")[" << TmpMaxX << "x" << MaxY << "]  Size=" << (ulong)(sizeof(FlowNode2D<double,NUM_COMPONENTS>)*TmpMaxX*MaxY)/(1024*1024) << " Mb\n"; 
                                
                CopyHostToDevice(TmpMatrixPtr,cudaSubmatrix,(sizeof(FlowNode2D<double,NUM_COMPONENTS>))*(TmpMaxX*MaxY));

                if(NumWallNodes > 0) {

                    
                    *o_stream << "CUDA threads: " << num_cuda_threads << endl;
                    *o_stream << "CUDA thread blocks : " << num_cuda_blocks << endl;
                    *o_stream << "\nParallel calc min distance to wall..." << endl;
                    *o_stream << "Run cuda_SetMinDistanceToWall2D kernel..." << flush;
                     
                    cuda_SetMinDistanceToWall2D<<<num_cuda_blocks,num_cuda_threads>>>(cudaSubmatrix,
                                                                                      TmpMaxX*MaxY,
                                                                                      cudaWallNodes,
                                                                                      NumWallNodes,
                                                                                      min(dx,dy),
                                                                                      max((x0+FlowNode2D<double,NUM_COMPONENTS>::dx*TmpMaxX), 
                                                                                             (FlowNode2D<double,NUM_COMPONENTS>::dy*MaxY)),
                                                                                      FlowNode2D<double,NUM_COMPONENTS>::dx,
                                                                                      FlowNode2D<double,NUM_COMPONENTS>::dy);       
                     
                     CUDA_BRRIER("cuda_SetMinDistanceToWall2D");
                     *o_stream << "OK" << endl;
                     
                      
                     *o_stream << "Run cuda_Recalc_y_plus kernel..." << flush;
                     cuda_Recalc_y_plus<<<num_cuda_blocks,num_cuda_threads>>>(cudaSubmatrix,
                                                                              TmpMaxX*MaxY,
                                                                              cudaWallNodes,
                                                                              NumWallNodes,
                                                                              min(dx,dy),
                                                                              max((x0+FlowNode2D<double,NUM_COMPONENTS>::dx*TmpMaxX), 
                                                                                     (FlowNode2D<double,NUM_COMPONENTS>::dy*MaxY)),
                                                                              FlowNode2D<double,NUM_COMPONENTS>::dx,
                                                                              FlowNode2D<double,NUM_COMPONENTS>::dy,
                                                                              MaxY);
                     
                     CUDA_BRRIER("cuda_Recalc_y_plus");
                     *o_stream << "OK" << endl;
                     

               }
                
                *o_stream << "Run cuda_SetInitBoundaryLayer kernel..." << flush;
                cuda_SetInitBoundaryLayer<<<num_cuda_blocks,num_cuda_threads>>>(cudaSubmatrix,
                                                                                TmpMaxX*MaxY, iX0, MaxY, 
                                                                                delta_bl,
                                                                                SigW,SigF,(TurbulenceExtendedModel)TurbExtModel, 
                                                                                FlowNode2D<double,NUM_COMPONENTS>::dx,
                                                                                FlowNode2D<double,NUM_COMPONENTS>::dy,
                                                                                cudaHu,
                                                                                FlowNode2D<double,NUM_COMPONENTS>::isSrcAdd,
                                                                                FlowNode2D<double,NUM_COMPONENTS>::FT);

                CUDA_BRRIER("cuda_SetInitBoundaryLayer");
                *o_stream << "OK" << endl;
                
                CopyDeviceToHost(cudaSubmatrix,TmpMatrixPtr,(sizeof(FlowNode2D<double,NUM_COMPONENTS>))*(TmpMaxX*MaxY));

                iX0 += TmpMaxX;
                TmpDim.SetXY(TmpMaxX,MaxY);
                cudaDimArray->AddElement(&TmpDim);
                o_stream->flush();
           }
            gettimeofday(&mark1,NULL);
            *o_stream << "OK\n" << "Time: " << (double)(mark1.tv_sec-mark2.tv_sec)+(double)(mark1.tv_usec-mark2.tv_usec)*1.e-6 << " sec." << endl; 
            
            DEEPS2D_Run((ofstream*)o_stream,        // ofstream* o_stream,                                 
                        J,                          // UMatrix2D< FlowNode2D< double,NUM_COMPONENTS > >*  
                        C,
                        cudaSubmatrix,              // FlowNode2D< double,NUM_COMPONENTS >*     
                        cudaCoreSubmatrix,          // FlowNodeCore2D< double,NUM_COMPONENTS >*  
                        cudaDimArray,               // UArray< XY<int> >*                                  
                        max_thread_block);          // int max_thread_block                            
            
            
            cudaState = cudaFree(cudaSubmatrix);
            
            if(cudaState != cudaSuccess ) {
               *o_stream << "\nError free Submatrix  from GPU memory." << endl;
               Exit_OpenHyperFLOW2D();
            }

            cudaState = cudaFree(cudaCoreSubmatrix);

            if(cudaState != cudaSuccess ) {
               *o_stream << "\nError free CoreSubmatrix from GPU memory." << endl;
               Exit_OpenHyperFLOW2D();
            }

            cudaState = cudaFree(cudaWallNodes);
            if(cudaState != cudaSuccess ) {
               *o_stream << "\nError free WallNodes array from GPU memory." << endl;
               Exit_OpenHyperFLOW2D();
            }
            
            cudaState = cudaFree(dt_global);
            if(cudaState != cudaSuccess ) {
               *o_stream << "\nError free dt_global from GPU memory." << endl;
               Exit_OpenHyperFLOW2D();
            }

            //>>>>>>>>>>>>>>>>>>>>>>
            //DataSnapshot(OutFileName,WM_REWRITE);
            //>>>>>>>>>>>>>>>>>>>>>>

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
        Exit_OpenHyperFLOW2D();
}
#endif //_CUDA_


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
*   deeps2d_core.cpp: CUDA kernels code.                                       *
*                                                                              *
*  last update: 14/04/2016                                                     *
********************************************************************************/

#ifdef _CUDA_
#define _PARALLEL_ONLY
#define _DEVICE_
#include "libDEEPS2D/deeps2d_core.hpp"
#undef _DEVICE_
#ifdef _CUDA_
__host__
#endif //_CUDA_ 
int LoadTable2GPU(Table* Src, Table*& Dst, int i_dev)
{
 Table* pTmpTable;
 pTmpTable = new Table(NULL,0); 
 pTmpTable->n = Src->n;
 
 if(cudaSetDevice(i_dev) != cudaSuccess ) {
    printf("\nError set CUDA device no: %d\n",i_dev);
    Exit_OpenHyperFLOW2D(1);
 }

 if(cudaMalloc( (void**)&pTmpTable->x, sizeof(FP)*Src->n ) == cudaErrorMemoryAllocation) {
    printf("\nError allocate GPU memory for Table %s\n",Src->GetName());
    Exit_OpenHyperFLOW2D(1);
 }

 if(cudaMalloc( (void**)&pTmpTable->y, sizeof(FP)*Src->n ) == cudaErrorMemoryAllocation) {
   printf("\nError allocate GPU memory for Table %s",Src->GetName());
   Exit_OpenHyperFLOW2D(1);
 }

 CopyHostToDevice(Src->x,pTmpTable->x,sizeof(FP)*Src->n); 
 CopyHostToDevice(Src->y,pTmpTable->y,sizeof(FP)*Src->n);
 
 if(cudaMalloc( (void**)&Dst, sizeof(Table) ) == cudaErrorMemoryAllocation) {
   printf("\nError allocate GPU memory for Table %s\n",Src->GetName());
   Exit_OpenHyperFLOW2D(1);
 }

 CopyHostToDevice(pTmpTable,Dst,sizeof(Table));
 
 return Src->n;
 
}

#ifdef _CUDA_
 __host__ __device__
#endif //_CUDA_ 
inline FP GetVal(register Table* t,
                 register FP _x ) {
    if ( !t )
        return 0.;
    
    register int  i, _n = t->n;
    
    register FP _y;

    if ( _n == 1 )
        return( t->y[0] );

    if ( _x <= t->x[0] ) {
        i = 1;
        goto EndGetVal;
    }

    if ( _x >= t->x[t->n-1] ) {
        i = _n - 1;
        goto EndGetVal;
    }
 
 #pragma unroll
    for ( i=1; i<_n; i++ ) {
        if ( (_x >= t->x[i-1]) && (_x < t->x[i]) )
            break;
    }

    EndGetVal:

    _y = t->y[i] + (t->y[i-1] - t->y[i])*(_x - t->x[i])/(t->x[i-1] - t->x[i]);

    return( _y );
}

#ifdef _CUDA_
__device__
#endif //_CUDA_ 
int inline cuda_CalcChemicalReactions(register FlowNode2D<FP,NUM_COMPONENTS>* CalcNode,
                                      register ChemicalReactionsModel cr_model, 
                                      register void* CRM_data,
                                      register SolverMode sm) {
    
    register ChemicalReactionsModelData2D* model_data = (ChemicalReactionsModelData2D*)CRM_data;
    register FP Yfu,Yox,Ycp,Yair;
    register FP rho_1=1.0/CalcNode->S[i2d_Rho];

    Yfu  = CalcNode->S[i2d_Yfu]*rho_1; // Fuel
    Yox  = CalcNode->S[i2d_Yox]*rho_1; // OX
    Ycp  = CalcNode->S[i2d_Ycp]*rho_1; // CP
    Yair = 1. - (Yfu+Yox+Ycp);         // Air (inert component)
        
    if(cr_model==CRM_ZELDOVICH) {
        if ( !CalcNode->isCond2D(CT_Y_CONST_2D) ) {
//--- chemical reactions (Zeldovich model) -------------------------------------------------->
            if ( !CalcNode->isCond2D(CT_Y_CONST_2D) ) {
                if ( CalcNode->Tg > CalcNode->Tf ) {
                      if ( Yox > Yfu*model_data->K0 * model_data->gamma ) { // Yox > Yfuel
                           Yox -= Yfu * model_data->K0 * model_data->gamma;
                           Ycp += Yfu * (1.0 + model_data->K0) * model_data->gamma;
                           Yfu -= Yfu * model_data->gamma;
                      } else {                                              // Yox < Yfuel
                           Yfu -= Yox/model_data->K0 * model_data->gamma;
                           Ycp += Yox * (1.0 + 1.0/model_data->K0) * model_data->gamma;
                           Yox -= Yox * model_data->gamma;
                      }
                   }
              }
            CalcNode->S[i2d_Yfu] = Yfu*CalcNode->S[i2d_Rho];
            CalcNode->S[i2d_Yox] = Yox*CalcNode->S[i2d_Rho];
            CalcNode->S[i2d_Ycp] = Ycp*CalcNode->S[i2d_Rho];
        }
//--- chemical reactions (Zeldovich model) -------------------------------------------------->
    } else if (cr_model==CRM_ARRENIUS) {
//--- chemical reactions (Arrenius model) -------------------------------------------------->
//--- chemical reactions (Arrenius model) -------------------------------------------------->
        CalcNode->S[i2d_Yfu] = Yfu*CalcNode->S[i2d_Rho];
        CalcNode->S[i2d_Yox] = Yox*CalcNode->S[i2d_Rho];
        CalcNode->S[i2d_Ycp] = Ycp*CalcNode->S[i2d_Rho];
    }
        
    CalcNode->R   = model_data->R_Fuel*Yfu+
                    model_data->R_OX*Yox+
                    model_data->R_cp*Ycp+
                    model_data->R_air*Yair;
    CalcNode->CP  = GetVal(model_data->Cp_Fuel,CalcNode->Tg)*Yfu+
                    GetVal(model_data->Cp_OX,CalcNode->Tg)*Yox+
                    GetVal(model_data->Cp_cp,CalcNode->Tg)*Ycp+
                    GetVal(model_data->Cp_air,CalcNode->Tg)*Yair;

    if ( sm == SM_NS ) {
        CalcNode->mu  = GetVal(model_data->mu_Fuel,CalcNode->Tg)*Yfu+
                        GetVal(model_data->mu_OX,CalcNode->Tg)*Yox+
                        GetVal(model_data->mu_cp,CalcNode->Tg)*Ycp+
                        GetVal(model_data->mu_air,CalcNode->Tg)*Yair;
        CalcNode->lam = GetVal(model_data->lam_Fuel,CalcNode->Tg)*Yfu+
                        GetVal(model_data->lam_OX,CalcNode->Tg)*Yox+
                        GetVal(model_data->lam_cp,CalcNode->Tg)*Ycp+
                        GetVal(model_data->lam_air,CalcNode->Tg)*Yair;
    }
    
    CalcNode->Y[h_fu]  = Yfu;
    CalcNode->Y[h_ox]  = Yox;
    CalcNode->Y[h_cp]  = Ycp;
    CalcNode->Y[h_air] = Yair;

    return 1;
}

void SetP2PAccess(int dev1, int dev2) {
 cudaError_t cudaState;
 int canAccess = 0;

 cudaState = cudaDeviceCanAccessPeer(&canAccess,dev1,dev2);

 if(cudaState != cudaSuccess) {
     printf("\nError set P2P access for devices %i<-->%i\n",dev1,dev2);
     printf("%s\n", cudaGetErrorString( cudaGetLastError() ) );
     Exit_OpenHyperFLOW2D(1);
 }

 if(!canAccess) {
     cudaSetDevice(dev1);
     cudaState = cudaDeviceEnablePeerAccess(dev2,0);

     if(cudaState != cudaSuccess) {
         printf("\nError set P2P access for devices %i<-->%i\n",dev1,dev2);
         printf("%s\n", cudaGetErrorString( cudaGetLastError() ) );
         Exit_OpenHyperFLOW2D(1);
     }
 }

 canAccess = 0;

 cudaSetDevice(dev2);

 cudaState = cudaDeviceCanAccessPeer(&canAccess,dev2,dev1);

 if(cudaState != cudaSuccess) {
    printf("\nError set P2P access for devices %i<-->%i\n",dev2,dev1);
    printf("%s\n", cudaGetErrorString( cudaGetLastError() ) );
    Exit_OpenHyperFLOW2D(1);
 }
 
 if(!canAccess) {
     cudaSetDevice(dev2);
     cudaState = cudaDeviceEnablePeerAccess(dev1,0);

     if(cudaState != cudaSuccess) {
         printf("\nError set P2P access for devices %i<-->%i\n",dev2,dev1);
         printf("%s\n", cudaGetErrorString( cudaGetLastError() ) );
         Exit_OpenHyperFLOW2D(1);
     }
 }
}


void DisableP2PAccess(int dev1, int dev2) {
 cudaError_t cudaState;
 int canAccess = 0;

 cudaSetDevice(dev1);

 cudaState = cudaDeviceCanAccessPeer(&canAccess,dev1,dev2);

 if(cudaState != cudaSuccess) {
     printf("\nError probe P2P access for devices %i-->%i\n",dev1,dev2);
     printf("%s\n", cudaGetErrorString( cudaGetLastError() ) );
 }

 if(!canAccess) {
     cudaState = cudaDeviceDisablePeerAccess(dev2);

     if(cudaState != cudaSuccess) {
         printf("\nError unset P2P access for devices %i-->%i\n",dev1,dev2);
         printf("%s\n", cudaGetErrorString( cudaGetLastError() ) );
         Exit_OpenHyperFLOW2D(1);
     }
 }

 canAccess = 1;

 cudaSetDevice(dev2);

 cudaState = cudaDeviceCanAccessPeer(&canAccess,dev2,dev1);

 if(cudaState != cudaSuccess) {
    printf("\nError probe P2P access for devices %i-->%i\n",dev2,dev1);
    printf("%s\n", cudaGetErrorString( cudaGetLastError() ) );
 }

 if(!canAccess) {
     cudaState = cudaDeviceDisablePeerAccess(dev1);

     if(cudaState != cudaSuccess) {
         printf("\nError unset P2P access for devices %i-->%i\n",dev2,dev1);
         printf("%s\n", cudaGetErrorString( cudaGetLastError() ) );
         Exit_OpenHyperFLOW2D(1);
     }
 }
}

void CUDA_BARRIER(char* KernelName) {
    cudaError_t cudaState = cudaDeviceSynchronize();
    if(cudaState != cudaSuccess) {
        printf("\nError in %s kernel...\n",KernelName);
        printf("%s\n", cudaGetErrorString( cudaGetLastError() ) );
        Exit_OpenHyperFLOW2D(1);
    }
}

void CopyDeviceToDeviceP2P(void*  src, int src_dev,
                           void*  dst, int dst_dev,
                           size_t length, 
                           void*  host_buff,
                           cudaStream_t cuda_stream) {
cudaError_t cudaState;
#ifdef _P2P_ACCESS_
    cudaState = cudaMemcpyPeerAsync(dst, dst_dev, src, src_dev, length, cuda_stream); 
#else
    cudaState = cudaSetDevice(src_dev);
    CopyDeviceToHost(src,host_buff,length,cuda_stream);
    cudaState = cudaSetDevice(dst_dev);
    CopyHostToDevice(host_buff,dst,length,cuda_stream);
#endif // _P2P_ACCESS_
    
    if(cudaState != cudaSuccess) {
     printf("\nError P2P copy device to device...\n");
     printf("%s\n", cudaGetErrorString( cudaGetLastError() ) );

       Exit_OpenHyperFLOW2D(1);
    }
}


void CopyDeviceToDevice(void* src, void* dst, size_t length, cudaStream_t stream) {
    cudaError_t cudaState = cudaMemcpyAsync(dst, src, length,cudaMemcpyDeviceToDevice, stream);
    if(cudaState != cudaSuccess) {
     printf("\nError copy device to device...\n");
     printf("%s\n", cudaGetErrorString( cudaGetLastError() ) );

       Exit_OpenHyperFLOW2D(1);
    }
}
#ifdef _CUDA_
__device__
#endif //_CUDA_ 
void  kernel_CopyDeviceToDevice(char* dst, char* src, size_t length)
{
#pragma unroll
    for (register size_t i=0; i < length; i++ ) {
         dst[i] = src[i]; 
    }
}

void CopyHostToDevice(void* src, void* dst, size_t length, cudaStream_t stream) {
    cudaError_t cudaState = cudaMemcpyAsync(dst, src, length,cudaMemcpyHostToDevice, stream);
    if(cudaState != cudaSuccess) {
     printf("\nError copy host to device...\n");
     printf("%s\n", cudaGetErrorString( cudaGetLastError() ) );
       Exit_OpenHyperFLOW2D(1);
    }
}

void SyncCopyHostToDevice(void* src, void* dst, size_t length) {
    cudaError_t cudaState = cudaMemcpy(dst, src, length,cudaMemcpyHostToDevice);
    if(cudaState != cudaSuccess) {
     printf("\nError copy host to device...\n");
     printf("%s\n", cudaGetErrorString( cudaGetLastError() ) );
       Exit_OpenHyperFLOW2D(1);
    }
}

void CopyDeviceToHost(void* src, void* dst, size_t length, cudaStream_t stream) {
    cudaError_t cudaState = cudaMemcpyAsync(dst, src, length,cudaMemcpyDeviceToHost, stream);
    if(cudaState != cudaSuccess) {
     printf("\nError copy device to host...\n");
     printf("%s\n", cudaGetErrorString( cudaGetLastError() ) );
       Exit_OpenHyperFLOW2D(1);
    }
}

void SyncCopyDeviceToHost(void* src, void* dst, size_t length) {
    cudaError_t cudaState = cudaMemcpy(dst, src, length,cudaMemcpyDeviceToHost);
    if(cudaState != cudaSuccess) {
     printf("\nError copy device to host...\n");
     printf("%s\n", cudaGetErrorString( cudaGetLastError() ) );
       Exit_OpenHyperFLOW2D(1);
    }
}

__global__ void 
cuda_SetInitBoundaryLayer(FlowNode2D<FP,NUM_COMPONENTS>* pJ2D,
                          unsigned long int index_limit,
                          int X0, int MAX_Y,
                          FP delta,
                          FP sig_w, 
                          FP sig_f,
                          TurbulenceExtendedModel etm,
                          FP _dx, FP _dy,
                          FP* _Hu,
                          int _isSrcAdd,
                          FlowType _FT,
                          SolverMode sm) {

    unsigned long int index = threadIdx.x + blockIdx.x * blockDim.x;

    if(index < index_limit) {

              register FlowNode2D< FP,NUM_COMPONENTS >* CurrentNode=&pJ2D[index];

              if(CurrentNode->CT != (ulong)(CT_SOLID_2D | CT_NODE_IS_SET_2D)) {

                  if(delta > 0.0  &&  CurrentNode->l_min <= delta) {
                     CurrentNode->S[i2d_RhoU] = CurrentNode->S[i2d_RhoU] * CurrentNode->l_min/delta;
                     CurrentNode->S[i2d_RhoV] = CurrentNode->S[i2d_RhoV] * CurrentNode->l_min/delta;
                     CurrentNode->FillNode2D(0,1,sig_w,sig_f,etm,delta,_dx,_dy,_Hu,_isSrcAdd,sm,_FT);
                  }

               if(CurrentNode->CT != (ulong)(NT_FC_2D)) {

                  register int  i = CurrentNode->ix - X0;
                  register int  j = CurrentNode->iy;

                  register int  n1 = CurrentNode->idXl;
                  register int  n2 = CurrentNode->idXr;
                  register int  n3 = CurrentNode->idYu;
                  register int  n4 = CurrentNode->idYd;

                  register int  N1 = i - n1;
                  register int  N2 = i + n2;
                  register int  N3 = j + n3;
                  register int  N4 = j - n4;

                  CurrentNode->UpNode    = &pJ2D[i*MAX_Y + N3];
                  CurrentNode->DownNode  = &pJ2D[i*MAX_Y + N4];
                  CurrentNode->RightNode = &pJ2D[N2*MAX_Y + j];
                  CurrentNode->LeftNode  = &pJ2D[N1*MAX_Y + j];
               }
          }
   }
}

__global__ void
cuda_SetMinDistanceToWall2D(FlowNode2D<FP,NUM_COMPONENTS>* pJ2D,
                            unsigned long int index_limit,
                            XY<int>* WallNodes2D, 
                            int NumWallNodes2D,
                            FP min_l_min,
                            FP max_l_min,
                            FP _dx, FP _dy,
                            FP  x0)   {

   unsigned long int index = threadIdx.x + blockIdx.x * blockDim.x;

   if(index < index_limit) {

       FlowNode2D<FP,NUM_COMPONENTS>* TmpNode = &pJ2D[index];

       if(TmpNode->CT != (ulong)(CT_SOLID_2D | CT_NODE_IS_SET_2D)) {   

           TmpNode->l_min = max_l_min;
 #pragma unroll
           for (int ii=0;ii<NumWallNodes2D;ii++) {

              register XY<int>* TmpWallNode = &WallNodes2D[ii]; 

              register FP L_x   = (TmpWallNode->X - TmpNode->ix)* _dx;// + x0;  
              register FP L_y   = (TmpWallNode->Y - TmpNode->iy)* _dy;          
              register FP l_min = sqrt(L_x*L_x + L_y*L_y);                      

              TmpNode->l_min = max(min(TmpNode->l_min,l_min),min_l_min);
              
              if (TmpNode->l_min == l_min) {
                  TmpNode->i_wall = TmpWallNode->X + (int)(x0/_dx);
                  TmpNode->j_wall = TmpWallNode->Y;
              }
            }
        }
   }
}
__global__ void 
cuda_Recalc_y_plus(FlowNode2D<FP,NUM_COMPONENTS>* pJ2D,
                   unsigned long int index_limit,
                   XY<int>* WallNodes2D,
                   int NumWallNodes2D,
                   FP min_l_min,
                   FP max_l_min,
                   FP _dx, 
                   FP _dy,
                   int max_y) {

    register unsigned long int index = threadIdx.x + blockIdx.x * blockDim.x;

    if(index < index_limit) {
        
        register FP tau_w;
        register FlowNode2D<FP,NUM_COMPONENTS>* TmpNode = &pJ2D[index];

        if(TmpNode->CT != (ulong)(CT_SOLID_2D | CT_NODE_IS_SET_2D)) {
           
           register unsigned long int wall_index = TmpNode->i_wall*max_y + TmpNode->j_wall;

           if(wall_index < index_limit) {

               register FlowNode2D<FP,NUM_COMPONENTS>* WallNode = &pJ2D[wall_index];   // x*nY + y
               
               tau_w = (fabs(WallNode->dUdy) + fabs(WallNode->dVdx)) * WallNode->mu;
               
               if(tau_w > 0.0) {
                   FP U_w   = sqrt(tau_w/WallNode->S[i2d_Rho]);
                   TmpNode->y_plus = U_w*TmpNode->l_min*TmpNode->S[i2d_Rho]/TmpNode->mu;
               } else {
                   TmpNode->y_plus = 0.0;
               }
           }
        }
     }
}
__global__  void
cuda_DEEPS2D_Stage1(register FlowNode2D<FP,NUM_COMPONENTS>*     pLJ,
                    register FlowNodeCore2D<FP,NUM_COMPONENTS>* pLC,
                    register unsigned long int index_limit,
                    register unsigned long r_limit,
                    register unsigned long l_limit,
                    register FP dxx, register FP dyy,
                    register FP dtdx, register FP dtdy,
                    register FP _dt,
                    register int _FT, register int Num_Eq,
                    register SolverMode sm,
                    register int isLocalTimeStep) {

    register size_t index = threadIdx.x + blockIdx.x * blockDim.x;

    if(index < index_limit) {
        FlowNode2D< FP,NUM_COMPONENTS >* CurrentNode=&pLJ[index];
          if(CurrentNode->CT != (ulong)(CT_SOLID_2D | CT_NODE_IS_SET_2D) &&
             CurrentNode->CT != (ulong)(NT_FC_2D) &&
             CurrentNode->ix <  r_limit &&
             CurrentNode->ix >= l_limit ) {
              register FlowNodeCore2D< FP,NUM_COMPONENTS >* NextNode=&pLC[index];
              register int  n1 = CurrentNode->idXl; 
              register int  n2 = CurrentNode->idXr;
              register int  n3 = CurrentNode->idYu;
              register int  n4 = CurrentNode->idYd;
              register FP   n_n_1; 
              register FP   m_m_1;
              register FP   dt;
              register FP   _dtdx;
              register FP   _dtdy;
              register FP   _y_1 = 1.0/(CurrentNode->iy+0.5);

              n_n_1 = 1./max(n1+n2,1); 
              m_m_1 = 1./max(n3+n4,1); 
              
              if (isLocalTimeStep) {
                  dt = CurrentNode->dt_local;
                  _dtdx = dt*dtdx;
                  _dtdy = dt*dtdy;
              } else {
                  dt    = _dt;
                  _dtdx = dtdx;
                  _dtdy = dtdy;
              }
              register FlowNode2D< FP,NUM_COMPONENTS >* UpNode    = CurrentNode->UpNode;
              register FlowNode2D< FP,NUM_COMPONENTS >* DownNode  = CurrentNode->DownNode;
              register FlowNode2D< FP,NUM_COMPONENTS >* RightNode = CurrentNode->RightNode;
              register FlowNode2D< FP,NUM_COMPONENTS >* LeftNode  = CurrentNode->LeftNode;
              
              // Scan equation system ... k - number of equation
#pragma unroll
              for (int k=0;k<Num_Eq;k++ ) {
                  register int c_flag, dx_flag, dx2_flag;
                  register int dy_flag, dy2_flag;
                  register FP  dXX,dYY;
                  
                  register FP beta = CurrentNode->beta[k];  
                  register FP _beta = 1. - beta;             

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
                    } else if (sm == SM_NS &&
                               ((CurrentNode->isTurbulenceCond2D(TCT_k_eps_Model_2D) ||
                                CurrentNode->isTurbulenceCond2D(TCT_Spalart_Allmaras_Model_2D)))) {
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
                    } else if(sm == SM_NS &&
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
                        if ( dx_flag) {
                            dXX = CurrentNode->dSdx[k] = (RightNode->A[k]-LeftNode->A[k])*n_n_1; //  
                        } else {
                            CurrentNode->S[k] = (LeftNode->S[k]*n2+RightNode->S[k]*n1)*n_n_1;
                            dXX = CurrentNode->dSdx[k] = 0.;                                     //    
                        }
                        if ( dy_flag ) {
                            dYY = CurrentNode->dSdy[k] = (UpNode->B[k]-DownNode->B[k])*m_m_1;    // 
                        } else {
                            CurrentNode->S[k] =  (UpNode->S[k]*n3+DownNode->S[k]*n4)*m_m_1;
                            dYY = CurrentNode->dSdy[k] = 0;                                      // 
                        }

                        // Cauchy BC
                        if ( dx2_flag ) {
                            dXX = (LeftNode->dSdx[k]+RightNode->dSdx[k])*0.5;
                        }
                        if ( dy2_flag ) {
                            dYY = (UpNode->dSdy[k]+DownNode->dSdy[k])*0.5;
                        }
                        
                        if ( _FT ) {
                            NextNode->S[k] = CurrentNode->S[k]*beta+_beta*(dxx*(LeftNode->S[k]+RightNode->S[k])+dyy*(UpNode->S[k]+DownNode->S[k]))*0.5
                                          - (_dtdx*dXX+_dtdy*(dYY+CurrentNode->F[k]*_y_1)) + (CurrentNode->Src[k])*dt+CurrentNode->SrcAdd[k];
                        } else {
                            NextNode->S[k] = CurrentNode->S[k]*beta+_beta*(dxx*(LeftNode->S[k]+RightNode->S[k])+dyy*(UpNode->S[k]+DownNode->S[k]))*0.5
                                          - (_dtdx*dXX+_dtdy*dYY) + (CurrentNode->Src[k])*dt+CurrentNode->SrcAdd[k];
                        }
                        
                    }
                }
             }
          }
 }

__global__  void 
cuda_DEEPS2D_Stage2(register FlowNode2D<FP,NUM_COMPONENTS>*     pLJ,
                    register FlowNodeCore2D<FP,NUM_COMPONENTS>* pLC,
                    register unsigned long int index_limit,
                    register unsigned long r_limit,
                    register unsigned long l_limit,
                    register FP beta_init, register FP  beta0, 
                    register int b_FF, register FP CFL, register FP nrbc_beta0,
                    register int Chemical_reactions_model,
                    register ChemicalReactionsModelData2D* pCRMD,
                    register int noTurbCond,
                    register FP SigW, register FP SigF, 
                    register FP dx_1, register FP dy_1, register FP delta_bl,
                    register FP dx, register FP dy,
                    register FlowType _FT, register int Num_Eq,
#ifdef _RMS_             
                    register int isAlternateRMS,
                    register RMS_pack* RMS_PACK,
                    register FP int2float_RMS_scale,    
#endif // _RMS_
                    register FP* _Hu,
                    register int _isSrcAdd,
                    register unsigned int* dt_min_device, register FP int2float_scale,
                    register TurbulenceExtendedModel TurbExtModel, 
                    register SolverMode sm,
                    register int isLocalTimeStep) {


    register unsigned long int index = threadIdx.x + blockIdx.x * blockDim.x;

    if(index < index_limit) {
        FlowNode2D< FP,NUM_COMPONENTS >* CurrentNode=&pLJ[index];

       if(CurrentNode->CT != (ulong)(CT_SOLID_2D | CT_NODE_IS_SET_2D) &&
          CurrentNode->CT != (ulong)(NT_FC_2D) &&
          CurrentNode->ix >= l_limit &&  
          CurrentNode->ix <  r_limit ) { 

              register  FP  beta_min;
              register FlowNodeCore2D< FP,NUM_COMPONENTS >* NextNode=&pLC[index];
              register int  n1 = CurrentNode->idXl; 
              register int  n2 = CurrentNode->idXr;
              register int  n3 = CurrentNode->idYu;
              register int  n4 = CurrentNode->idYd;
              
              register FP n_n_1 = 1./max(n1+n2,1);  
              register FP m_m_1 = 1./max(n3+n4,1);  

              register FlowNode2D< FP,NUM_COMPONENTS >* UpNode    = CurrentNode->UpNode;
              register FlowNode2D< FP,NUM_COMPONENTS >* DownNode  = CurrentNode->DownNode;
              register FlowNode2D< FP,NUM_COMPONENTS >* RightNode = CurrentNode->RightNode;
              register FlowNode2D< FP,NUM_COMPONENTS >* LeftNode  = CurrentNode->LeftNode;

              register FP dx_1xn_n_1;
              register FP dy_1xm_m_1;

              dx_1xn_n_1=dx_1*n_n_1; 
              dy_1xm_m_1=dy_1*m_m_1; 
              
              // Scan equation system ... k - number of equation
#pragma unroll
              for (int k=0;k<Num_Eq;k++ ) {

                  register int    c_flag = 0;

                  if ( k < 4 ) // Make bit flags for future test for current equation 
                      c_flag  = CT_Rho_CONST_2D   << k;
                  else if (k<(4+NUM_COMPONENTS))  // 7 ?
                      c_flag  = CT_Y_CONST_2D;
                  else if(sm == SM_NS &&
                          (CurrentNode->isTurbulenceCond2D(TCT_k_eps_Model_2D) ||
                           CurrentNode->isTurbulenceCond2D(TCT_Spalart_Allmaras_Model_2D) )) {
                      c_flag  = TCT_k_CONST_2D << (k-4-NUM_COMPONENTS);
                  }
                                    
                  if ( !CurrentNode->isCond2D((CondType2D)c_flag) && 
                        CurrentNode->S[k] != 0. ) {
                        
                        register FP sqrt_RES;
                        register FP absDD;
                        register FP DD_local;
                        register FP Tmp;

                        if(k == i2d_RhoU && k == i2d_RhoV ) {
                            Tmp = max(fabs(CurrentNode->S[i2d_RhoU]),fabs(CurrentNode->S[i2d_RhoV]));// max Flux
                        } else {
                            Tmp = CurrentNode->S[k];
                        }
                        
                        absDD = NextNode->S[k]-CurrentNode->S[k];

                        if(Tmp != 0.0 && absDD !=0  &&
                           !isnan(Tmp) && !isnan(absDD)) {
                            DD_local = fabs(absDD/Tmp);
                            sqrt_RES = sqrt(DD_local);
                        } else {
                            DD_local = 0.0;
                            sqrt_RES = 0.0;
                            Tmp      = 0.0;
                            absDD    = 0.0;
                        }
                        
                        if(CurrentNode->isCond2D(CT_NONREFLECTED_2D)) {
                           beta_min = nrbc_beta0;
                        } else {
                           beta_min = min(beta0,beta_init); 
                        }
                        
                        if( b_FF == BFF_L) {
                         //LINEAR locally adopted blending factor function  (LLABFF)
                           CurrentNode->beta[k] = min(beta_min,(beta_min*beta_min)/(beta_min+DD_local));
                         } else if( b_FF == BFF_LR) {
                         //LINEAR locally adopted blending factor function with relaxation (LLABFFR)
                           CurrentNode->beta[k] = min((beta_min+CurrentNode->beta[k])*0.5,(beta_min*beta_min)/(beta_min+DD_local));
                         } else if( b_FF == BFF_S) {
                         //SQUARE locally adopted blending factor function (SLABF)
                           CurrentNode->beta[k] = min(beta_min,(beta_min*beta_min)/(beta_min+DD_local*DD_local));
                         } else if (b_FF == BFF_SR) {
                         //SQUARE locally adopted blending factor function with relaxation (SLABFFR)
                           CurrentNode->beta[k] = min((beta_min+CurrentNode->beta[k])*0.5,(beta_min*beta_min)/(beta_min+DD_local*DD_local));
                         } else if( b_FF == BFF_SQR) {
                         //SQRT() locally adopted blending factor function (SQRLABF)
                           CurrentNode->beta[k] = min(beta_min,(beta_min*beta_min)/(beta_min+sqrt_RES));
                         } else if( b_FF == BFF_SQRR) {
                         //SQRT() locally adopted blending factor function with relaxation (SQRLABFFR)
                           CurrentNode->beta[k] = min((beta_min+CurrentNode->beta[k])*0.5,(beta_min*beta_min)/(beta_min+sqrt_RES)); 
                         } else {
                           // Default->SQRLABF */
                           CurrentNode->beta[k] = min(beta_min,(beta_min*beta_min)/(beta_min+sqrt_RES));
                         }
#ifdef _RMS_
                         if (isAlternateRMS) {
                              atomicAdd(&RMS_PACK->sumDiv[k],(float)(Tmp*Tmp));
                              atomicAdd(&RMS_PACK->RMS[k],(float)(absDD*absDD));
                          } else {
                              atomicAdd(&RMS_PACK->RMS[k],(float)(DD_local*DD_local));
                              atomicAdd(&RMS_PACK->sum_iRMS[k],(unsigned long long int)(1L));
                          }
#endif // RMS
                  }
                  
                  if (k<(4+NUM_COMPONENTS)) {
                      if ( !CurrentNode->isCond2D((CondType2D)c_flag) )
                            CurrentNode->S[k]   = NextNode->S[k];
                  } else if (sm == SM_NS &&
                             (CurrentNode->isTurbulenceCond2D(TCT_k_eps_Model_2D) ||
                              CurrentNode->isTurbulenceCond2D(TCT_Spalart_Allmaras_Model_2D)) ){
                      if ( !CurrentNode->isTurbulenceCond2D((TurbulenceCondType2D)c_flag) )
                            CurrentNode->S[k]   =  NextNode->S[k];
                  }
              }

              CurrentNode->beta[i2d_RhoV] = CurrentNode->beta[i2d_RhoU] = max(CurrentNode->beta[i2d_RhoU],CurrentNode->beta[i2d_RhoV]);  // for symmetry keeping

              if(sm == SM_NS) {

                  CurrentNode->droYdx[NUM_COMPONENTS]=CurrentNode->droYdy[NUM_COMPONENTS]=0.;

                  register FP  rhoY_air_Right = RightNode->S[i2d_Rho];
                  register FP  rhoY_air_Left  = LeftNode->S[i2d_Rho]; 
                  register FP  rhoY_air_Up    = UpNode->S[i2d_Rho];   
                  register FP  rhoY_air_Down  = DownNode->S[i2d_Rho]; 
                  register FP  rho_1          = 1./CurrentNode->S[i2d_Rho];
                  
                  CurrentNode->droYdx[NUM_COMPONENTS]=CurrentNode->droYdy[NUM_COMPONENTS]=0.;
#pragma unroll
                  for (int k=4;k<FlowNode2D<FP,NUM_COMPONENTS>::NumEq-2;k++ ) {
                      if ( !CurrentNode->isCond2D(CT_dYdx_NULL_2D) ) {
                          CurrentNode->droYdx[k-4]=(RightNode->S[k]-LeftNode->S[k])*dx_1xn_n_1;
                          rhoY_air_Right -= RightNode->S[k];
                          rhoY_air_Left  -= LeftNode->S[k];
                      }
                      if ( !CurrentNode->isCond2D(CT_dYdy_NULL_2D) ) {
                            CurrentNode->droYdy[k-4]=(UpNode->S[k]-DownNode->S[k])*dy_1xm_m_1;
                            rhoY_air_Up    -= UpNode->S[k];
                            rhoY_air_Down  -= DownNode->S[k];
                      }
                  }

                  if ( !CurrentNode->isCond2D(CT_dYdx_NULL_2D) ) {
                      CurrentNode->droYdx[NUM_COMPONENTS]=(rhoY_air_Right - rhoY_air_Left)*dx_1xn_n_1;
                  }

                  if ( !CurrentNode->isCond2D(CT_dYdy_NULL_2D) ) {
                      CurrentNode->droYdy[NUM_COMPONENTS]=(rhoY_air_Up - rhoY_air_Down)*dy_1xm_m_1;
                  }

                  if (CurrentNode->isCond2D(CT_WALL_NO_SLIP_2D) || CurrentNode->isCond2D(CT_WALL_LAW_2D) )  {
                      CurrentNode->dUdx=(RightNode->U*n1-LeftNode->U*n2)*dx_1xn_n_1;
                      CurrentNode->dVdx=(RightNode->V*n1-LeftNode->V*n2)*dx_1xn_n_1;

                      CurrentNode->dUdy=(UpNode->U*n3-DownNode->U*n4)*dy_1xm_m_1;
                      CurrentNode->dVdy=(UpNode->V*n3-DownNode->V*n4)*dy_1xm_m_1;

                      if(CurrentNode->isTurbulenceCond2D(TCT_k_eps_Model_2D)){
                        CurrentNode->dkdx   =(RightNode->S[i2d_k]*n1-LeftNode->S[i2d_k]*n2)*dx_1xn_n_1*rho_1;
                        CurrentNode->depsdx =(RightNode->S[i2d_eps]*n1-LeftNode->S[i2d_eps]*n2)*dx_1xn_n_1*rho_1;

                        CurrentNode->dkdy   =(UpNode->S[i2d_k]*n3-DownNode->S[i2d_k]*n4)*dy_1xm_m_1*rho_1;
                        CurrentNode->depsdy =(UpNode->S[i2d_eps]*n3-DownNode->S[i2d_eps]*n4)*dy_1xm_m_1*rho_1;
                      } else if (CurrentNode->isTurbulenceCond2D(TCT_Spalart_Allmaras_Model_2D)) {
                                 CurrentNode->dkdx   =(RightNode->S[i2d_k]*n1-LeftNode->S[i2d_k]*n2)*dx_1xn_n_1*rho_1;
                                 CurrentNode->dkdy   =(UpNode->S[i2d_k]*n3-DownNode->S[i2d_k]*n4)*dy_1xm_m_1*rho_1;
                      }
                  } else {
                      CurrentNode->dUdx   =(RightNode->U-LeftNode->U)*dx_1xn_n_1;
                      CurrentNode->dVdx   =(RightNode->V-LeftNode->V)*dx_1xn_n_1;

                      CurrentNode->dUdy   =(UpNode->U-DownNode->U)*dy_1xm_m_1;
                      CurrentNode->dVdy   =(UpNode->V-DownNode->V)*dy_1xm_m_1;
                      if(CurrentNode->isTurbulenceCond2D(TCT_k_eps_Model_2D)){
                        CurrentNode->dkdx   =(RightNode->S[i2d_k]-LeftNode->S[i2d_k])*dx_1xn_n_1*rho_1;
                        CurrentNode->depsdx =(RightNode->S[i2d_eps]-LeftNode->S[i2d_eps])*dx_1xn_n_1*rho_1;

                        CurrentNode->dkdy   =(UpNode->S[i2d_k]-DownNode->S[i2d_k])*dy_1xm_m_1*rho_1;
                        CurrentNode->depsdy =(UpNode->S[i2d_eps]-DownNode->S[i2d_eps])*dy_1xm_m_1*rho_1;
                      } else if (CurrentNode->isTurbulenceCond2D(TCT_Spalart_Allmaras_Model_2D)) {
                                 CurrentNode->dkdx   =(RightNode->S[i2d_k]-LeftNode->S[i2d_k])*dx_1xn_n_1*rho_1;
                                 CurrentNode->dkdy   =(UpNode->S[i2d_k]-DownNode->S[i2d_k])*dy_1xm_m_1*rho_1;
                      }
                  }

                  CurrentNode->dTdx=(RightNode->Tg-LeftNode->Tg)*dx_1xn_n_1;
                  CurrentNode->dTdy=(UpNode->Tg-DownNode->Tg)*dy_1xm_m_1;
              }

              cuda_CalcChemicalReactions(CurrentNode,(ChemicalReactionsModel)Chemical_reactions_model, (void*)(pCRMD),sm);
                
              CurrentNode->FillNode2D(sm,noTurbCond,SigW,SigF,TurbExtModel,delta_bl,dx,dy,_Hu,_isSrcAdd,sm,_FT);

              if( CurrentNode->Tg < 0. || isnan(CurrentNode->Tg) ) {
                  *dt_min_device = 0;  // Computational instability
              }  else {
                  register FP AAA          = sqrt(CurrentNode->k*CurrentNode->R*CurrentNode->Tg); 
                  register FP dt_min_local = CFL*min(1.0/(dx_1*(AAA+fabs(CurrentNode->U))),1.0/(dy_1*(AAA+fabs(CurrentNode->V))));
                  if (isLocalTimeStep) {
                      CurrentNode->dt_local = dt_min_local;
                      *dt_min_device = 1;
                  } else {
#if FP == double
#warning use FP64
                    atomicMin(dt_min_device,double2uint(int2float_scale*dt_min_local));
#else
#warning use FP32                  
                    atomicMin(dt_min_device,float2uint(int2float_scale*dt_min_local));
#endif
                  }
              }
         } else if (CurrentNode->CT == (ulong)(NT_FC_2D) &&
                    CurrentNode->ix >= l_limit &&
                    CurrentNode->ix <  r_limit ) {
#define _DEVICE_                    
                    CurrentNode->FillNode2D(1,0,SigW,SigF,TurbExtModel,delta_bl,dx,dy,_Hu,_isSrcAdd,sm,_FT);
#undef  _DEVICE_
         }
    }
}

__global__  void 
cuda_CalcHeatOnWallSources(FlowNode2D<FP,NUM_COMPONENTS>*  pJ,
                           unsigned long int index_limit,
                           int MAX_X, int MAX_Y,
                           unsigned long r_limit,
                           unsigned long l_limit,
                           FP dx,FP dy,FP _dt, int isLocalTimeStep) {
    
    register unsigned long int index = threadIdx.x + blockIdx.x * blockDim.x;

    if(index < index_limit) {

       FlowNode2D< FP,NUM_COMPONENTS >* CurrentNode=&pJ[index];

       if(CurrentNode->CT != (ulong)(CT_SOLID_2D | CT_NODE_IS_SET_2D) &&
          CurrentNode->ix >= l_limit && 
          CurrentNode->ix <  r_limit) { 

              register FP dx_local, dy_local;
              register FlowNode2D< FP,NUM_COMPONENTS >* UpNode    = CurrentNode->UpNode;
              register FlowNode2D< FP,NUM_COMPONENTS >* DownNode  = CurrentNode->DownNode;
              register FlowNode2D< FP,NUM_COMPONENTS >* RightNode = CurrentNode->RightNode;
              register FlowNode2D< FP,NUM_COMPONENTS >* LeftNode  = CurrentNode->LeftNode;

              if ( CurrentNode->isCond2D(CT_WALL_LAW_2D) || 
                   CurrentNode->isCond2D(CT_WALL_NO_SLIP_2D)) { 

                  register FP lam_eff;
                  register FP Q_conv;
                  register FP dt; 
                  register int num_near_nodes = 0;

                  lam_eff = Q_conv = 0;
                  
                  dx_local = dx;
                  dy_local = dy;
                  
                  if (isLocalTimeStep) {
                      dt = CurrentNode->dt_local;
                  } else {
                      dt = _dt;
                  }
                  
                  if (DownNode->isCond2D(CT_SOLID_2D) ) {

                      num_near_nodes = 1;
                      lam_eff = CurrentNode->lam+CurrentNode->lam_t; 

                      if (CurrentNode->UpNode) {
                          lam_eff +=CurrentNode->UpNode->lam + CurrentNode->UpNode->lam_t;
                          num_near_nodes++;
                      }

                      lam_eff = lam_eff/num_near_nodes;

                      if(Q_conv > 0.)
                         Q_conv = (Q_conv - lam_eff*(DownNode->Tg - CurrentNode->Tg)/dy_local)*0.5;
                      else
                         Q_conv = -lam_eff*(DownNode->Tg - CurrentNode->Tg)/dy_local;

                      CurrentNode->SrcAdd[i2d_RhoE] = -dt*Q_conv/dy;
                  }

                  if (UpNode->isCond2D(CT_SOLID_2D) ) {

                      num_near_nodes = 1;
                      lam_eff = CurrentNode->lam+CurrentNode->lam_t; 

                      if (CurrentNode->DownNode) {
                          lam_eff +=CurrentNode->DownNode->lam + CurrentNode->DownNode->lam_t;
                          num_near_nodes++;
                      }

                      lam_eff = lam_eff/num_near_nodes;

                      if(Q_conv > 0.)
                         Q_conv = (Q_conv - lam_eff*(UpNode->Tg - CurrentNode->Tg)/dy_local)*0.5;
                      else
                         Q_conv = -lam_eff*(UpNode->Tg - CurrentNode->Tg)/dy_local;

                      CurrentNode->SrcAdd[i2d_RhoE] = -dt*Q_conv/dy;
                  }

                  if (LeftNode->isCond2D(CT_SOLID_2D) ) {

                      num_near_nodes = 1;
                      lam_eff = CurrentNode->lam+CurrentNode->lam_t; 

                      if (CurrentNode->RightNode) {
                          lam_eff +=CurrentNode->RightNode->lam + CurrentNode->RightNode->lam_t;
                          num_near_nodes++;
                      }

                      lam_eff = lam_eff/num_near_nodes;

                      if(Q_conv > 0.)
                         Q_conv = (Q_conv - lam_eff*(LeftNode->Tg - CurrentNode->Tg)/dx_local)*0.5;
                      else
                         Q_conv = -lam_eff*(LeftNode->Tg - CurrentNode->Tg)/dx_local;

                      CurrentNode->SrcAdd[i2d_RhoE] = -dt*Q_conv/dx;
                  }

                  if (RightNode->isCond2D(CT_SOLID_2D) ) {

                      num_near_nodes = 1;
                      lam_eff = CurrentNode->lam+CurrentNode->lam_t; 

                      if (CurrentNode->LeftNode) {
                          lam_eff +=CurrentNode->LeftNode->lam + CurrentNode->LeftNode->lam_t;
                          num_near_nodes++;
                      }

                      lam_eff = lam_eff/num_near_nodes;

                      if(Q_conv > 0.)
                         Q_conv = (Q_conv - lam_eff*(RightNode->Tg - CurrentNode->Tg)/dx_local)*0.5;
                      else
                         Q_conv = -lam_eff*(RightNode->Tg - CurrentNode->Tg)/dx_local;

                      CurrentNode->SrcAdd[i2d_RhoE] = -dt*Q_conv/dx;
                  }
              }
          }
       }
}
#endif // _CUDA_


/*******************************************************************************
*   OpenHyperFLOW2D                                                            *
*                                                                              *
*   Version  1.0.0                                                             *
*   Copyright (C)  1995-2013 by Serge A. Suchkov                               *
*   Copyright policy: LGPL V3                                                  *
*                                                                              *
*  hf2d_start.cpp: OpenHyperFLOW2D solver init code....                        *
*                                                                              *
*  last update: 15/12/2013                                                     *
********************************************************************************/
#include "libDEEPS2D/deeps2d_core.hpp"

#ifdef _MPI
int rank;
int last_rank;
UArray< UMatrix2D< FlowNode2D<double,NUM_COMPONENTS> >* >* ArraySubmatrix  = NULL;
double  x0;
#endif // _MPI

double  dx;
double  dy;

UArray< double >*     WallNodesUw_2D = NULL;
int                   NumWallNodes;

int main( int argc, char **argv )
{
    const  float  ver=_VER;
    static char    inFile[256];
#ifdef _MPI
    FlowNode2D<double,NUM_COMPONENTS>* TmpMatrixPtr;
    int TmpMaxX;
    UMatrix2D< FlowNode2D<double,NUM_COMPONENTS> >*            TmpSubmatrix    = NULL;
    UMatrix2D< FlowNodeCore2D<double,NUM_COMPONENTS> >*        TmpCoreSubmatrix= NULL;
#endif // _MPI
    ostream*     o_stream = &cout;
#ifndef   _WIN32
#ifdef    _DEBUG_0
    __ExceptLib.SetSystemException(SIGINT, AT_HANDLED);
    __ExceptLib.SetSystemException(SIGFPE, AT_HANDLED);
    __ExceptLib.SetSystemException(SIGSEGV,AT_HANDLED);
    __ExceptLib.SetSystemException(SIGBUS, AT_HANDLED);
#endif //_DEBUG_0
#endif //_WIN32

#ifdef _DEBUG_0
#ifndef _WIN32
    ___try {
#endif  // _WIN32
#endif  // _DEBUG_0
        if (argc < 2) {
            printf("OpenHyperFLOW2D/DEEPS solver v %'.2f ",ver);
#ifndef _MPI
#ifndef _OPEN_MP
            printf(" (serial version)\n");
#else
            printf(" (parallel OpenMP version)\n");
#endif // OPEN_MP
#else
            printf("(parallel MPI version)\n");
#endif // _MPI
            printf("Copyright (C) 1995-2013 by Serge A. Suchkov\nCopyright policy: LGPL V3\nUsage: %s [{input_data_file}]\n",argv[0]);

            printf("\n\t* Density-based 2D-Navier-Stokes solver for ");
#ifdef _UNIFORM_MESH_
            printf("uniform cartesian mesh");
#else
            printf("non-uniform mesh");
#endif //_UNIFORM_MESH_
            printf("\n");
            exit(0);
        } else {
#ifdef _MPI
            MPI::Init(argc, argv);
rank      = MPI::COMM_WORLD.Get_rank();
#endif // _MPI
            sprintf(inFile,"%s",argv[1]);
            Data = new InputData(inFile,DS_FILE,o_stream,0,10
#ifdef _MPI
                                 ,rank
#endif // _MPI
                                 );
            if (Data->GetDataError()!=0) {
#ifdef _MPI
                if( rank == 0 ) {
#endif // _MPI
                    *o_stream << "\nInput data error.\n" ;
                    o_stream->flush();
#ifdef _MPI

                }
                Exit_OpenHyperFLOW2D();
#else
                exit(0);
#endif // _MPI
            }
        }
#ifdef _MPI
//------------------------- MPI version ------------------------------------
	    last_rank = MPI::COMM_WORLD.Get_size()-1;
            if(last_rank == 0) {
	       printf("At least 2 CPU should be used.\n Check mpiexec or mpirun parameters\n");
	       Exit_OpenHyperFLOW2D();
	    }
	        InitSharedData(Data,&chemical_reactions,rank);        // Init shared data
            MPI::COMM_WORLD.Barrier();
            if (rank == 0) {
               ArraySubmatrix = new UArray< UMatrix2D< FlowNode2D<double,NUM_COMPONENTS> >* >();
               InitDEEPS2D((void*)o_stream);                       // Init solver (only for rank=0)
               // Scan area for seek wall nodes      
               WallNodes = GetWallNodes((ofstream*)o_stream,J,Data->GetIntVal((char*)"isVerboseOutput")); 
#ifndef     _PARALLEL_ONLY
               *o_stream << "\nSerial calc min distance to wall..." << flush;
               SetMinDistanceToWall2D(J,WallNodes);
               *o_stream << "OK\n" << flush;
               Recalc_y_plus(J,WallNodes);                        // Calculate initial y+ value         
               SetInitBoundaryLayer(J,delta_bl);                  // Set Initial boundary layer profile
#else
               *o_stream << "\nParallel calc min distance to wall...OK" << endl;
#endif //   _PARALLEL_ONLY_
               NumWallNodes = WallNodes->GetNumElements();
               
               TmpMatrixPtr=J->GetMatrixPtr();
               int SubStartIndex, SubMaxX,
                            r_Overlap=0,
                            l_Overlap=0;
               *o_stream << "Allocate SubMatrix:\n";
               SubStartIndex = 0;
               
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
               TmpMatrixPtr = (FlowNode2D<double,NUM_COMPONENTS>*)((ulong)J->GetMatrixPtr()+(ulong)(sizeof(FlowNode2D<double,NUM_COMPONENTS>)*(SubStartIndex)*MaxY));                      
               
               x0 = SubStartIndex*FlowNode2D<double,NUM_COMPONENTS>::dx;
               
               *o_stream << "SubMatrix("<<i<<")[" << TmpMaxX << "x" << MaxY << "]  Size=" << (ulong)(sizeof(FlowNode2D<double,NUM_COMPONENTS>)*TmpMaxX*MaxY)/(1024*1024) << " Mb\n"; 
               TmpSubmatrix = new UMatrix2D< FlowNode2D<double,NUM_COMPONENTS> >(TmpMatrixPtr,TmpMaxX,MaxY);

#ifdef     _PARALLEL_ONLY_
               SetMinDistanceToWall2D(TmpSubmatrix,WallNodes,x0);
#endif //   _PARALLEL_ONLY_

               ArraySubmatrix->AddElement(&TmpSubmatrix);
               o_stream->flush();
               if(i>0) {
                   MPI::COMM_WORLD.Send(&MaxY,1,MPI::INT,i,tag_MaxY);
                   MPI::COMM_WORLD.Send(&TmpMaxX,1,MPI::INT,i,tag_MaxX);
#ifdef _IMPI_
                   LongMatrixSend(i, TmpSubmatrix->GetMatrixPtr(), TmpSubmatrix->GetMatrixSize());// Low Mem Send subdomain
#else
                   MPI::COMM_WORLD.Send(TmpSubmatrix->GetMatrixPtr(),
                                        TmpSubmatrix->GetMatrixSize(),
                                        MPI::BYTE,i,tag_Matrix); 
#endif // _IMPI_
               
                   MPI::COMM_WORLD.Send(&NumWallNodes,1,MPI::INT,i,tag_NumWallNodes);              // Send wall nodes array size
                   MPI::COMM_WORLD.Send(WallNodes->GetArrayPtr(),NumWallNodes*sizeof(XY<int>),     // Send wall nodes array
                                        MPI::BYTE,i,tag_WallNodesArray);
                   MPI::COMM_WORLD.Send(&x0,1,MPI::DOUBLE,i,tag_X0);                               // Send x0 for submatrix
               }
           }
           TmpMaxX      = ArraySubmatrix->GetElement(0)->GetX();
           MaxY         = ArraySubmatrix->GetElement(0)->GetY();
           TmpSubmatrix = ArraySubmatrix->GetElement(0);
      } else {
           MPI::COMM_WORLD.Recv(&MaxY,1,MPI::INT,0,tag_MaxY);
           MPI::COMM_WORLD.Recv(&TmpMaxX,1,MPI::INT,0,tag_MaxX);
           TmpSubmatrix = new UMatrix2D< FlowNode2D<double,NUM_COMPONENTS> >(TmpMaxX,MaxY);
#ifdef _IMPI_
             LongMatrixRecv(0,TmpSubmatrix->GetMatrixPtr(),TmpSubmatrix->GetMatrixSize());
#else
             MPI::COMM_WORLD.Recv(TmpSubmatrix->GetMatrixPtr(),
                                  TmpSubmatrix->GetMatrixSize(),
                                  MPI::BYTE,0,tag_Matrix);
#endif // _IMPI_
           MPI::COMM_WORLD.Recv(&NumWallNodes,1,MPI::INT,0,tag_NumWallNodes);                      // Recive wall nodes array size
           WallNodes = new UArray< XY<int> >(NumWallNodes,-1);                                     // Create wall nodes array
           MPI::COMM_WORLD.Recv(WallNodes->GetArrayPtr(),NumWallNodes*sizeof(XY<int>),             // Recive wall nodes array
                                MPI::BYTE,0,tag_WallNodesArray);
           WallNodesUw_2D = new UArray<double>(NumWallNodes,-1);                                   // Create friction velosity array
           MPI::COMM_WORLD.Recv(&x0,1,MPI::DOUBLE,0,tag_X0);                                       // Recive x0 for submatrix
       }
       MPI::COMM_WORLD.Barrier();
       TmpCoreSubmatrix = new UMatrix2D< FlowNodeCore2D<double,NUM_COMPONENTS> >(TmpMaxX,MaxY);
       MPI::COMM_WORLD.Barrier();
     if(rank == 0) {
           *o_stream << "\nStart computation...\n" << flush;
     }
       DEEPS2D_Run((ofstream*)o_stream, 
                   TmpSubmatrix,
                   TmpCoreSubmatrix,
                   rank, last_rank);
//------------------------- MPI version ------------------------------------
#else
//---------------------- OpenMP/Single thread version ----------------------
       InitSharedData(Data,&chemical_reactions);          // Init shared data
       InitDEEPS2D((void*)o_stream);                      // Init solver 
       // Scan area for seek wall nodes
       WallNodes = GetWallNodes((ofstream*)o_stream,J,Data->GetIntVal((char*)"isVerboseOutput")); 
       SetMinDistanceToWall2D(J,WallNodes);
       Recalc_y_plus(J,WallNodes);                        // Calculate initial y+ value
       SetInitBoundaryLayer(J,delta_bl);                  // Set Initial boundary layer profile

       UMatrix2D< FlowNode2D<double,NUM_COMPONENTS> >* TmpSubmatrix=NULL;
       UMatrix2D< FlowNodeCore2D<double,NUM_COMPONENTS> >* TmpCoreSubmatrix=NULL;
       FlowNode2D<double,NUM_COMPONENTS>* TmpMatrixPtr=J->GetMatrixPtr();
       int SubStartIndex, Overlap, SubMaxX;
       *o_stream << "Allocate SubMatrix:\n";
       for (unsigned int i=0;i<GlobalSubmatrix->GetNumElements();i++) {
            SubStartIndex = GlobalSubmatrix->GetElementPtr(i)->GetX();  
            SubMaxX = GlobalSubmatrix->GetElementPtr(i)->GetY();
            if(i == GlobalSubmatrix->GetNumElements()-1)
               Overlap = 0;
            else
               Overlap = 2;
            TmpMatrixPtr=(FlowNode2D<double,NUM_COMPONENTS>*)((ulong)J->GetMatrixPtr()+(ulong)(sizeof(FlowNode2D<double,NUM_COMPONENTS>)*SubStartIndex*MaxY));
            TmpSubmatrix     = new UMatrix2D< FlowNode2D<double,NUM_COMPONENTS> >(TmpMatrixPtr,(SubMaxX-SubStartIndex)+Overlap,MaxY);
            TmpCoreSubmatrix = new UMatrix2D< FlowNodeCore2D<double,NUM_COMPONENTS> >((SubMaxX-SubStartIndex)+Overlap,MaxY);
            SubmatrixArray->AddElement(&TmpSubmatrix);
            CoreSubmatrixArray->AddElement(&TmpCoreSubmatrix);
           }
       SubStartIndex = 0;
       for(unsigned int i=0;i<SubmatrixArray->GetNumElements();i++) {
          *o_stream << "SubMatrix(" << i << ")[" << SubmatrixArray->GetElement(i)->GetX() << "x" << \
          SubmatrixArray->GetElement(i)->GetY() << "] Size= " << SubmatrixArray->GetElement(i)->GetMatrixSize()/(1024*1024) << "+" \
          << CoreSubmatrixArray->GetElement(i)->GetMatrixSize()/(1024*1024) << " Mb\n"; 
           o_stream->flush();
          }
       *o_stream << "\nStart computation...\n" << flush;
       o_stream->flush();
       DEEPS2D_Run((ofstream*)o_stream);
//---------------------- OpenMP version ------------------------------------
#endif // _MPI

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
#ifdef _MPI
    if(rank==0)  {
#endif // _MPI
        *o_stream << "Computation stopped.\n";
        o_stream->flush();
#ifdef _MPI
    }
#endif // _MPI
    Exit_OpenHyperFLOW2D();
}


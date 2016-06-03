/*******************************************************************************
*   OpenHyperFLOW2D-CUDA                                                       *
*                                                                              *
*   Version  2.0.1                                                             *
*   Copyright (C)  1995-2016 by Serge A. Suchkov                               *
*   Copyright policy: LGPL V3                                                  *
*   http://github.com/sergeas67/openhyperflow2d                                *
*                                                                              *
*   hf2d_start.cpp: OpenHyperFLOW2D solver init code....                       *
*                                                                              *
*  last update: 14/04/2016                                                     *
********************************************************************************/
#include "libDEEPS2D/deeps2d_core.hpp"
int main( int argc, char **argv )
{
    ostream*     o_stream = &cout;
    const  FP     ver=_VER;

#ifdef    _DEBUG_0
    __ExceptLib.SetSystemException(SIGINT, AT_HANDLED);
    __ExceptLib.SetSystemException(SIGFPE, AT_HANDLED);
    __ExceptLib.SetSystemException(SIGSEGV,AT_HANDLED);
    __ExceptLib.SetSystemException(SIGBUS, AT_HANDLED);
#endif //_DEBUG_0
    if (argc < 2) {
        printf("OpenHyperFLOW2D/DEEPS/FP%u solver v%'.2f ",(unsigned int)(sizeof(FP)*8),ver);
        printf(" (parallel CUDA version)");
        printf("\nCopyright (C) 1995-2016 by Serge A. Suchkov\nCopyright policy: LGPL V3\n\nUsage: %s [{input_data_file}]\n",argv[0]);
        printf("\n\t* Density-based 2D-Navier-Stokes solver for ");
#ifdef _UNIFORM_MESH_
        printf("uniform cartesian mesh");
#else
        printf("non-uniform mesh");
#endif //_UNIFORM_MESH_

        printf("\n\n\tFlowNode2D size = %u bytes\n\n",(unsigned int)(sizeof(FlowNode2D<FP,NUM_COMPONENTS>)));
#ifdef  _P2P_ACCESS_
        printf("\tDirect P2P MultiGPU exchange\n\n");
#else
        printf("\tMultiGPU exchange from host memory\n\n");
#endif  //_P2P_ACCESS_

#ifdef  _DEVICE_MMAP_
        printf("\tMMap device to host for global time step (dt)\n\n");
#endif  //_DEVICE_MMAP_

#ifdef  _PARALLEL_RECALC_Y_PLUS_
        printf("\tParallel recalc y+\n\n");
#endif  //_PARALLEL_RECALC_Y_PLUS_

        Exit_OpenHyperFLOW2D();
    } else {
//------------------------- CUDA version ------------------------------------
             start_OpenHyperFLOW2D_GPU(argc,argv);
//------------------------- CUDA version ------------------------------------
             Exit_OpenHyperFLOW2D(0);
    }
 return 0;
}


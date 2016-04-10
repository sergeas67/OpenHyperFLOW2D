/*******************************************************************************
*   OpenHyperFLOW2D                                                            *
*                                                                              *
*   Version  1.0.3                                                             *
*   Copyright (C)  1995-2016 by Serge A. Suchkov                               *
*   Copyright policy: LGPL V3                                                  *
*   http://openhyperflow2d.googlecode.com                                      *
*                                                                              *
*   last update: 04/07/2016                                                    *
*******************************************************************************/
#ifdef _MPI
#include <mpi.h>
#endif // _MPI
#include "libExcept/except.hpp"
#include "libOpenHyperFLOW2D/hyper_flow_area.hpp"
// <------------- 2D --------------->
void Abort_OpenHyperFLOW2D() {
#ifdef _MPI
    printf("\nTask[%d] aborted.\n",MPI::COMM_WORLD.Get_rank());
    MPI::COMM_WORLD.Abort(0);
    MPI::Finalize();
#endif // _MPI
    exit(0);
}

void Exit_OpenHyperFLOW2D() {
#ifdef _MPI
    printf("\nTask[%d] stopped.\n",MPI::COMM_WORLD.Get_rank());
    MPI::COMM_WORLD.Abort(0);
    MPI::Finalize();
#endif // _MPI
    exit(0);
};
// Area constructor
Area2D::Area2D(char* name, UMatrix2D< FlowNode2D<FP,NUM_COMPONENTS> >* J):AreaName(name),pMFN(J) {
    ;
}

// Fill area by Flow
void Area2D::FillArea2D(unsigned int X,
                        unsigned int Y,
                        ulong     bnt,
                        Flow*   pf,
                        FP* p_Y,
                        ulong     att,
                        int MaterialID) {
    Flow2D F2D(*pf);
    FillArea2D((unsigned int)X,(unsigned int)Y,
               bnt | CT_NODE_IS_SET_2D,&F2D,p_Y,att,MaterialID);
}
         
// Fill area by custom type nodes
void Area2D::FillArea2D(unsigned int X,unsigned int Y, 
                        ulong bnt,ulong att,
                        int MaterialID) {
    FillArea2D((unsigned int)X,(unsigned int)Y,bnt | CT_NODE_IS_SET_2D,(Flow2D*)NULL,NULL,att,MaterialID);
}

// Fill area by custom type nodes
void Area2D::FillArea2D(FP X, FP Y, ulong bnt,ulong att,
                        int MaterialID) {
    FillArea2D((unsigned int)(X/pMFN->GetValue(0,0).dx),(unsigned int)(Y/pMFN->GetValue(0,0).dy),bnt | CT_NODE_IS_SET_2D,(Flow2D*)NULL,NULL,att,MaterialID);
}

// Fill area by Flow2D
void Area2D::FillArea2D(unsigned int X,
                        unsigned int Y,
                        ulong    bnt,
                        Flow2D*  pf2d,
                        FP*  p_Y,
                        ulong    att,
                        int MaterialID) {

    static unsigned int      ii,i,XMax,YMax,tX,tY;
    static XY<unsigned int>  TmpXY;
    ulong                    TestNT;

    static UArray<XY <unsigned int> > FNA;
    static UArray<XY <unsigned int> > BNA;

    StartX  =  X;
    StartY  =  Y;

    YMax    =  pMFN->GetY();
    XMax    =  pMFN->GetX();

    if ( XMax>X && YMax>Y ) {
        TestNT = pMFN->GetValue(X,Y).CT;
        if (!pMFN->GetValue(X,Y).isCond2D(CT_NODE_IS_SET_2D)) {
            ANT     = bnt | CT_NODE_IS_SET_2D;
            ATT     = att;
            TmpXY.SetX(X);
            TmpXY.SetY(Y);
            pMFN->GetValue(X,Y).CT = ANT;
            pMFN->GetValue(X,Y).TurbType = ATT;
            BNA.AddElement(&TmpXY);
            do {
                for (i=0;i<BNA.GetNumElements();i++) {
                    TmpXY = BNA.GetElement(i);

                    tX = TmpXY.GetX();
                    tY = TmpXY.GetY();

                    if (NUM_COMPONENTS  && p_Y)
                        for (ii=0;ii<NUM_COMPONENTS+1;ii++)
                            (pMFN->GetValue(tX,tY)).Y[ii]=p_Y[ii];

                    if (pf2d != NULL)
                        pMFN->GetValue(tX,tY) = *pf2d;

                    pMFN->GetValue(tX,tY).BGX =1.;
                    pMFN->GetValue(tX,tY).BGY =1.;
                    pMFN->GetValue(tX,tY).NGX =1;
                    pMFN->GetValue(tX,tY).NGY =1;
                    pMFN->GetValue(tX,tY).idXl=1;
                    pMFN->GetValue(tX,tY).idYu=1;
                    pMFN->GetValue(tX,tY).idXr=1;
                    pMFN->GetValue(tX,tY).idYd=1;

                    if (tX > 0) {
                        if (!pMFN->GetValue(tX-1,tY).isCond2D(CT_NODE_IS_SET_2D)) {
                            TmpXY.SetX(tX-1);
                            TmpXY.SetY(tY);
                            pMFN->GetValue(TmpXY.GetX(),TmpXY.GetY()).CT = ANT;
                            pMFN->GetValue(TmpXY.GetX(),TmpXY.GetY()).TurbType = ATT;
                            FNA.AddElement(&TmpXY);
                        } else if (!pMFN->GetValue(tX-1,tY).isCond2D(CT_SOLID_2D) &&
                                   pMFN->GetValue(tX,tY).isCond2D(CT_SOLID_2D)) {
                            pMFN->GetValue(tX-1,tY).NGX = 0;
                            pMFN->GetValue(tX-1,tY).idXr= 0;
                        }
                    }

                    if (tX < XMax-1) {
                        if (!pMFN->GetValue(tX+1,tY).isCond2D(CT_NODE_IS_SET_2D)) {
                            TmpXY.SetX(tX+1);
                            TmpXY.SetY(tY);
                            pMFN->GetValue(TmpXY.GetX(),TmpXY.GetY()).CT = ANT;
                            pMFN->GetValue(TmpXY.GetX(),TmpXY.GetY()).TurbType = ATT; 
                            FNA.AddElement(&TmpXY);
                        } else if (!pMFN->GetValue(tX+1,tY).isCond2D(CT_SOLID_2D) &&
                                   pMFN->GetValue(tX,tY).isCond2D(CT_SOLID_2D)) {
                            pMFN->GetValue(tX+1,tY).NGX = 0;
                            pMFN->GetValue(tX+1,tY).idXl= 0;
                        }
                    }
                    if (tY > 0) {
                        if (!pMFN->GetValue(tX,tY-1).isCond2D(CT_NODE_IS_SET_2D)) {
                            TmpXY.SetX(tX);
                            TmpXY.SetY(tY-1);
                            pMFN->GetValue(TmpXY.GetX(),TmpXY.GetY()).CT = ANT;
                            pMFN->GetValue(TmpXY.GetX(),TmpXY.GetY()).TurbType = ATT; 
                            FNA.AddElement(&TmpXY);
                        } else if (!pMFN->GetValue(tX,tY-1).isCond2D(CT_SOLID_2D) &&
                                   pMFN->GetValue(tX,tY).isCond2D(CT_SOLID_2D)) {
                            pMFN->GetValue(tX,tY-1).NGY = 0;
                            pMFN->GetValue(tX,tY-1).idYu= 0;
                        }
                    }
                    if (tY < YMax-1) {
                        if (!pMFN->GetValue(tX,tY+1).isCond2D(CT_NODE_IS_SET_2D)) {
                            TmpXY.SetX(tX);
                            TmpXY.SetY(tY+1);
                            pMFN->GetValue(TmpXY.GetX(),TmpXY.GetY()).CT = ANT;
                            pMFN->GetValue(TmpXY.GetX(),TmpXY.GetY()).TurbType = ATT; 
                            FNA.AddElement(&TmpXY);
                        } else if (!pMFN->GetValue(tX,tY+1).isCond2D(CT_SOLID_2D) &&
                                   pMFN->GetValue(tX,tY).isCond2D(CT_SOLID_2D)) {
                            pMFN->GetValue(tX,tY+1).NGY = 0;
                            pMFN->GetValue(tX,tY+1).idYd= 0;
                        }
                    }
                  //pMFN->GetValue(tX,tY).NodeID = MaterialID;
                  pMFN->GetValue(tX,tY).FillNode2D(1);
                }
                BNA = FNA;
                FNA.CleanArray();
            }while (BNA.GetNumElements()>0);
        } else {
            as = AS_ERR_INIT_POINT; throw(this);
        }
    } else {
        as = AS_ERR_OUT_OF_RANGE;throw(this);
    }
    as = AS_OK;
}

void Area2D::FillArea2D(FP  X,
                        FP  Y,
                        ulong   bnt,
                        Flow2D*  pf2d,
                        FP*  p_Y,
                        ulong    att,
                        int MaterialID) {

  FillArea2D((unsigned int)(X/pMFN->GetValue(0,0).dx),
             (unsigned int)(Y/pMFN->GetValue(0,0).dy),
             bnt,
             pf2d,
             p_Y,
             att,
             MaterialID);
}

//Area destructor
Area2D::~Area2D() {
    ;
}
// <------------- 2D --------------->           


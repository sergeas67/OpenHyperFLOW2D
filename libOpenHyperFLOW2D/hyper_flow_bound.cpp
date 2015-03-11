/*******************************************************************************
*   OpenHyperFLOW2D                                                            *
*                                                                              *
*   Version  1.0.2                                                             *
*   Copyright (C)  1995-2015 by Serge A. Suchkov                               *
*   Copyright policy: LGPL V3                                                  *
*   http://openhyperflow2d.googlecode.com                                      *
*                                                                              *
*   last update: 01/02/2015                                                    *
*******************************************************************************/
#include "libExcept/except.hpp"
#include "libOpenHyperFLOW2D/hyper_flow_bound.hpp"

// <------------- 2D --------------->           

typedef FlowNode2D< FP, NUM_COMPONENTS>*  NODE2D ; 

// Init bound 
void Bound2D::InitBound(char* name,  
                        UMatrix2D< FlowNode2D< FP, NUM_COMPONENTS> >* J,
                        int     bt,
                        FP* p_Y,
                        int     btc) {
    BoundName = name;
    bs        = BND_INACTIVE;
    BNT       = bt;
    BTC       = btc;
    pMFN      = J;
    if (NUM_COMPONENTS != 0 && p_Y != 0) {
        pY = new FP[NUM_COMPONENTS+1];
        for (unsigned int i=0;i<NUM_COMPONENTS+1;i++)
            pY[i]=p_Y[i];          // Y[a] - Y1,Y2,Y3,...,Ya
    }

}
// Bound2D constructor (Flow,uint,uint,uint,uint)
Bound2D::Bound2D(char* name,   
                 UMatrix2D< FlowNode2D<FP, NUM_COMPONENTS> >* J,
                 unsigned int X1, unsigned int Y1,
                 unsigned int X2, unsigned int Y2,
                 int      bt,
                 Flow*    pInFlow,
                 FP*  p_Y,
                 int     btc) {
    
    pBoundFlow   = pInFlow;
    pBoundFlow2D = NULL; 

    InitBound(name,J,bt,p_Y,btc);
    
    Start.SetX(X1);
    Start.SetY(Y1);

    End.SetX(X2);
    End.SetY(Y2);

    fStart.SetX((FP)X1);
    fStart.SetY((FP)Y1);

    fEnd.SetX((FP)X2);
    fEnd.SetY((FP)Y2);
}


// Bound2D constructor (Flow2D,uint,uint,uint,uint)
Bound2D::Bound2D(char* name,   
                 UMatrix2D< FlowNode2D<FP, NUM_COMPONENTS> >* J,
                 unsigned int X1, unsigned int Y1,
                 unsigned int X2, unsigned int Y2,
                 int      bt,
                 Flow2D*  pInFlow2D,
                 FP*  p_Y,
                 int     btc) {
    pBoundFlow2D = pInFlow2D;

    pBoundFlow   = NULL;

    InitBound(name,J,bt,p_Y,btc);
    
    Start.SetX(X1);
    Start.SetY(Y1);

    End.SetX(X2);
    End.SetY(Y2);

    fStart.SetX((FP)X1);
    fStart.SetY((FP)Y1);

    fEnd.SetX((FP)X2);
    fEnd.SetY((FP)Y2);
}

// Bound2D constructor (Flow,XY<uint,uint>,XY<uint,uint>)
Bound2D::Bound2D(char* name,   
                 UMatrix2D< FlowNode2D<FP, NUM_COMPONENTS> >* J,
                 XY<unsigned int>*  p_start,
                 XY<unsigned int>*  p_end,
                 int                bt,
                 Flow*              pInFlow,
                 FP*            p_Y,
                 int                btc) {
    
    pBoundFlow   = pInFlow;

    pBoundFlow2D = NULL; 

    InitBound(name,J,bt,p_Y,btc);
    
    Start.SetXY(p_start);
    End.SetXY(p_end);

    fStart.SetX(static_cast<FP>(Start.GetX()));
    fStart.SetY(static_cast<FP>(Start.GetY()));

    fEnd.SetX(static_cast<FP>(End.GetX()));
    fEnd.SetY(static_cast<FP>(End.GetY()));
}

// Bound2D constructor (Flow2D,XY<uint,uint>,XY<uint,uint>)
Bound2D::Bound2D(char*  name,   
                 UMatrix2D< FlowNode2D<FP, NUM_COMPONENTS> >* J,
                 XY<unsigned int>*  p_start,
                 XY<unsigned int>*  p_end,
                 int                bt,
                 Flow2D*            pInFlow2D,
                 FP*            p_Y,
                 int                btc) {

    pBoundFlow2D = pInFlow2D;

    pBoundFlow   = NULL;

    InitBound(name,J,bt,p_Y,btc);

    Start.SetXY(p_start);
    End.SetXY(p_end);

    fStart.SetX(static_cast<FP>(Start.GetX()));
    fStart.SetY(static_cast<FP>(Start.GetY()));

    fEnd.SetX(static_cast<FP>(End.GetX()));
    fEnd.SetY(static_cast<FP>(End.GetY()));
}

// Bound constructor (Flow,FP,FP,FP,FP)
Bound2D::Bound2D(char* name,   
                 UMatrix2D< FlowNode2D<FP, NUM_COMPONENTS> >* J,
                 FP  x1,
                 FP  y1,
                 FP  x2,
                 FP  y2,
                 int     bt,
                 Flow*   pInFlow,
                 FP* p_Y,
                 int     btc) {

    pBoundFlow   = pInFlow;

    pBoundFlow2D = NULL; 

    InitBound(name, J,
              bt,p_Y,btc);

    Start.SetX(static_cast<unsigned int>(x1/dx));
    Start.SetY(static_cast<unsigned int>(y1/dy));

    End.SetX(static_cast<unsigned int>(x2/dx));
    End.SetY(static_cast<unsigned int>(y2/dy));
    
    fStart.SetX(x1);
    fStart.SetY(y1);

    fEnd.SetX(x2);
    fEnd.SetY(y2);
}

// Bound constructor (Flow2D,FP,FP,FP,FP)
Bound2D::Bound2D(char* name,   
                 UMatrix2D< FlowNode2D<FP, NUM_COMPONENTS> >* J,
                 FP  x1,
                 FP  y1,
                 FP  x2,
                 FP  y2,
                 int     bt,
                 Flow2D*           pInFlow2D,
                 FP* p_Y,
                 int     btc) {

    pBoundFlow2D = pInFlow2D;

    pBoundFlow   = NULL;

    InitBound(name, J, bt,p_Y,btc);
    
    Start.SetX(static_cast<unsigned int>(x1/dx));
    Start.SetY(static_cast<unsigned int>(y1/dy));

    End.SetX(static_cast<unsigned int>(x2/dx));
    End.SetY(static_cast<unsigned int>(y2/dy));
    
    fStart.SetX(x1);
    fStart.SetY(y1);

    fEnd.SetX(x2);
    fEnd.SetY(y2);
}

// Bound2D constructor (Flow,XY<uint,uint>,XY<uint,uint>)
Bound2D::Bound2D(char* name,   
                 UMatrix2D< FlowNode2D<FP, NUM_COMPONENTS> >* J,
                 XY<FP>*  p_start,
                 XY<FP>*  p_end,
                 int                bt,
                 Flow*              pInFlow,
                 FP*            p_Y,
                 int                btc) {

    pBoundFlow   = pInFlow;

    pBoundFlow2D = NULL; 

    InitBound(name, J, bt,p_Y,btc);
    fStart.SetXY(p_start);
    fEnd.SetXY(p_end);
    
    Start.SetX(static_cast<unsigned int>(fStart.GetX()));
    Start.SetY(static_cast<unsigned int>(fStart.GetY()));

    End.SetX(static_cast<unsigned int>(fEnd.GetX()));
    End.SetY(static_cast<unsigned int>(fEnd.GetY()));
}

// Bound2D constructor (Flow2D,XY<uint,uint>,XY<uint,uint>)
Bound2D::Bound2D(char* name,   
                 UMatrix2D< FlowNode2D<FP, NUM_COMPONENTS> >* J,
                 XY<FP>*  p_start,
                 XY<FP>*  p_end,
                 int      bt,
                 Flow2D*  pInFlow2D,
                 FP*      p_Y,
                 int      btc) {

    pBoundFlow2D = pInFlow2D;

    pBoundFlow   = NULL;

    InitBound(name, J, bt,p_Y,btc);
    fStart.SetXY(p_start);
    fEnd.SetXY(p_end);
    Start.SetX(static_cast<unsigned int>(fStart.GetX()));
    Start.SetY(static_cast<unsigned int>(fStart.GetY()));

    End.SetX(static_cast<unsigned int>(fEnd.GetX()));
    End.SetY(static_cast<unsigned int>(fEnd.GetY()));
}

// Set bound
BoundState Bound2D::SetBound(int MaterialID) {
    static  FP  Alpha,K;
    static  unsigned int   ii,i,j1,j2,k1,j=0;
    static  FP  DX,DY;
    XY<FP>      fCurrentPoint;

    if (Start.GetX() > pMFN->GetX()||
        Start.GetY() > pMFN->GetY()||
        End.GetX()   > pMFN->GetX()||
        End.GetY()   > pMFN->GetY()) {
        bs = BND_ERR;
        return bs;
    }

    if (Start.GetX() == pMFN->GetX()) Start.SetX(pMFN->GetX()-1);
    if (Start.GetY() == pMFN->GetY()) Start.SetY(pMFN->GetY()-1);
    if (End.GetX()   == pMFN->GetX()) End.SetX(pMFN->GetX()-1);
    if (End.GetY()   == pMFN->GetY()) End.SetY(pMFN->GetY()-1);

    DX=fStart.GetX()-fEnd.GetX();
    DY=fStart.GetY()-fEnd.GetY();

    if (DX != 0) {
        K=DY/DX;Alpha = atan(K);
    } else {
        Alpha = pi/2.;K=0;
    }

    if (fabs(DX)>fabs(DY)) {
        j1 = min(Start.GetX(),End.GetX());

        if (j1 == Start.GetX()) k1 = Start.GetY();
        else                    k1 = End.GetY();
        
        j2 = max(Start.GetX(),End.GetX());
        
        for (i=j1;i<=j2;i++) {
            j = k1 +(int)((FP)(i-j1)*tan(Alpha));
            pMFN->GetValue(i,j).CT = pMFN->GetValue(i,j).CT | BNT | CT_NODE_IS_SET_2D;
            pMFN->GetValue(i,j).TurbType  = BTC;
            pMFN->GetValue(i,j).NGX = (3-pMFN->GetValue(i,j).idXr-pMFN->GetValue(i,j).idXl);
            pMFN->GetValue(i,j).NGY = (3-pMFN->GetValue(i,j).idYu-pMFN->GetValue(i,j).idYd);
            pMFN->GetValue(i,j).BGX = cos(Alpha);
            pMFN->GetValue(i,j).BGY = sin(Alpha);
            if (NUM_COMPONENTS != 0  && pY != 0)
                for (ii=0;ii<NUM_COMPONENTS+1;ii++)
                    (pMFN->GetValue(i,j)).Y[ii]=pY[ii];

            //pMFN->GetValue(i,j).NodeID   = MaterialID;

            if (pBoundFlow)
                pMFN->GetValue(i,j)     = *pBoundFlow;
            else if (pBoundFlow2D)
                pMFN->GetValue(i,j)     = *pBoundFlow2D;
            else
                pMFN->GetValue(i,j).FillNode2D();
        }
    } else {
        j1 = min(Start.GetY(),End.GetY());

        if (j1 == Start.GetY()) k1 = Start.GetX();
        else              k1 = End.GetX();
        
        j2 = max(Start.GetY(),End.GetY());
        
        for (i=j1;i<=j2;i++) {

            if (tan(Alpha) != 0.)  j = k1 + (int)((FP)(i-j1)/tan(Alpha));
            else                   j = k1;
            
            pMFN->GetValue(j,i).CT  = pMFN->GetValue(j,i).CT | BNT | CT_NODE_IS_SET_2D;
            pMFN->GetValue(j,i).TurbType  = BTC;
            pMFN->GetValue(j,i).NGX = 3-pMFN->GetValue(j,i).idXr-pMFN->GetValue(j,i).idXl;
            pMFN->GetValue(j,i).NGY = 3-pMFN->GetValue(j,i).idYd-pMFN->GetValue(j,i).idYu;
            pMFN->GetValue(j,i).BGX = cos(Alpha);
            pMFN->GetValue(j,i).BGY = sin(Alpha);
            
            if (NUM_COMPONENTS != 0  && pY != 0)
                for (ii=0;ii<NUM_COMPONENTS+1;ii++)
                    (pMFN->GetValue(j,i)).Y[ii]=pY[ii];

            //pMFN->GetValue(j,i).NodeID   = MaterialID;

            if (pBoundFlow)
                pMFN->GetValue(j,i)     = *pBoundFlow;
            else if (pBoundFlow2D)
                pMFN->GetValue(j,i)     = *pBoundFlow2D;
            else
                pMFN->GetValue(j,i).FillNode2D();
        }
    }
    bs = BND_OK;
    return bs;
}
// Set bound. Push result to array
BoundState Bound2D::SetBound(UArray< FlowNode2D<FP, NUM_COMPONENTS>* >* node_array, int MaterialID) {
    static  FP  Alpha,K;
    static  unsigned int   ii,i,j1,j2,k1,j=0;
    static  FP  DX,DY;
    XY<FP>      fCurrentPoint;
    FlowNode2D<FP, NUM_COMPONENTS>* TmpNodePtr;

    if(node_array==NULL) {
        bs = BND_ERR;
        return bs;
    }

    if (Start.GetX() > pMFN->GetX()||
        Start.GetY() > pMFN->GetY()||
        End.GetX()   > pMFN->GetX()||
        End.GetY()   > pMFN->GetY()) {
        bs = BND_ERR;
        return bs;
    }

    if (Start.GetX() == pMFN->GetX()) Start.SetX(pMFN->GetX()-1);
    if (Start.GetY() == pMFN->GetY()) Start.SetY(pMFN->GetY()-1);
    if (End.GetX()   == pMFN->GetX()) End.SetX(pMFN->GetX()-1);
    if (End.GetY()   == pMFN->GetY()) End.SetY(pMFN->GetY()-1);

    DX=fStart.GetX()-fEnd.GetX();
    DY=fStart.GetY()-fEnd.GetY();

    if (DX != 0) {
        K=DY/DX;Alpha = atan(K);
    } else {
        Alpha = pi/2.;K=0;
    }

    if (fabs(DX)>fabs(DY)) {
        j1 = min(Start.GetX(),End.GetX());
        if (j1 == Start.GetX()) k1 = Start.GetY();
        else              k1 = End.GetY();
        j2 = max(Start.GetX(),End.GetX());
        for (i=j1;i<=j2;i++) {
            j = k1 +(int)((FP)(i-j1)*tan(Alpha));
            pMFN->GetValue(i,j).CT  = pMFN->GetValue(i,j).CT |  BNT | CT_NODE_IS_SET_2D;
            pMFN->GetValue(i,j).TurbType  = BTC;
            pMFN->GetValue(i,j).NGX = (3-pMFN->GetValue(i,j).idXr-pMFN->GetValue(i,j).idXl);
            pMFN->GetValue(i,j).NGY = (3-pMFN->GetValue(i,j).idYu-pMFN->GetValue(i,j).idYd);
            pMFN->GetValue(i,j).BGX = cos(Alpha);
            pMFN->GetValue(i,j).BGY = sin(Alpha);
            if (NUM_COMPONENTS != 0  && pY != 0)
                for (ii=0;ii<NUM_COMPONENTS+1;ii++)
                    (pMFN->GetValue(i,j)).Y[ii]=pY[ii];

            //pMFN->GetValue(i,j).NodeID   = MaterialID;

            if (pBoundFlow)
                pMFN->GetValue(i,j)     = *pBoundFlow;
            else if (pBoundFlow2D)
                pMFN->GetValue(i,j)     = *pBoundFlow2D;
            else
                pMFN->GetValue(i,j).FillNode2D();

            TmpNodePtr = &(pMFN->GetValue(i,j));  

            node_array->AddElement(&TmpNodePtr);
        }
    } else {
        j1 = min(Start.GetY(),End.GetY());
        if (j1 == Start.GetY()) k1 = Start.GetX();
        else              k1 = End.GetX();
        j2 = max(Start.GetY(),End.GetY());
        for (i=j1;i<=j2;i++) {
            if (tan(Alpha) != 0.)  j = k1 + (int)((FP)(i-j1)/tan(Alpha));
            else                   j = k1;
            pMFN->GetValue(j,i).CT  = pMFN->GetValue(j,i).CT | BNT | CT_NODE_IS_SET_2D;
            pMFN->GetValue(j,i).TurbType  = BTC;
            pMFN->GetValue(j,i).NGX = 3-pMFN->GetValue(j,i).idXr-pMFN->GetValue(j,i).idXl;
            pMFN->GetValue(j,i).NGY = 3-pMFN->GetValue(j,i).idYd-pMFN->GetValue(j,i).idYu;
            pMFN->GetValue(j,i).BGX = cos(Alpha);
            pMFN->GetValue(j,i).BGY = sin(Alpha);
            if (NUM_COMPONENTS != 0  && pY != 0)
                for (ii=0;ii<NUM_COMPONENTS+1;ii++)
                    (pMFN->GetValue(j,i)).Y[ii]=pY[ii];

            //pMFN->GetValue(i,j).NodeID   = MaterialID;

            if (pBoundFlow)
                pMFN->GetValue(j,i)     = *pBoundFlow;
            else if (pBoundFlow2D)
                pMFN->GetValue(j,i)     = *pBoundFlow2D;
            else
                pMFN->GetValue(j,i).FillNode2D();
            
            TmpNodePtr = &(pMFN->GetValue(j,i));  
            
            node_array->AddElement(&TmpNodePtr);
        }
    }
    bs = BND_OK;
    return bs;
}

// Bound destructor
Bound2D::~Bound2D() {
    if (pY != 0)
        delete[] pY;
}

FP* Bound2D::GetYArray() {
    return pY;
}

// Get num components
unsigned int Bound2D::GetNumComponents() {
    return NUM_COMPONENTS+1;
}
// Get bound state
BoundState   Bound2D::GetBoundState() {
    return bs;
}

// Get start X bound coord (nodes)
unsigned int Bound2D::GetStartX() {
    return Start.GetX();
}
// Get start Y bound coord (nodes)
unsigned int Bound2D::GetStartY() {
    return Start.GetY();
}
// Get end X bound coord (nodes)
unsigned int Bound2D::GetEndX() {
    return End.GetX();
}
// Get end Y bound coord (nodes)
unsigned int Bound2D::GetEndY() {
    return End.GetY();
}
// Set start X bound coord (nodes)
void  Bound2D::SetStartX(unsigned int x) {
    Start.SetX(x);
    fStart.SetX(static_cast<FP>(x*dx));
}
// Set start Y bound coord (nodes)
void  Bound2D::SetStartY(unsigned int y) {
    Start.SetY(y);
    fStart.SetY(static_cast<FP>(y*dy));
}
// Set end X bound coord (nodes)
void  Bound2D::SetEndX(unsigned int x) {
    End.SetX(x);
    fEnd.SetX(static_cast<FP>(x*dx));
}
// Set end Y bound coord (nodes)
void   Bound2D::SetEndY(unsigned int y) {
    End.SetY(y);
    fEnd.SetY(static_cast<FP>(y*dy));
}

// Get start X bound coord (m)
FP  Bound2D::GetStartFX() {
    return fStart.GetY();
}
// Get start Y bound coord (m)
FP  Bound2D::GetStartFY() {
    return fStart.GetY();
}
// Get end X bound coord (m)
FP  Bound2D::GetEndFX() {
    return fEnd.GetX();
}
// Get end Y bound coord (m)
FP  Bound2D::GetEndFY() {
    return fEnd.GetY();
}
// Set start X bound coord (m)
void    Bound2D::SetStartFX(FP x) {
    fStart.SetX(x);
    Start.SetX(static_cast<unsigned int>(x/dx));
}

// Set start Y bound coord (m)
void    Bound2D::SetStartFY(FP y) {
    fStart.SetY(y);
    Start.SetY(static_cast<unsigned int>(y/dy));
}

// Set end X bound coord (m)
void    Bound2D::SetEndFX(FP x) {
    fEnd.SetX(x);
    End.SetX(static_cast<unsigned int>(x/dx));
}

// Set end Y bound coord (m)
void    Bound2D::SetEndFY(FP y) {
    fEnd.SetY(y);
    End.SetY(static_cast<unsigned int>(y/dy));
}

// Set start XY bound coord (nodes)
void  Bound2D::SetStartXY(XY<unsigned int>* xy) {
    Start.SetXY(xy);
    fStart.SetX(static_cast<FP>(Start.GetX()));
    fStart.SetY(static_cast<FP>(Start.GetY()));
}

// Set start XY bound coord (m)
void  Bound2D::SetStartFXY(XY<FP>* xy) {
    fStart.SetXY(xy);
    Start.SetX(static_cast<unsigned int>(fStart.GetX()));
    Start.SetY(static_cast<unsigned int>(fStart.GetY()));
}

// Set end XY bound coord (nodes)
void   Bound2D::SetEndXY(XY<unsigned int>* xy) {
    End.SetXY(xy);
    fEnd.SetX(static_cast<FP>(End.GetX()));
    fEnd.SetY(static_cast<FP>(End.GetY()));
}

// Set end XY bound coord (m)
void  Bound2D::SetEndFXY(XY<FP>* xy) {
    fEnd.SetXY(xy);
    End.SetX(static_cast<unsigned int>(fEnd.GetX()));
    End.SetY(static_cast<unsigned int>(fEnd.GetY()));
}


// Rotate bound
// Bound can be rotated _only_ in INACTIVE state (before SetBound() call )!!!
int Bound2D::RotateBound2D(FP x0, FP y0, FP angle) {
    if (bs != BND_INACTIVE) return 0;
    FP dxs = fStart.GetX()-x0,
       dys = fStart.GetY()-y0,
       fi1 = atan2(dxs,dys),
       dxe = fEnd.GetX()-x0,
       dye = fEnd.GetY()-y0,
       fi2 = atan2(dxe,dye),
       r1  = sqrt(dxs*dxs+dys*dys+1.e-30),
       r2  = sqrt(dxe*dxe+dye*dye+1.e-30),
       xsn,ysn,xen,yen;

    xsn = x0 + (r1*sin(fi1+angle)); 
    ysn = y0 + (r1*cos(fi1+angle));
    xen = x0 + (r2*sin(fi2+angle)); 
    yen = y0 + (r2*cos(fi2+angle)); 

    if (xsn < 0. || ysn < 0. || xen < 0.|| yen < 0.) return 0;
    
    if (pMFN->GetX()<static_cast<unsigned int>(xsn/dx))    return 0;
    if (pMFN->GetY()<static_cast<unsigned int>(ysn/dy))    return 0;
    if (pMFN->GetX()<static_cast<unsigned int>(xen/dx))    return 0;
    if (pMFN->GetY()<static_cast<unsigned int>(yen/dy))    return 0;
    
    SetStartFX(xsn);
    SetStartFY(ysn);

    SetEndFX(xen);
    SetEndFY(yen);

    return 1;
}
// Test bound rotation
int Bound2D::TestRotateBound2D(FP x0, FP y0, FP angle) {
    if (bs != BND_INACTIVE) return 0;
    FP dxs = fStart.GetX()-x0,
           dys = fStart.GetY()-y0,
           fi1 = atan2(dxs,dys),
           dxe = fEnd.GetX()-x0,
           dye = fEnd.GetY()-y0,
           fi2 = atan2(dxe,dye),
           r1  = sqrt(dxs*dxs+dys*dys+1.e-30),
           r2  = sqrt(dxe*dxe+dye*dye+1.e-30),
           xsn,ysn,xen,yen;

    xsn = x0 + (r1*sin(fi1+angle)); 
    ysn = y0 + (r1*cos(fi1+angle));
    xen = x0 + (r2*sin(fi2+angle)); 
    yen = y0 + (r2*cos(fi2+angle)); 

    if (xsn < 0. || ysn < 0. || xen < 0.|| yen < 0.)    return 0;
    
    if (pMFN->GetX()<static_cast<unsigned int>(xsn/dx))    return 0;
    if (pMFN->GetY()<static_cast<unsigned int>(ysn/dy))    return 0;
    if (pMFN->GetX()<static_cast<unsigned int>(xen/dx))    return 0;
    if (pMFN->GetY()<static_cast<unsigned int>(yen/dy))    return 0;

    return 1;
}
// <------------- 2D --------------->           



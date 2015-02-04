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
#include "libOpenHyperFLOW2D/hyper_flow_bound_contour.hpp"
// <------------- 2D --------------->           
// Bound contour constructor
BoundContour2D::BoundContour2D(char* Name, 
                               UMatrix2D< FlowNode2D< FP, NUM_COMPONENTS> >* fnm,
                               unsigned int x,
                               unsigned int y
                               ) {
    current_x=first_x=x;
    current_y=first_y=y;
    FlowNodeMatrixPtr=fnm;
    isContourClosed=isActivateContour=0;
    BoundContourName = Name; 
}

// Clear bound contour
void BoundContour2D::ClearBoundContour() {
    for (int i=0;i<GetNumBounds();i++) {
        delete GetBound(i);
    }
    CleanArray();
}
// Bound Contour destructor
BoundContour2D::~BoundContour2D() {
    ClearBoundContour();
}

// Add bound in bound contour
int BoundContour2D::AddBound2D(char* name,
                               unsigned int x,
                               unsigned int y,
                               ulong        bt,
                               Flow*        pInFlow,
                               Flow2D*      pInFlow2D,
                               FP*      Y,
                               ulong        btt) {
    Bound2D* TmpBound=NULL;
    if (isActivateContour ||
        isContourClosed)     return -1;
    if (pInFlow)
        TmpBound = new Bound2D(name,FlowNodeMatrixPtr,
                               current_x,current_y,
                               x,y,bt,pInFlow,Y,btt);
    else if (pInFlow2D)
        TmpBound = new Bound2D(name,FlowNodeMatrixPtr,
                               current_x,current_y,
                               x,y,bt,pInFlow2D,Y,btt);
    else
        TmpBound = new Bound2D(name,FlowNodeMatrixPtr,
                               current_x,current_y,
                               x,y,bt,pInFlow,Y,btt);
    current_x = x;
    current_y = y;
    return AddElement(&TmpBound);
}
// To close a contour
int BoundContour2D::CloseContour2D(char*    name,
                                   ulong    bt,
                                   Flow*    pInFlow,
                                   Flow2D*  pInFlow2D,
                                   FP*  Y,
                                   ulong    btt) {
    Bound2D* TmpBound;
    if (isActivateContour ||
        isContourClosed)     return -1;
    if (GetNumBounds()<2) return -1;

    if (pInFlow)
        TmpBound = new Bound2D(name,
                              FlowNodeMatrixPtr,
                             current_x,current_y,
                             first_x,first_y,
                             bt,pInFlow,Y,btt);
    else if (pInFlow2D)
      TmpBound = new Bound2D(name,
                             FlowNodeMatrixPtr,
                             current_x,current_y,
                             first_x,first_y,
                             bt,pInFlow2D,Y,btt);
    else
        TmpBound = new Bound2D(name,
                               FlowNodeMatrixPtr,
                               current_x,current_y,
                               first_x,first_y,
                               bt,pInFlow,Y,btt);
    current_x = first_x;
    current_y = first_y;
    isContourClosed=1;
    return AddElement(&TmpBound);
}
// Insert bound in contour (closed contour too)
int BoundContour2D::InsertBound2D(char* name,
                                  unsigned int nb,
                                  unsigned int x,
                                  unsigned int y,
                                  ulong    bt,
                                  Flow*    pInFlow,
                                  Flow2D*  pInFlow2D,
                                  FP*   Y,
                                  ulong btt) {
    return -1; // function is not implemented now...
}
// Delete bound from contour
int BoundContour2D::DelBound2D(int nb) {
    Bound2D* RemovedBound;
    Bound2D* NextBound;
    int      nextIndex=0;

    if (isActivateContour)     return -1;
    if (nb < 0 || nb >= GetNumBounds()) return -1;

    RemovedBound=GetBound(nb);

    if (nb+1 == GetNumBounds())
        if (isContourClosed) nextIndex=0;
        else return -1;
    else nextIndex=nb+1;

    NextBound=GetBound(nextIndex);

    NextBound->SetStartX(RemovedBound->GetStartX());
    NextBound->SetStartY(RemovedBound->GetStartY());

    delete RemovedBound;
    return DelElement(nb);
}
// Set bounds in contour
int BoundContour2D::SetBounds(int MaterialID) {
    int i;
    if (!isContourClosed) return -1;
    for (i=0;i<GetNumBounds();i++) {
        GetBound(i)->SetBound(MaterialID);
        if (GetBound(i)->GetBoundState() == BND_ERR) {
            throw(i);
        }
    }
    isActivateContour=1;
    return GetNumBounds();  
}
// Set bounds in contour and push result in to the array
int BoundContour2D::SetBounds(UArray<FlowNode2D< FP, NUM_COMPONENTS>* >* node_array,int MaterialID) {
    int i;
    if (!isContourClosed) return -1;
    for (i=0;i<GetNumBounds();i++) {
        GetBound(i)->SetBound(node_array,MaterialID);
        if (GetBound(i)->GetBoundState() == BND_ERR) {
            throw(i);
        }
    }
    isActivateContour=1;
    return GetNumBounds();  
}

// Get current X coord (nodes)
int    BoundContour2D::GetCurrentX() {
    return current_x;
}
// Get current Y coord (nodes)
int    BoundContour2D::GetCurrentY() {
    return current_y;
}
// Set first bound X coord (nodes)
void   BoundContour2D::SetFirstX(int x) {
    first_x=x;
}
// Set first bound Y coord (nodes)
void   BoundContour2D::SetFirstY(int y) {
    first_y=y;
}
// Set current bound X coord (nodes)
void   BoundContour2D::SetCurrentX(int x) {
    current_x=x;
}
// Set current bound Y coord (nodes)
void   BoundContour2D::SetCurrentY(int y) {
    current_y=y;
}
// Get current bound X coord (m)
FP BoundContour2D::GetCurrentFX() {
    return f_current_x;
}
// Get current bound Y coord (m)
FP BoundContour2D::GetCurrentFY() {
    return f_current_y;
}
// Set current bound X coord (m)
void   BoundContour2D::SetCurrentFX(FP x) {
    f_current_x=x;
}
// Set current bound Y coord (m)
void   BoundContour2D::SetCurrentFY(FP y) {
    f_current_y=y;
}
// Set first bound X coord (m)
void   BoundContour2D::SetFirstFX(FP x) {
    f_first_x=x;
}
// Set first bound Y coord (m)
void   BoundContour2D::SetFirstFY(FP y) {
    f_first_y=y;
}
// Get num bounds in contour
int    BoundContour2D::GetNumBounds() {
    return(int)GetNumElements();
}
// Get bound by index
Bound2D* BoundContour2D::GetBound(int nb) {
    return GetElement(nb);
}
// Is contour Closed
int    BoundContour2D::IsContourClosed() {
    return isContourClosed;
}
// Is contour activated
int    BoundContour2D::IsContourActivate() {
    return isActivateContour;
}
// Rotate bound contour
int    BoundContour2D::RotateBoundContour2D(FP x0, FP y0, FP angle) {
    int i,a=1;
    if (IsContourActivate()) return 0;

    for (i=0;i<GetNumBounds();i++) {
        a *= GetBound(i)->TestRotateBound2D(x0,y0,angle);
    }
    if (!a) return a;

    for (i=0;i<GetNumBounds();i++) {
        GetBound(i)->RotateBound2D(x0,y0,angle);
    }
    return a;
}
// <------------- 2D --------------->           


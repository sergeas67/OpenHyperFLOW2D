/*******************************************************************************
*   OpenHyperFLOW2D                                                            *
*                                                                              *
*   Version  2.0.1                                                             *
*   Copyright (C)  1995-2016 by Serge A. Suchkov                               *
*   Copyright policy: LGPL V3                                                  *
*   http://github.com/sergeas67/openhyperflow2d                                *
*                                                                              *
*   Shape2D object implementation                                              *
*                                                                              *
*   last update: 14/04/2016                                                    *
*******************************************************************************/
#include "libExcept/except.hpp"
#include "libOpenHyperFLOW2D/hyper_flow_shape2d.hpp"
// <------------- 2D --------------->           
// Bound contour constructor
Shape2D::Shape2D(char* Name, 
                 UMatrix2D< FlowNode2D< FP, NUM_COMPONENTS> >* fnm,
                 unsigned int x,
                 unsigned int y ) {
    
    current_x=first_x=x;
    current_y=first_y=y;
    FlowNodeMatrixPtr=fnm;
    FlowNodeMatrixPtr=fnm;
    isActivateShape=0;
    ShapeName = Name; 
}

// Clean Shape2D
void Shape2D::CleanShape2D() {
    for (int i=0;i<GetNumBounds();i++) {
        delete GetBound(i);
    }
    CleanArray();
}

// Shape2D destructor
Shape2D::~Shape2D() {
    CleanShape2D();
}

// Add bound in shape
int Shape2D::AddBound2D(char* name,
                        unsigned int x,
                        unsigned int y,
                        ulong        bt,
                        Flow*        pInFlow,
                        Flow2D*      pInFlow2D,
                        FP*      Y,
                        ulong        btt) {
    Bound2D* TmpBound=NULL;
    if (isActivateShape)     return -1;
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

// Add bound in shape
int Shape2D::AddBound2D(char* name,
                        FP           x,
                        FP           y,
                        ulong        bt,
                        Flow*        pInFlow,
                        Flow2D*      pInFlow2D,
                        FP*          Y,
                        ulong        btt) {
    Bound2D* TmpBound=NULL;
    if (isActivateShape)     return -1;
    if (pInFlow)
        TmpBound = new Bound2D(name,FlowNodeMatrixPtr,
                               f_current_x,f_current_y,
                               x,y,bt,pInFlow,Y,btt);
    else if (pInFlow2D)
        TmpBound = new Bound2D(name,FlowNodeMatrixPtr,
                               f_current_x,f_current_y,
                               x,y,bt,pInFlow2D,Y,btt);
    else
        TmpBound = new Bound2D(name,FlowNodeMatrixPtr,
                               f_current_x,f_current_y,
                               x,y,bt,pInFlow,Y,btt);
    f_current_x = x;
    f_current_y = y;
    return AddElement(&TmpBound);
}


// Insert bound in shape 
int Shape2D::InsertBound2D(char* name,
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


// Insert bound in shape 
int Shape2D::InsertBound2D(char* name,
                           unsigned int nb,
                           FP       x,
                           FP       y,
                           ulong    bt,
                           Flow*    pInFlow,
                           Flow2D*  pInFlow2D,
                           FP*   Y,
                           ulong btt) {
    return -1; // function is not implemented now...
}

// Delete bound from shape
int Shape2D::DelBound2D(int nb) {
    Bound2D* RemovedBound;
    Bound2D* NextBound;
    int      nextIndex=0;

    if (isActivateShape)     return -1;
    if (nb < 0 || nb >= GetNumBounds()) return -1;

    RemovedBound=GetBound(nb);

    if (nb+1 == GetNumBounds()) {
        delete RemovedBound;
        return DelElement(nb);
    } else { 
        nextIndex=nb+1;
    }

    NextBound=GetBound(nextIndex);

    NextBound->SetStartX(RemovedBound->GetStartX());
    NextBound->SetStartY(RemovedBound->GetStartY());

    delete RemovedBound;
    return DelElement(nb);
}

// Set bounds in shape
int Shape2D::SetBounds(int MaterialID) {
    int i;
    for (i=0;i<GetNumBounds();i++) {
        GetBound(i)->SetBound(MaterialID);
        if (GetBound(i)->GetBoundState() == BND_ERR) {
            throw(i);
        }
    }
    isActivateShape=1;
    return GetNumBounds();  
}

// Set bounds in shape and push result in to the array
int Shape2D::SetBounds(UArray<FlowNode2D< FP, NUM_COMPONENTS>* >* node_array,int MaterialID) {
    int i;
    for (i=0;i<GetNumBounds();i++) {
        GetBound(i)->SetBound(node_array,MaterialID);
        if (GetBound(i)->GetBoundState() == BND_ERR) {
            throw(i);
        }
    }
    isActivateShape=1;
    return GetNumBounds();  
}

// Get current X coord (nodes)
int    Shape2D::GetCurrentX() {
    return current_x;
}
// Get current Y coord (nodes)
int    Shape2D::GetCurrentY() {
    return current_y;
}
// Set first bound X coord (nodes)
void   Shape2D::SetFirstX(int x) {
    first_x=x;
}
// Set first bound Y coord (nodes)
void   Shape2D::SetFirstY(int y) {
    first_y=y;
}
// Set current bound X coord (nodes)
void   Shape2D::SetCurrentX(int x) {
    current_x=x;
}
// Set current bound Y coord (nodes)
void   Shape2D::SetCurrentY(int y) {
    current_y=y;
}
// Get current bound X coord (m)
FP Shape2D::GetCurrentFX() {
    return f_current_x;
}
// Get current bound Y coord (m)
FP Shape2D::GetCurrentFY() {
    return f_current_y;
}
// Set current bound X coord (m)
void   Shape2D::SetCurrentFX(FP x) {
    f_current_x=x;
}
// Set current bound Y coord (m)
void   Shape2D::SetCurrentFY(FP y) {
    f_current_y=y;
}
// Set first bound X coord (m)
void   Shape2D::SetFirstFX(FP x) {
    f_first_x=x;
}
// Set first bound Y coord (m)
void   Shape2D::SetFirstFY(FP y) {
    f_first_y=y;
}
// Get num bounds in contour
int    Shape2D::GetNumBounds() {
    return(int)GetNumElements();
}
// Get bound by index
Bound2D* Shape2D::GetBound(int nb) {
    return GetElement(nb);
}

// Is contour activated
int    Shape2D::IsShapeActivate() {
    return isActivateShape;
}
// Rotate shape
int    Shape2D::RotateShape2D(FP x0, FP y0, FP angle) {
    int i,a=1;
    if (IsShapeActivate()) return 0;

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


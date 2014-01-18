/*******************************************************************************
*   OpenHyperFLOW2D                                                            *
*                                                                              *
*   Version  1.0.1                                                             *
*   Copyright (C)  1995-2014 by Serge A. Suchkov                               *
*   Copyright policy: LGPL V3                                                  *
*                                                                              *
*   last update: 16/01/2014                                                    *
*******************************************************************************/
#include "libExcept/except.hpp"
#include "libOpenHyperFLOW2D/hyper_flow_field.hpp"

// FlowField2D constructor
FlowField2D::FlowField2D(char* filename, InputData*  data):
UMatrix2D< FlowNode2D< double, NUM_COMPONENTS> >(
              (FlowNode2D< double, NUM_COMPONENTS>*)
               LoadSwapFile2D(filename,
               data->GetIntVal((char*)"MaxX"),
               data->GetIntVal((char*)"MaxY"),
               sizeof(FlowNode2D< double, NUM_COMPONENTS>),
               &p,
               &fd,
               (ofstream*)(data->GetMessageStream())),
               data->GetIntVal((char*)"MaxX"),
               data->GetIntVal((char*)"MaxY"))
                                                                                                                           
{
  FlowFieldName = data->GetStringVal((char*)"ProjectName");;
}
// Field2D constructor
FlowField2D::FlowField2D(char* name, int x, int y):UMatrix2D< FlowNode2D< double, NUM_COMPONENTS> >(x,y) {
}
// FlowField2D constructor
FlowField2D::FlowField2D(FlowField2D* f2d):UMatrix2D< FlowNode2D< double, NUM_COMPONENTS> >(f2d->GetMatrixPtr(),f2d->GetX(),f2d->GetY()){
    if(f2d->GetFlowFieldName()) {
      FlowFieldName = new char[strlen(f2d->GetFlowFieldName())+1];
      strcpy(FlowFieldName,f2d->GetFlowFieldName());
    }
}
// FlowField2D constructor
FlowField2D::FlowField2D(char* name, UMatrix2D< FlowNode2D< double, NUM_COMPONENTS> >* m2d):UMatrix2D< FlowNode2D< double, NUM_COMPONENTS> >(m2d->GetMatrixPtr(),m2d->GetX(),m2d->GetY()){
    if(name) {
      FlowFieldName = new char[strlen(name)+1];
      strcpy(FlowFieldName,name);
    }
}
// FlowField2D constructor
FlowField2D::~FlowField2D() {
  delete FlowFieldName;
  //delete (UMatrix2D< FlowNode2D< double, NUM_COMPONENTS> >*)(this)
}
// FlowField2D constructor
int FlowField2D::SaveFlowField2D(char* filename){
    int   fd_g,pf;
    ostream*     o_stream = &cout;
    ofstream*    f_stream = (ofstream*)(o_stream);
    
    ssize_t max_write = 1024L*1024L*1024L;                                                                                                      
    ssize_t one_write = 0L;                                                                                                                     
    ssize_t len  = 0L;                                                                                                                        
    off_t  off = 0L;                                                                                                                           
    LoadSwapFile2D(filename,GetX(),GetY(),sizeof(FlowNode2D< double, NUM_COMPONENTS>),&pf,&fd_g,f_stream);
    char*  TmpPtr=(char*)GetMatrixPtr();                                                                                                             
    if(GetMatrixSize() > max_write) {                                                                                                                 
        for(off = 0L,one_write = max_write; len < GetMatrixSize(); off += max_write) {                                                                
           len += pwrite64(fd_g,TmpPtr+off,one_write,off);                                                                                       
           if(GetMatrixSize() - len < max_write)                                                                                                      
             one_write = GetMatrixSize() - len;                                                                                                       
           }                                                                                                                                   
          if(len != GetMatrixSize())                                                                                                                    
          *f_stream << "Error: len(" << len << ") != MatrixSize(" << GetMatrixSize() << ") " << endl << flush;                                               
    } else {                                                                                                                                   
        len = pwrite64(fd_g,GetMatrixPtr(),GetMatrixSize(),0L);                                                                                               
    }
    close(fd_g);
    return len;
}



/*******************************************************************************
*   OpenHyperFLOW2D                                                            *
*                                                                              *
*   Version  1.0.2                                                             *
*   Copyright (C)  1995-2015 by Serge A. Suchkov                               *
*   Copyright policy: LGPL V3                                                  *
*                                                                              *
*   last update: 11/01/2015                                                    *
*******************************************************************************/
#ifndef _flow_list_hpp_
#define _flow_list_hpp_

#include "utl/uarray.hpp"
#include "obj_data/obj_data.hpp"
#include "libOpenHyperFLOW2D/hyper_flow_node.hpp"

class FlowList : public UArray<Flow*> {
InputData* data;
public:
 FlowList(InputData*);
 ~FlowList();
};
class Flow2DList : public UArray<Flow2D*> {
InputData* data;
public:
 Flow2DList(InputData*);
 ~Flow2DList();
};
#endif // _flow_list_hpp_

/***********************************************************************
*  Exceptlib v0.2beta (System Exceptions C++ Library)                  *
*                                                                      *
*  Copyright (C)  1995-2016 by Serge A. Suchkov                        *
*  Copyright policy: LGPL V3                                           *
*  http://github.com/sergeas67/openhyperflow2d                         *
***********************************************************************/

#ifndef  _except_hpp_
#define  _except_hpp_

#ifdef __ICC
#include "intel_compatibility.h"
#endif //__ICC
 
#include <stdlib.h>
#include <stdio.h>
#include <setjmp.h>
#include <signal.h>
#include <iostream>

using namespace std;

struct ExceptContext {jmp_buf EC;};

#define  _SAFE_ACCESS_       // pthreads effective only
#include "utl/uarray.hpp"
#include "utl/ustack.hpp"

/*
Usage:

___try {
...
throw({type1})
...
throw({type2})
...
throw({type2})
...
}__except(SysException) {   // System errors (SIGNALS)
...
}__except({type1}) {
...
}__except({type2}) {
...
}__except({type3}) {
...
}__end_except;

*/

typedef int SysException;

enum ActionType
     {
      AT_DEFAULT,
      AT_IGNORE,
      AT_HANDLED
     };


// ExceptLib
class ExceptLib : public UStack<ExceptContext>
{
 ActionType   signal2exception[32];

public:

 ExceptLib();
 ~ExceptLib();
 int  SetSystemException(int sig, ActionType at);
 void SignalHandlerOn(int sig=0);
 void SignalHandlerOff();
 int  GetHandlerState();  // if 0 -> _off_ else _on_
 ActionType GetHandlerAction(int sig);
};

// ContextControl
class ContextControl 
{
 ExceptLib*    except_lib_ptr;
 ExceptContext tmp_env;
public:
 ContextControl(ExceptLib*);
 ~ContextControl();
 ExceptContext& GetCurrentContext();
 unsigned int    GetStackSize();
};

// Signal handler
void GlobalErrorHandler(int sig);

extern ExceptLib   __ExceptLib;    // ExceptiLib instance

#define ___try try {\
     ContextControl CC(&__ExceptLib);\
     __ExceptLib.SignalHandlerOn();\
     SysException e=(SysException)setjmp(CC.GetCurrentContext().EC);\
     if(e!=0) throw(e);\
        else

#define  __except(a)   } catch(a)  {
//
//  WARNING!
//  Defenition "__except(a)" parse uncorrect ("comma problem"), when "a" is a
//  multiparameters template, same as SomeClass<a,b>.
//  In this case use:
//   - typedef MaskedClass  SomeClass<a,b>
//     "...}__except(MaskedClass){...",
//   - "}}catch(SomeClass<a,b>){{" expression
//   - dummy  subclassing: "class SubClass : public SomeClass<a,b>{;};"
//     "...}__except(SubClass){..."

#define  __end_except  } catch(SysException e){cout <<  "\n<<LibExcept>>: Throw unhandled system exception (signal " << e << " )\n" << flush;\
         }catch(...){cout <<  "\n<<LibExcept>>: Throw unhandled exception.\n" << flush;}

#endif // _except_hpp_


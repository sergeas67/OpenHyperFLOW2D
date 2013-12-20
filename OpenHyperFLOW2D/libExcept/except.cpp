/***********************************************************************
*  Exceptlib v0.2beta (System Exceptions C++ Library)                  *
*                                                                      *
*  Copyright (C)  1995-2013 by Serge A. Suchkov                        *
*  Copyright policy: LGPL V3                                           *
***********************************************************************/

#include "libExcept/except.hpp"

ExceptLib   __ExceptLib;    // ExceptiLib instance

// +Constructor+
ExceptLib::ExceptLib():UStack<ExceptContext>() {
    for (int i=0;i<32;i++)
        signal2exception[i]=AT_DEFAULT;
}

// -Destructor-
ExceptLib::~ExceptLib() {
    ::exit(0);
}

// Get context stack size
int ExceptLib::GetHandlerState() {
    return GetStackSize();
}

// Get action type from signal num.
ActionType ExceptLib::GetHandlerAction(int sig) {
    if (sig<0 || sig > 31) return AT_DEFAULT;
    else                  return signal2exception[sig];
}

// Signal handler
void GlobalErrorHandler(int sig) {
    if (__ExceptLib.GetStackSize()>0) {
        __ExceptLib.SignalHandlerOn(sig);
        longjmp(__ExceptLib.Top().EC,sig);   // jump to exception handler
    }
}

// Set exception handler
int ExceptLib::SetSystemException(int sig, ActionType at) {
    if (sig<0 || sig > 31) return -1;
    else signal2exception[sig]=at;
    SignalHandlerOn();
    return 0;
}

// Set all signals handlers
void ExceptLib::SignalHandlerOn(int sig) {
    if (sig==0 ) {
        for (int i=0;i<32;i++) {
            if (signal2exception[i]==AT_DEFAULT)
#ifndef _MIT_POSIX_THREADS
                signal(i,SIG_DFL);
#else
            {
                ;
            }
#endif
            else if (signal2exception[i]==AT_IGNORE)  signal(i,SIG_IGN);
            else if (signal2exception[i]==AT_HANDLED) signal(i,GlobalErrorHandler);
        }
    } else if ( sig > 0 && sig < 32) {
        if (signal2exception[sig]==AT_DEFAULT)
#ifndef _MIT_POSIX_THREADS
            signal(sig,SIG_DFL);
#else
        {
            ;
        }
#endif
        else if (signal2exception[sig]==AT_IGNORE)  signal(sig,SIG_IGN);
        else if (signal2exception[sig]==AT_HANDLED) signal(sig,GlobalErrorHandler);
    }

}

// Remove all signals handlers
void ExceptLib::SignalHandlerOff() {
    if (GetStackSize()==0)
        for (int i=0;i<32;i++)
#ifndef _MIT_POSIX_THREADS
            signal(i,SIG_DFL);
#else
        {
            ;
        }
#endif
}

// +Constructor+
ContextControl::ContextControl(ExceptLib* stack_ptr) {
    except_lib_ptr=stack_ptr;
    except_lib_ptr->Push(&tmp_env);
}

// -Destructor-
ContextControl::~ContextControl() {
    except_lib_ptr->Pop();
}

// Get current context
ExceptContext& ContextControl::GetCurrentContext() {
    return except_lib_ptr->Top();
}

// Get context stack size
unsigned int ContextControl::GetStackSize() {
    return except_lib_ptr->GetStackSize();
}

/************************************************************
*  UNKNOWN TEMPLATES LIBRARY (UTL) v1.3L                    *
*  Copyright (C) 1994-2013 Serge A. Suchkov                 *
*  Copyright policy: LGPL V3                                *
*  Threads Locker base class.                               *
*                                                           *
*  Please report all bugs and problems to "sergeas67@gmail.com". *
************************************************************/

#ifndef _locker_hpp_
#define _locker_hpp_

#ifdef __ICC
#include "intel_compatibility.h"
#endif //__ICC

#include <stdlib.h>
#include <stdio.h>
#include <signal.h>
#include <string.h>

#ifdef _PTHREAD_
#ifndef __USE_GNU
    #define  __USE_GNU
#endif // __USE_GNU
    #include <pthread.h>
#endif // _PTHREAD_

#ifdef _WIN32
    #include <windows.h>
#else
    #include <unistd.h>
#endif // _WIN32

class Locker {
    #ifdef _WIN32
    CRITICAL_SECTION CS;
    #endif // _WIN32

    #ifdef _PTHREAD_
    pthread_mutex_t  PM;
    #endif //  _PTHREAD_

    int   LockCount;

public:

// Locker default  constructor
    Locker()
    {
    #ifdef _WIN32
        LockCount=0;InitializeCriticalSection(&CS);
    #endif // _WIN32
    #ifdef _PTHREAD_
        pthread_mutex_t pm=PTHREAD_RECURSIVE_MUTEX_INITIALIZER_NP;
        LockCount=0;
        PM = pm;
    #endif // _PTHREAD_
    }

// Locker  destructor
    virtual  ~Locker()
    {
        Free();
    #ifdef _WIN32
        DeleteCriticalSection(&CS);
    #endif // _WIN32
    #ifdef _PTHREAD_
        pthread_mutex_destroy(&PM);
    #endif //_PTHREAD_
    }

// lock function
    virtual  int     Lock()
    {LockCount++;
    #ifdef _WIN32
        EnterCriticalSection(&CS);
    #endif // _WIN32
    #ifdef _PTHREAD_
        pthread_mutex_lock(&PM);
    #endif // _PTHREAD_
        return LockCount;}

// unlock function
    virtual  int     UnLock()
    {LockCount--;
    #ifdef _WIN32
        LeaveCriticalSection(&CS);
    #endif // _WIN32
    #ifdef _PTHREAD_
        pthread_mutex_unlock(&PM);
    #endif // __PTHREAD_
        return LockCount;}
// Free all locks
    virtual  int Free()
    {
        int i;
        for (i=0;i<LockCount;i++)
    #ifdef _WIN32
            LeaveCriticalSection(&CS);
    #endif // _WIN32
    #ifdef _PTHREAD_
        pthread_mutex_unlock(&PM);
    #endif // __PTHREAD_
        LockCount=0;
        return i;
    }

// Get Locks count function
    virtual  int     GetLockCount() {return LockCount;}

// Get Locker object
    Locker* GetLocker()             {return  this;}
};

// Check Locker object
class CheckLocker {
    Locker*  lk;
public:
    CheckLocker(Locker* l) {lk=l;l->Lock();}
    ~CheckLocker()         {lk->UnLock();}
};

#endif //_locker_hpp_

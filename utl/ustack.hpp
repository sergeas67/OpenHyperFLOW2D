/****************************************************************
*  UNKNOWN TEMPLATES LIBRARY (UTL) v1.3                         *
*  Copyright (C) 1994-2010 Serge A. Suchkov                     *
*  Copyright policy: LGPL V2.0                                  *
*  Stack template (LIFO object).                                *
*  Please report all bugs and problems to "sergeas67@gmail.com".*
*****************************************************************/
#ifndef _ustack_hpp_
    #define _ustack_hpp_

    #include  "utl/uarray.hpp"

typedef  ArrayState StackState;

template <class T>
class UStack  : private UArray<T> {

    T   TopElement;

public:

    UStack(int MaxSize=-1);
    UStack(UStack& Stack);
    UStack(UArray<T>& Array);
    ~UStack();

    inline int          GetMaxStackSize() {return UArray<T>::GetMaxNumElements();}
    inline unsigned int GetStackSize()    {return UArray<T>::GetNumElements();}
    inline StackState   GetStackState()   {return UArray<T>::GetArrayState();}
    inline void         ClearStack()      {UArray<T>::CleanArray();}

    inline unsigned int Push(T*);
    inline T&           Pop(T* t=NULL);
    inline T&           Top();
    inline T&           Peek(unsigned int sp=0);
    inline UStack<T>&   operator = (UStack<T>& NewT);
    inline UStack<T>&   operator = (UArray<T>& NewT);
    inline UStack<T>&   operator << (T&);
    inline UStack<T>&   operator >> (T&);
};

// +++ Stack constructors +++
template <class T>
UStack<T>::UStack(int MaxSize):UArray<T>(MaxSize){;}
template <class T>
UStack<T>::UStack(UStack& Stack):UArray<T>(Stack){;}
template <class T>
UStack<T>::UStack(UArray<T>& Array):UArray<T>(Array){;}

template <class T>
UStack<T>::~UStack() {;}

// --- Stack destructor --- 
// Push
template <class T>
unsigned int UStack<T>::Push(T* pt)
{
    UArray<T>::AddElement(pt);
    return GetStackSize();
}

// Pop
template <class T>
T& UStack<T>::Pop(T* pt) {
    if (pt==NULL) {
        memcpy(&TopElement,UArray<T>::GetElementPtr(GetStackSize()-1),sizeof(T));
        UArray<T>::DelElement(GetStackSize()-1);
        return TopElement;
    } else {
        memcpy(pt,UArray<T>::GetElementPtr(GetStackSize()-1),sizeof(T));
        UArray<T>::DelElement(GetStackSize()-1);
        return  *pt;
    }
}

// Pop
template <class T> 
T& UStack<T>::Top() {
    return(*this)[GetStackSize()-1];
}

// Peek
template <class T>
T& UStack<T>::Peek(unsigned int sp) {
    return(*this)[sp];
}

// operator  '='  (UStack<T>=UStack<T>)
template <class T>
UStack<T>&   UStack<T>::operator = (UStack<T>& NewT) {
    (*this)=NewT;
    return *this;
}

// operator  '='  (UStack<T>=UArray<T>)
template <class T>
UStack<T>&   UStack<T>::operator = (UArray<T>& NewT) {
    *this=NewT;
    return *this;
}

// stream operator <<
template <class T>
UStack<T>&  UStack<T>::operator << (T& t) {
    Push(&t);
    return *this;
}

// stream operator >>
template <class T>
UStack<T>&  UStack<T>::operator >> (T& t) {
    t=Pop();
    return *this;
}
#endif  //_ustack_hpp_

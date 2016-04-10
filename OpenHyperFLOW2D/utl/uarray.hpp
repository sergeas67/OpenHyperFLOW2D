/****************************************************************
*  UNKNOWN TEMPLATES LIBRARY (UTL) v1.3                         *
*  Copyright (C) 1994-2010 Serge A. Suchkov                     *
*  Copyright policy: LGPL V2.0                                  *
*  Array  template.                                             *
*  http://github.com/sergeas67/openhyperflow2d                        *
*  Please report all bugs and problems to "sergeas67@gmail.com".*
*****************************************************************/

#ifndef _uarray_hpp_
#define _uarray_hpp_

#include "utl/locker.hpp"

enum ArrayState
    {
     AS_NOERROR,         // no errors
     AS_ERR_MEM,         // memory allocation error during operation
     AS_ERR_OUT_OF_INDEX,// attempt out of index range
     AS_ERR_MAX_SIZE,    // attempt of oversizing
     AS_ERR_VALUE        // attempt use bad value
    };

template <class T>
class UArray : public Locker
{
    T*             ArrayPtr;
    int            MaxNumElements;  // -1 - unlimited...
    unsigned int   NumElements;
    ArrayState     as;

public:
    UArray(unsigned int Num,int  MaxNum); // Full specified array constructor
    UArray(unsigned int Num);             // Size specified array constructor (MaxSize=unlimited)
    UArray(int  MaxNum=-1);               // MaxSize specified array constructor (Size=0)
    UArray(UArray& Array);                // Copy array constructor

    virtual  ~UArray();
    virtual  int      AddElement(T*);
    virtual  int      DelElement(unsigned int);
    virtual  T*       GetElementPtr(unsigned int i);
    virtual  T*       GetArrayPtr();
    virtual  T&       GetElement(unsigned int i);
    virtual  int      FindElement(T* Ptr,unsigned int start_i=0);
    virtual  void     SetElement(unsigned int i, T*);
    virtual  void     SetMaxNumElements(int NumMax);
    UArray<T>& operator = (UArray<T>& NewT);
    T& operator [] (unsigned int i) {return *GetElementPtr(i);}
    UArray<T>& operator += (T& ET);
    UArray<T>& operator += (UArray<T>& AT);

    virtual  void     CleanArray();
    virtual  unsigned int GetNumElements() {return NumElements;}
    virtual  int      GetMaxNumElements()  {return MaxNumElements;}
    virtual  unsigned int GetElementSize() {return sizeof(T);}
    ArrayState GetArrayState()             {return as;}

};

// +++ Array Constructors +++
// Full specified array constructor
template <class T>
UArray<T>::UArray(unsigned int Num,int  MaxNum):Locker(),ArrayPtr(NULL)
{
#ifdef _SAFE_ACCESS_
    CheckLocker cl(GetLocker());
#endif //_SAFE_ACCESS_ 

 if(MaxNum==0)     
     {as=AS_ERR_VALUE;throw(this);}
 if(MaxNum != -1 &&
    (int)Num >  MaxNum) 
     {as=AS_ERR_MAX_SIZE;throw(this);}
 MaxNumElements = MaxNum;
 if(Num>0)
   {
    ArrayPtr = new T[Num];
    if(ArrayPtr == NULL) {as=AS_ERR_MEM;throw(this);}
    else  NumElements = Num;
   }
 else NumElements = 0;
 as=AS_NOERROR;
}

// Size specified array constructor (MaxSize=unlimited)
template <class T>
UArray<T>::UArray(unsigned int Num):Locker(),ArrayPtr(NULL)
{
#ifdef _SAFE_ACCESS_
    CheckLocker cl(GetLocker());
#endif //_SAFE_ACCESS_ 

 MaxNumElements = -1;   
 if(Num>0)
   {
    ArrayPtr = new T[Num];
    if(ArrayPtr == NULL) {as=AS_ERR_MEM;throw(this);}
    else  NumElements = Num;
   }
 else NumElements = 0;
 as=AS_NOERROR;
}

// MaxSize specified array constructor (Size=0)
template <class T>
UArray<T>::UArray(int MaxNum):Locker(),ArrayPtr(NULL)
{
#ifdef _SAFE_ACCESS_
 CheckLocker cl(GetLocker());
#endif //_SAFE_ACCESS_ 
 if(MaxNum==0) {as=AS_ERR_VALUE;throw(this);}
 MaxNumElements = MaxNum;   
 NumElements    = 0;
 as=AS_NOERROR;
}

// Copy array constructor
template <class T>
UArray<T>::UArray(UArray& Array):Locker(),ArrayPtr(NULL)
{
 #ifdef _SAFE_ACCESS_
    CheckLocker cl1(GetLocker());
    CheckLocker cl2(Array.GetLocker());
 #endif //_SAFE_ACCESS_ 

 *this=Array;
}

// ---Array destructor---
template <class T>
UArray<T>::~UArray()
 {
  CleanArray();
 }


// Find Element
template <class T>
int UArray<T>::FindElement(T* Ptr,unsigned int start_i) 
{
#ifdef _SAFE_ACCESS_
  CheckLocker cl(GetLocker());
#endif //_SAFE_ACCESS_ 
  
 if(NumElements==0 ) {
      as=AS_NOERROR;
      return -1;
  }
  
 if(start_i>NumElements) {
    as=AS_ERR_OUT_OF_INDEX;
    throw(this);
   }
 
 for(register unsigned int i=start_i;i<NumElements;i++)
    {
     if(memcmp(Ptr,&ArrayPtr[i],sizeof(T))==0)
       {
        as=AS_NOERROR;
        return (int)i;
       }
    }
 as=AS_NOERROR;
 return -1;
}

// Add Element
template <class T>
int UArray<T>::AddElement(T* Ptr)  
{
    T* TmpTPtr;
#ifdef _SAFE_ACCESS_
    CheckLocker cl(GetLocker());
#endif //_SAFE_ACCESS_ 
    if((int)NumElements ==  MaxNumElements) {as=AS_ERR_MAX_SIZE;throw(this);}
    TmpTPtr = new T[NumElements + 1];
    if(TmpTPtr == NULL) {as=AS_ERR_MEM;throw(this);}
    if(NumElements > 0)
      {
       memcpy(TmpTPtr, ArrayPtr, sizeof(T)*NumElements);
       delete[] ArrayPtr;
      }
    memcpy((TmpTPtr+NumElements), Ptr, sizeof(T));
    ArrayPtr = TmpTPtr;
    NumElements++;
    as=AS_NOERROR;
    return NumElements-1;
}

// Del Element
template <class T>
int UArray<T>::DelElement(unsigned int iDel)  
{
    register unsigned int i = 0;
    register unsigned int j = 0;
    T*   TmpTPtr=(T*)NULL;
#ifdef _SAFE_ACCESS_
    CheckLocker cl(GetLocker());
#endif //_SAFE_ACCESS_ 
    if(iDel >= NumElements) {as=AS_ERR_OUT_OF_INDEX;throw(this);}
    else if(NumElements == 1)
    {
     delete[] ArrayPtr;
     NumElements = 0;
     ArrayPtr = (T*)NULL;
     as=AS_NOERROR;
     return 0;
    }
    TmpTPtr = new T[NumElements-1];
    if (TmpTPtr == NULL) {as=AS_ERR_MEM;throw(this);}
    while(j<NumElements)
         {
          if(j != iDel)
            {
             memcpy(TmpTPtr+i, ArrayPtr+j,sizeof(T));
             i++;
            }
           j++;
         }
    NumElements--;
    delete[] ArrayPtr;
    ArrayPtr = TmpTPtr;
    as=AS_NOERROR;
    return 0;
}

// Get Element
template <class T>
T* UArray<T>::GetElementPtr(unsigned int iGet)  
{
#ifdef _SAFE_ACCESS_
    CheckLocker cl(GetLocker());
#endif //_SAFE_ACCESS_ 
    if(NumElements==0  || iGet >= NumElements) {as=AS_ERR_OUT_OF_INDEX;throw(this);}
    as=AS_NOERROR;
    return  &ArrayPtr[iGet];
}

template <class T>
T&  UArray<T>::GetElement(unsigned int i)
{
 return *GetElementPtr(i);
}

// Set Array Element
template <class T>
void UArray<T>::SetElement(unsigned int iSet, T* TPtr)  
{
#ifdef _SAFE_ACCESS_
    CheckLocker cl(GetLocker());
#endif //_SAFE_ACCESS_ 
    if(NumElements==0  || iSet > NumElements) {as=AS_ERR_OUT_OF_INDEX;throw(this);}
    else memcpy(&ArrayPtr[iSet], TPtr, sizeof(T));
    as=AS_NOERROR;
    return;
}

// Operator '='
template <class T>
UArray<T>& UArray<T>::operator = (UArray<T>& NewT)  
{
#ifdef _SAFE_ACCESS_
    CheckLocker cl1(GetLocker());
    CheckLocker cl2(NewT.GetLocker());
#endif //_SAFE_ACCESS_ 
    CleanArray();
    MaxNumElements = NewT.GetMaxNumElements();
    if(NewT.GetNumElements() > 0)
      {
        ArrayPtr = new T[NewT.NumElements];
        if(!ArrayPtr) {as=AS_ERR_MEM;throw(this);}
        else memcpy(ArrayPtr, NewT.GetElementPtr(0), sizeof(T)*NewT.GetNumElements());
        NumElements = NewT.GetNumElements();
      }
    else NumElements=0;
    as=AS_NOERROR;
    return *this;
}

// Operator '+=' (UArray<T>+T)
template <class T>
UArray<T>& UArray<T>::operator += (T& ET)    
{
#ifdef _SAFE_ACCESS_
    CheckLocker cl(GetLocker());
#endif //_SAFE_ACCESS_ 
    AddElement(&ET);
    return *this;
}

// Operator '+=' (UArray<T>+UArray<T>)
template <class T>
UArray<T>& UArray<T>::operator += (UArray<T>& AT)   
{
    T   *tmpPtr;
#ifdef _SAFE_ACCESS_
    CheckLocker cl1(GetLocker());
    CheckLocker cl2(AT.GetLocker());
#endif //_SAFE_ACCESS_ 
    if(GetNumElements()+AT.GetNumElements() > GetMaxNumElements()) {as=AS_ERR_MAX_SIZE;throw(this);} 
    tmpPtr = new T[GetNumElements() + AT.GetNumElements()];
    if(tmpPtr==NULL) {as=AS_ERR_MAX_SIZE;throw(this);}
    memcpy(tmpPtr,ArrayPtr,sizeof(T)*GetNumElements());
    memcpy(tmpPtr+GetNumElements(),AT.ArrayPtr,sizeof(T)*AT.GetNumElements());
    NumElements = GetNumElements()+AT.GetNumElements();
    delete[] ArrayPtr;
    ArrayPtr = tmpPtr;
    as=AS_NOERROR;
    return *this;
}

// Clear Array
template <class T>
void UArray<T>::CleanArray()
{
#ifdef _SAFE_ACCESS_
    CheckLocker cl(GetLocker());
#endif //_SAFE_ACCESS_ 
    if (NumElements > 0) {
        if (ArrayPtr)
            delete [] ArrayPtr;
        ArrayPtr    = (T*)NULL;
        NumElements = 0;
    }
    as          = AS_NOERROR;
    return;
}
// Get Array Ptr
template <class T>
T* UArray<T>::GetArrayPtr() {
 return ArrayPtr;
}
// Set Max array size
template <class T>
void UArray<T>::SetMaxNumElements(int NumMax) 
{
#ifdef _SAFE_ACCESS_
  CheckLocker cl(GetLocker());
#endif //_SAFE_ACCESS_ 
  if(NumMax < (int)NumElements) {as=AS_ERR_VALUE;throw(this);}
  MaxNumElements=NumMax;
  as=AS_NOERROR;
 }
#endif //  _uarray_hpp_


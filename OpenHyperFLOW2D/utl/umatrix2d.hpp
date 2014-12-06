/****************************************************************
*  UNKNOWN TEMPLATES LIBRARY (UTL) v1.3                         *
*  Copyright (C) 1994-2010 Serge A. Suchkov                     *
*  Copyright policy: LGPL V2.0                                  *
*  Matrix2D template.                                           *
*  Please report all bugs and problems to "sergeas67@gmail.com".*
*****************************************************************/
#ifndef _umatrix2d_hpp_
#define _umatrix2d_hpp_

#include "utl/locker.hpp"
#define UTL_PAGE_SIZE 1024*1024*1024
#undef _SAFE_ACCESS_
// Common 2D/3D types
// Matrix type
enum MatrixType {
    MXT_MEM=1000,   // Need memory allocation
    MXT_MAP,        // Use preallocate memory area
    MXT_PAGE_MEMMAP // Use preallocate paged memory area
};

// Matrix state
enum MatrixState {
    MXS_OK=2000,         // no errors
    MXS_ERR_OUT_OF_INDEX,// attempt out of index range
    MXS_ERR_MEM,         // memory allocation error during operation
    MXS_ERR_MAP,         // memory mapping error
    MXS_ERR_PAGE         // page error
};

ssize_t  _XY(ssize_t x,ssize_t y, ssize_t atf_N);
ssize_t  _YX(ssize_t x,ssize_t y, ssize_t atf_N);
// Store date format
enum MatrixStorageOrder2D {
    MSO_XY=4000,
    MSO_YX,
};

typedef ssize_t (*access_type_func_2d)(ssize_t x,
                                      ssize_t y,
                                      ssize_t n);
template<class T>
class XY
{
    T    X;
    T    Y;
public:
    XY()          {;}
    XY(T fX, T fY) {
        SetXY(fX,fY);
    }
    inline void SetXY(T fX, T fY) {
        X = static_cast<T>(fX);
        Y = static_cast<T>(fY);
    }
    inline XY&    GetXY()                 { return *this;}
    inline void   SetXY(XY<T>* xy)        { memcpy(this,xy,sizeof(XY<T>));}
    inline XY<T>& operator =  (XY<T>& xy) { memcpy(this,&xy,sizeof(XY<T>));return *this;}
    inline int    operator == (XY<T>& xy) { return memcmp(this,&xy,sizeof(XY<T>));}
    inline T      GetX()                  { return X;}
    inline T      GetY()                  { return Y;}
    inline void   SetX(T x)               { X=x;}
    inline void   SetY(T y)               { Y=y;}
};


template <class T>
class UMatrix2D :  public Locker
{
    unsigned int         nX;
    unsigned int         nY;
    T*   Ptr;
    MatrixType           mt;
    MatrixState          ms;
    MatrixStorageOrder2D mso;
    access_type_func_2d  ATF;
    int                  atfN;  

public:
    UMatrix2D(unsigned int,unsigned int,MatrixStorageOrder2D m=MSO_YX);
    UMatrix2D(T*,unsigned int,unsigned int,MatrixStorageOrder2D m=MSO_YX);
    ~UMatrix2D();
    unsigned int  GetX()             {return nX;}
    unsigned int  GetY()             {return nY;}
    MatrixState   GetMatrixState()   {return ms;}
    MatrixType    GetMatrixType()    {return mt;}
    T*            GetMatrixPtr()     {return Ptr;}
    void          SetMatrixPtr(T* x) {Ptr=x;}
    T& operator   () (unsigned int,unsigned int);
    UMatrix2D<T>&   operator  = (UMatrix2D<T>&);
    T&            GetValue(unsigned int,unsigned int);
    ssize_t       GetMatrixSize()    {return nX*nY*sizeof(T);}

    ssize_t  GetRowSize() {
        return nX*sizeof(T);
    }

    ssize_t  GetColSize() {
        return nY*sizeof(T);
    }

};

template <class T>
inline UMatrix2D<T>& UMatrix2D<T>::operator  = (UMatrix2D<T>& M)
                                        {
#ifdef _SAFE_ACCESS_
    CheckLocker cl1(GetLocker());
    CheckLocker cl2(M.GetLocker());
#endif //_SAFE_ACCESS_
    nX  = M.GetX();
    nY  = M.GetY();
    if(mt == MXT_MEM) {
        if(Ptr !=  NULL) 
#ifdef __ICC
            _mm_free(Ptr);
#else
            free(Ptr);
#endif //__ICC        
#ifdef __ICC
            Ptr = (T*)(_mm_malloc(sizeof(T)*nX*nY,_ALIGN));
#else         
            Ptr = (T*)(malloc(sizeof(T)*nX*nY));
#endif         
        memcpy(Ptr,M.GetMatrixPtr(),sizeof(T)*nX*nY);
    } else {
        Ptr = M.GetMatrixPtr();
    }

    ms = M.GetMatrixState();

    return *this;
}

template <class T>
UMatrix2D<T>::UMatrix2D(T* ptr, unsigned int x, unsigned int y,  MatrixStorageOrder2D m) :Locker()
{
#ifdef _SAFE_ACCESS_
    CheckLocker cl(GetLocker());
#endif //_SAFE_ACCESS_
    mt = MXT_MAP;
    mso = m;
    Ptr = ptr;
    nX  = x;
    nY  = y;
    if(Ptr == NULL) {
        ms = MXS_ERR_MAP;throw(this);
    } else ms = MXS_OK;
    if(mso==MSO_XY) {
        ATF  = _XY; 
        atfN = nY;
    } else {
        ATF  = _YX;
        atfN = nX;
    }
}

template <class T>
UMatrix2D<T>::UMatrix2D(unsigned int x, unsigned int y,  MatrixStorageOrder2D m)  :Locker()
{
#ifdef _SAFE_ACCESS_
    CheckLocker cl(GetLocker());
#endif //_SAFE_ACCESS_
    mt = MXT_MEM;
    mso = m;
    Ptr = NULL;
    nX  = x;
    nY  = y;
#ifdef __ICC
            Ptr = (T*)(_mm_malloc(sizeof(T)*nX*nY,_ALIGN));
#else         
            Ptr = (T*)(malloc(sizeof(T)*nX*nY));
#endif         
    if(Ptr == NULL) {
        ms = MXS_ERR_MEM;throw(this);
    } else ms = MXS_OK;
    memset(Ptr,0,sizeof(T)*nX*nY);
    if(mso==MSO_XY) {
        ATF  = _XY; 
        atfN = nY;
    } else {
        ATF  = _YX;
        atfN = nX;
    }
}

template <class T>
UMatrix2D<T>::~UMatrix2D()
{
#ifdef _SAFE_ACCESS_
    CheckLocker cl(GetLocker());
#endif //_SAFE_ACCESS_
    if(mt == MXT_MEM)
        if(Ptr) {
#ifdef __ICC
            _mm_free(Ptr);
#else
            free(Ptr);
#endif //__ICC
        }
    Ptr =  NULL;
}

template <class T>
inline T& UMatrix2D<T>::operator () (unsigned int x,unsigned int y)
{
#ifdef _SAFE_ACCESS_
    CheckLocker cl(GetLocker());
#endif //_SAFE_ACCESS_
    if(x>=nX) {
        ms=MXS_ERR_OUT_OF_INDEX;throw(this);
    } else if(y>=nY) {
        ms=MXS_ERR_OUT_OF_INDEX;throw(this);
    } else ms = MXS_OK;

    if ( mso == MSO_XY )
        return Ptr[y*nX + x];
    else
        return Ptr[x*nY + y];
}

template <class T>
inline T& UMatrix2D<T>::GetValue(unsigned int x,unsigned int y)
{
#ifdef _SAFE_ACCESS_
    CheckLocker cl(GetLocker());
#endif //_SAFE_ACCESS_
    if(x>=nX) {
        ms=MXS_ERR_OUT_OF_INDEX;throw(this);
    } else if(y>=nY) {
        ms=MXS_ERR_OUT_OF_INDEX;throw(this);
    } else ms = MXS_OK;

    if(mso == MSO_XY)
        return Ptr[y*nX + x];
    else
        return Ptr[x*nY + y];
    //  return Ptr[ATF(x,y,atfN)]; 
}
#endif  // _umatrix2d_hpp_

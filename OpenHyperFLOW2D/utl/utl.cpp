/****************************************************************
*  UNKNOWN TEMPLATES LIBRARY (UTL) v1.3                         *
*  Copyright (C) 1994-2010 Serge A. Suchkov                     *
*  Copyright policy: LGPL V2.0                                  *
*  Matrix template.                                             *
*  Please report all bugs and problems to "sergeas67@gmail.com".*
*****************************************************************/
#include "utl/umatrix2d.hpp"
#ifdef _3D
#include "utl/umatrix3d.hpp"
ssize_t _XYZ(ssize_t x,ssize_t y, ssize_t z, ssize_t atf_N, ssize_t atf_M) {
     return (ssize_t)z*(ssize_t)atf_M + (ssize_t)y*(ssize_t)atf_N + (ssize_t)x;
}

ssize_t _ZXY(ssize_t x,ssize_t y, ssize_t z, ssize_t atf_N, ssize_t atf_M) {
     return (ssize_t)y*(ssize_t)atf_M + (ssize_t)x*(ssize_t)atf_N + (ssize_t)z;
}

ssize_t _YZX(ssize_t x,ssize_t y, ssize_t z, ssize_t atf_N, ssize_t atf_M) {
    return (ssize_t)x*(ssize_t)atf_M + (ssize_t)z*(ssize_t)atf_N + (ssize_t)y;
}
#endif // _3D

ssize_t  _XY(ssize_t  x,
         ssize_t  y,
         ssize_t  atf_N) {
 return (ssize_t)x*(ssize_t)atf_N + (ssize_t)y;
}

ssize_t  _YX(ssize_t  x,
         ssize_t  y,
         ssize_t  atf_N) {
 return (ssize_t)y*(ssize_t)atf_N + (ssize_t)x;
}


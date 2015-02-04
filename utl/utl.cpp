/****************************************************************
*  UNKNOWN TEMPLATES LIBRARY (UTL) v1.3                         *
*  Copyright (C) 1994-2010 Serge A. Suchkov                     *
*  Copyright policy: LGPL V2.0                                  *
*  Matrix template.                                             *
*  http://openhyperflow2d.googlecode.com                        *
*  Please report all bugs and problems to "sergeas67@gmail.com".*
*****************************************************************/
#include "utl/umatrix2d.hpp"

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


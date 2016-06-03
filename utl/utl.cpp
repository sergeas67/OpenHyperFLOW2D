/****************************************************************
*  UNKNOWN TEMPLATES LIBRARY (UTL) v1.3                         *
*  Copyright (C) 1994-2010 Serge A. Suchkov                     *
*  Copyright policy: LGPL V2.0                                  *
*  Matrix template.                                             *
*  http://github.com/sergeas67/openhyperflow2d                  *
*  Please report all bugs and problems to "sergeas67@gmail.com".*
*****************************************************************/
#include "utl/umatrix2d.hpp"

ssize_t  _XY(register ssize_t  x,
             register ssize_t  y,
             register ssize_t  atf_N) {
 return (ssize_t)x*(ssize_t)atf_N + (ssize_t)y;
}

ssize_t  _YX(register ssize_t  x,
             register ssize_t  y,
             register ssize_t  atf_N) {
 return (ssize_t)y*(ssize_t)atf_N + (ssize_t)x;
}


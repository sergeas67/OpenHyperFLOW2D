#  UNKNOWN TEMPLATES LIBRARY (UTL) v1.3
#  Copyright (C) 1994-2011 Serge A. Suchkov
#  Copyright policy: LGPL V2.0
#  http://openhyperflow2d.googlecode.com 

include   ../.compiler
include   ../.parallel
include   ../.models

ifeq     ("$(COMPILER)","Intel")
include   ../icc.compiler
else
include   ../gcc.compiler
endif

DEBUGLEVEL   = 0

CFLAGS       = $(PICFLAGS) $(OPTIONS) $(FP_OPTS)

TARGET_2D    = libutl2d.a
SOURCES      = utl.cpp
OBJECTS      = utl.o
INCLUDES_2D  = uarray.hpp ustack.hpp umatrix2d.hpp
TARGET       = $(TARGET_2D)

all: $(TARGET)
	ln -sf `pwd`/*.a ../lib
	rm -f *.o

.SUFFIXES: .cpp
	
.cpp.o:
	@if [  ! -a $(TARGET_2D)  ];\
        then $(CXXC) -c $(INCPATH) $(CXXOPTIONS) $(CFLAGS) $<; \
        ar -r $(TARGET_2D) $(OBJECTS);\
        fi
$(TARGET_2D): $(OBJECTS)

clean:
	rm -f *.o $(TARGET) ../lib/$(TARGET_2D)
 
#  UNKNOWN TEMPLATES LIBRARY (UTL) v1.3
#  Copyright (C) 1994-2011 Serge A. Suchkov
#  Copyright policy: LGPL V2.0
include   ../.compiler
include   ../.parallel
include   ../.models

ifeq     ("$(COMPILER)","Intel")
include   ../icc.compiler
else
include   ../gcc.compiler
ifneq ("$CUDA","")
include  cuda.compiler
CXXC=$(NVCC)
endif
endif

DEBUGLEVEL   = 0

CFLAGS       = $(PICFLAGS) $(OPTIONS)

TARGET_2D    = libutl2d.a
SOURCES      = utl.cu
OBJECTS      = utl.o
INCLUDES_2D  = uarray.cuh ustack.cuh umatrix2d.cuh
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
 
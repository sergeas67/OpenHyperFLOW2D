# DEEPS2D
#   Version  1.0.0                                                             
#   Copyright (C)  1995-2013 by Serge A. Suchkov                               
#   Copyright policy: LGPL V3                                                  

include   ../.version
include   ../.compiler
include   ../.parallel
include   ../.models

ifeq     ("$(COMPILER)","Intel")
include   ../icc.compiler
else
include   ../gcc.compiler
endif


DEBUGLEVEL  = 0

CFLAGS    = $(PICFLAGS) $(OPTIONS) -D_STREAM_COMPAT  \
           -D_REMOVE_SWAPFILE  -D_GNUPLOT_ -D_TECPLOT_ -D_WRITE_LARGE_FILE_ \
           -DNUM_COMPONENTS=3  $(MODELS)

TARGET_LIBS_DEEPS2D    = libDEEPS2D.a
SOURCES_LIBS_DEEPS2D   = deeps2d_core.cpp
OBJECTS_LIBS_DEEPS2D   = deeps2d_core.o
INCLUDES               =
INCPATH                = -I ../

TARGET    +=  $(TARGET_LIBS_DEEPS2D)

LLIBS     += -L ../lib -lflow2d -lobj_data -lutl2d -lexcept -lhf2d

#include ../.make
all:  $(TARGET)
	ln -sf `pwd`/*.a ../lib
.SUFFIXES: .cpp  .c
	
.cpp.o:
	$(CXXC) -c  $(INCPATH) $(CXXOPTIONS) $(CFLAGS) $<
$(TARGET_LIBS_DEEPS2D): $(OBJECTS_LIBS_DEEPS2D)
	ar -r $(TARGET_LIBS_DEEPS2D) $(OBJECTS_LIBS_DEEPS2D)
clean:
	rm -f $(TARGET_LIBS_DEEPS2D) *.o *.a
	rm -f ../lib/$(TARGET_LIBS_DEEPS2D)

rebuild: clean all
	

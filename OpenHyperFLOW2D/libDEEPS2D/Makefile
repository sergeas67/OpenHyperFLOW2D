#   DEEPS2D
#   Version  1.0.2
#   Copyright (C)  1995-2015 by Serge A. Suchkov
#   Copyright policy: LGPL V3
#   http://openhyperflow2d.googlecode.com

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
ASM_LIBS_DEEPS2D       = deeps2d_core.S
INCLUDES               =
INCPATH                = -I ../

TARGET    +=  $(TARGET_LIBS_DEEPS2D)

LLIBS     += -L ../lib -lhf2d -lflow2d -lobj_data -lutl2d -lexcept

#include ../.make
all:  $(TARGET)
	ln -sf `pwd`/*.a ../lib
.SUFFIXES: .cpp  .c  .S
	

.cpp.S:
	$(CXXC) -S -c  $(INCPATH) $(CXXOPTIONS) $(CFLAGS) $<
	

.cpp.o:
	$(CXXC) -c  $(INCPATH) $(CXXOPTIONS) $(CFLAGS) $<
	
$(TARGET_LIBS_DEEPS2D): $(ASM_LIBS_DEEPS2D) $(OBJECTS_LIBS_DEEPS2D)
	ar -r $(TARGET_LIBS_DEEPS2D) $(OBJECTS_LIBS_DEEPS2D)
clean:
	rm -f $(TARGET_LIBS_DEEPS2D) *.o *.a
	rm -f ../lib/$(TARGET_LIBS_DEEPS2D)

rebuild: clean all
	

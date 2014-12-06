# DEEPS2D
#   Version  1.0.2                                                             
#   Copyright (C)  1995-2014 by Serge A. Suchkov                               
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

ifneq ("$CUDA","")
include  cuda.compiler
endif

DEBUGLEVEL  = 0

CFLAGS    +=  -D_STREAM_COMPAT  \
           -D_REMOVE_SWAPFILE  -D_GNUPLOT_ -D_TECPLOT_ -D_WRITE_LARGE_FILE_ \
           -DNUM_COMPONENTS=3 -D_UNIFORM_MESH_ $(MODELS) $(FP_OPTS)

TARGET_LIBS_DEEPS2D    = libDEEPS2D.a
SOURCES_LIBS_DEEPS2D   = deeps2d_core.cu init_deeps2d.cpp cuda_deeps2d_kernels.cu
OBJECTS_LIBS_DEEPS2D   = deeps2d_core.o init_deeps2d.o  cuda_deeps2d_kernels.o
INCLUDES               =
INCPATH                = -I ../

TARGET    +=  $(TARGET_LIBS_DEEPS2D)

LLIBS     += -L ../lib -lflow2d -lobj_data -lutl2d -lexcept -lhf2d

#include ../.make
all:  $(TARGET)
	ln -sf `pwd`/*.a ../lib
.SUFFIXES: .cpp  .cu
	
.cpp.o:
	$(CXXC) -c  $(INCPATH) $(CXXOPTIONS)  $(OPTIONS) $(CFLAGS) $<
$(TARGET_LIBS_DEEPS2D): $(OBJECTS_LIBS_DEEPS2D)
	ar -r $(TARGET_LIBS_DEEPS2D) $(OBJECTS_LIBS_DEEPS2D)
.cu.o:
	$(NVCC) -c  $(INCPATH)  $(CUDA_OPTIONS) $(CFLAGS) $<
$(TARGET_LIBS_DEEPS2D): $(OBJECTS_LIBS_DEEPS2D)
	ar -r $(TARGET_LIBS_DEEPS2D) $(OBJECTS_LIBS_DEEPS2D)
clean:
	rm -f $(TARGET_LIBS_DEEPS2D) *.o *.a
	rm -f ../lib/$(TARGET_LIBS_DEEPS2D)

rebuild: clean all
	

#   OpenHyperFLOW2D                                                            
#   Version  1.0.1                                                             
#   Copyright (C)  1995-2014 by Serge A. Suchkov                               
#   Copyright policy: LGPL V3                                                  

include   ../.compiler
include   ../.parallel
include   ../.models

ifeq     ("$(COMPILER)","Intel")
include   ../icc.compiler
else
include   ../gcc.compiler
endif

DEBUGLEVEL=0

CFLAGS    =  $(PICFLAGS) $(OPTIONS)  \
            -D_REMOVE_SWAPFILE_ -D_GNUPLOT_  -D_STREAM_COMPAT  -Wno-deprecated $(MODELS) $(FP_OPTS)

TARGET    = libOutCFD.a

# File list
SOURCES   = out_cfd_param.cpp
OBJECTS   = out_cfd_param.o
INCLUDES  = out_cfd_param.hpp

all: $(TARGET)
	ln -sf `pwd`/*.a ../lib
	rm -f *.o

.SUFFIXES: .cpp  .c

.cpp.s:
	@rm -f $(TARGET)
	$(CXXC) -S $(INCPATH) $(CFLAGS) $<
.cpp.o:
	@rm -f $(TARGET)
	$(CXXC) -c $(INCPATH) $(CFLAGS)  $<
.c.s:
	@rm -f $(TARGET)
	$(CC) -S $(INCPATH) $(CFLAGS) $<
.c.o:
	@rm -f $(TARGET)
	$(CC) -c $(INCPATH) $(CFLAGS) $<
$(TARGET):$(OBJECTS)
	ar -r $(TARGET) $(OBJECTS)
clean:
	rm -f *.o $(TARGET)
	rm -f ../lib/$(TARGET)

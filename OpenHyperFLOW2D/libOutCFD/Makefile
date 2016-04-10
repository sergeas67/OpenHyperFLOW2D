#   OpenHyperFLOW2D                                                            
#   Version  1.0.3                                                             
#   Copyright (C)  1995-2016 by Serge A. Suchkov                               
#   Copyright policy: LGPL V3                                                  
#   http://openhyperflow2d.googlecode.com

include   ../.compiler
include   ../.parallel
include   ../.models

ifeq     ("$(COMPILER)","Intel")
include   ../icc.compiler
else
include   ../gcc.compiler
endif

DEBUGLEVEL=0
#-D_REMOVE_SWAPFILE_
CFLAGS    =  $(PICFLAGS) $(OPTIONS)  -D_GNUPLOT_  -D_STREAM_COMPAT  -Wno-deprecated $(MODELS) 

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

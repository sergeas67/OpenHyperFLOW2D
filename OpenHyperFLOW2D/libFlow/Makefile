# Libflow v2.2
#  Copyright (C)  1995-2015 by Serge A. Suchkov
#  Copyright policy: LGPL V3
#  http://openhyperflow2d.googlecode.com

include   ../.compiler
include   ../.parallel
include   ../.models


ifeq     ("$(COMPILER)","Intel")
include   ../icc.compiler
else
include   ../gcc.compiler
endif

CFLAGS    = $(PICFLAGS) $(OPTIONS) $(FP_OPTS)   

TARGET_2D = libflow2d.a

# File list
SOURCES_2D  = flow.cpp flow2d.cpp
OBJECTS_2D  = flow.o flow2d.o 
INCLUDES_2D = flow.hpp flow2d.hpp

TARGET = $(TARGET_2D)

all: $(TARGET)
	ln -sf `pwd`/*.a ../lib
.SUFFIXES: .cpp  .c
	
.cpp.s:
	$(CXXC) -S $(INCPATH) $(CFLAGS) $<
.cpp.o:
	$(CXXC) -c $(INCPATH) $(CFLAGS)  $<
.c.s:
	$(CC) -S $(INCPATH) $(CFLAGS) $<
.c.o:
	$(CC) -c $(INCPATH) $(CFLAGS) $<

$(TARGET_2D): $(OBJECTS_2D)
	ar -r $(TARGET_2D) $(OBJECTS_2D)
	rm -f *.o

clean:
	rm -f *.o $(TARGET)
	rm -f ../lib/$(TARGET_2D)

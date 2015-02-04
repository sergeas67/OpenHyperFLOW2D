# Makefile for Exceptlib v0.1beta (System Exceptions C++ Library)
#
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

CFLAGS    = $(OPTIONS) $(PICFLAGS) $(FP_OPTS)

TARGET    = libexcept.a

SOURCES   = except.cpp
OBJECTS   = except.o
INCLUDES  = except.hpp

all: $(TARGET)
	ln -sf `pwd`/*.a ../lib
	rm -f *.o
.SUFFIXES: .cpp  .c
	
.cpp.s:
	rm -f $(TARGET)
	$(CXXC) -S $(INCPATH) $(CXXOPTIONS) $(CFLAGS) $<
.cpp.o:
	rm -f $(TARGET)
	$(CXXC) -c $(INCPATH) $(CXXOPTIONS) $(CFLAGS)  $<
.c.s:
	rm -f $(TARGET)
	$(CC) -S $(INCPATH) $(CFLAGS) $<
.c.o:
	rm -f $(TARGET)
	$(CC) -c $(INCPATH) $(CFLAGS) $<

$(TARGET):$(OBJECTS)
	ar -r $(TARGET) $(OBJECTS)
clean:
	rm -f *.o $(TARGET)
	rm -f ../lib/$(TARGET)
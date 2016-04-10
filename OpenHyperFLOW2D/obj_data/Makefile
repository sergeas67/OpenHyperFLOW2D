# Object data loader (ver.1.5.2)
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
CFLAGS    =  $(PICFLAGS) $(MODELS) $(OPTIONS) -D_STREAM_COMPAT  -Wno-deprecated \
            -DBUFF_SIZE=16384 

TARGET    = libobj_data.a

# File list
SOURCES   = obj_data.cpp
OBJECTS   = obj_data.o
INCLUDES  = obj_data.hpp

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

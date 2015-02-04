# OpenOpenHyperFLOW2D project (Version  1.0.2)
# Written by  Serge A. Suchkov (sergeas67@gmail.com)
# http://openhyperflow2d.googlecode.com

include   ../.version
include   ../.compiler
include   ../.parallel
include   ../.models

# for utils only
#OPENMP=-D_OPEN_MP

ifeq     ("$(COMPILER)","Intel")
include   ../icc.compiler
else
include   ../gcc.compiler
endif

DEBUGLEVEL= 0
CFLAGS    = $(PICFLAGS) $(OPTIONS) -D_STREAM_COMPAT -DUTILS_OPENMP  \
            -D_CLEAN_WALL_NODES_ -DNUM_COMPONENTS=3  $(MODELS) -D_WRITE_LARGE_FILE_  $(FP_OPTS)

TARGET_2D   = libhf2d.a
TARGET      = $(TARGET_2D)

SOURCES_2D  = hyper_flow_bound.cpp \
	      hyper_flow_bound_contour.cpp \
	      hyper_flow_solid_bound_rect.cpp \
	      hyper_flow_bound_circle.cpp \
	      hyper_flow_airfoil.cpp \
	      hyper_flow_area.cpp \
	      hyper_flow_field.cpp \
	      hyper_flow_list.cpp \
	      hyper_flow_source.cpp
OBJECTS_2D  = hyper_flow_bound.o \
	      hyper_flow_bound_contour.o \
	      hyper_flow_solid_bound_rect.o \
	      hyper_flow_bound_circle.o \
	      hyper_flow_airfoil.o \
	      hyper_flow_area.o \
	      hyper_flow_field.o \
	      hyper_flow_list.o \
	      hyper_flow_source.o
INCLUDES_2D = hyper_flow2d.hpp \
	      hyper_flow_node.hpp \
	      hyper_flow_bound.hpp \
	      hyper_flow_bound_contour.hpp \
	      hyper_flow_solid_bound_rect.hpp \
	      hyper_flow_bound_circle.hpp \
	      hyper_flow_airfoil.hpp \
	      hyper_flow_area.hpp \
	      hyper_flow_field.hpp \
	      hyper_flow_list.hpp \
	      hyper_flow_source.hpp
all: $(TARGET)
	ln -sf `pwd`/*.a ../lib
	rm -f *.o
.SUFFIXES: .cpp  .c
	
.cpp.o:
	@if [  ! -a $(TARGET_2D)  ];\
	then echo $(CXXC) -c $(INCPATH) $(CXXOPTIONS) $(CFLAGS) $< -o $@;\
        $(CXXC) -c $(INCPATH) $(CXXOPTIONS) $(CFLAGS) $< -o $@; ar -r $(TARGET_2D) $@;\
	fi ;
	
$(TARGET_2D): $(OBJECTS_2D)
	
clean2d:
	rm -f $(OBJECTS_2D)
clean: clean2d
	rm -f $(TARGET) 
	rm -f ../lib/$(TARGET_2D)

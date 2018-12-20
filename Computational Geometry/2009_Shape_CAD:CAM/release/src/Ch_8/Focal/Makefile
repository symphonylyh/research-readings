# Copyright (C) 2008 Massachusetts Institute of Technology, Cambridge, MA
# All rights reserved
#
# make cplus
# will make focalc and focals

BIN = .
INC = ../../../include
LIB = ../../../lib

FOCALC      = $(BIN)/focalc
FOCALS      = $(BIN)/focals

OBJC = focalc.o
OBJS = focals.o

CFLAGS = $(DEBUG) -I$(INC)
FFLAGS = $(DEBUG)
LIB3 = $(LIB)/libbspl.a $(LIB)/libgen.a $(LIB)/libeditor.a 
LIB2 = $(LIB3)
LIB1 = $(LIB2)
LIBS = $(LIB)/libiges.a $(LIB1)
SYSLB2 = -lm

# this is to link the C++ programs
LIBI = $(LIB)/libIdist.a $(LIB)/libIsimp.a $(LIB)/libInt.a

# this is to compile the C++ programs
.SUFFIXES: .cc
.cc.o:
	CC -c -DUSE_INTERVAL $(DEBUG) -I$(INC) +p +w $<


cplus: $(FOCALC) $(FOCALS)


# link focalc focals
OBJC = focalc.o
OBJS = focals.o

$(FOCALC): $(OBJC) $(LIBS)
	cc -o $(FOCALC) $(OBJC) $(LIBS) $(SYSLB2) 

$(FOCALS): $(OBJS) $(LIBS)
	cc -o $(FOCALS) $(OBJS) $(LIBS) $(SYSLB2) 


clean:
	rm -f *~ *.o \#*

remake:
	rm -f *.o


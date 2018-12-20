# Copyright (C) 2008 Massachusetts Institute of Technology, Cambridge, MA
# All rights reserved
#
# make cplus
# will make orthoc orthos

BIN = .
INC = ../../../include
LIB = ../../../lib

ORTHOC      = $(BIN)/orthoc
ORTHOS      = $(BIN)/orthos

OBJC = orthoc.o
OBJS = orthos.o

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


cplus: $(ORTHOC) $(ORTHOS)


# link orthoc orthos
OBJC = orthoc.o
OBJS = orthos.o

$(ORTHOC): $(OBJC) $(LIBS)
	cc -o $(ORTHOC) $(OBJC) $(LIBS) $(SYSLB2) 

$(ORTHOS): $(OBJS) $(LIBS)
	cc -o $(ORTHOS) $(OBJS) $(LIBS) $(SYSLB2) 


clean:
	rm -f *~ *.o \#*

remake:
	rm -f *.o


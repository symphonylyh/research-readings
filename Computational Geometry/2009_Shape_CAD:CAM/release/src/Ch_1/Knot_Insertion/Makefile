# Copyright (C) 2008 Massachusetts Institute of Technology, Cambridge, MA
# All rights reserved
#
# make cplus
# will make knotc and knots

BIN = .
INC = ../../../include
LIB = ../../../lib

KNOTC      = $(BIN)/knotc
KNOTS      = $(BIN)/knots

OBJC = knotc.o
OBJS = knots.o

CFLAGS = $(DEBUG) -I$(INC)
FFLAGS = $(DEBUG)
LIB3 = $(LIB)/libbspl.a $(LIB)/libgen.a
LIB2 = $(LIB3)
LIB1 = $(LIB2)
SYSLB2 = -lm

# this is to link the C++ programs
LIBI = $(LIB)/libIdist.a $(LIB)/libIsimp.a $(LIB)/libInt.a

# this is to compile the C++ programs
.SUFFIXES: .cc
.cc.o:
	CC -c -DUSE_INTERVAL $(DEBUG) -I$(INC) +p +w $<


cplus: $(KNOTC) $(KNOTS)

# link knotc knots
OBJC = knotc.o
OBJS = knots.o

$(KNOTC): $(OBJC) $(LIB1)
	cc -o $(KNOTC) $(OBJC) $(LIB1) $(SYSLB2)

$(KNOTS): $(OBJS) $(LIB1)
	cc -o $(KNOTS) $(OBJS) $(LIB1) $(SYSLB2)


clean:
	rm -f *~ *.o \#*

remake:
	rm -f *.o


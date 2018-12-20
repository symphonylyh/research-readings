# Copyright (C) 2008 Massachusetts Institute of Technology, Cambridge, MA
# All rights reserved
#
# make cplus
# will make evalc and evals

BIN = .
INC = ../../../include
LIB = ../../../lib

EVALC      = $(BIN)/evalc
EVALS      = $(BIN)/evals


OBJC = evalc.o
OBJS = evals.o


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


cplus: $(EVALC) $(EVALS)


# link evalc evals
OBJC = evalc.o
OBJS = evals.o

$(EVALC): $(OBJC) $(LIB1)
	cc -o $(EVALC) $(OBJC) $(LIB1) $(SYSLB2) 

$(EVALS): $(OBJS) $(LIB1)
	cc -o $(EVALS) $(OBJS) $(LIB1) $(SYSLB2) 

clean:
	rm -f *~ *.o \#*

remake:
	rm -f *.o


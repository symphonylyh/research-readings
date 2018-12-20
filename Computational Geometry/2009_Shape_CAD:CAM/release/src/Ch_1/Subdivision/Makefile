# Copyright (C) 2008 Massachusetts Institute of Technology, Cambridge, MA
# All rights reserved
#
# make cplus
# will make subc and subs

BIN = .
INC = ../../../include
LIB = ../../../lib


SUBC      = $(BIN)/subc
SUBS      = $(BIN)/subs


OBJC = subc.o
OBJS = subs.o


CFLAGS = $(DEBUG) -I$(INC)
FFLAGS = $(DEBUG)
LIB3 = $(LIB)/libbspl.a $(LIB)/libgen.a
LIB2 = $(LIB3)
LIB1 = $(LIB)/libeditor.a $(LIB2)
SYSLB2 = -lm

# this is to link the C++ programs
LIBI = $(LIB)/libIdist.a $(LIB)/libIsimp.a $(LIB)/libInt.a

# this is to compile the C++ programs
.SUFFIXES: .cc
.cc.o:
	CC -c -DUSE_INTERVAL $(DEBUG) -I$(INC) +p +w $<


cplus: $(SUBC) $(SUBS)

# link subc subs
OBJC = subc.o
OBJS = subs.o

$(SUBC): $(OBJC) $(LIB1)
	cc -o $(SUBC) $(OBJC) $(LIB1) $(SYSLB2)

$(SUBS): $(OBJS) $(LIB1)
	cc -o $(SUBS) $(OBJS) $(LIB1) $(SYSLB2)


clean:
	rm -f *~ *.o \#*

remake:
	rm -f *.o


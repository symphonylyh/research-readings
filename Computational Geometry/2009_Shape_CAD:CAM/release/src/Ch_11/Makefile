# Copyright (C) 2008 Massachusetts Institute of Technology, Cambridge, MA
# All rights reserved
#
# make cplus
# will make offsetc and offsets

BIN = .
INC = ../../include
LIB = ../../lib

OFFSETC      = $(BIN)/offsetc
OFFSETS      = $(BIN)/offsets

OBJC = offsetc.o
OBJS = offsets.o

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


cplus: $(OFFSETC) $(OFFSETS)


# link offsetc offsets
OBJC = offsetc.o
OBJS = offsets.o

$(OFFSETC): $(OBJC) $(LIBS)
	cc -o $(OFFSETC) $(OBJC) $(LIBS) $(SYSLB2) 

$(OFFSETS): $(OBJS) $(LIBS)
	cc -o $(OFFSETS) $(OBJS) $(LIBS) $(SYSLB2) 


clean:
	rm -f *~ *.o \#*

remake:
	rm -f *.o


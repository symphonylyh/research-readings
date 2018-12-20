# Copyright (C) 2008 Massachusetts Institute of Technology, Cambridge, MA
# All rights reserved
#
# make cplus
# will make mdist-...

BIN = .
INC = ../../include
LIB = ../../lib

MDISTP2C = $(BIN)/mdist-p2c
MDISTP2S = $(BIN)/mdist-p2s
MDISTC2C = $(BIN)/mdist-c2c
MDISTC2S = $(BIN)/mdist-c2s
MDISTS2S = $(BIN)/mdist-s2s

OBJM = mdist-p2c.o
OBJN = mdist-p2s.o
OBJO = mdist-c2c.o
OBJQ = mdist-c2s.o
OBJP = mdist-s2s.o

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
	g++ -c -DUSE_INTERVAL $(DEBUG) -I$(INC) $<

# originally for linking the C++ stand-alone minimum distance program
cplus: $(MDISTP2C) $(MDISTP2S) $(MDISTC2C) $(MDISTC2S) $(MDISTS2S)


# link
OBJM = mdist-p2c.o
OBJN = mdist-p2s.o
OBJO = mdist-c2c.o
OBJQ = mdist-c2s.o
OBJP = mdist-s2s.o

$(MDISTP2C): $(OBJM) $(LIB1)
	g++ -o $(MDISTP2C) $(OBJM) $(LIBI) $(LIB1) $(SYSLB2) 

$(MDISTP2S): $(OBJN) $(LIB1)
	g++ -o $(MDISTP2S) $(OBJN) $(LIBI) $(LIB1) $(SYSLB2) 

$(MDISTC2C): $(OBJO) $(LIB1)
	g++ -o $(MDISTC2C) $(OBJO) $(LIBI) $(LIB1) $(SYSLB2) 

$(MDISTC2S): $(OBJQ) $(LIB1)
	g++ -o $(MDISTC2S) $(OBJQ) $(LIBI) $(LIB1) $(SYSLB2) 

$(MDISTS2S): $(OBJP) $(LIB1)
	g++ -o $(MDISTS2S) $(OBJP) $(LIBI) $(LIB1) $(SYSLB2) 


clean:
	rm -f *~ *.o \#*

remake:
	rm -f *.o


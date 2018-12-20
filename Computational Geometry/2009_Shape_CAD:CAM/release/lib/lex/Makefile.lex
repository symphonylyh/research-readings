# Copyright (C) Massachusetts Institute of Technology, 1995-2008
# All rights reseved

# Makefile for the lexical library of the projected-polyhedron
# nonlinear system solver

ROOT = ../..
EINC = $(ROOT)/include
LEX = ../lib$(FI)lex.a

.SUFFIXES: .cc
.cc.o:
	g++ -c $(DUSE) $(DEBUG) -I$(EINC) $<

OBJS = lexConv.o monoToBern.o bernToMono.o consolidate.o lex.pp.o 

$(LEX): Makefile $(OBJS)
	rm -f $(LEX)
	ar ruv $@ $(OBJS)

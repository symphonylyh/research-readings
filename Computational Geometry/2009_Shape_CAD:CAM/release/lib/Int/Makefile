# Copyright (C) Massachusetts Institute of Technology, 1995-2008
# All rights reseved

# Makefile for the lexical library of the projected-polyhedron
# nonlinear system solver
# make both the floating point and interval versions


all: ria 

ria:
	touch *.cc
	make -f Makefile.int DEBUG=$(DEBUG) FI=I DUSE=-DUSE_INTERVAL

# force complete recompile and relink
remake:
	touch *.cc

# delete backups (*~)
clean:
	rm -f *~ \#*

# delete backups (*~) and objects (*.o)
empty:
	rm -f *.o *~ \#*

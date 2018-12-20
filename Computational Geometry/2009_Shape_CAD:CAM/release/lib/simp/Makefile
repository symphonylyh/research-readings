# Copyright (C) Massachusetts Institute of Technology, 1995-2008
# All rights reserved

# projected-polyhedron nonlinear equation system solver

all: fpa ria


# compile and link the floating point version
fpa:
	touch *.cc
	rm -f *.o
	make -f Makefile.simp DEBUG=$(DEBUG) FI=f \
		DUSE= SEPARATE= WRAP_SOLVER=wrap_solver-no-separate.o

# compile and link the interval version
ria:
	touch *.cc
	rm -f *.o
	make -f Makefile.simp DEBUG=$(DEBUG) FI=I \
		DUSE=-DUSE_INTERVAL SEPARATE=separate_root.o WRAP_SOLVER=wrap_solver-no-separate.o

# force recompile and relink
remake:
	rm -f *.o

# delete backups (*~) and objects (*.o)
clean:
	rm -f *.o *~ \#*

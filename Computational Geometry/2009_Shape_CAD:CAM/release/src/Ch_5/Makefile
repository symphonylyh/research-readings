# Copyright (C) 2008 Massachusetts Institute of Technology, Cambridge, MA
# All rights reserved

INCDIR =../../include
LIBDIR =../../lib

LIBS = -lm

all: ria


ria:	
	touch *.cc
	rm -f *.o
	make -f Makefile.ria
	

# force recompile and relink
remake:
	rm -f *.o

# delete backups (*~) and objects (*.o)
clean:
	rm -f *.o *~ \#*
	rm solbern-ria solpow-ria 





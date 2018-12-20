# Copyright (C) 2008 Massachusetts Institute of Technology, Cambridge, MA
# All rights reserved

INCDIR =../../include
LIBDIR =../../lib

LIBS = -lfl -lm

all: fpa ria


fpa:
	touch *.cc
	rm -f *.o
	make -f Makefile.float

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
	rm sollex-fpa sollex-ria solbern-fpa solbern-ria \
		solpow-fpa solpow-ria 





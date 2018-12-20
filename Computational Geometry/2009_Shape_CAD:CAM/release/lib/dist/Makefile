# Copyright (C) Massachusetts Institute of Technology, 1995-2008
# All rights reserved

# robust distance library

all: fpa ria

# compile and link the floating point version
fpa:
	touch *.cc
	rm -f *.o
	make -f Makefile.dist DEBUG=$(DEBUG) FI=f \
		DUSE=

# compile and link the interval version
ria:
	touch *.cc
	rm -f *.o
	make -f Makefile.dist DEBUG=$(DEBUG) FI=I \
		DUSE=-DUSE_INTERVAL

# stand-alone linking of interval distance library
dist:
	make -f Makefile.dist DEBUG=-g FI=I DUSE=-DUSE_INTERVAL

# force recompile and relink
remake:
	rm -f *.o

# delete backups (*~) and objects (*.o)
clean:
	rm -f *.o *~ \#*

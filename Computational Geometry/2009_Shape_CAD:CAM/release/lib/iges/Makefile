# Copyright (C) 1995-2008 Massachusetts Institute of Technology, Cambridge, MA
# All rights reserved

INCL = ../../include
CFLAGS = $(DEBUG) -I$(INCL)
IGES = ../libiges.a

OBJS = check.o copious.o decode.o encode.o pgeom.o powcurv.o \
	powsurf.o print.o read.o sgeom.o tmatrix.o utils.o \
	write.o

$(IGES): Makefile $(OBJS)
	rm -f $(IGES)
	ar ruv $@ $(OBJS)

# force recompile and relink
remake:
	rm -f *.o

# delete backups (*~) and objects (*.o)
clean:
	rm -f *.o *~ \#*

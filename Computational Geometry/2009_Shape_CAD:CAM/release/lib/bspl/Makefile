# Copyright (C) 1995-2008 Massachusetts Institute of Technology, Cambridge, MA
# All rights reserved

INCL = ../../include 
CFLAGS = $(DEBUG) -I$(INCL)
FFLAGS = $(DEBUG)
BSPL = ../libbspl.a

OBJS = bernrs.o berntomono.o boehm.o boehms.o boundbox.o bsplutil.o \
	convexbox.o ctrans.o curvature.o cxdeb.o evalbsp.o \
	evalsurf.o evalsurf1.o evalsurf_per.o find.o hodog.o merge.o \
	miscel.o monotobern.o oslo1.o oslo3.o percurv.o persurf.o \
	prim.o raise.o setop.o solve.o subbez.o subdiv.o subdivs.o \
	util.o 

$(BSPL): Makefile $(OBJS)
	rm -f $(BSPL)
	ar ruv $@ $(OBJS)

# force recompile and relink
remake:
	rm -f *.o

# delete backups (*~) and objects (*.o)
clean:
	rm -f *.o *~ \#*

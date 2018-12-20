# Copyright (C) 1995-2008 Massachusetts Institute of Technology, Cambridge, MA
# All rights reserved 

INCL = ../../include
CFLAGS = $(DEBUG) -I$(INCL)
GEN = ../libgen.a

OBJS = alloc.o avec02b.o bit_array.o copy.o dbl_array.o dispatch.o \
	distance.o dynamic_memory.o egeom.o error.o fgeom.o fileopen.o \
	flt_array.o gen_array.o int_array.o matrix.o nag_funs.o \
	ptr_array.o read.o readper.o sht_array.o string.o \
	vec_array.o vect_trans.o vectarith.o vectarith2.o write.o \
	writeper.o

$(GEN): Makefile $(OBJS)
	rm -f $(GEN)
	ar ruv $@ $(OBJS)

# force recompile and relink
remake:
	rm -f *.o

# delete backups (*~) and objects (*.o)
clean:
	rm -f *.o *~ \#*

# Copyright (C) 1995-2008 Massachusetts Institute of Technology, Cambridge, MA
# All rights reserved

# Makefile for lib

all: libs

# link all C libraries
libs:
	cd bspl;   make -k DEBUG=$(DEBUG)
	cd editor; make -k DEBUG=$(DEBUG)
	cd gen;    make -k DEBUG=$(DEBUG)
	cd iges;   make -k DEBUG=$(DEBUG)

# link all C++ libraries
cplus:
	cd dist;   make -k DEBUG=$(DEBUG) all
	cd Int;    make -k DEBUG=$(DEBUG)
	cd simp;   make -k DEBUG=$(DEBUG) all
	cd lex;    make -k DEBUG=$(DEBUG) all

# force complete recompile relink of libraries
remake:
	cd bspl;   make -k remake
	cd dist;   make -k remake
	cd editor; make -k remake
	cd gen;    make -k remake
	cd iges;   make -k remake
	cd Int;    make -k remake
	cd simp;   make -k remake
	cd lex;    make -k remake

# delete backups (*~) and objects (*.o)
clean:
	cd bspl;   make -k clean
	cd dist;   make -k clean
	cd editor; make -k clean
	cd gen;    make -k clean
	cd iges;   make -k clean
	cd Int;    make -k clean
	cd simp;   make -k clean
	cd lex;    make -k clean
	rm -f *~ \#*

# delete library archives (lib*.a)
nolib:
	rm lib*.a

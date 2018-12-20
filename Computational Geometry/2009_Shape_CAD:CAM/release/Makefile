# Copyright (C) 2008 Massachusetts Institute of Technology, Cambridge, MA
# All rights reserved

# Makefile for hyperbook codes

COMPILER = cplus

# directories
ROOT = "/usr/local"
BIN = $(ROOT)/bin
HOME = "/usr/hyperbook/release"

SRC = $(HOME)/src
CH1 = $(SRC)/Ch_1
CH2 = $(SRC)/Ch_2
CH3 = $(SRC)/Ch_3
EVALUATION = $(CH1)/Evaluation
SUBDIVISION = $(CH1)/Subdivision
KNOT_INSERTION = $(CH1)/Knot_Insertion
CH4 = $(SRC)/Ch_4
CH5 = $(SRC)/Ch_5
CH7 = $(SRC)/Ch_7
FOCAL = $(SRC)/Ch_8/Focal
RADIAL = $(SRC)/Ch_8/Radial
ORTHOTOMICS = $(SRC)/Ch_8/Orthotomics
CH9 = $(SRC)/Ch_9
CH10 = $(SRC)/Ch_10
CH11 = $(SRC)/Ch_11


# debug flag for hyperbook source
SRCDBG =
#SRCDBG = -g

# debug flag for libraries
LIBDBG =
#LIBDBG = -g


# soft link 
LN = ln

all: $(COMPILER)


# compile and link the C++ versions
cplus: ebook_cplus


# update libraries and hyperbook executables
ebook_cplus:
	cd lib; make -k PARENT=$(HOME) DEBUG=$(LIBDBG)
	cd lib; make -k PARENT=$(HOME) DEBUG=$(LIBDBG) cplus
	cd src; make -k PARENT=$(HOME) DEBUG=$(SRCDBG) cplus
	cd $(BIN); rm -f evalc evals subc subs knotc knots geomc geoms \
	                 solbern-fpa solbern-ria sollex-fpa sollex-ria solpow-fpa solpow-ria \
				ex.5.6 ex.5.9 \
				mdist-p2c mdist-p2s mdist-c2c mdist-c2s mdist-s2s \
				focalc radialc orthoc focals orthos \
				loc geod \
				offsetc offsets
	

	$(LN) -s $(EVALUATION)/evalc $(BIN)/evalc
	$(LN) -s $(EVALUATION)/evals $(BIN)/evals
	$(LN) -s $(SUBDIVISION)/subc $(BIN)/subc
	$(LN) -s $(SUBDIVISION)/subs $(BIN)/subs
	$(LN) -s $(KNOT_INSERTION)/knotc $(BIN)/knotc
	$(LN) -s $(KNOT_INSERTION)/knots $(BIN)/knots
	$(LN) -s $(CH2)/geomc $(BIN)/geomc
	$(LN) -s $(CH3)/geoms $(BIN)/geoms
	$(LN) -s $(CH4)/solbern-fpa $(BIN)/solbern-fpa
	$(LN) -s $(CH4)/solbern-ria $(BIN)/solbern-ria
	$(LN) -s $(CH4)/sollex-fpa $(BIN)/sollex-fpa
	$(LN) -s $(CH4)/sollex-ria $(BIN)/sollex-ria
	$(LN) -s $(CH4)/solpow-fpa $(BIN)/solpow-fpa
	$(LN) -s $(CH4)/solpow-ria $(BIN)/solpow-ria
	$(LN) -s $(CH5)/ex.5.6 $(BIN)/ex.5.6
	$(LN) -s $(CH5)/ex.5.9 $(BIN)/ex.5.9
	$(LN) -s $(CH7)/mdist-p2c $(BIN)/mdist-p2c
	$(LN) -s $(CH7)/mdist-p2s $(BIN)/mdist-p2s
	$(LN) -s $(CH7)/mdist-c2c $(BIN)/mdist-c2c
	$(LN) -s $(CH7)/mdist-c2s $(BIN)/mdist-c2s
	$(LN) -s $(CH7)/mdist-s2s $(BIN)/mdist-s2s
	$(LN) -s $(FOCAL)/focalc $(BIN)/focalc
	$(LN) -s $(FOCAL)/focals $(BIN)/focals
	$(LN) -s $(RADIAL)/radialc $(BIN)/radialc
	$(LN) -s $(ORTHOTOMICS)/orthoc $(BIN)/orthoc
	$(LN) -s $(ORTHOTOMICS)/orthos $(BIN)/orthos
	$(LN) -s $(CH9)/loc $(BIN)/loc
	$(LN) -s $(CH10)/geod $(BIN)/geod
	$(LN) -s $(CH11)/offsetc $(BIN)/offsetc
	$(LN) -s $(CH11)/offsets $(BIN)/offsets


# force complete recompile and relink of libraries and hyperbook executables
remake:
	cd lib; make -k remake
	cd src; make -k remake

# delete all backups (*~) and object files (*.o)
clean:
	cd lib; make -k clean
	cd src; make -k clean
	rm *~


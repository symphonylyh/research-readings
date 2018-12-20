# Copyright (C) 2008 Massachusetts Institute of Technology, Cambridge, MA
# All rights reserved

# Makefile for hyperbook executables

all: cplus

cplus:
	cd ./Ch_1/Evaluation;   make -k DEBUG=$(DEBUG) cplus
	cd ./Ch_1/Subdivision;   make -k DEBUG=$(DEBUG) cplus
	cd ./Ch_1/Knot_Insertion;   make -k DEBUG=$(DEBUG) cplus
	cd ./Ch_2;   make -k DEBUG=$(DEBUG) cplus
	cd ./Ch_3;   make -k DEBUG=$(DEBUG) cplus
	cd ./Ch_4;   make -k DEBUG=$(DEBUG) 
	cd ./Ch_5;   make -k DEBUG=$(DEBUG) 
	cd ./Ch_7;   make -k DEBUG=$(DEBUG) 
	cd ./Ch_8/Focal;   make -k DEBUG=$(DEBUG) cplus
	cd ./Ch_8/Radial;   make -k DEBUG=$(DEBUG) cplus
	cd ./Ch_8/Orthotomics;   make -k DEBUG=$(DEBUG) cplus
	cd ./Ch_9;   make -k DEBUG=$(DEBUG)
	cd ./Ch_10;   make -k DEBUG=$(DEBUG)
	cd ./Ch_11;   make -k DEBUG=$(DEBUG) 
      
# force complete recompile relink of libraries
remake:
	cd ./Ch_1/Evaluation;   make -k remake
	cd ./Ch_1/Subdivision;   make -k remake
	cd ./Ch_1/Knot_Insertion;   make -k remake
	cd ./Ch_2;   make -k remake
	cd ./Ch_3;   make -k remake
	cd ./Ch_4;   make -k remake
	cd ./Ch_5;   make -k remake
	cd ./Ch_7;   make -k remake
	cd ./Ch_8/Focal;   make -k remake
	cd ./Ch_8/Radial;   make -k remake
	cd ./Ch_8/Orthotomics;   make -k remake
	cd ./Ch_9;   make -k remake
	cd ./Ch_10;   make -k remake
	cd ./Ch_11;   make -k remake
	
# delete backups (*~) and objects (*.o)
clean:
	cd ./Ch_1/Evaluation;   make -k clean
	cd ./Ch_1/Subdivision;   make -k clean
	cd ./Ch_1/Knot_Insertion;   make -k clean
	cd ./Ch_2;   make -k clean
	cd ./Ch_3;   make -k clean
	cd ./Ch_4;   make -k clean
	cd ./Ch_5;   make -k clean
	cd ./Ch_7;   make -k clean
	cd ./Ch_8/Focal;   make -k clean
	cd ./Ch_8/Radial;   make -k clean
	cd ./Ch_8/Orthotomics;   make -k clean
	cd ./Ch_9;   make -k clean
	cd ./Ch_10;   make -k clean
	cd ./Ch_11;   make -k clean
	rm -f *~ \#*


/****************************************************************
 Copyright (C) 2008 Massachusetts Institute of Technology, Cambridge, MA
 All rights reserved

					     knotc.c

  Insert new knots into the knot vectors of an open NURBS curve at u 

  make knotc
  knotc -i input_curve_file_name -u u_value -k #_of_knots_to_be_inserted [-o output_curve_file_name]

  e.g.: Example 1.5.1 in pages 30--31 
        (insert 3 knots at theta = 45 deg, i.e. at u = 0.414213562)
        knotc -i c.curv  -u 0.414213562  -k 3  -o c1.curv
  ***************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include "gen.h"
#include "bspl.h"

void print_usage()
{
	fprintf(stderr, "usage: knotc -i input_curve_file_name -u u_value  -k #_of_knots_to_be_inserted  [-o output_curve_file_name]\n");

	return;
}


int main(int argc, char* argv[])
{
	ParCurv *egeom; /* input and output curve */
	double **inpoly, **outpoly, *t; /* control pts of an input curve and a new set of control pts
	                                   and new augmented knot vector */ 
	double param; /* insert knots at param */
	short i, j, nk, nKnots, order; /* nKnots is the # of knots to be inserted */

	FILE *fi, *fo;  
	char c, *input, *para, *num_k, *output;
	double u_min, u_max;

	int errflg;

	input = para = num_k = output = 0;
	errflg = 0;

	if (argc > 1) {
		while ((c = getopt(argc, argv, "i:u:k:o:")) != -1)
			switch (c) {
			case 'i':
				input = optarg;
				break;

			case 'u':
				para = optarg;
				break;

			case 'k':
				num_k = optarg;
				break;

			case 'o':
				output = optarg;
				break;

			default:
				errflg = 1;
		}
	}
	else
		errflg = 1;

	if (errflg) {
		print_usage();
		exit(-1);
	}


	if(fi = fopen(input, "r")) { 
		egeom = ReadParCurv(fi, NULL); /* read an input curve */
		fclose(fi);

		if(para) { /* read a parameter value where the knots are inserted */
			param = atof(para); /* string to double */
			u_min = egeom->knots[0];
			u_max = egeom->knots[egeom->ncontpts + egeom->order -1];
			if(param < u_min || param > u_max) { /* check param range */
				fprintf(stderr, "Paramer value should be within [%lf,%lf]\n", u_min, u_max);
				print_usage();
			    exit(-1);
			}
			
			nk = egeom->ncontpts + egeom->order; /* # of knots */
	
			for (i=j=0; i<nk; i++)
				if(egeom->knots[i] == param)
					j++;

		    if (j == egeom->order) {
				fprintf(stderr, "Already %d knots at %g", j, param);
				exit(-1);
			}

			
          if(num_k) { /* read the # of knots to be inserted */
		    	nKnots = atoi(num_k); /* string to int */
			    if(nKnots < 0 || nKnots > egeom->order - j) { /* check nKnots range [0, order-j] */
				   fprintf(stderr, "# of knots to be inserted should be within [0,%d]\n", egeom->order-j);
				   print_usage();
			       exit(-1);
				}

			inpoly = dbl_array2(egeom->ncontpts, 4);
	
			for(i=0; i<egeom->ncontpts; i++) {
				inpoly[i][0] = egeom->contpts[i]->x;
				inpoly[i][1] = egeom->contpts[i]->y;
				inpoly[i][2] = egeom->contpts[i]->z;
				inpoly[i][3] = egeom->contpts[i]->w;
			}

			nk = egeom->ncontpts + egeom->order + nKnots;
			t = dbl_array1(nk);
			t[0] = egeom->knots[0];
			for (i=1,j=0; i < egeom->ncontpts + egeom->order; i++) {
				if (egeom->knots[i] >= param && egeom->knots[i-1] < param)
					for(j=0; j<nKnots; j++)
						t[i+j] = param;
				t[i+j] = egeom->knots[i];
			}
			outpoly = dbl_array2(nk-egeom->order, 4);
			curve_oslo1(egeom->order, egeom->ncontpts, nk,
					egeom->knots, t, inpoly, outpoly);
	
			order = egeom->order;
			free_egeom(egeom);
			egeom = egeomalloc1(order, nk - order);

			for(i=0; i < nk - egeom->order; i++) {
				egeom->contpts[i]->x = outpoly[i][0];
				egeom->contpts[i]->y = outpoly[i][1];
				egeom->contpts[i]->z = outpoly[i][2];
				egeom->contpts[i]->w = outpoly[i][3];
			}

			for (i=0; i<nk; i++)
				egeom->knots[i] = t[i];

			free_darray2(inpoly);
			free_darray2(outpoly);
			free_darray1(t);

		  }/* if(num_k) */
		  else {
			fprintf(stderr, "Can't read the # of knots to be inserted\n");
			print_usage();
			exit(-1);
		  }

		} /* if(para) { */
		else {
			fprintf(stderr, "Can't read a parametric point param\n");
			print_usage();
			exit(-1);
		}

		if(fo = fopen(output, "w")) {
			WriteParCurv(fo, egeom);
			fclose(fo);
		}

		WriteParCurv(stdout, egeom);
	
		free_egeom(egeom);

	}/* if(fi = fopen(input, "r")) { */
	else {
		fprintf(stderr, "Can't open a file %s\n", input);
		print_usage();
		exit(-1);
	}

	return 0;
}


	

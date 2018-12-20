/****************************************************************
 Copyright (C) 2008 Massachusetts Institute of Technology, Cambridge, MA
 All rights reserved

					     subc.c

  Subdivide an open NURBS curve at u

  make subc
  subc -i input_curve_file  -u u_value  [-o output_curve_file_1,output_curve_file_2]

  Note: No spaces before or after the "comma" in -o command line option.

  e.g.: Example 1.5.1 in pages 30--31 
        (subdivide the curve at theta = 45 deg, i.e. at u = 0.414213562)
		subc -i c.curv  -u 0.414213562  -o c1.curv,c2.curv
  ***************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include "gen.h"
#include "bspl.h"

void print_usage()
{
	fprintf(stderr, "usage: subc -i input_curve_file  -u u_value  [-o output_curve_file_1,output_curve_file_2]\n");

	return;
}


int main(int argc, char* argv[])
{
	ParCurv *egeom, *egeom1, *egeom2; /* input curve and its subdivided curves */
	double **inpoly, **outpoly, *t; /* control pts of an input curve and a new set of control pts
	                                   and new augmented knot vector */ 
	double param; /* subdivide the input curve at param */
	short i, indx = 1, j, k, nk, nKnots, order;

	FILE *fi, *fo1, *fo2;  
	char c, *input, *para, output1[80], output2[80];
	double u_min, u_max;

	input = para = 0;
	short errflg = 0;
	output1[0] = output2[0] = '\0';

	if (argc > 1) {
		while ((c = getopt(argc, argv, "i:u:o:")) != -1)
			switch (c) {
			case 'i':
				input = optarg;
				break;

			case 'u':
				para = optarg;
				break;

			case 'o':
				for (i=0; i<strlen(optarg); i++)
					if (optarg[i] == ',')
						optarg[i] = ' ';
				sscanf(optarg, "%s %s", output1, output2);
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

		if(para) { /* read a parameter value for subdivision */
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

			nKnots = egeom->order;
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


			for (i=egeom->order;  i<egeom->ncontpts; i++) {
				/* allow the curve to have multiple internal knots, other than at param */
				if (egeom->knots[i] == param && egeom->knots[i] == egeom->knots[i+1]) {
					indx = indx + 1;
					if (indx == egeom->order) {
						indx = i + 2;

						egeom1 = egeomalloc1(egeom->order, indx-egeom->order);
						egeom2 = egeomalloc1(egeom->order, egeom->ncontpts-indx+egeom->order);
	      
						for (j=0; j<indx; j++)
							egeom1->knots[j] = egeom->knots[j];
						k=0;
						for (j=indx-egeom->order; j<egeom->order+egeom->ncontpts; j++) {
							egeom2->knots[k] = egeom->knots[j];
							k++;
						}
						for (j=0; j < indx - egeom->order; j++)
							copyvector(egeom->contpts[j], egeom1->contpts[j]);
						k = egeom->ncontpts - indx + egeom->order;
						for(j = egeom->ncontpts - 1; j >= indx - egeom->order; j--) {
							k--;
							copyvector(egeom->contpts[j], egeom2->contpts[k]);
						}
				
						knot_normalize(egeom1->knots, egeom1->kmem);
						knot_normalize(egeom2->knots, egeom2->kmem);

					}
				}
			}
		} /* if(para) { */
		else {
			fprintf(stderr, "Can't read a parametric point u\n");
			print_usage();
			exit(-1);
		}

		free_egeom(egeom);

	}/* if(fi = fopen(input, "r")) { */
	else {
		fprintf(stderr, "Can't open a file %s\n", input);
		print_usage();
		exit(-1);
	}


	if(fo1 = fopen(output1, "w")) {
		WriteParCurv(fo1, egeom1);
		fclose(fo1);
	}

	if(fo2 = fopen(output2, "w")) {
		WriteParCurv(fo2, egeom2);
		fclose(fo2);
	}

	WriteParCurv(stdout, egeom1);
	WriteParCurv(stdout, egeom2);

	free_egeom(egeom1);
	free_egeom(egeom2);

	return 0;
}


	

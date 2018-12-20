/****************************************************************
 Copyright (C) 2008 Massachusetts Institute of Technology, Cambridge, MA
 All rights reserved

					     evalc.c

  Evaluate i-th derivative of an open NURBS curve 
  at a parametric point u where i = 0,1,2...
  So, if i = 0, then it will evaluate the 0-th derivative, i.e. a position at u.

  make evalc
  evalc -i input_curve_file_name -u u_value -d #_of_derivative [-o output_file_name]

  e.g.: Example 1.5.1 in pages 30--31 
        (position vector (i.e. i=0) for theta=45 deg i.e. at u=0.414213562)
        evalc  -i c.curv  -u 0.414213562  -d 0  -o c.out
  ***************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include "gen.h"
#include "bspl.h"

void print_usage()
{
	fprintf(stderr, "usage: evalc -i input_curve_file_name -u u_value -d #_of_derivative [-o output_file_name]\n");

	return;
}


int main(int argc, char* argv[])
{
	FILE *fi, *fo;  
	ParCurv* cv;
	double u;
	char c, *input, *param, *deriv, *output;
	vector* pt;
	double w;
	int ideriv, errflg;
	double u_min, u_max;

	input = param = deriv = output = 0;
	errflg = 0;

	if (argc > 1) {
		while ((c = getopt(argc, argv, "i:u:d:o:")) != -1)
			switch (c) {
			case 'i':
				input = optarg;
				break;

			case 'u':
				param = optarg;
				break;

			case 'd':
				deriv = optarg;
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


	if(fi = fopen(input, "r")) { /* input is the file name storing the input curve */
		cv = ReadParCurv(fi, NULL);

		if(param) { /* read a parameter value */
			u = atof(param); /* convert a string to double */
			u_min = cv->knots[0];
			u_max = cv->knots[cv->ncontpts + cv->order -1];
			if(u < u_min || u > u_max) { /* check param range */
				fprintf(stderr, "Paramer u_value should be within [%lf,%lf]\n", u_min, u_max);
				print_usage();
			    exit(-1);
			}
		 if(deriv) { /* read which derivative */
		   if((ideriv = atoi(deriv)) >= 0) { /* ideriv = 0,1,2,... */
			   fprintf(stdout, "At u = %+.16le w/ # of derivative = %d:\n", u, ideriv);
			if(pt = rbspeval(cv, u, ideriv)) { /* evaluate a curve */
				w = pt->w;
				if(output) {
					if(fo = fopen(output, "w")) {
						fprintf(fo, "At u = %+.16le w/ # of derivative = %d:\n", u, ideriv);
						fprintf(fo, "In homogeneous coordinates:\n");
						Writevector(fo, pt); /* write a result in homogeneous coords */
						fprintf(fo, "So, in 3D coordinates:\n"); /* and in 3D coords */
						fprintf(fo,"%+.16le %+.16le %+.16le \n",pt->x/w, pt->y/w, pt->z/w);

						fclose(fo);
					}
				}
				
				fprintf(stdout, "In homogeneous coordinates:\n");
				Writevector(stdout, pt);
                fprintf(stdout, "So, in 3D coordinates:\n");
			    fprintf(stdout,"%+.16le %+.16le %+.16le \n",pt->x/w, pt->y/w, pt->z/w);
				
				vectfree(pt);
			} /* if(pt = rbspeval(cv, u, 0)) { */
		   } /* if((ideriv = atoi(deriv)) >= 0) { */
		   else { 
             fprintf(stderr, "#_of_derivative should be a non-negative integer\n");
			 print_usage();
			 exit(-1);
		   }
          } /* if(deriv) { */
		  else {
		   fprintf(stderr, "Can't read #_of_derivative\n");
		   print_usage();
		   exit(-1);
		  }
		} /* if(param) { */
		else {
			fprintf(stderr, "Can't read parameter u_value\n");
			print_usage();
			exit(-1);
		}

		fclose(fi);
		free_egeom(cv);
	}

	else {
		fprintf(stderr, "Can't open a file %s\n", input);
		print_usage();
		exit(-1);
	}


	return 0;

}

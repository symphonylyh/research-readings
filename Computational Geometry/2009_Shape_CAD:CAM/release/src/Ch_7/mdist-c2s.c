/****************************************************************
 Copyright (C) 2008 Massachusetts Institute of Technology, Cambridge, MA
 All rights reserved

					     mdist-c2s.c

  Compute the minimum distance from a NURBS curve to a NURBS surface
  using the IPP solver

  How to make: 
  prompt> make mdist-c2s
  How to run:
  prompt> mdist-c2s -f curve_file_name -t surface_file_name -e root_tolerance_for_IPP_solver  [-o output_file_name]

  e.g.: mdist-c2s  -f ex.7.4.curv -t ex.7.4.surf  -e 1.e-8  -o ex.7.4.out

  ***************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include "gen.h"
#include "bspl.h"


void print_usage()
{
	fprintf(stderr, "usage: mdist-c2s -f curve_file_name -t surface_file_name -e root_tolerance_for_IPP_solver  [-o output_file_name]\n");

	return;
}


int main(int argc, char* argv[])
{
	FILE *cv_fi, *sf_fi, *fo;  
	ParCurv *cv;	ParSurf *sf;
	
	char c, *cv_input, *sf_input, *root_tol, *output;
	vector *pt1, *pt2;
	double w1, w2;
	int errflg;
	double tol, u1, u2,v2, d0;

	cv_input = sf_input = root_tol = output = 0;
	errflg = 0;

	if (argc > 1) {
		while ((c = getopt(argc, argv, "f:t:e:o:")) != -1)
			switch (c) {
			case 'f':
				cv_input = optarg;
				break;

			case 't':
				sf_input = optarg;
				break;

			case 'e':
				root_tol = optarg;
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


	if(cv_fi = fopen(cv_input, "r")) { /* cv_input is the file name storing the input curve */
		cv = ReadParCurv(cv_fi, NULL);

		if(sf_fi = fopen(sf_input, "r")) { /* sf_input is the file name storing the input surface */
			sf = ReadParSurf(sf_fi, NULL);

			if(root_tol) { /* read a root tolerance for the IPP solver */
				tol = atof(root_tol); /* convert a string to double */
				if(tol <= 0) { /* check tol range */
					fprintf(stderr, "Root tolerance for the IPP solver, %lf should be a positive number\n", tol);
					print_usage();
					exit(-1);
				}

				/* Compute the minimum distance from the input curve cv to the input surface sf */
				/* The minimum distance is d0 and the corresponding parameter values for cv, sf are u1, (u2,v2) */
				curvToSurf_n(cv, sf, tol, &u1, &u2,&v2, &d0);
		 	 
				if(output) { /* optional output */
					if(fo = fopen(output, "w")) {
						fprintf(fo, "* Minimum distance from the input NURBS curve:\n\n");
						WriteParCurv(fo, cv);
						fprintf(fo, "\n  to the NURBS surface:\n\n");
						WriteParSurf(fo, sf);
						fprintf(fo, "\n");
						fprintf(fo, "* (Signed) minimum distance = %lf is found\n", d0);
						fprintf(fo, "  from the input curve at the parameter value u1 = %lf\n", u1);
						fprintf(fo, "  to the input surface at the parameter values (u2,v2) = (%lf,%lf), and\n", u2,v2);
						pt1 = rbspeval(cv, u1, 0);		pt2 = revalderivsurf(sf, u2,v2, 0, 0);
						w1 = pt1->w;					w2 = pt2->w;
						fprintf(fo, "  the corresponding points on the curve and surface are (%lf, %lf, %lf),\n", pt1->x/w1, pt1->y/w1, pt1->z/w1);
						fprintf(fo, "  and (%lf, %lf, %lf), respectively.\n", pt2->x/w2, pt2->y/w2, pt2->z/w2);

						fclose(fo);
					}
				}
				
				/* output */
				fprintf(stdout, "* Minimum distance from the input NURBS curve:\n\n");
				WriteParCurv(stdout, cv);
				fprintf(stdout, "\n  to the NURBS surface:\n\n");
				WriteParSurf(stdout, sf);
				fprintf(stdout, "\n");
				fprintf(stdout, "* (Signed) minimum distance = %lf is found\n", d0);
				fprintf(stdout, "  from the input curve at the parameter value u1 = %lf\n", u1);
				fprintf(stdout, "  to the input surface at the parameter values (u2,v2) = (%lf,%lf), and\n", u2,v2);
				pt1 = rbspeval(cv, u1, 0);		pt2 = revalderivsurf(sf, u2,v2, 0, 0);
				w1 = pt1->w;					w2 = pt2->w;
				fprintf(stdout, "  the corresponding points on the curve and surface are (%lf, %lf, %lf),\n", pt1->x/w1, pt1->y/w1, pt1->z/w1);
				fprintf(stdout, "  and (%lf, %lf, %lf), respectively.\n", pt2->x/w2, pt2->y/w2, pt2->z/w2);
				
				vectfree(pt1);		vectfree(pt2);
			
			} /* if(root_tol) { */
			else {
				fprintf(stderr, "Can't read a root tolerance for the IPP solver\n");
				print_usage();
				exit(-1);
			}

			fclose(cv_fi);		fclose(sf_fi);
			free_egeom(cv);	free_fgeom(sf);
		} /* if(sf_fi = fopen(sf_input, "r")) { */

		else {
			fprintf(stderr, "Can't open a file %s\n", sf_input);
			print_usage();
			exit(-1);
		}

	} /* if(cv_fi = fopen(cv_input, "r")) { */
	else {
		fprintf(stderr, "Can't open a file %s\n", cv_input);
		print_usage();
		exit(-1);
	}

	

	return 0;

}

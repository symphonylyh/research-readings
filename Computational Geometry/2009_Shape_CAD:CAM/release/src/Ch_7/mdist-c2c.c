/****************************************************************
 Copyright (C) 2008 Massachusetts Institute of Technology, Cambridge, MA
 All rights reserved

					     mdist-c2c.c

  Compute the minimum distance from a NURBS curve to a NURBS curve
  using the IPP solver

  How to make: 
  prompt> make mdist-c2c
  How to run:
  prompt> mdist-c2c -f curve1_file_name -t curve2_file_name -e root_tolerance_for_IPP_solver  [-o output_file_name]

  e.g.: mdist-c2c  -f ex.7.3.a.curv -t ex.7.3.b.curv  -e 1.e-4  -o ex.7.3.out

  ***************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include "gen.h"
#include "bspl.h"


void print_usage()
{
	fprintf(stderr, "usage: mdist-c2c -f curve1_file_name -t curve2_file_name -e root_tolerance_for_IPP_solver  [-o output_file_name]\n");

	return;
}


int main(int argc, char* argv[])
{
	FILE *cv1_fi, *cv2_fi, *fo;  
	ParCurv *cv1, *cv2;
	
	char c, *cv1_input, *cv2_input, *root_tol, *output;
	vector *pt1, *pt2;
	double w1, w2;
	int errflg;
	double tol, u1, u2, d0;

	cv1_input = cv2_input = root_tol = output = 0;
	errflg = 0;

	if (argc > 1) {
		while ((c = getopt(argc, argv, "f:t:e:o:")) != -1)
			switch (c) {
			case 'f':
				cv1_input = optarg;
				break;

			case 't':
				cv2_input = optarg;
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


	if(cv1_fi = fopen(cv1_input, "r")) { /* cv1_input is the file name storing the input curve #1 */
		cv1 = ReadParCurv(cv1_fi, NULL);

		if(cv2_fi = fopen(cv2_input, "r")) { /* cv2_input is the file name storing the input curve #2*/
			cv2 = ReadParCurv(cv2_fi, NULL);

			if(root_tol) { /* read a root tolerance for the IPP solver */
				tol = atof(root_tol); /* convert a string to double */
				if(tol <= 0) { /* check tol range */
					fprintf(stderr, "Root tolerance for the IPP solver, %lf should be a positive number\n", tol);
					print_usage();
					exit(-1);
				}

				/* Compute the minimum distance from the input curve cv1 to the input curve cv2 */
				/* The minimum distance is d0 and the corresponding parameter values for cv1, cv2 are u1, u2 */
				curvToCurv_n(cv1, cv2, tol, &u1, &u2, &d0);
		 	 
				if(output) { /* optional output */
					if(fo = fopen(output, "w")) {
						fprintf(fo, "* Minimum distance from the 1st input NURBS curve:\n\n");
						WriteParCurv(fo, cv1);
						fprintf(fo, "\n  to the 2nd input NURBS curve:\n\n");
						WriteParCurv(fo, cv2);
						fprintf(fo, "\n");
						fprintf(fo, "* (Signed) minimum distance = %+.16le is found\n", d0);
						fprintf(fo, "  from the 1st input curve at the parameter value u1 = %+.16le\n", u1);
						fprintf(fo, "  to the 2nd input curve at the parameter value u2 = %+.16le, and\n", u2);
						pt1 = rbspeval(cv1, u1, 0);		pt2 = rbspeval(cv2, u2, 0);
						w1 = pt1->w;					w2 = pt2->w;
						fprintf(fo, "  the corresponding points on the curves are (%+.16le, %+.16le, %+.16le),\n", pt1->x/w1, pt1->y/w1, pt1->z/w1);
						fprintf(fo, "  and (%+.16le, %+.16le, %+.16le), respectively.\n", pt2->x/w2, pt2->y/w2, pt2->z/w2);

						fclose(fo);
					}
				}
				
				/* output */
				fprintf(stdout, "* Minimum distance from the 1st input NURBS curve:\n\n");
				WriteParCurv(stdout, cv1);
				fprintf(stdout, "\n  to the 2nd input NURBS curve:\n\n");
				WriteParCurv(stdout, cv2);
				fprintf(stdout, "\n");
				fprintf(stdout, "* (Signed) minimum distance = %lf is found\n", d0);
				fprintf(stdout, "  from the 1st input curve at the parameter value u1 = %lf\n", u1);
				fprintf(stdout, "  to the 2nd input curve at the parameter value u2 = %lf, and\n", u2);
				pt1 = rbspeval(cv1, u1, 0);		pt2 = rbspeval(cv2, u2, 0);
				w1 = pt1->w;					w2 = pt2->w;
				fprintf(stdout, "  the corresponding points on the curves are (%lf, %lf, %lf),\n", pt1->x/w1, pt1->y/w1, pt1->z/w1);
				fprintf(stdout, "  and (%lf, %lf, %lf), respectively.\n", pt2->x/w2, pt2->y/w2, pt2->z/w2);
				
				vectfree(pt1);		vectfree(pt2);
			
			} /* if(root_tol) { */
			else {
				fprintf(stderr, "Can't read a root tolerance for the IPP solver\n");
				print_usage();
				exit(-1);
			}

			fclose(cv1_fi);		fclose(cv2_fi);
			free_egeom(cv1);	free_egeom(cv2);
		} /* if(cv2_fi = fopen(cv2_input, "r")) { */

		else {
			fprintf(stderr, "Can't open a file %s\n", cv2_input);
			print_usage();
			exit(-1);
		}

	} /* if(cv1_fi = fopen(cv1_input, "r")) { */
	else {
		fprintf(stderr, "Can't open a file %s\n", cv1_input);
		print_usage();
		exit(-1);
	}

	

	return 0;

}

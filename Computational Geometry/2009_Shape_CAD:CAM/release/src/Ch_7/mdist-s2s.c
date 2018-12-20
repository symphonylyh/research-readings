/****************************************************************
 Copyright (C) 2008 Massachusetts Institute of Technology, Cambridge, MA
 All rights reserved

					     mdist-s2s.c

  Compute the minimum distance from a NURBS surface to a NURBS surface
  using the IPP solver

  How to make: 
  prompt> make mdist-s2s
  How to run:
  prompt> mdist-s2s -f surface1_file_name -t surface2_file_name -e root_tolerance_for_IPP_solver  [-o output_file_name]

  e.g.: mdist-s2s  -f ex.7.5.a.surf -t ex.7.5.b.surf  -e 1.e-8  -o ex.7.5.out

  ***************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include "gen.h"
#include "bspl.h"


void print_usage()
{
	fprintf(stderr, "usage: mdist-s2s -f surface1_file_name -t surface2_file_name -e root_tolerance_for_IPP_solver  [-o output_file_name]\n");

	return;
}


int main(int argc, char* argv[])
{
	FILE *sf1_fi, *sf2_fi, *fo;  
	ParSurf *sf1, *sf2;
	char c, *sf1_input, *sf2_input, *root_tol, *output;
	vector *pt1, *pt2;
	double w1, w2;
	int errflg;
	double tol, u1,v1, u2,v2, d0;

	sf1_input = sf2_input = root_tol = output = 0;
	errflg = 0;

	if (argc > 1) {
		while ((c = getopt(argc, argv, "f:t:e:o:")) != -1)
			switch (c) {
			case 'f':
				sf1_input = optarg;
				break;

			case 't':
				sf2_input = optarg;
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


	if(sf1_fi = fopen(sf1_input, "r")) { /* sf1_input is the file name storing the input surface #1 */
		sf1 = ReadParSurf(sf1_fi, NULL);

		if(sf2_fi = fopen(sf2_input, "r")) { /* sf2_input is the file name storing the input surface #2*/
			sf2 = ReadParSurf(sf2_fi, NULL);

			if(root_tol) { /* read a root tolerance for the IPP solver */
				tol = atof(root_tol); /* convert a string to double */
				if(tol <= 0) { /* check tol range */
					fprintf(stderr, "Root tolerance for the IPP solver, %lf should be a positive number\n", tol);
					print_usage();
					exit(-1);
				}

				/* Compute the minimum distance from the input surface sf1 to the input surface sf2 */
				/* The minimum distance is d0 and the corresponding parameter values for sf1, sf2 are (u1,v1) and (u2,v2) */
				surfToSurf_n(sf1, sf2, tol, &u1,&v1, &u2,&v2, &d0);
	 	 
				if(output) { /* optional output */
					if(fo = fopen(output, "w")) {
						fprintf(fo, "* Minimum distance from the 1st input NURBS surface:\n\n");
						WriteParSurf(fo, sf1);
						fprintf(fo, "\n  to the 2nd input NURBS surface:\n\n");
						WriteParSurf(fo, sf2);
						fprintf(fo, "\n");
						fprintf(fo, "* (Signed) minimum distance = %lf is found\n", d0);
						fprintf(fo, "  from the 1st input surface at the parameter value (u1,v1) = (%lf,%lf)\n", u1,v1);
						fprintf(fo, "  to the 2nd input surface at the parameter value (u2,v2) = (%lf,%lf), and\n", u2,v2);
						pt1 = revalderivsurf(sf1, u1, v1, 0, 0);		pt2 = revalderivsurf(sf2, u2, v2, 0, 0);		
						w1 = pt1->w;					w2 = pt2->w;
						fprintf(fo, "  the corresponding points on the surfaces are (%lf, %lf, %lf),\n", pt1->x/w1, pt1->y/w1, pt1->z/w1);
						fprintf(fo, "  and (%lf, %lf, %lf), respectively.\n", pt2->x/w2, pt2->y/w2, pt2->z/w2);

						fclose(fo);
					}
				}
				
				/* output */
				fprintf(stdout, "* Minimum distance from the 1st input NURBS surface:\n\n");
				WriteParSurf(stdout, sf1);
				fprintf(stdout, "\n  to the 2nd input NURBS surface:\n\n");
				WriteParSurf(stdout, sf2);
				fprintf(stdout, "\n");
				fprintf(stdout, "* (Signed) minimum distance = %lf is found\n", d0);
				fprintf(stdout, "  from the 1st input surface at the parameter value (u1,v1) = (%lf,%lf)\n", u1,v1);
				fprintf(stdout, "  to the 2nd input surface at the parameter value (u2,v2) = (%lf,%lf), and\n", u2,v2);
				pt1 = revalderivsurf(sf1, u1, v1, 0, 0);		pt2 = revalderivsurf(sf2, u2, v2, 0, 0);		
				w1 = pt1->w;					w2 = pt2->w;
				fprintf(stdout, "  the corresponding points on the surfaces are (%lf, %lf, %lf),\n", pt1->x/w1, pt1->y/w1, pt1->z/w1);
				fprintf(stdout, "  and (%lf, %lf, %lf), respectively.\n", pt2->x/w2, pt2->y/w2, pt2->z/w2);
				
				vectfree(pt1);		vectfree(pt2);
			
			} /* if(root_tol) { */
			else {
				fprintf(stderr, "Can't read a root tolerance for the IPP solver\n");
				print_usage();
				exit(-1);
			}

			fclose(sf1_fi);		fclose(sf2_fi);
			free_fgeom(sf1);	free_fgeom(sf2);
		} /* if(sf2_fi = fopen(sf2_input, "r")) { */

		else {
			fprintf(stderr, "Can't open a file %s\n", sf2_input);
			print_usage();
			exit(-1);
		}

	} /* if(sf1_fi = fopen(sf1_input, "r")) { */
	else {
		fprintf(stderr, "Can't open a file %s\n", sf1_input);
		print_usage();
		exit(-1);
	}

	

	return 0;

}

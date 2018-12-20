/****************************************************************
 Copyright (C) 2008 Massachusetts Institute of Technology, Cambridge, MA
 All rights reserved

					     mdist-p2s.c

  Compute the minimum distance from a point to a NURBS surface 
  using the IPP solver

  How to make: 
  prompt> make mdist-p2s
  How to run:
  prompt> mdist-p2s -f point_file_name -t surface_file_name -e root_tolerance_for_IPP_solver  [-o output_file_name]

  e.g.: mdist-p2s  -f ex.7.2.pt -t ex.7.2.surf  -e 1.e-8  -o ex.7.2.out

  ***************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include "gen.h"
#include "bspl.h"


void print_usage()
{
	fprintf(stderr, "usage: mdist-p2s -f point_file_name -t surface_file_name -e root_tolerance_for_IPP_solver  [-o output_file_name]\n");

	return;
}


int main(int argc, char* argv[])
{
	FILE *pt_fi, *sf_fi, *fo;  
	ParSurf* sf;
	
	char c, *pt_input, *sf_input, *root_tol, *output;
	double x, y, z;
	vector* pt;
	double w;
	int errflg;
	double tol, u0, v0, d0;

	pt_input = sf_input = root_tol = output = 0;
	errflg = 0;

	if (argc > 1) {
		while ((c = getopt(argc, argv, "f:t:e:o:")) != -1)
			switch (c) {
			case 'f':
				pt_input = optarg;
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


	if(pt_fi = fopen(pt_input, "r")) { /* pt_input is the file name storing the input point */
		fscanf(pt_fi, "%lf %lf %lf\n", &x, &y, &z); /* read x, y, z */

		if(sf_fi = fopen(sf_input, "r")) { /* sf_input is the file name storing the input surface */
			sf = ReadParSurf(sf_fi, NULL);
 
			if(root_tol) { /* read a root tolerance for the IPP solver */
				tol = atof(root_tol); /* convert a string to double */
				if(tol <= 0) { /* check tol range */
					fprintf(stderr, "Root tolerance for the IPP solver, %lf should be a positive number\n", tol);
					print_usage();
					exit(-1);
				}

				/* Compute the minimum distance from the input point (x,y,z) to the input surface sf */
				/* The minimum distance is d0 and the corresponding parameter value for sf is (u0,v0) */
				pointToSurf_n(x,y,z, sf, tol, &u0, &v0, &d0);
			 
				if(output) { /* optional output */
					if(fo = fopen(output, "w")) {
						fprintf(fo, "* Minimum distance from the input point (%lf,%lf,%lf)\n", x,y,z);
						fprintf(fo, "  to the input NURBS surface:\n\n");
						WriteParSurf(fo, sf);
						fprintf(fo, "\n");
						fprintf(fo, "* (Signed) minimum distance = %lf is found\n", d0);
						fprintf(fo, "  to the input surface at the parameter value (u,v) = (%lf,%lf)\n", u0,v0);
						pt = revalderivsurf(sf, u0, v0, 0, 0);
						w = pt->w;
						fprintf(fo, "  i.e. (%lf, %lf, %lf) in 3D coordinates ",pt->x/w, pt->y/w, pt->z/w);
						fprintf(fo, "on the surface\n"); 	


						fclose(fo);
					}
				}
				
				/* output */
				fprintf(stdout, "* Minimum distance from the input point (%lf,%lf,%lf)\n", x,y,z);
				fprintf(stdout, "  to the input NURBS surface:\n\n");
				WriteParSurf(stdout, sf);
				fprintf(stdout, "\n");
				fprintf(stdout, "* (Signed) minimum distance = %lf is found\n", d0);
				fprintf(stdout, "  to the input surface at the parameter value (u,v) = (%lf,%lf)\n", u0,v0);
				pt = revalderivsurf(sf, u0, v0, 0, 0);
				w = pt->w;
				fprintf(stdout, "  i.e. (%lf, %lf, %lf) in 3D coordinates ",pt->x/w, pt->y/w, pt->z/w);
				fprintf(stdout, "on the surface\n"); 	

				
				vectfree(pt);
			
			} /* if(root_tol) { */
			else {
				fprintf(stderr, "Can't read a root tolerance for the IPP solver\n");
				print_usage();
				exit(-1);
			}

		fclose(sf_fi);
		free_fgeom(sf);
	}

	else {
		fprintf(stderr, "Can't open a file %s\n", sf_input);
		print_usage();
		exit(-1);
	}

	}
	else {
		fprintf(stderr, "Can't open a file %s\n", pt_input);
		print_usage();
		exit(-1);
	}

	

	return 0;

}

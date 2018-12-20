/****************************************************************
 Copyright (C) 2008 Massachusetts Institute of Technology, Cambridge, MA
 All rights reserved

					     evals.c

  Evaluate i,j-th derivatives of an open NURBS surface at (u,v)
  where i,j = 0,1,2...
  So, if i = j = 0, then it will evaluate the 0-th derivative in u,v i.e. a position at (u,v).

  make evals
  evals -i input_surface_file_name -p u_value,v_value -d #_of_u-deriv,#_of_v-deriv [-o output_file_name]
  
  Note: No spaces before or after the "comma" in -p and -d  command line options.

  e.g.: Example 1.5.2 in pages 32--33 (for theta = phi = 45 deg)
	`	evals -i s.surf  -p 0.414213562,0.414213562  -d 0,0  -o s.out
  ***************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include "gen.h"
#include "bspl.h"

void print_usage()
{
	fprintf(stderr, "usage: evals -i input_surface_file_name -p u_value,v_value -d #_of_u-deriv,#_of_v-deriv [-o output_file_name]\n");

	return;
}


int main(int argc, char* argv[])
{
	FILE *fi, *fo;  
	ParSurf* sf;
	double u,v;
	char c, *input, *output;
	vector* pt;
	double w;
	int i, iu_deriv, iv_deriv, errflg;
	double u_min, u_max, v_min, v_max;

	input = output = 0;
	errflg = 0;
	u = v = -1.0;	
	iu_deriv = iv_deriv = -1;

	if (argc > 1) {
		while ((c = getopt(argc, argv, "i:p:d:o:")) != -1)
			switch (c) {
			case 'i':
				input = optarg;
				break;

			case 'p':
				for (i=0; i<strlen(optarg); i++)
					if (optarg[i] == ',')
						optarg[i] = ' ';
				if(sscanf(optarg, "%lf %lf", &u, &v) == 0) {
					fprintf(stderr, "Can't read parameter values (u,v)\n");
					print_usage();
					exit(-1);
				}
				break;

			case 'd':
				for (i=0; i<strlen(optarg); i++)
					if (optarg[i] == ',')
						optarg[i] = ' ';
				if(sscanf(optarg, "%d %d", &iu_deriv, &iv_deriv) == 0) {
					fprintf(stderr, "Can't read # of derivatives in (u,v)\n");
					print_usage();
					exit(-1);
				}
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


	if(fi = fopen(input, "r")) { /* input is the file name storing the input surface */
		sf = ReadParSurf(fi, NULL);

		u_min = sf->uknots[0];
		u_max = sf->uknots[sf->ucontpts + sf->uorder -1];
		if(u < u_min || u > u_max) { /* check param range */
			fprintf(stderr, "Paramer u_value should be within [%lf,%lf]\n", u_min, u_max);
			print_usage();
			exit(-1);
		}
			
		v_min = sf->vknots[0];
		v_max = sf->vknots[sf->vcontpts + sf->vorder -1];
		if(v < v_min || v > v_max) { /* check param range */
			fprintf(stderr, "Paramer v_value should be within [%lf,%lf]\n", v_min, v_max);
			print_usage();
			exit(-1);
		}

			 
		if(iu_deriv >= 0) { /* iu_deriv = 0,1,2,... */  		 
			if(iv_deriv >= 0) { /* iv_deriv = 0,1,2,... */ 
				
				if(pt = revalderivsurf(sf, u, v, iu_deriv, iv_deriv)) { /* evaluate a surface */
					w = pt->w;
					if(output) {
						if(fo = fopen(output, "w")) {
							fprintf(fo, "At (u,v) = (%+.16le,%+.16le) w/ # of derivatives in u,v = (%d,%d):\n", u,v, iu_deriv,iv_deriv);
							fprintf(fo, "In homogeneous coordinates:\n");
							Writevector(fo, pt); /* write a result in homogeneous coords */
							fprintf(fo, "So, in 3D coordinates:\n"); /* and in 3D coords */
							fprintf(fo,"%+.16le %+.16le %+.16le \n",pt->x/w, pt->y/w, pt->z/w);

							fclose(fo);
						}
					}
				
					fprintf(stdout, "At (u,v) = (%+.16le,%+.16le) w/ # of derivatives in u,v = (%d,%d):\n", u,v, iu_deriv,iv_deriv);
					fprintf(stdout, "In homogeneous coordinates:\n");
					Writevector(stdout, pt);
                    fprintf(stdout, "So, in 3D coordinates:\n");
					fprintf(stdout,"%+.16le %+.16le %+.16le \n",pt->x/w, pt->y/w, pt->z/w);
				
					vectfree(pt);
				}/* if(pt = revalderivsurf(sf, u, v, 0, 0)) */

			}/* if(iv_deriv >= 0) { */
			else {
					fprintf(stderr, "You should enter a non-negative integer for #_of_v-deriv\n");
					print_usage();
					exit(-1);
			}
		}/* if(iu_deriv >= 0) */
		else {
			fprintf(stderr, "You should enter a non-negative integer for #_of_u-deriv\n");
			print_usage();
			exit(-1);
		}
		
		fclose(fi);
		free_fgeom(sf);

	}/* if(fi = fopen(input, "r")) */

	else {
		fprintf(stderr, "Can't open a file %s\n", input);
		print_usage();
		exit(-1);
	}

	return 0;
}

/****************************************************************
 Copyright (C) 2008 Massachusetts Institute of Technology, Cambridge, MA
 All rights reserved

					     subs.c

  Subdivide an open NURBS surface along u (or v)

  make subs
  subs -i input_surface_file  -u u_value (or -v v_value)  [-o output_surface_file_1,output_surface_file_2]

  Note: No spaces before or after the "comma" in -o command line options.

  e.g.: Example 1.5.2 in pages 32--33 

        (subdivide the surface along theta = 45 deg, i.e. along u = 0.414213562)
        subs  -i s.surf  -u 0.414213562  -o s1.surf,s2.surf

		(subdivide the surface along phi = 45 deg, i.e. along v = 0.414213562)
        subs  -i s.surf  -v 0.414213562  -o s1.surf,s2.surf
  ***************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include "gen.h"
#include "bspl.h"
#include "editor.h"

void print_usage()
{
  fprintf(stderr, "usage: subs -i input_surface_file  -u u_value (or -v v_value)  [-o output_surface_file_1,output_surface_file_2]\n");

  return;
}


int main(int argc, char* argv[])
{
  double uParam, vParam; /* subdivide the input surface at uParam or vParam */
  short i, indx = 1, j, k, ku, kv, n, nSurfs, uKnots = 0, vKnots = 0;

  ParSurf *fgeom, *fgeom1, *fgeom2; /* input surface and its subdivided surfaces */

  FILE *fi, *fo1, *fo2;  
  char c, *input, *para, output1[80], output2[80];
  double u_min, u_max, v_min, v_max;

  short u_sub = 0, v_sub = 0;
  input = para = 0;
  short errflg = 0;
  output1[0] = output2[0] = '\0';

  if (argc > 1) {
		while ((c = getopt(argc, argv, "i:u:v:o:")) != -1)
			switch (c) {
			case 'i':
				input = optarg;
				break;

			case 'u':
				para = optarg;
				u_sub = 1;
				break;

			case 'v':
				para = optarg;
				v_sub = 1;
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
	  fgeom = ReadParSurf(fi, NULL); /* read an input surface */
	  fclose(fi);

  if (u_sub == 0 && v_sub == 1) {
	  uKnots = 0;
	  vKnots = fgeom->vorder;
  }
  else if (u_sub == 1 && v_sub == 0) {
	  uKnots = fgeom->uorder;
	  vKnots = 0;
  }
  else {
	  fprintf(stderr, "Enter a subdivision direction (i.e. along u or v).\n");
	  print_usage();
	  exit(-1);
  }

  ku = fgeom->ucontpts + fgeom->uorder;

  if(para) { /* read a parameter value for subdivision */		

	  if (uKnots) {
		  uParam = atof(para); /* string to double */
		  u_min = fgeom->uknots[0];
		  u_max = fgeom->uknots[fgeom->ucontpts + fgeom->uorder -1];
		  if(uParam < u_min || uParam > u_max) { /* check param range */
			  fprintf(stderr, "Paramer u-value should be within [%lf,%lf]\n", u_min, u_max);
			  print_usage();
			  exit(-1);
		  }

		  for (i=j=0; i<ku; i++)
			  if(fgeom->uknots[i] == uParam)
				  j++;
	
		  if (j == fgeom->uorder) {
			  fprintf(stderr, "Already %d knots at u = %g", j, uParam);
			  exit(-1);
		  }

	  }/* if (uKnots) { */


	  kv = fgeom->vcontpts + fgeom->vorder;
	  if (vKnots) {
		  vParam = atof(para); /* string to double */
		  v_min = fgeom->vknots[0];
		  v_max = fgeom->vknots[fgeom->vcontpts + fgeom->vorder -1];
		  if(vParam < v_min || vParam > v_max) { /* check param range */
			  fprintf(stderr, "Paramer v-value should be within [%lf,%lf]\n", v_min, v_max);
			  print_usage();
			  exit(-1);
		  }

		  for (i=j=0; i<kv; i++)
			  if(fgeom->vknots[i] == vParam)
				  j++;
    
		  if (j == fgeom->vorder) {
			  fprintf(stderr, "Already %d knots at u = %g", j, vParam);
			  exit(-1);
		  }
	  }/* if (vKnots) { */

	  SplitSurf(fgeom, uKnots, uKnots ? uParam : vParam, &fgeom1, &fgeom2); 

	  if(fo1 = fopen(output1, "w")) {
		  WriteParSurf(fo1, fgeom1);
		  fclose(fo1);
	  }

	  if(fo2 = fopen(output2, "w")) {
		  WriteParSurf(fo2, fgeom2);
		  fclose(fo2);
	  }

	  WriteParSurf(stdout, fgeom1);
	  WriteParSurf(stdout, fgeom2);

	  free_fgeom(fgeom1);
	  free_fgeom(fgeom2);

  }/* if(para) { */
  else {
	  fprintf(stderr, "Can't read a parametric point alonh which the surface is divided\n");
	  print_usage();
	  exit(-1);
  }
  
 }/* if(fi = fopen(input, "r")) { */
 else {
	 fprintf(stderr, "Can't open a file %s\n", input);
	 print_usage();
	 exit(-1);
 }

 return 0;
}

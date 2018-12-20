/****************************************************************
 Copyright (C) 2008 Massachusetts Institute of Technology, Cambridge, MA
 All rights reserved

					     knots.c

  Insert new knots into the knot vectors of an open NURBS surface at (u,v) 

  make knots
  knots -i input_surface_file_name  -p u_value,v_value  -k #_of_u-knots_inserted,#_of_v-knots_inserted  [-o output_surface_file_name]

  Note: No spaces before or after the "comma" in -p and -k command line options.

  e.g.: Example 1.5.2 in pages 32--33 
        (insert 3 knots at u = 0.414213562 (i.e. corresponding to theta = 45 deg) and
		 also   3 knots at v = 0.414213562 (i.e. corresponding to phi = 45 deg))
        knots -i s.surf  -p 0.414213562,0.414213562  -k 3,3  -o s1.surf
  ***************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include "gen.h"
#include "bspl.h"

void print_usage()
{
	fprintf(stderr, "usage: knots -i input_surface_file_name  -p u_value,v_value  -k #_of_u-knots_inserted,#_of_v-knots_inserted  [-o output_surface_file_name]\n");

	return;
}


int main(int argc, char* argv[])
{
  double uParam, vParam; /* insert new knots at u=uParam, v=vParam */
  short i, j, ku, kv, uKnots = -1, vKnots = -1, errflg = 0;

  ParSurf *fgeom; /* input and output surface */
  ParSurf *temp;

  FILE *fi, *fo;  
  char c, *input, *output;
  double u_min, u_max, v_min, v_max;

  input = output = 0;
  uParam = vParam = -1.0;	

  if (argc > 1) {
	  while ((c = getopt(argc, argv, "i:p:k:o:")) != -1)
		  switch (c) {
			case 'i':
				input = optarg;
				break;

			case 'p':
				for (i=0; i<strlen(optarg); i++)
					if (optarg[i] == ',')
						optarg[i] = ' ';
				if(sscanf(optarg, "%lf %lf", &uParam, &vParam) == 0) {
					fprintf(stderr, "Can't read parameter values (u,v)\n");
					print_usage();
					exit(-1);
				}
				break;

			case 'k':
				for (i=0; i<strlen(optarg); i++)
					if (optarg[i] == ',')
						optarg[i] = ' ';
				if(sscanf(optarg, "%hd %hd", &uKnots, &vKnots) == 0) {
					fprintf(stderr, "Can't read # of knots to be inserted in (u,v)\n");
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


  if(fi = fopen(input, "r")) { 
	  fgeom = ReadParSurf(fi, NULL); /* read an input surface */
	  fclose(fi);

	  ku = fgeom->ucontpts + fgeom->uorder;

	  u_min = fgeom->uknots[0];
	  u_max = fgeom->uknots[fgeom->ucontpts + fgeom->uorder -1];
	  if(uParam < u_min || uParam > u_max) { /* check param range */
		  fprintf(stderr, "You must enter a paramer u-value within [%lf,%lf]\n", u_min, u_max);
		  print_usage();	  
		  exit(-1);
	  }

	  v_min = fgeom->vknots[0];
	  v_max = fgeom->vknots[fgeom->vcontpts + fgeom->vorder -1];
	  if(vParam < v_min || vParam > v_max) { /* check param range */
		  fprintf(stderr, "You must enter a paramer v-value within [%lf,%lf]\n", v_min, v_max);
		  print_usage();
		  exit(-1);
	  }

	  for (i=j=0; i<ku; i++)
		  if(fgeom->uknots[i] == uParam)
			  j++;
		  

	  if (j == fgeom->uorder) {
		  fprintf(stderr, "Already %d knots at u = %g\n", j, uParam);
		  exit(-1);
	  }

	  if(uKnots < 0 || uKnots > fgeom->uorder -j) { /* check uKnots range */
		  fprintf(stderr, "You must enter # of u-knots to be inserted within [0,%d]\n", fgeom->uorder -j);
	      print_usage();	  
		  exit(-1);
	  }

	  kv = fgeom->vcontpts + fgeom->vorder;

	  for (i=j=0; i<kv; i++)
		  if(fgeom->vknots[i] == vParam)
			  j++;
    
	  if (j == fgeom->vorder) {
		  fprintf(stderr, "Already %d knots at v = %g\n", j, vParam);
		  exit(-1);
	  }

      if(vKnots < 0 || vKnots > fgeom->vorder -j) { /* check vKnots range */
		  fprintf(stderr, "You must enter # of v-knots to be inserted within [0,%d]\n", fgeom->vorder -j);
		  print_usage();
		  exit(-1);
	  }

	  if(uKnots == 0 && vKnots == 0) /* return the original surface unchanged */
		  goto RETURN;

	  temp = copyfgeom(fgeom, 0);
	  free_fgeom(fgeom);
	  fgeom = fgeomalloc1(temp->uorder, temp->vorder, temp->ucontpts + uKnots, temp->vcontpts + vKnots);
	  copyfgeom(temp, fgeom);
	  fgeom->ucontpts = temp->ucontpts + uKnots;
	  fgeom->vcontpts = temp->vcontpts + vKnots;

	  if (uKnots) {
		  fgeom->uknots[0] = temp->uknots[0];
		  for(i = 1; i < temp->ucontpts + fgeom->uorder; i++) {
			  if(temp->uknots[i] >= uParam && temp->uknots[i-1] < uParam)
				  for(j=0; j<uKnots; j++)
					  fgeom->uknots[i+j] = uParam;
				  fgeom->uknots[i+j] = temp->uknots[i];
		  }
	  }

	  if (vKnots) {
		  fgeom->vknots[0] = temp->vknots[0];
		  for(i = 1; i < temp->vcontpts + fgeom->vorder; i++) {
			  if(temp->vknots[i] >= vParam && temp->vknots[i-1] < vParam)
				  for(j=0; j<vKnots; j++)
					  fgeom->vknots[i+j] = vParam;
				  fgeom->vknots[i+j] = temp->vknots[i];
		  }
	  }

	  surfoslo3(temp, fgeom, 4);

	  free_fgeom(temp);

RETURN:
	  if(fo = fopen(output, "w")) {
		  WriteParSurf(fo, fgeom);
		  fclose(fo);
	  }

	  WriteParSurf(stdout, fgeom);

	  free_fgeom(fgeom);

	}/* if(fi = fopen(input, "r")) { */
	else {
		fprintf(stderr, "Can't open a file %s\n", input);
		print_usage();
		exit(-1);
	}

  return 0;
}

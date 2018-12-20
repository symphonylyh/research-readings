/***********************************************************************************
 Copyright (C) 2008 Massachusetts Institute of Technology, Cambridge, MA
 All rights reserved

						geoms.c

  Evaluate a position, normalized tangent directions in u,v, normal direction, 
  max/min principal curvatures, and normal u,v curvatures
  of an open NURBS surface at (u,v)
  
	Make:
	make geoms
	Run:
	geoms  -i input_surface_file  -p u_value,v_value  [-o output_file_name]
    Output: 
    Evaluation results

	Note: No spaces before or after the "comma" in -p command line option.
 
   e.g.: Example 1.5.2 in pages 32--33 (for theta = phi = 45 deg)
         geoms  -i s.surf  -p 0.414213562,0.414213562  -o s.out 
  *********************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include "gen.h"
#include "bspl.h"
#include "praxiteles.h" 
#include "iges.h" 


static double curvatureZero = 1.0e-10;
static double distanceZero  = 1.0e-10;
static double knotZero      = 1.0e-10;

void print_usage()
{
	fprintf(stderr, "usage: geoms  -i input_surface_file  -p u_value,v_value  [-o output_file_name]\n");

	return;
}

double GetZero(int type)
{
  double zero = 0.0;

  switch (type) {
  case PRAX_ZERO_CURVATURE:
    zero = curvatureZero;
    break;
  case PRAX_ZERO_DISTANCE:
    zero = distanceZero;
    break;
  case PRAX_ZERO_KNOT:
    zero = knotZero;
    break;
  }

  return zero;
}



/* Evaluate a position, normalized tangent directions in u,v, normal direction, 
   max/min principal curvatures, and normal u,v curvatures
   of an open NURBS surface at (u,v)     */
vector *GetSurfValues(ParSurf *fgeom, float u, float v, vector **tangu, vector **tangv, 
					  vector **norm, double *kmax, double *kmin, double *normu, double *normv)
{
  struct vector *r[6], *ruv, *ruu, *rvv;
  double curvatureZero, e, f, g, l, m, n, H, K;
  double2 tempmod, Gdet, Ddet;

  evalrsurf(fgeom, u, v, 5, r);  /* pts, 1st, and 2nd derivatives */

  if (mag(r[1]) < PRAX_ZERO10)   /* normal vector */
    (*norm) = cross(r[3], r[2]);	
  else if(mag(r[2]) < PRAX_ZERO10)
    (*norm) = cross(r[1], r[3]);
  else
    (*norm) = cross(r[1], r[2]);
  unitvector1((*norm), (*norm));

  curvatureZero = GetZero(PRAX_ZERO_CURVATURE);


    e = dot(r[1], r[1]);   /* fundamental coefficients */
    f = dot(r[1], r[2]);
    g = dot(r[2], r[2]);

    l = dot((*norm), r[4]);
    m = dot((*norm), r[3]);
    n = dot((*norm), r[5]);

    *normu = -l/e;   /* normal u curvature */
    *normv = -n/g;   /* normal v curvature */

    if (fabs(*normu) < curvatureZero) *normu = 0.0; /* if curvature is small */
    if (fabs(*normv) < curvatureZero) *normv = 0.0; /* enough, consider it 0 */

    Gdet = e*g - f*f;
    Ddet = l*n - m*m;
  
    H = ((f*m + f*m) - (e*n + g*l))/(2.0*Gdet);   /* mean */
    K = Ddet/Gdet;                                /* gauss */

    tempmod = fabs(H*H - K);

      *kmax = H + sqrt(tempmod);   /* max principal */
      *kmin = H - sqrt(tempmod);   /* min principal */


  unitvector1(r[1], r[1]);
  unitvector1(r[2], r[2]);
  (*tangu) = r[1];
  (*tangv) = r[2];
  vectfree(r[3]);
  vectfree(r[4]);
  vectfree(r[5]);

  return (r[0]);
}


int main(int argc, char* argv[])
{
	vector *norm, *r, *tangu, *tangv;
    	double kn, kx, nu, nv;

	FILE *fi, *fo;  
	ParSurf* fgeom;
	double u,v;
	char c, *input, *output;
	double u_min, u_max, v_min, v_max;

	short i, errflg = 0;
	input = output = 0;
	u = v = -1.0;

	if (argc > 1) {
		while ((c = getopt(argc, argv, "i:p:o:")) != -1)
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
		fgeom = ReadParSurf(fi, NULL);
		
		u_min = fgeom->uknots[0];
		u_max = fgeom->uknots[fgeom->ucontpts + fgeom->uorder -1];
		if(u < u_min || u > u_max) { /* check param range */
			fprintf(stderr, "Paramer u-value should be within [%lf,%lf]\n", u_min, u_max);
			print_usage();
			exit(-1);
		}

		v_min = fgeom->vknots[0];
		v_max = fgeom->vknots[fgeom->vcontpts + fgeom->vorder -1];
		if(v < v_min || v > v_max) { /* check param range */
			fprintf(stderr, "Paramer v-value should be within [%lf,%lf]\n", v_min, v_max);
			print_usage();
			exit(-1);
		}

		if(r = GetSurfValues(fgeom, u, v, &tangu, &tangv, &norm, &kx, &kn, &nu, &nv)) {

			if(output) {
				if(fo = fopen(output, "w")) {
					fprintf(fo, "# X,Y,Z = %+.15e %+.15e %+.15e\n", r->x/r->w, r->y/r->w, r->z/r->w);
					fprintf(fo, "# Tangent U = %+.15e %+.15e %+.15e\n", tangu->x/tangu->w, tangu->y/tangu->w, tangu->z/tangu->w);
					fprintf(fo, "# Tangent V = %+.15e %+.15e %+.15e\n", tangv->x/tangv->w, tangv->y/tangv->w, tangv->z/tangv->w);
					fprintf(fo, "# Normal    = %+.15e %+.15e %+.15e\n", norm->x/norm->w, norm->y/norm->w, norm->z/norm->w);
					fprintf(fo, "# Kmin, Kmax = %+.15e %+.15e\n", kn, kx);
					fprintf(fo, "# Normal U,V = %+.15e %+.15e\n", nu, nv);

					fclose(fo);
				}/* if(fo = fopen(output, "w")) { */
			}/* if(output) */

			fprintf(stdout, "# X,Y,Z = %9.6g %9.6g %9.6g\n", r->x/r->w, r->y/r->w, r->z/r->w);
			fprintf(stdout, "# Tangent U = %9.6g %9.6g %9.6g\n", tangu->x/tangu->w, tangu->y/tangu->w, tangu->z/tangu->w);
			fprintf(stdout, "# Tangent V = %9.6g %9.6g %9.6g\n", tangv->x/tangv->w, tangv->y/tangv->w, tangv->z/tangv->w);
			fprintf(stdout, "# Normal    = %9.6g %9.6g %9.6g\n", norm->x/norm->w, norm->y/norm->w, norm->z/norm->w);
			fprintf(stdout, "# Kmin, Kmax = %9.6g %9.6g\n", kn, kx);
			fprintf(stdout, "# Normal U,V = %9.6g %9.6g\n", nu, nv);

			free(r);
			free(tangu);
			free(tangv);
			free(norm);							

		} /* if(pt = GetSurfValues(fgeom, u, v, &tangu, &tangv, &norm, &kx, &kn, &nu, &nv)) */	

	fclose(fi);
	free_fgeom(fgeom);

	} /* if(fi = fopen(input, "r")) { */
	else {
		fprintf(stderr, "Can't open a file %s\n", input);
		print_usage();
		exit(-1);
	}

	return 0;
}


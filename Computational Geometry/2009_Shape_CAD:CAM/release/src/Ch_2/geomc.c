/***********************************************************************************
 Copyright (C) 2008 Massachusetts Institute of Technology, Cambridge, MA
 All rights reserved

						geomc.c

  Evaluate a position, tangent, normal, curvature, and torsion of an open NURBS curve at u
  
	Make:
	make geomc
	Run:
	geomc  -i input_curve_file  -u u_value  [-o output_file]
    Output: 
    Evaluation results
 
    example: quadrant of an ellipse in Example 1.5.1 in pages 30--31
             for theta = 45 deg
             prompt> geomc -i c.curv  -u 0.414213562  -o c.out
  *********************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include "gen.h"
#include "bspl.h"
#include "praxiteles.h" 
#include "iges.h" 

#define IS_OPEN 1

static double curvatureZero = 1.0e-10;
static double distanceZero  = 1.0e-10;
static double knotZero      = 1.0e-10;

void print_usage()
{
	fprintf(stderr, "usage: geomc  -i input_curve_file  -u u_value  [-o output_file]\n");

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


/* Evaluate a position, tangent, normal, curvature, and torsion of an open NURBS curve */
vector *GetCurvValues(ParCurv *egeom, double t, int is_open, int is_planar,
		      vector **tang, vector **norm, double *kcurv, double *s)
{
  vector *c, *r, *r2, *r3, *u;
  double a, b, k;

  r = evalbsp(egeom, t);
  if (egeom->order > 2) {
	  (*tang) = rbspeval(egeom, t, 1);
      r2 = rbspeval(egeom, t, 2);
      r3 = rbspeval(egeom, t, 3);
  }
  
  if (egeom->order > 2) {

    if (is_planar) {
      a = (*tang)->x*(*tang)->x + (*tang)->y*(*tang)->y;
      b = sqrt(a);

      (*kcurv) = ((*tang)->x*r2->y - (*tang)->y*r2->x)/sqrt(a*a*a);

      *s = 0.0;                /* planar curve has no torsion */

      (*norm) = vectalloc();
      (*norm)->x = -(*tang)->y/b;
      (*norm)->y =  (*tang)->x/b;
      (*norm)->z = 0.0;
      (*norm)->w = 1.0;
    }
    else {
      a = mag((*tang));
      u = cross((*tang), r2);
      b = mag(u);
      (*kcurv) = b/(a*a*a);

      if (fabs(*kcurv) < GetZero(PRAX_ZERO_CURVATURE)) /* if curvature is */
	     (*kcurv) = 0.0;                 /* small enough, consider it zero */

      *s = -dot(u, r3)/(b*b);

      a = dot(r2, (*tang));
      b = (mag((*tang))*mag((*tang)));
      scale_vect1(a/b, (*tang), u);
      (*norm) = sub_vect(r2, u);
      vectfree(u);
    }
    vectfree(r2);
    vectfree(r3);
    unitvector1((*tang), (*tang));
    if (*kcurv != 0.0)
      unitvector1((*norm), (*norm));
  }
  else {
    (*kcurv) = 0.0;
    *s = 0.0;
    *tang = vectalloc();
    (*tang)->x = (*tang)->y = (*tang)->z = 0.0;
    (*tang)->w = 1.0;
    *norm = vectalloc();
    (*norm)->x = (*norm)->y = (*norm)->z = 0.0;
    (*norm)->w = 1.0;
  }

  return (r);
}


int main(int argc, char* argv[])
{
	int is_planar;
	double kc, s;
	vector *nvect, *pt, *tang;
	double x, y, z;

	FILE *fi, *fo;  
	ParCurv* egeom;
	double u;
	char c, *input, *param, *output;
	double u_min, u_max;

	short errflg = 0;
	input = param = output = 0;

	if (argc > 1) {
		while ((c = getopt(argc, argv, "i:u:o:")) != -1)
			switch (c) {
			case 'i':
				input = optarg;
				break;

			case 'u':
				param = optarg;
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
		egeom = ReadParCurv(fi, NULL);

		if(param) { /* read a parameter value */
			u = atof(param); /* convert a string to double */
			u_min = egeom->knots[0];
			u_max = egeom->knots[egeom->ncontpts + egeom->order -1];
			if(u < u_min || u > u_max) { /* check param range */
				fprintf(stderr, "Paramer value should be within [%lf,%lf]\n", u_min, u_max);
				print_usage();
			    exit(-1);
			}

			is_planar = CheckCurvPlanar(egeom, &x,&y,&z);

			if(pt = GetCurvValues(egeom, u, IS_OPEN, is_planar, &tang, &nvect, &kc, &s)) {

				if(output) {
					if(fo = fopen(output, "w")) {
  
						fprintf(fo, "# X,Y,Z   = %+.16le %+.16le %+.16le\n", pt->x/pt->w, pt->y/pt->w,
									pt->z/pt->w);
						fprintf(fo, "# Tangent = %+.16le %+.16le %+.16le\n", tang->x/tang->w,
									tang->y/tang->w, tang->z/tang->w);
						fprintf(fo, "# Normal  = %+.16le %+.16le %+.16le\n", nvect->x/nvect->w,
									nvect->y/nvect->w, nvect->z/nvect->w);
						fprintf(fo, "# Curvature = %+.16le\n", kc);
						fprintf(fo, "# Torsion   = %+.16le\n", s);

						fclose(fo);
					}/* if(fo = fopen(output, "w")) { */
				}/* if(output) */

						fprintf(stdout, "# X,Y,Z   = %9.6g %9.6g %9.6g\n", pt->x/pt->w, pt->y/pt->w,
									pt->z/pt->w);
						fprintf(stdout, "# Tangent = %9.6g %9.6g %9.6g\n", tang->x/tang->w,
									tang->y/tang->w, tang->z/tang->w);
						fprintf(stdout, "# Normal  = %9.6g %9.6g %9.6g\n", nvect->x/nvect->w,
									nvect->y/nvect->w, nvect->z/nvect->w);
						fprintf(stdout, "# Curvature = %9.6g\n", kc);
						fprintf(stdout, "# Torsion   = %9.6g\n", s);

 
						free(pt);
						free(tang);
						free(nvect);

			} /* if(pt = GetCurvValues(egeom, u, 1, is_planar, &tang, &nvect, &kc, &s)) */

		} /* if(param) */
		else {
			fprintf(stderr, "Can't read a parametric point u\n");
			print_usage();
			exit(-1);
		}

		fclose(fi);
		free_egeom(egeom);

	} /* if(fi = fopen(input, "r")) { */
	else {
		fprintf(stderr, "Can't open a file %s\n", input);
		print_usage();
		exit(-1);
	}

	return 0;
}


/***********************************************************************************
 Copyright (C) 2008 Massachusetts Institute of Technology, Cambridge, MA
 All rights reserved

							focalc.c

  Find a focal curve for an input open NURBS curve

	Make:
	make focalc

	Run:
	focalc -i input_curve_file  -n number_of_segments_per_knot_span  -s focal_scale_factor  -o output_file

      Output: 
       Resulting focal curve in .VECT format
 
      Example:
           prompt> focalc -i c.CURV  -n 50  -s 1.0  -o focalc.VECT
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
	fprintf(stderr, "usage: focalc  -i input_curve_file	-n number_of_segments_per_knot_span	-s focal_scale_factor	-o output_file\n");

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


void StoreCurvValues(ParCurv *egeom, int nsegspan, int is_open, int is_planar,
		     float ***pts, float ***norms, double **kcurv,
		     short *nsegs)
/* A user-input "nsegspan" is the number of straight line segments per knot span 
   used to draw each internal knot span of the curve */
{
  vector *nvect, *r, *tang;
  double t, t0, t1, dt, kc, s;
  int i, j, k, nspans;

  nspans = egeom->ncontpts - egeom->order + 1; 
 
  *nsegs = nspans*nsegspan + 1;

  (*pts)   = flt_array2(*nsegs, 3);
  (*norms) = flt_array2(*nsegs, 3);
  (*kcurv) = dbl_array1(*nsegs);

  for (j=k=0; j<nspans; j++) {
    
    t0 = egeom->knots[egeom->order + j - 1]; 
    t1 = egeom->knots[egeom->order + j];    
    
    dt = t1 - t0;
    for (i = (j ? 1 : 0); i<=nsegspan; i++, k++) {
      t = t0 + dt*(double)i/(double)nsegspan;
      r = GetCurvValues(egeom, t, is_open, is_planar, &tang, &nvect, &kc, &s);

      (*pts)[k][0] = r->x/r->w;
      (*pts)[k][1] = r->y/r->w;
      (*pts)[k][2] = r->z/r->w;
      vectfree(r);
      vectfree(tang);

      (*norms)[k][0] = nvect->x/nvect->w;
      (*norms)[k][1] = nvect->y/nvect->w;
      (*norms)[k][2] = nvect->z/nvect->w;
      vectfree(nvect);

      (*kcurv)[k] = kc;
    }
  }
}


void OutputCurvFocal(short nsegs, float scal, float **pts, float **norms,
		   double *kcurv, FILE *fo)
{
  float v[3], k;
  short i;

  fprintf(fo, "M ");
  for (i=0; i<nsegs; i++) {
    if (kcurv[i] != 0.0)
      k = scal/kcurv[i];
    else
      k = 1.0e30;
    v[0] = pts[i][0] + k*norms[i][0];
    v[1] = pts[i][1] + k*norms[i][1];
    v[2] = pts[i][2] + k*norms[i][2];
	
	if(i == 0)
		fprintf(fo, "%+.16le %+.16le %+.16le\n", v[0],v[1],v[2]);
	else
		fprintf(fo, "D %+.16le %+.16le %+.16le\n", v[0],v[1],v[2]);

  }
  fprintf(fo, "E\n");
}



int main(int argc, char* argv[])
{
	int is_planar;
	double x, y, z;

	int nsegspan;
	float **pts, **norms;	double *kcurv;
	float focalScale;
	short nsegs;

	FILE *fi, *fo;  
	ParCurv* egeom;
	char c, *input, *num_seg_span, *focal_scale, *output;

	short errflg = 0;
	input = num_seg_span = focal_scale = output = 0;

	if (argc > 1) {
		while ((c = getopt(argc, argv, "i:n:s:o:")) != -1)
			switch (c) {
			case 'i':
				input = optarg;
				break;

			case 'n':
				num_seg_span = optarg;
				break;

			case 's':
				focal_scale = optarg;
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

		if(num_seg_span) { /* read the # of segments per knot span */
			nsegspan = atoi(num_seg_span); 
			if(nsegspan < 1) { /* check nsegspan range */
				fprintf(stderr, "Number of segments per knot span should be greater than zero\n");
				print_usage();
			    exit(-1);
			}

			if(focal_scale) {
				focalScale = atof(focal_scale);

				is_planar = CheckCurvPlanar(egeom, &x,&y,&z);

				StoreCurvValues(egeom, nsegspan, IS_OPEN, is_planar, &pts, &norms, &kcurv, &nsegs);

				if(output) {
					if(fo = fopen(output, "w")) {

						OutputCurvFocal(nsegs, focalScale, pts, norms, kcurv, fo);

						fclose(fo);
					}/* if(fo = fopen(output, "w")) { */
				}/* if(output) */


			} /* if(focal_scale) */
			else {
				fprintf(stderr, "Can't read the focal scale\n");
				print_usage();
				exit(-1);
			}


		} /* if(num_seg_span) */
		else {
			fprintf(stderr, "Can't read the number of segments per knot span\n");
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


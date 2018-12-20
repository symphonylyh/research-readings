/***********************************************************************************
 Copyright (C) 2008 Massachusetts Institute of Technology, Cambridge, MA
 All rights reserved

							orthos.c

  Find an orthotomic surface for an input open NURBS surface w.r.t. a point P

	Make:
	make orthos

	Run:
	orthos -i input_surface_file   -x x_coord_of_P -y y_coord_of_P -z z_coord_of_P
		   -m number_of_segments_per_u-knot_span -n number_of_segments_per_v-knot_span
	       -s scale_factor    -o output_file

    Output: 
    Resulting orthotomic surface in .VECT format
 
    Example:
           prompt> orthos -i s.SURF  -x 0.0 -y 0.0 -z 0.0  -m 20 -n 20  -s 2.0  -o orthos.VECT
  *********************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include "gen.h"
#include "bspl.h"
#include "praxiteles.h" 
#include "iges.h" 
#include "editor.h"

#define IS_OPEN 1
#define OPEN 0

static double curvatureZero = 1.0e-10;
static double distanceZero  = 1.0e-10;
static double knotZero      = 1.0e-10;

void print_usage()
{
	fprintf(stderr, "usage: orthos -i input_surface_file  -x x_coord_of_P -y y_coord_of_P -z z_coord_of_P  -m number_of_segments_per_u-knot_span -n number_of_segments_per_v-knot_span  -s scale_factor  -o output_file\n");

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


/* Evaluate position, tangents, normals, curvatures of an open NURBS surface */
vector *GetSurfValues(ParSurf *fgeom, float u, float v, vector **tangu,
		      vector **tangv, vector **norm, double *kmax,
		      double *kmin, double *normu, double *normv)
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


void StoreSurfValues(ParSurf *fgeom, short nsegspanu, short nsegspanv,
		     float ****pts, float ****norms, double ***kmax,
		     double ***kmin, double ***normu, double ***normv,
		     short *nsegu, short *nsegv, int flag)
/* A user-input "nsegspanu(v)" is the number of straight line segments per u(v)_knot span 
   used to draw each internal knot span of the surface */
{
  vector *r, *norm, *tangu, *tangv;
  double du, dv, u, v, u0, u1, v0, v1;
  double kx, kn, nu, nv;
  int i, j, k, m, nspanu, nspanv, s, t;

  if (flag) {
    *nsegu = nsegspanu;
    *nsegv = nsegspanv;
  }
  else {
    nspanu = fgeom->ucontpts - fgeom->uorder + 1;
    nspanv = fgeom->vcontpts - fgeom->vorder + 1;
    *nsegu = nspanu*(nsegspanu-1) + 1;
    *nsegv = nspanv*(nsegspanv-1) + 1;
  }

  (*pts)   = flt_array3(*nsegu, *nsegv, 3);
  (*norms) = flt_array3(*nsegu, *nsegv, 3);
  (*kmin)  = dbl_array2(*nsegu, *nsegv);
  (*kmax)  = dbl_array2(*nsegu, *nsegv);
  (*normu) = dbl_array2(*nsegu, *nsegv);
  (*normv) = dbl_array2(*nsegu, *nsegv);

  if (flag) {
    u0 = fgeom->uknots[fgeom->uorder-1];
    u1 = fgeom->uknots[fgeom->ucontpts];
    v0 = fgeom->vknots[fgeom->vorder-1];
    v1 = fgeom->vknots[fgeom->vcontpts];
    du = u1 - u0;
    dv = v1 - v0;

    for (i=0; i<*nsegu; i++) {
      u = u0 + du*(double)i/(double)(*nsegu-1);
      for (j=0; j<*nsegv; j++) {
	v = v0 + dv*(double)j/(double)(*nsegv-1);
	r = GetSurfValues(fgeom, u, v, &tangu, &tangv, &norm, &kx, &kn,
			  &nu, &nv);
	(*pts)[i][j][0] = r->x/r->w;
	(*pts)[i][j][1] = r->y/r->w;
	(*pts)[i][j][2] = r->z/r->w;
	vectfree(r);
	vectfree(tangu);
	vectfree(tangv);

	(*norms)[i][j][0] = norm->x/norm->w;
	(*norms)[i][j][1] = norm->y/norm->w;
	(*norms)[i][j][2] = norm->z/norm->w;
	vectfree(norm);

	(*kmin)[i][j] = kn;
	(*kmax)[i][j] = kx;
	(*normu)[i][j] = nu;
	(*normv)[i][j] = nv;
      }
    }
  }
  else {
    for (i=s=0; i<nspanu; i++) {
      u0 = fgeom->uknots[fgeom->uorder + i - 1];
      u1 = fgeom->uknots[fgeom->uorder + i];
      du = u1 - u0;
      for (k = (i ? 1 : 0); k<nsegspanu; k++, s++) {
	u = u0 + du*(double)k/(double)(nsegspanu-1);

	for (j=t=0; j<nspanv; j++) {
	  v0 = fgeom->vknots[fgeom->vorder + j - 1];
	  v1 = fgeom->vknots[fgeom->vorder + j];
	  dv = v1 - v0;
	  for (m = (j ? 1 : 0); m<nsegspanv; m++, t++) {
	    v = v0 + dv*(double)m/(double)(nsegspanv-1);
	    r = GetSurfValues(fgeom, u, v, &tangu, &tangv, &norm, &kx, &kn,
			      &nu, &nv);
	    (*pts)[s][t][0] = r->x/r->w;
	    (*pts)[s][t][1] = r->y/r->w;
	    (*pts)[s][t][2] = r->z/r->w;
	    vectfree(r);
	    vectfree(tangu);
	    vectfree(tangv);

	    (*norms)[s][t][0] = norm->x/norm->w;
	    (*norms)[s][t][1] = norm->y/norm->w;
	    (*norms)[s][t][2] = norm->z/norm->w;
	    vectfree(norm);

	    (*kmin)[s][t] = kn;
	    (*kmax)[s][t] = kx;
	    (*normu)[s][t] = nu;
	    (*normv)[s][t] = nv;
	  }
	}
      }
    }
  }
}
 

void OutputSurfOrtho(ParSurf *fgeom, short nsegu, short nsegv, short nsubu,
		   short nsubv, float scal, float *pt, float ***pts,
		   float ***norms, FILE *fo)
{
  vector *r, *norm, *tangu, *tangv;
  double du, dv, u, v, u0, u1, v0, v1;
  double kn, kx, nu, nv;
  float d, np[3], p[3], xp[3];
  short i, j, k;

  u0 = fgeom->uknots[fgeom->uorder-1];
  u1 = fgeom->uknots[fgeom->ucontpts];
  v0 = fgeom->vknots[fgeom->uorder-1];
  v1 = fgeom->vknots[fgeom->vcontpts];
  du = u1 - u0;
  dv = v1 - v0;
  for (i=0; i<nsegu; i++) {
    u = u0 + du*(double)i/(double)(nsegu-1);

	fprintf(fo, "M ");

    for (j=0; j<nsegv; j++) {
      Sub(3, pts[i][j], pt, xp);
      d = Dot(xp, norms[i][j]);
      Scale(3, norms[i][j], scal*d, p);
      Add(3, pt, p, p);

	  if(j == 0)
		fprintf(fo, "%+.16le %+.16le %+.16le\n", p[0],p[1],p[2]);
	  else
		fprintf(fo, "D %+.16le %+.16le %+.16le\n", p[0],p[1],p[2]);

      if (j<nsegv-1)
	for (k=1; k<nsubv; k++) {
	  v = v0 + dv*(double)(j*nsubv+k)/(double)((nsegv-1)*nsubv);
	  r = GetSurfValues(fgeom, u, v, &tangu, &tangv, &norm, &kx, &kn,
			    &nu, &nv);
	  p[0] = r->x/r->w;
	  p[1] = r->y/r->w;
	  p[2] = r->z/r->w;
	  Sub(3, p, pt, xp);
	  np[0] = norm->x/norm->w;
	  np[1] = norm->y/norm->w;
	  np[2] = norm->z/norm->w;
	  d = Dot(xp, np);
	  Scale(3, np, scal*d, p);
	  Add(3, pt, p, p);

	  fprintf(fo, "D %+.16le %+.16le %+.16le\n", p[0],p[1],p[2]);

	  free(r);
	  free(tangu);
	  free(tangv);
	  free(norm);
	}
    }

  }
  for (j=0; j<nsegv; j++) {
    v = v0 + dv*(double)j/(double)(nsegv-1);

	fprintf(fo, "M ");
    for (i=0; i<nsegu; i++) {
      Sub(3, pts[i][j], pt, xp);
      d = Dot(xp, norms[i][j]);
      Scale(3, norms[i][j], scal*d, p);
      Add(3, pt, p, p);

	  if(i == 0)
		fprintf(fo, "%+.16le %+.16le %+.16le\n", p[0],p[1],p[2]);
	  else
		fprintf(fo, "D %+.16le %+.16le %+.16le\n", p[0],p[1],p[2]);

      if (i<nsegu-1)
	for (k=1; k<nsubu; k++) {
	  u = u0 + du*(double)(i*nsubu+k)/(double)((nsegu-1)*nsubu);
	  r = GetSurfValues(fgeom, u, v, &tangu, &tangv, &norm, &kx, &kn,
			    &nu, &nv);
	  p[0] = r->x/r->w;
	  p[1] = r->y/r->w;
	  p[2] = r->z/r->w;
	  Sub(3, p, pt, xp);
	  np[0] = norm->x/norm->w;
	  np[1] = norm->y/norm->w;
	  np[2] = norm->z/norm->w;
	  d = Dot(xp, np);
	  Scale(3, np, scal*d, p);
	  Add(3, pt, p, p);

	  fprintf(fo, "D %+.16le %+.16le %+.16le\n", p[0],p[1],p[2]);

	  free(r);
	  free(tangu);
	  free(tangv);
	  free(norm);
	}
    }

  }
  fprintf(fo, "E\n");
}



int main(int argc, char* argv[])
{
	int nsegspanu, nsegspanv;
	float pt[3];
	float ***pts, ***norms;	double **kmax, **kmin, **normu, **normv;
	float orthoScale;
	short nsegu, nsegv,  nsubu, nsubv;

	FILE *fi, *fo;  
	ParSurf* egeom;
	char c, *input, *px, *py, *pz, *num_seg_u, *num_seg_v, *ortho_scale, *output;

	short errflg = 0;
	input = px = py = pz = num_seg_u = num_seg_v = ortho_scale = output = 0;

	if (argc > 1) {
		while ((c = getopt(argc, argv, "i:x:y:z:m:n:s:o:")) != -1)

			switch (c) {
			case 'i':
				input = optarg;
				break;

			case 'x':
				px = optarg;
				break;

			case 'y':
				py = optarg;
				break;

			case 'z':
				pz = optarg;
				break;

			case 'm':
				num_seg_u = optarg;
				break;

			case 'n':
				num_seg_v = optarg;
				break;

			case 's':
				ortho_scale = optarg;
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
		egeom = ReadParSurf(fi, NULL);

		if(px && py && pz) {
			pt[0] = atof(px);	pt[1] = atof(py);	pt[2] = atof(pz); /* input constant point P */

		if(num_seg_u) { /* read the # of segments per u-knot span */
			nsegspanu = atoi(num_seg_u); 
			if(nsegspanu < 1) { /* check nsegspanu range */
				fprintf(stderr, "Number of segments per knot span should be greater than zero\n");
				print_usage();
			    exit(-1);
			}

		if(num_seg_v) { /* read the # of segments per v-knot span */
			nsegspanv = atoi(num_seg_v); 
			if(nsegspanv < 1) { /* check nsegspanv range */
				fprintf(stderr, "Number of segments per knot span should be greater than zero\n");
				print_usage();
			    exit(-1);
			}

			if(ortho_scale) {
				orthoScale = atof(ortho_scale);

				StoreSurfValues(egeom, nsegspanu, nsegspanv, &pts, &norms, &kmax, &kmin,
					            &normu, &normv, &nsegu, &nsegv, OPEN);


				if(output) {
					if(fo = fopen(output, "w")) {
						nsubu = nsubv = 1;
						OutputSurfOrtho(egeom, nsegu, nsegv, nsubu, nsubv, orthoScale, pt, 
							         pts, norms, fo);

						fclose(fo);
					}/* if(fo = fopen(output, "w")) { */
				}/* if(output) */


			} /* if(ortho_scale) */
			else {
				fprintf(stderr, "Can't read the orthotomic scale\n");
				print_usage();
				exit(-1);
			}


		} /* if(num_seg_v) */
		else {
			fprintf(stderr, "Can't read the number of segments per v-knot span\n");
			print_usage();
			exit(-1);
		}

		} /* if(num_seg_u) */
		else {
			fprintf(stderr, "Can't read the number of segments per u-knot span\n");
			print_usage();
			exit(-1);
		}

		} /* if(px && py && pz) */
		else {
			fprintf(stderr, "Can't read x, y, z coordinates of the point P\n");
			print_usage();
			exit(-1);
		}

		fclose(fi);
		free_fgeom(egeom);

	} /* if(fi = fopen(input, "r")) { */
	else {
		fprintf(stderr, "Can't open a file %s\n", input);
		print_usage();
		exit(-1);
	}

	return 0;
}


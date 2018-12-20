/************************************************************************
 *									*
                       Copyright (C) 1992 by
	Massachusetts Institute of Technology, Cambridge, MA
			 All rights reserved

 *									*
 ************************************************************************/
/*   Filename:  sample_blend1.c

     Written by:  Peter C. Filkins
     Date:        25 February 1991

 ***********************************************************************
     Description:  A modification to sample_blend.c to add a capability
                   to sample the blend surface for curvature continuity
		   as well as tangent continuity.  Introduces some of
		   the blend surface input control
 ***********************************************************************/
#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include "gen.h"
#include "bspl.h"
#include "appr.h"
#include "editor.h"

/********************
#include <stdlib.h>
********************/
void qsort(void *, size_t, size_t, int (*)(double *, double *));
/*******************/

#define M_PI 3.14159265358979323846
#define RAD ((double) M_PI/180.0)
#define KNOTEPS 1.0e-05
#define CDIFF 1.0e-02
#define MIN_CURV 1.0e-5

short sample_blend(ParCurv *(*fn)(double, ParSurf **, ParCurv **, short,
				  short, double, double), ParSurf *fgb,
		   double *knots, double *eps, short *n, short nn,
		   short bc1, short bc2, ParSurf **fgm, ParCurv **egm,
		   double bias1, double bias2)
{
  double u,t,*k,delta,angle,*tau,dist[6];
  int i,nold,numaddpts = 0,ncontpts;
  ParCurv *iso, *test;
  vector *v, *v1, *norm;
  double *f, k_err, max_dist = 0.0, max_curv = 0.0;
  double dudv[2], max_pos = 0.0;

  if ((*n) > nn)
    return (0);

  angle = RAD*eps[1];
  k_err = CDIFF*eps[2];

  k = dbl_array1(*n);
  tau = dbl_array1(*n);
  nodes_bsp(knots, *n, 4, tau) ;
  dudv[0] = 0.0;
  dudv[1] = 1.0;

  for (i=0; i < (*n)-1; i++) {
    t = (tau [i+1] + tau [i]) / 2;
    test = fn(t, fgm, egm, bc1, bc2, bias1, bias2);

               /* test for tangent continuity first.  Don't test for
		  curvature continuity until tangent is sat....*/

    if (bc1 > 0) {
      v = rbspeval(test, 0.0, 1);
      dist[0] = test_tangent(fgb, v, t, 0.0);
      vectfree(v);
    }
    if (bc2 > 0) {
      v = rbspeval(test, 1.0, 1);
      dist[1] = test_tangent(fgb, v, t, 1.0);
      vectfree(v);
    }
    dist[0] = MAX(dist[0], dist[1]);

    if (dist[0] > max_dist)
      max_dist = dist[0];

    if (dist[0] > angle) {
      numaddpts++;
      k[numaddpts] = t;
      printf("Distance = %f, Angle = %f\n", dist[0], angle);
    }
    free_egeom(test);
  }

  printf("Maximum mid-point angle = %+e\n", max_dist);

  if ((bc1 == 2 || bc2 == 2) && !numaddpts) {
    iso = egeomalloc1(fgb->uorder, fgb->ucontpts);
    for (i=0; i < (*n)-1; i++) {
      t = (tau [i+1] + tau [i]) / 2;
      test = fn(t, fgm, egm, bc1, bc2, bias1, bias2);
      ParCurv_iso(fgb, t, 0, iso);
      if (bc1 == 2) dist[2] = test_curvature0(test, iso);
      if (bc2 == 2) dist[3] = test_curvature1(test, iso);
      dist[2] = MAX(dist[2], dist[3]);
      if (dist[2] > max_curv) max_curv = dist[2];
      if (dist[2] > k_err) {
	numaddpts++;
	k[numaddpts] = t;
	printf("Distance = %f, Percent = %f\n", dist[2], k_err);
      }
      free_egeom(test);
    } 
    printf("Maximum curvature error = %+e\n", max_curv);
    free_egeom(iso);
  }

  if (numaddpts) {     /* We have not met tolerance */
    nold = (*n);
    *n = numaddpts+nold;
    if ((*n) > nn) {
      free_darray1(k);
      return (0);
    }
    for (i=nold; i < (*n); i++)
      knots[i+4] = k[i-nold+1];
    qsort(knots, (*n)+4, sizeof(double), compar);
    for (i=4; i< (*n) - 1; i++)
      if (knots[i+1] - knots[i] < KNOTEPS) {
	knots[i] = knots[i-1] + (knots[i+2] - knots[i-1]) /3.0;
	knots[i+1] = knots[i-1] + 2.0*(knots[i+2] - knots[i-1]) /3.0;
      }
  }

  printf("Number of cross-link curves added = %d\n", numaddpts);
  fflush(stdout);
  free_darray1(k);
  free_darray1(tau);

  if (!numaddpts)
    return (0);

  return(1);
}

double test_tangent(ParSurf *fgm, vector *v0, double u, double v)
{
  vector *v1;
  double angle;
  
  v1 = revalderivsurf(fgm, u, v, 0, 1);
  
  angle = vec_angle(v0, v1);
  vectfree(v1);
  
  if (angle > M_PI)
    angle = 0.0;
  else if (angle > M_PI/2.0)
    angle = M_PI - angle;
  
  return (angle);
}

double test_curvature0(ParCurv *clink, ParCurv *iso)
{
  double ks, kb, diff;

  ks = iso_curv0(clink);
  kb = iso_curv0(iso);
  
  if ((diff = fabs(ks - kb)) > MIN_CURV)
    diff = diff / (fabs(kb));

  return (diff);
}

double test_curvature1(ParCurv *clink, ParCurv *iso)
{
  double ks, kb, diff;

  ks = iso_curv1(clink);
  kb = iso_curv1(iso);
  
  if ( (diff = fabs(ks - kb)) > MIN_CURV)
    diff = diff / (fabs(kb));

  return (diff);
}

double iso_curv0(ParCurv *iso)
{
  vector *r1,*r2;
  double diff, *f, r13, temp, k;

  f = knot_factor0(iso->order, iso->knots);
  r1 = sub_vect(iso->contpts[1], iso->contpts[0]);
  scale_vect1(f[0], r1, r1);
  r2 = scale_vect(f[1], r1);
  sub_vect1(iso->contpts[2], r2, r2);
  sub_vect1(r2,iso->contpts[0], r2);
  scale_vect1((1.0/f[2]), r2, r2);
  temp = mag(r1);
  r13 = temp*temp*temp;
  k = mag(cross(r1, r2)) / temp;

  free(f);
  vectfree(r1);
  vectfree(r2);

  return (k);
}

double iso_curv1(ParCurv *iso)
{
  vector *r1, *r2;
  double *f, k, temp, r13; 
  short m;

  m = iso->ncontpts; 
  f = knot_factor1(iso->order, iso->ncontpts, iso->knots);
  r1 = sub_vect(iso->contpts[m-2], iso->contpts[m-1]);
  scale_vect1(f[0], r1, r1);
  r2 = scale_vect(f[1], r1);
  sub_vect1(iso->contpts[m-3], r2, r2);
  sub_vect1(r2,iso->contpts[m-1], r2);
  scale_vect1((1.0/f[2]), r2, r2);
  temp = mag(r1);
  r13 = temp*temp*temp;
  k = mag(cross(r1, r2)) / temp;

  free(f);
  vectfree(r1);
  vectfree(r2);

  return (k);
}

double position_err(ParCurv *test, ParSurf *fgb, double u, double v)
{
  double diff;
  vector *v1, *v2;

  v1 = rbspeval(test, v, 0);
  v2 = revalderivsurf(fgb, u, v, 0, 0);

  diff = distance(v1, v2);
  vectfree(v1);
  vectfree(v2);

  return (diff);
}

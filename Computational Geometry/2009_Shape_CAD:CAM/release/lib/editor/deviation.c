/* Copyright (C) 1997 Massachusetts Institute of Technology, Cambridge, MA
   All rights reserved
*/

/* deviation.c */

/* CurvDeviation()
*/

#include <stdlib.h>
#include <math.h>
#include "gen.h"
#include "bspl.h"
#include "editor.h"

void e02baf_(int *, int *, double *, double *, double *, double *, double *,
	     double *, double *, double *, int *);

ParCurv *CurvDeviation(ParCurv *egeom1, ParCurv *egeom2, int kSample,
		       int nSegs, int is_planar, int cubic)
{
  ParCurv *egeom, *egm[2];
  vector *norm, *r1, *r2, *tang, *v;
  double a, b, dt, epsilon, *knots, **lim, *pts, t, t0, t1;
  int index[3], ncontpts, nknots, nPts, nSpans, order;
  short i, j, k;

  /* make sure that both curves have the same order */

  epsilon = 1.0e-15;
  order = MAX(egeom1->order, egeom2->order);
  ncontpts = MAX(egeom1->ncontpts, egeom2->ncontpts);
  nknots = order + ncontpts;
  knots = dbl_array1(nknots*order);

  egm[0] = egeomalloc1(order, order*(order+ncontpts));
  egm[1] = egeomalloc1(order, order*(order+ncontpts));
  
  copyegeom(egeom1, egm[0]);
  copyegeom(egeom2, egm[1]);

  /* raise both curves to same order */

  if (order > egm[0]->order) {
    lim = dbl_array2((order+ncontpts)*order, ncontpts*order);
    for (i=0; i < order - egeom1->order; i++) {
      k = compute_lim(egm[0], lim, knots, epsilon);
      raise_byone(egm[0], lim, knots, k);
    }
    free_darray2(lim);
  }

  if (order > egm[1]->order) {
    lim =dbl_array2((order+ncontpts)*order, ncontpts*order);
    for(i=0; i < order - egeom2->order; i++) {
      k = compute_lim(egm[1], lim, knots, epsilon);
      raise_byone(egm[1], lim, knots, k);
    }
    free_darray2(lim);
  }

  /* merge knot vectors */

  nknots = merge_knotv(egm, 0, 2, knots, index, epsilon);
  addpoints(egm, 0, 2, knots, nknots);
  free_darray1(knots);

  if (kSample) {            /* uniform sample over each internal knot span */
    nSpans = egm[0]->ncontpts - egm[0]->order + 1;
    nPts = nSpans*(nSegs-1) + 1;
  }
  else                              /* uniform sample over parameter space */
    nPts = nSegs;

  pts = dbl_array1(nPts);           /* allocate array to hold sample points */
  knots = dbl_array1(nPts);         /* array to hold sample locations */

  /* sample the points on both curves and find positional deviation */

  if (kSample) {   /* uniform sample over each internal knot span */
    for (i=k=0; i<nSpans; i++) {
      t0 = egm[0]->knots[egm[0]->order + i - 1];
      t1 = egm[0]->knots[egm[0]->order + i];
      dt = t1 - t0;
      for (j=(i ? 1 : 0); j<nSegs; j++,k++) {
	t = t0 + dt*(double)j/(double)(nSegs-1);
	knots[k] = t;
	r1 = evalbsp(egm[0], t);
	r2 = evalbsp(egm[1], t);
	pts[k] = distance(r1, r2);
	if (is_planar) {
	  tang = rbspeval(egm[0], t, 1);

	  a = tang->x*tang->x + tang->y*tang->y;
	  b = sqrt(a);

	  norm = vectalloc();      /* unit vector on curve 1 */
	  norm->x = -tang->y/b;
	  norm->y =  tang->x/b;
	  vectfree(tang);
	  unitvector1(norm, norm);

	  v = sub_vect(r2, r1);     /* deviation vector between curves 1 & 2 */
	  unitvector1(v, v);

	  if (dot(norm, v) < 0.0)   /* if deviation vector is on opposite */
	    pts[k] *= -1.0;         /* side of normal, deviation is negative */

	  vectfree(norm);
	  vectfree(v);
	}
	vectfree(r1);
	vectfree(r2);
      }
    }
  }
  else {           /* uniform sample over parameter space */
    t0 = egm[0]->knots[egm[0]->order-1];
    t1 = egm[0]->knots[egm[0]->ncontpts];
    dt = t1 - t0;
    for (i=0; i<nPts; i++) {
      t = t0 + dt*(double)i/(double)(nPts-1);
      knots[i] = t;
      r1 = evalbsp(egm[0], t);
      r2 = evalbsp(egm[1], t);
      pts[i] = distance(r1, r2);
      if (is_planar) {
	tang = rbspeval(egm[0], t, 1);
	
	a = tang->x*tang->x + tang->y*tang->y;
	b = sqrt(a);

	norm = vectalloc();         /* unit vector on curve 1 */
	norm->x = -tang->y/b;
	norm->y =  tang->x/b;
	vectfree(tang);
	unitvector1(norm, norm);

	v = sub_vect(r2, r1);       /* deviation vector between curves 1 & 2 */
	unitvector1(v, v);

	if (dot(norm, v) < 0.0)     /* if deviation vector is on opposite */
	  pts[k] *= -1.0;           /* side of normal, deviation is negative */

	vectfree(norm);
	vectfree(v);
      }
      vectfree(r1);
      vectfree(r2);
    }
  }
  free_egeom(egm[0]);
  free_egeom(egm[1]);

  /* fit B-spline through points */
  egeom = FitFunction(nPts, pts, NULL, NULL, knots, 0);
  free_darray1(pts);
  free_darray1(knots);

  return egeom;
}

ParCurv *FitFunction(int m, double *t1, double *t2, double *t3, double *knots,
		     int cubic)
{
  ParCurv *egeom = NULL;
  double *x, *y, *w, *lamda, *work1, *work2, s, s0, ratio, err, max_error;
  double *c1, *c2 = NULL, *c3 = NULL;
  double accuracy = 0.01, RATIO_LIMIT = 0.6;
  int i, ifail, ncap7;

  if (!cubic) {        /* fit linear B-spline through points */
    egeom = egeomalloc1(2, m);
    for (i=0; i<m; i++) {
      egeom->contpts[i]->x = t1[i];
      if (t2) {
	egeom->contpts[i]->y = t2[i];
	if (t3)
	  egeom->contpts[i]->z = t3[i];
	else
	  egeom->contpts[i]->z = (double)i/(double)(m-1);
      }
      else {
	egeom->contpts[i]->y = (double)i/(double)(m-1);
	egeom->contpts[i]->z = 0.0;
      }
      egeom->contpts[i]->w = 1.0;
      
      if (!knots)
	if (i) {                        /* accumulate total arc length */
	  s0 += distance(egeom->contpts[i], egeom->contpts[i-1]);
	}
	else
	  s0 = 0.0;
    }

    egeom->knots[0] = 0.0;            /* build knot vector */
    for (i=0; i<m; i++) {
      if (knots) {
	egeom->knots[i+1] = knots[i];
      }
      else {
	if (i) {                        /* accumulate total arc length */
	  s += distance(egeom->contpts[i], egeom->contpts[i-1]);
	}
	else
	  s = 0.0;
	egeom->knots[i+1] = s/s0;
      }
      egeom->knots[m+1] = 1.0;
    }
  }
  else {    /* fit cubic B-spline through points */
    x        = dbl_array1(m);         /* holds x-coordinate (uniform) */
    y        = dbl_array1(m);         /* holds y-value (copy from t) */
    w        = dbl_array1(m);         /* weights */
    lamda    = dbl_array1(m+7);       /* knot vector (uniform) */
    work1    = dbl_array1(m);         /* workbench */
    work2    = dbl_array1(4*(m+7));   /* another workbench?!? */
    c1       = dbl_array1(m+7);       /* holds return values */
    if (t2) {
      c2     = dbl_array1(m+7);       /* holds return values */
      if (t3)
	c3   = dbl_array1(m+7);       /* holds return values */
    }

    if (m<20)      /* ratio depends on num. Points */
      ratio=0.5;
    else
      ratio=0.25;

    for (i=0; i<m; i++) {  /* set x-axis uniformly and no weights */
      x[i] = i/(m-1.);
      w[i] = (-4*i/(m-1.0)*(1-i/(m-1.0))+1.0)*4+1;  /* parabola 5..1..5 */
    }

    do {
      ncap7     = m*ratio + 7;               /* set num. of points */
      if (egeom)
	free_egeom(egeom);                   /* every turn new egeom */
      egeom     = egeomalloc1(4, ncap7-4);
      ifail     = -1;

      for (i=0; i< ncap7-6; i++)             /* uniform knot vector */
	lamda[i+3] = i/(ncap7-7.0);

      for (i=0; i<m; i++)
	y[i] = t1[i];
      e02baf_(&m, &ncap7, x, y, w, lamda, work1, work2, c1, &s, &ifail);
      max_error = sqrt(s/ncap7);

      if (t2) {
	for (i=0; i<m; i++)
	  y[i] = t2[i];
	e02baf_(&m, &ncap7, x, y, w, lamda, work1, work2, c2, &s, &ifail);
	err = sqrt(s/ncap7);
	if (err > max_error) max_error = err;

	if (t3) {
	  for (i=0; i<m; i++)
	    y[i] = t3[i];
	  e02baf_(&m, &ncap7, x, y, w, lamda, work1, work2, c3, &s, &ifail);
	  err = sqrt(s/ncap7);
	  if (err > max_error) max_error = err;
	}
      }
    
      if (max_error*5 <= accuracy)
	break;
      ratio += (1-ratio)*0.2;                /* increment ratio slowly */

      if (ratio>RATIO_LIMIT)
	break;
    } while (1);

    for (i=0; i< ncap7-4; i++) {
      egeom->contpts[i]->x = c1[i];
      if (t2) {
	egeom->contpts[i]->y = c2[i];
	if (t3)
	  egeom->contpts[i]->z = c3[i];
	else
	  egeom->contpts[i]->z = (double)i/(double)(ncap7-5);
      }
      else {
	egeom->contpts[i]->y = (double)i/(double)(ncap7-5);
	egeom->contpts[i]->z = 0.0;
      }
      egeom->contpts[i]->w = 1.0;
    }

    for (i=0; i<ncap7; i++)
      egeom->knots[i] = lamda[i];

    free_darray1(x);                         /* clean up (the mess) */
    free_darray1(y);
    free_darray1(w);
    free_darray1(lamda);
    free_darray1(work1);
    free_darray1(work2);
    free_darray1(c1);
    if (t2) {
      free_darray1(c2);
      if (t3)
	free_darray1(c3);
    }
  }

  return egeom;
}

void SurfDeviation(ParSurf *fgeom1, ParSurf *fgeom2, int kSample, int nu,
		   int nv, int *nPts, double ***pts, int **iEnd, int *nAdj,
		   int **iAdj)
{
  ParSurf *fg[2];
  ParCurv *eg[2];
  vector *norm, *r, *r1, *r2;
  double dis, du, dv, *kn[2], u, u0, u1, v, v0, v1;
  int el, i, is, j, k, n, np2, nps, nspan2, nspans, p, q, s, t, *tn;
  int nmax = 500;

  kn[0] = dbl_array1(nmax);
  kn[1] = dbl_array1(nmax);
  tn = int_array1(nmax);

  for(j=0; j<2; j++) {   /* make copy of original surfaces */
    fg[j] = fgeomalloc2(fgeom1->uorder, fgeom1->vorder, nmax, nmax);
    eg[j] = egeomalloc2(fgeom1->uorder, nmax);
    for(i=0; i<fgeom2->ucontpts; i++)
      eg[j]->contpts[i] = vectalloc();
  }

  is = fgeom1->uorder-fgeom2->uorder;
  el = fgeom1->vorder-fgeom2->vorder;
  if (is > -1 && el > -1) {
    copyfgeom(fgeom2, fg[0]);
    copyfgeom(fgeom1, fg[1]);
  }
  else {
    copyfgeom(fgeom2, fg[1]);
    copyfgeom(fgeom1, fg[0]);
    is = abs(is);
    el = abs(el);
  }

  raise_surf(fg[0], ZERO, is, el);       /* raise surfaces to same degree */
  merge_fgeom(fg[0], fg[1], eg, kn[0], kn[1], tn, ZERO);   /* merge knots */
  free_egeom(eg[0]);
  free_egeom(eg[1]);

  free_darray1(kn[0]);
  free_darray1(kn[1]);

  free_iarray1(tn);

  if (kSample) {
    nspans = fg[0]->ucontpts - fg[0]->uorder + 1;
    nspan2 = fg[0]->vcontpts - fg[0]->vorder + 1;
    nps = nspans*(nu-1) + 1;
    np2 = nspan2*(nv-1) + 1;
    *nPts = nps*np2;
  }
  else
    *nPts = nu*nv;

  *pts = dbl_array2(*nPts, 3);

  if (kSample) {   /* uniform sample over internal knot spans */
    for (j=n=q=0; j<nspans; j++) {
      u0 = fg[0]->uknots[fg[0]->uorder + j - 1];
      u1 = fg[0]->uknots[fg[0]->uorder + j];
      du = u1 - u0;
      for (k=(j ? 1 : 0); k<nu; k++,n++) {
	u = u0 + du*(double)k/(double)(nu-1);
	for (s=p=0; s<nspan2; s++) {
	  v0 = fg[0]->vknots[fg[0]->vorder + s - 1];
	  v1 = fg[0]->vknots[fg[0]->vorder + s];
	  dv = v1 - v0;
	  for (t=(s ? 1 : 0); t<nv; t++,p++,q++) {
	    v = v0 + dv*(double)t/(double)(nv-1);
	    r1 = evalsurf(fg[0], u, v);
	    r2 = evalsurf(fg[1], u, v);
	    (*pts)[q][0] = u;
	    (*pts)[q][1] = v;
	    (*pts)[q][2] = distance(r1, r2);

	    norm = normalsurf(fg[0], u, v); /* normal of surface 1 */
	    unitvector1(norm, norm);
	    r = sub_vect(r2, r1);    /* deviation vector between */
	    unitvector1(r, r);       /* surfaces 1 and 2 */

	    if (dot(norm, r) < 0.0)  /* if deviation vector is on opposite */
	      (*pts)[q][2] *= -1.0;  /* side of norm, deviation is negative */

	    vectfree(norm);
	    vectfree(r);
	    vectfree(r1);
	    vectfree(r2);
	  }
	}
      }
    }
  }
  else {           /* uniform sample over parameter space */
    u0 = fg[0]->uknots[fg[0]->uorder-1];
    u1 = fg[0]->uknots[fg[0]->ucontpts];
    v0 = fg[0]->vknots[fg[0]->vorder-1];
    v1 = fg[0]->vknots[fg[0]->vcontpts];
    du = u1 - u0;
    dv = v1 - v0;
    for (j=n=0; j<nu; j++) {
      u = u0 + du*(double)j/(double)(nu-1);
      for (k=0; k<nv; k++,n++) {
	v = v0 + dv*(double)k/(double)(nv-1);
	r1 = evalsurf(fg[0], u, v);
	r2 = evalsurf(fg[1], u, v);
	(*pts)[n][0] = u;
	(*pts)[n][1] = v;
	(*pts)[n][2] = distance(r1, r2);

	norm = normalsurf(fg[0], u, v);     /* normal of surface 1 */
	unitvector1(norm, norm);
	r = sub_vect(r2, r1);        /* deviation vector between */
	unitvector1(r, r);           /* surfaces 1 and 2 */

	if (dot(norm, r) < 0.0)      /* if deviation vector is on opposite */
	  (*pts)[n][2] *= -1.0;      /* side of norm, deviation is negative */

	vectfree(norm);
	vectfree(r);
	vectfree(r1);
	vectfree(r2);
      }
    }
  }
  free_fgeom(fg[0]);
  free_fgeom(fg[1]);

  *iEnd = int_array1(*nPts);     /* triangulate the points */
  Delaunay(*nPts, *pts, *iEnd, nAdj, iAdj, 0);
}

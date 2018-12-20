/* Copyright (C) 1993 Massachusetts Institute of Technology, Cambridge, MA
   All rights reserved
*/

/* intersect_curve.c */

/* CurveOrientation()
   functcv()
   functcv_sqr()
   intersect_curve()
   min_bez_curve_to_axis_intrsct()
   partl_u_functcv_sqr()
   RecursiveInt()
*/

#include <math.h>
#include <malloc.h>
#include "gen.h"
#include "bspl.h"
#include "editor.h"

static ParCurv *global_egeom;

void e04bbf_(void(double *, double *, double *), double *, double *,
	     double *, double *, int *, double *, double *, double *, int *);

int RecursiveInt(ParCurv *egeom, double vmin, double vmax, int dir)
{
  static count = 0;
  ParCurv *egm;
  double ct;
  int orient;

  if (++count > 4) {
    count--;
    return (-1);
  }

  ct = (vmin + vmax) / 2.0;

  egm = copyegeom(egeom, NULL);
  orient = CurveOrientation(egm, ct, dir);
  free_egeom(egm);

  if (orient > -1) {
    count--;
    return (orient);
  }

  if ((orient = RecursiveInt(egeom, vmin, ct, dir)) > -1) {
    count--;
    return (orient);
  }

  orient = RecursiveInt(egeom, ct, vmax, dir);
  count--;
  return (orient);
}

int CurveOrientation(ParCurv *egeom, double ct, int dir)
{
  ParCurv *egm;
  vector *r;
  double *knots, **inpoly, **outpoly, *pts, tang, umin;
  int nKnots, nPts;
  short i, intersects = 0, j, nc, order;

  order = egeom->order;
  nc    = egeom->ncontpts;

  nKnots = checknot(egeom);

  knots   = dbl_array1(nKnots);
  prodknot(egeom, knots);

  inpoly  = dbl_array2(nKnots, 4);
  outpoly = dbl_array2(nKnots, 4);

  for(i=0; i<nc; i++) {
    inpoly[i][0] = egeom->contpts[i]->x;
    inpoly[i][1] = egeom->contpts[i]->y;
    inpoly[i][2] = egeom->contpts[i]->z;
    inpoly[i][3] = egeom->contpts[i]->w;
  }

  curve_oslo1(egeom->order, egeom->ncontpts, nKnots, egeom->knots, knots,
	      inpoly, outpoly);
  free_darray2(inpoly);

  umin = 1.0e30;
  intersects = 0;
  for (j=0; j*order < nKnots-order; j++) {
    egm = egeomalloc1(order, order);
    egm->order    = order;
    egm->ncontpts = order;

    for(i=0; i<order; i++)
      egm->knots[i] = 0.0;

    for(i=order; i<2*order; i++)
      egm->knots[i] = 1.0;

    for(i=0; i<egm->order; i++) {
      egm->contpts[i]->x = outpoly[i + j*order][0];
      egm->contpts[i]->y = outpoly[i + j*order][1];
      egm->contpts[i]->z = outpoly[i + j*order][2];
      egm->contpts[i]->w = outpoly[i + j*order][3];
    }

    global_egeom = copyegeom(egm, NULL);
    for (i=0; i<global_egeom->ncontpts; i++) {
      global_egeom->contpts[i]->z = global_egeom->contpts[i]->y - ct;
      global_egeom->contpts[i]->y = 0.0;
      global_egeom->contpts[i]->x = (double)i /
	                            (double)(global_egeom->ncontpts-1);
    }
      
    normalize_egeom(global_egeom);

    pts = dbl_array1(EDITOR_MAX_INTERSECT_POINT);
    nPts = 0;
    intersect_curve_to_axis(global_egeom, dir, &nPts, pts);
    if (nPts) {
      intersects = 1;
      for (i=0; i<nPts; i++) {
	r = evalbsp(egm, pts[i]);
	if (r->x/r->w < umin) {
	  umin = r->x/r->w;
	  vectfree(r);
	  r = rbspeval(egm, pts[i], 1);
	  tang = r->y/r->w;
	}
	vectfree(r);
      }
    }
    free_darray1(pts);    

    free_egeom(global_egeom);
    free_egeom(egm);
  }

  free_darray2(outpoly);

  if (intersects)
    return (tang < 0.0 ? 1 : 0);

  return (-1);
}


int intersect_curve_to_axis(ParCurv *egeom, int dir, int *n, double *pts)
{
  ParCurv *kid_1, *kid_2, *clipped_curve;
  double min, max, mid;
  int clip, total_inte;

  clip = bezier_clip_decision_egeom(egeom, dir, &min, &max);

  if (!clip)
    total_inte = 0;
  else if (clip == 1)
    total_inte = min_bez_curve_to_axis_intrsct(min, max, pts, n);
  else if (clip == 2) {
    mid = (egeom->knots[0] + egeom->knots[egeom->kmem-1]) / 2.0;
    decasteljau_curve(egeom, mid, &kid_1, &kid_2);

    total_inte  = intersect_curve_to_axis(kid_1, dir, n, pts);
    total_inte += intersect_curve_to_axis(kid_2, dir, n, pts);

    free_egeom(kid_1);
    free_egeom(kid_2);
  }
  else if (clip == 3) {
    clipped_curve = clip_bezier_curve(egeom, min, max);

    total_inte = intersect_curve_to_axis(clipped_curve, dir, n, pts);
	  
    free_egeom(clipped_curve);
  }

  return (total_inte);
}

/* This routine searches for a minimum, in a given finite interval, of a 
   continuous function of a single variable, using function and first
   derivative values.  NAG routine E04BBF is used.
*/

int min_bez_curve_to_axis_intrsct(double min, double max, double *point_inter,
				  int *num_point_inter)
{
  double a, b, eps, f, g, t, x;
  int ifail , maxcal;
  int x_min_comp, x_max_comp, a_b_comp;

  a_b_comp = lf2comp(max, min, EDITOR_TOL);
  if (a_b_comp <= 0) {
    a =  min - EDITOR_TOL/2.;
    b =  max + EDITOR_TOL/2.;
  }
  else {
    a = min;
    b = max;
  }

  maxcal =  30;
  ifail = 1;
  t = 0.0;
  eps = 0.0; 

     /*eps and t are set to zero so that e04bbf will reset them
       to default values */

  e04bbf_(functcv ,&eps, &t, &a, &b, &maxcal, &x, &f, &g, &ifail);   

  x_min_comp = lf2comp(x, min, EDITOR_ZERO_MIN);
  x_max_comp = lf2comp(x, max, EDITOR_ZERO_MIN);
     
  if ( f <= EDITOR_ZERO_MIN_SOL && x_min_comp >= 0 && x_max_comp <= 0) {  
    point_inter[*num_point_inter] = x;   /* u value */
    (*num_point_inter)++;
      
    return(EDITOR_SURF_INTRSCT);
  }
  else
    return(EDITOR_SURF_NT_INTRSCT);
}

/*-----------------------------------------------------------------------
  Calculates the values of the objective function F(x) and its first
  derivatives dF/dx at any point x
 ----------------------------------------------------------------------*/

void functcv(double *xc, double *fc, double *gc)
{
  *fc = functcv_sqr(*xc);
  *gc = partl_u_functcv_sqr(*xc);
}

double functcv_sqr(double xc)
{
  double f, value;
  vector *vect;
     
  vect = eval_curve_bounded(global_egeom, xc, 0);
  f = vect->z / vect->w;
  vectfree(vect);

  value = f * f;
  
  return value;
}

double partl_u_functcv_sqr(double xc)
{
  double f, fu, value;
  vector *vect;

  vect = eval_curve_bounded(global_egeom, xc, 0);
  f = vect->z / vect->w;
  vectfree(vect);
     
  vect = eval_curve_bounded(global_egeom, xc, 1);
  fu = vect->z / vect->w;
  vectfree(vect);
     
  value = 2. * f * fu;
     
  return value;
}

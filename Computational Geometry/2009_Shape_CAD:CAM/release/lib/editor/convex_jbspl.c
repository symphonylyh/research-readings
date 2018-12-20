/* Copyright (C) 1993 Massachusetts Institute of Technology, Cambridge, MA
   All rights reserved
*/

/* convex_jbspl.c */

/* ConvexHullSurf()
   ConvexHullTest()
*/

#include <stdlib.h>
#include "gen.h"
#include "bspl.h"
#include "editor.h"

/* convex hull test for two curves */

double ConvexHullTest(ParCurv *egeom, ParCurv *polyg)
{
  ParCurv *egm[2];
  double diff, dist, *knot, *node, **lim, epsilon;
  int index[3];
  short i, k, order, ncontpts, nknots;

  epsilon = 1.0e-15;
  order = MAX(egeom->order, polyg->order);
  ncontpts = MAX(egeom->ncontpts, polyg->ncontpts);
  nknots = order + ncontpts;
  knot = dbl_array1(nknots*order);

  egm[0] = egeomalloc1(order, order*(order+ncontpts));
  egm[1] = egeomalloc1(order, order*(order+ncontpts));
  
  copyegeom(egeom, egm[0]);
  copyegeom(polyg, egm[1]);

  /* raise both curves to same order */

  if (order > egm[0]->order) {
    lim = dbl_array2((order+ncontpts)*order, ncontpts*order);
    for (i=0; i < order - egeom->order; i++) {
      k = compute_lim(egm[0], lim, knot, epsilon);
      raise_byone(egm[0], lim, knot, k);
    }
    free_darray2(lim);
  }

  if (order > egm[1]->order) {
    lim =dbl_array2((order+ncontpts)*order, ncontpts*order);
    for(i=0; i < order - polyg->order; i++) {
      k = compute_lim(egm[1], lim, knot, epsilon);
      raise_byone(egm[1], lim, knot, k);
    }
    free_darray2(lim);
  }

  /* merge knot vectors */

  nknots = merge_knotv(egm, 0, 2, knot, index, epsilon);
  addpoints(egm, 0, 2, knot, nknots);
  free_darray1(knot);

/* check convex hull distances */

  dist = -1.0e+30;
  for (i=0; i<egm[0]->ncontpts; i++) {
    diff = distance(egm[0]->contpts[i], egm[1]->contpts[i]);

    if (diff > dist) dist = diff; 
  }
  free_egeom(egm[0]);
  free_egeom(egm[1]);

  return dist;
}

/* convex hull test for two surfaces */

/* fgeom and ogeom are the two surfaces to be compared */
/* nmax is the maximum number knots allowed when raising the surfaces to
   the same degree - a good value is 500 */
/* nknots is the number of knots to add to the surfaces before the
   convex hull check - to add no knots set equal to 0 */

double ConvexHullSurf(ParSurf *fgeom, ParSurf *ogeom, int nmax, int nknots)
{
  ParSurf *fg[2];
  ParCurv *eg[2];
  double dis, *kn[2], maxErr = 0.0;
  int *tn;
  short i, j, is, el;

  kn[0] = dbl_array1(nmax);
  kn[1] = dbl_array1(nmax);
  tn = int_array1(nmax);

  for(j=0; j<2; j++) {   /* make copy of original surfaces */
    fg[j] = fgeomalloc2(fgeom->uorder, fgeom->vorder, nmax, nmax);
    eg[j] = egeomalloc2(fgeom->uorder, nmax);
    for(i=0; i<ogeom->ucontpts; i++)
      eg[j]->contpts[i] = vectalloc();
  }

  is = fgeom->uorder-ogeom->uorder;
  el = fgeom->vorder-ogeom->vorder;
  if (is > -1 && el > -1) {
    copyfgeom(ogeom, fg[0]);
    copyfgeom(fgeom, fg[1]);
  }
  else {
    copyfgeom(ogeom, fg[1]);
    copyfgeom(fgeom, fg[0]);
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
  
  /*** we used to use fg[1]->ucontpts and vcontpts as the loop counters ***/
  /*** we should use fg[0] instead, in case extra control points have ***/
  /*** been added to the new bezier patch (S. L. Abrams, 10/31/96) ***/
  for (i=0; i<fg[0]->ucontpts; i++)      /* convex hull test */
    for (j=0; j<fg[0]->vcontpts; j++) {
      dis = distance(fg[1]->contpts[i][j], fg[0]->contpts[i][j]);
      if (dis > maxErr) maxErr = dis;
    }

  free_fgeom(fg[0]);
  free_fgeom(fg[1]);

  return maxErr;
}

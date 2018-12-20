/* Copyright (C) 1992 Massachusetts Institute of Technology, Cambridge, MA
   All rights reserved
*/

/* max_curvature.c */

/* max_curve_curvature()
*/

#include <malloc.h>
#include "gen.h"
#include "bspl.h"
#include "editor.h"

double max_curve_curvature(ParCurv *egeom, short nsegs)
{
  vector *(*evlrbsp)(ParCurv *, double, int);
  vector *v, *v0, *v1, *v2;
  double a,b,du,minimum_u,maximum_u,k,maxk,u;
  short i;

  minimum_u = egeom->knots[0];
  if (egeom->type == PCurvePer) {
    evlrbsp = rbspeval_per;
    maximum_u = egeom->knots[egeom->ncontpts];
  }
  else {
    evlrbsp = rbspeval;
    maximum_u = egeom->knots[egeom->order+egeom->ncontpts-1];
  }

  u = 0.0;
  du = (maximum_u-minimum_u)/(double)nsegs;
  maxk = EDITOR_EMIN;
  for (i=1; i<nsegs; i++) {
    u = (double)i*du;  
    v0 = (*evlrbsp)(egeom, u, 0);
    v0->x = v0->x/v0->w;
    v0->y = v0->y/v0->w;
    v0->z = v0->z/v0->w;
    v0->w = 1.0;
    v1 = (*evlrbsp)(egeom, u, 1);
    v2 = (*evlrbsp)(egeom, u, 2);
    
    a = mag(v1);
    v = cross(v1, v2);
    b = mag(v);
    k = b/(a*a*a);
    if (k > maxk) maxk = k;

    vectfree(v);
    vectfree(v0);
    vectfree(v1);
    vectfree(v2);
  }
  return (maxk);
}

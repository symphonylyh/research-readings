/* Copyright (C) 1992 Massachusetts Institute of Technology, Cambridge, MA
   All rights reserved
*/

/* evaluate.c */

/* EvaluateCurv()
   EvaluateSurf()
*/

#include <math.h>
#include <malloc.h>
#include "gen.h"
#include "bspl.h"
#include "editor.h"

vector *EvaluateCurv(ParCurv *egeom, float t, short is_open, vector **tang,
		     vector **norm, double *kcurv)
{
  vector *r, *r2, *u;
  double a, b, k;

  if (is_open) {
    r = evalbsp(egeom, t);
    (*tang) = rbspeval(egeom, t, 1);
    r2 = rbspeval(egeom, t, 2);
  }
  else {
    r = evalbsp_per(egeom, t);
    (*tang) = rbspeval_per(egeom, t, 1);
    r2 = rbspeval_per(egeom, t, 2);
  }
  a = mag((*tang));
  u = cross((*tang), r2);
  b = mag(u);
  (*kcurv) = b/(a*a*a);            /* curvature */
  
  a = dot(r2, (*tang));
  b = (mag((*tang))*mag((*tang)));
  scale_vect1(a/b, (*tang), u);
  (*norm) = sub_vect(r2, u);

  vectfree(r2);
  vectfree(u);

  unitvector1((*tang), (*tang));   /* tangent */
  unitvector1((*norm), (*norm));   /* normal */

  return (r);
}

vector *EvaluateSurf(ParSurf *fgeom, float u, float v, vector **tangu,
		     vector **tangv, vector **norm, double *kmax,
		     double *kmin, double *normu, double *normv)
{
  struct vector *r[6], *ruv, *ruu, *rvv;
  double e, f, g, l, m, n, H, K;
  double2 tempmod, Gdet, Ddet;
 
  evalrsurf(fgeom, u, v, 5, r);    /* pts, 1st, and 2nd derivatives */

  if (mag(r[1]) < EDITOR_ZERO_10)  /* normal vector */
    (*norm) = cross(r[3], r[2]);	
  else if(mag(r[2]) < EDITOR_ZERO_10)
    (*norm) = cross(r[1], r[3]);
  else
    (*norm) = cross(r[1], r[2]);
  unitvector1((*norm), (*norm));

  e = dot(r[1], r[1]);   /* fundamental coefficients */
  f = dot(r[1], r[2]);
  g = dot(r[2], r[2]);

  l = dot((*norm), r[4]);
  m = dot((*norm), r[3]);
  n = dot((*norm), r[5]);

  *normu = -l/e;         /* normal u */
  *normv = -n/g;         /* normal v */

  Gdet = e*g - f*f;
  Ddet = l*n - m*m;
  
  H = ((f*m + f*m) - (e*n + g*l))/(2.0*Gdet);   /* mean */
  K = Ddet/Gdet;                                /* gauss */

  tempmod = fabs(H*H - K);
  if (tempmod < EDITOR_ZERO_6) {
    *kmax = H;                   /* max principal */
    *kmin = H;                   /* min principal */
  }
  else {
    *kmax = H + sqrt(tempmod);
    *kmin = H - sqrt(tempmod);
  }

  unitvector1(r[1], r[1]);
  unitvector1(r[2], r[2]);
  (*tangu) = r[1];               /* tangent in U */
  (*tangv) = r[2];               /* tangent in V */

  vectfree(r[3]);
  vectfree(r[4]);
  vectfree(r[5]);

  return (r[0]);
}

vector *EvaluateCos(ParCurv *egeom, float t, short is_open, ParSurf *fgeom,
		    vector **tangu, vector **tangv, vector **nvect,
		    double *kmax, double *kmin, double *normu, double *normv)
{
  vector *r;
  double u, v;
  short i;

  if (is_open)
    r = evalbsp(egeom, t);
  else
    r = evalbsp_per(egeom, t);

  u = r->x/r->w;
  v = r->y/r->w;
  vectfree(r);

  r = EvaluateSurf(fgeom, u, v, tangu, tangv, nvect, kmax, kmin, normu, normv);

  return (r);
}

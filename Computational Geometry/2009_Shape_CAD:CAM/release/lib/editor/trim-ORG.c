/* Copyright (C) 1993 Massachusetts Institute of Technology, Cambridge, MA
   All rights reserved

   Last edited: April 11, 2003 for v.11.1

   v.11.1: (1) free_trim(...) changed.
*/

/* trim.c */

/* Circumcenter()
   Circumradius()
   CopyTrimSurf()
   CurveConvexHull()
   Det()
   free_trim()
   GetIntersect()
   InsertTriangle()
   InvertTrimSurf()
   IsoLoopIntersect()
   MakeTriangle()
   NoIntersects()
   Overlap()
   ScramblePoints()
   SegmentMinMax()
   ToLeft()
   TrimCompare()
   TrimCompareXY()
   TrimPoints()
*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "gen.h"
#include "bspl.h"
#include "editor.h"
char *tempnam(const char *, const char *);

void IsoLoopIntersect(int flag, double w, int nLoops, TrimLoop *loops,
		      int nsegspan, double *in, int *ni, int max)
{
  vector *r;
  double dt, **pts, t, t0, t1, umin, umax, vmin, vmax;
  int i, j, k, m, n, nsegs, nspans, order;

  for (i=0; i<nLoops; i++) {
    for (j=0; j<loops[i].nCurves; j++) {
      CurveConvexHull(loops[i].egeoms[j], &umin, &umax, &vmin, &vmax);

      /* umin,umax,vmin,vmax is bounding box of this curve */

      if ((flag == EDITOR_TRIM_U && (w < umin || w > umax)) ||
	  (flag == EDITOR_TRIM_V && (w < vmin || w > vmax)))
	continue;   /* iso line doesn't intersect this curve, go to next */

      order = loops[i].egeoms[j]->order;
      nspans = loops[i].egeoms[j]->ncontpts - order + 1;
      nsegs = nspans*nsegspan + 1;

      pts = dbl_array2(nsegs, 2);   /* linear approximation of curve */

      for (m=n=0; m<nspans; m++) {
	t0 = loops[i].egeoms[j]->knots[order + m - 1];
	t1 = loops[i].egeoms[j]->knots[order + m];
	dt = t1 - t0;
	for (k = (m ? 1 : 0); k<=nsegspan; k++, n++) {
	  t = t0 + dt*(double)k/(double)nsegspan;
	  r = evalbsp(loops[i].egeoms[j], t);
	  pts[n][0] = r->x/r->w;
	  pts[n][1] = r->y/r->w;
	  vectfree(r);
	}
      }
      for (k=0; k<nsegs-1 && *ni < max; k++) {
	SegmentMinMax(pts[k], pts[k+1], &umin, &umax, &vmin, &vmax);
      /* umin,umax,vmin,vmax is bounding box of this linear segment of curve */
	if (flag == EDITOR_TRIM_U && umin <= w && w < umax)
	  in[(*ni)++] = pts[k][1] +
	                 (w - pts[k][0])*(pts[k+1][1] - pts[k][1])/
	                                 (pts[k+1][0] - pts[k][0]);
	else if (flag == EDITOR_TRIM_V && vmin <= w && w < vmax)
	  in[(*ni)++] = pts[k][0] +
	                 (w - pts[k][1])*(pts[k+1][0] - pts[k][0])/
	                                 (pts[k+1][1] - pts[k][1]);
      }
      free_darray2(pts);
    }
  }                                              /* sort intersections in */
  qsort(in, *ni, sizeof(double), TrimCompare);   /* ascending order */
}

void CurveConvexHull(ParCurv *egeom, double *umin, double *umax,
		     double *vmin, double *vmax)
{
  double u, v;
  int i;

  *umin = *vmin = EDITOR_EMAX;
  *umax = *vmax = EDITOR_EMIN;
  for (i=0; i<egeom->ncontpts; i++) {
    u = egeom->contpts[i]->x/egeom->contpts[i]->w;
    v = egeom->contpts[i]->y/egeom->contpts[i]->w;
    if (u < *umin) *umin = u;
    if (u > *umax) *umax = u;
    if (v < *vmin) *vmin = v;
    if (v > *vmax) *vmax = v;
  }
}

void SegmentMinMax(double *p0, double *p1, double *umin, double *umax,
		   double *vmin, double *vmax)
{
  if (p0[0] < p1[0]) {
    *umin = p0[0];
    *umax = p1[0];
  }
  else {
    *umin = p1[0];
    *umax = p0[0];
  }

  if (p0[1] < p1[1]) {
    *vmin = p0[1];
    *vmax = p1[1];
  }
  else {
    *vmin = p1[1];
    *vmax = p0[1];
  }
}

int TrimCompare(const void *v1, const void *v2)
{
  if (*(double *)v1 < *(double *)v2)
    return (-1);
  else if (*(double *)v1 > *(double *)v2)
    return (1);

  return (0);
}

TrimSurf *CopyTrimSurf(TrimSurf *from, TrimSurf *to)
{
  int i, j;

  if (!to)
    to = (TrimSurf *)gen_array1(1, sizeof(TrimSurf));

  to->fgeom = copyfgeom(from->fgeom, NULL);
  to->bLoop = from->bLoop;
  to->nLoops = from->nLoops;
  to->loops = (TrimLoop *)gen_array1(to->nLoops, sizeof(TrimLoop));
  for (i=0; i<to->nLoops; i++) {
    to->loops[i].nCurves = from->loops[i].nCurves;
    to->loops[i].egeoms = (ParCurv **)gen_array1(to->loops[i].nCurves,
						 sizeof(ParCurv *));
    for (j=0; j<to->loops[i].nCurves; j++)
      to->loops[i].egeoms[j] = copyegeom(from->loops[i].egeoms[j], NULL);
  }

  return (to);
}

int TrimCompareXY(const void *v1, const void *v2)
{
  double *a1, *a2;

  a1 = (double *)v1;
  a2 = (double *)v2;

  if (a1[0] < a2[0])
    return (-1);
  else if (a1[0] > a2[0])
    return (1);
  else {
    if (a1[1] < a2[1])
      return (-1);
    else if (a1[1] > a2[1])
      return (1);
  }
  return (0);
}

void TrimPoints(TrimSurf *trim, int nsegu, int nsegv, int nsegspan,
		char *tmpfil, int *lcc, int *nNodes, int boundary)
{
  FILE *fp;
  vector *r;
  double dt, du, dv, t, told, t0, t1, u, v, u0, v0, u1, v1;
  double in[EDITOR_TRIM_MAX], uinc, vinc;
  int flag, i, ic = 0, j, k, m, ni, np, nspans, order;
  int ue, us, ve, vs;

  u0 = trim->fgeom->uknots[trim->fgeom->uorder-1];
  v0 = trim->fgeom->vknots[trim->fgeom->vorder-1];
  u1 = trim->fgeom->uknots[trim->fgeom->ucontpts];
  v1 = trim->fgeom->vknots[trim->fgeom->vcontpts];
  du = u1 - u0;
  dv = v1 - v0;
  uinc = du/(double)(nsegu-1);
  vinc = dv/(double)(nsegv-1);
  
  np = 0;
  
  if (boundary) {
    us = vs = 1;
    ue = nsegu - 1;
    ve = nsegv - 1;
  }
  else {
    us = vs = 0;
    ue = nsegu;
    ve = nsegv;
  }

  fp = fopen(tmpfil, "w");
  for (i=us; i<ue; i++) {
    u = u0 + i*uinc;
    
    if (trim->bLoop) {
      in[0] = v0;
      in[1] = v1;
      ni = 2;
    }
    else
      ni = 0;

    /* "in" holds the "ni" break points, the parametric values where */
    /* the trimming loops intersect the isoparameter curve */
    
    IsoLoopIntersect(EDITOR_TRIM_U, u, trim->nLoops, trim->loops, nsegspan,
		     in, &ni, EDITOR_TRIM_MAX);

    v = v0;                     /* add points in interior of trimmed surface */
    for (j=0; j<ni; j+= 2) {        /* number of break points is always even */
      while (j < ni-1 && v < in[j+1]) {
	v += vinc;

	if (in[j] < v && v < 0.99*in[j+1]) {   /* segment ends between */
	  fprintf(fp, "%9.7f %9.7f\n", u, v);        /* break points */
	  np++;
	}
      }
    }
  }
  *nNodes  = np;

  if (trim->bLoop) {                 /* add points on boundary of surface */
    lcc[ic++] = *nNodes + 1;         /* lcc is offset-1, not offset-0 */
    *nNodes += 2*(nsegu+nsegv-2);
    
    u = u0;                          /* note that the points are defined */
    for (i=1; i<nsegv; i++) {        /* in clockwise orientation */
      v = v0 + i*vinc;
      fprintf(fp, "%9.7f %9.7f\n", u, v);
    }
    for (i=1; i<nsegu; i++) {
      u = u0 + i*uinc;
      fprintf(fp, "%9.7f %9.7f\n", u, v);
    }
    for (i=1; i<nsegv; i++) {
      v = v1 - i*vinc;
      fprintf(fp, "%9.7f %9.7f\n", u, v);
    }
    for (i=1; i<nsegu; i++) {
      u = u1 - i*uinc;
      fprintf(fp, "%9.7f %9.7f\n", u, v);
    }
  }

  for (i=0; i<trim->nLoops; i++) {   /* add points on trimming loops */
    lcc[ic++] = *nNodes + 1;         /* lcc is offset-1, not offset-0 */

    /* the TrimSuf file format has the outmost trimming curve
     * counter-clockwise
     * the new TRIPACK triangulation routines want the outmost curve
     * to be clockwise
     * so we have to define the points on the curves in the opposite
     * direction
     * start with the last span of the last curve and work towards the
     * first span of the first curve
     */

    told = -999.999;         /* initialize last t value */
/***************************************************************************/
/***** the original lines (commented out) assummed the other orientation ***/
/***** they start with the first span of the first curve, and so on ********/
/** for (j=0; j<trim->loops[i].nCurves; j++) { **/
/************************************************/
    for (j=trim->loops[i].nCurves-1; j >= 0; j--) {
      order = trim->loops[i].egeoms[j]->order;
      nspans = trim->loops[i].egeoms[j]->ncontpts - order + 1;
      
/**   for (m=0; m<nspans; m++) { **/
      for (m=nspans-1; m >= 0; m--) {
	t0 = trim->loops[i].egeoms[j]->knots[order + m - 1];
	t1 = trim->loops[i].egeoms[j]->knots[order + m];
	dt = t1 - t0;
/**	for (k=1; k<=nsegspan; k++) { **/
	for (k=nsegspan; k > 0; k--) {
	  t = t0 + dt*(double)k/(double)nsegspan;
	  if (t != told) {   /* make sure t is different from last value */
	    r = evalbsp(trim->loops[i].egeoms[j], t);
	    u = r->x/r->w;
	    v = r->y/r->w;
	    vectfree(r);
	    fprintf(fp, "%9.7f %9.7f\n", u, v);

	    told = t;
	    (*nNodes)++;
	  }
	}
      }
    }
  }
  fclose(fp);
}

int MakeTriangle(double **pts, int *lambda, int *sigma, int *nNodes,
		   double uinc, double vinc, double sinc, int *actual,
		   AdjType *adjs, int *nTriang)
{
  double ab, cc[2], r, rk, umin, umax, vmin, vmax;
  int a = -1, b, c, i, j, k, *kappa, kc, kmax, nKappa, outside;
  int cross_lambda, cross_over = 0, new_lambda, try_again = 1;

  for (i=0; i < *nNodes; i++)
    if (lambda[i] == 1 || lambda[i] == 2) {
      a = i;
      break;
    }
  if (a < 0)
    return (0);

  while (try_again) {   /* sigma[a] points to the next node on boundary */
    b = sigma[a];       /* in this case, b */

    if (sigma[c = sigma[b]] == a) {   /* case 4 */
      InsertTriangle(pts, adjs, actual[a], actual[b], actual[c]);
      InsertTriangle(pts, adjs, actual[b], actual[c], actual[a]);
      InsertTriangle(pts, adjs, actual[c], actual[a], actual[b]);
      (*nTriang)++;

      lambda[a] = 0;    /* all three nodes, a, b, c, have been triangulated */
      lambda[b] = 0;
      lambda[c] = 0;
    
      return (1);
    }

    kappa = int_array1(*nNodes);

    nKappa = 0;
    for (j=0; nKappa == 0 && j<3; j++)
      for (k=0; k < *nNodes; k++)
	if (lambda[k] && k != a && k != b) {
	  SegmentMinMax(pts[a], pts[b], &umin, &umax, &vmin, &vmax);
	  if ((j == 0 &&
	       umin - uinc <= pts[k][0] && pts[k][0] <= umax + uinc &&
	       vmin - vinc <= pts[k][1] && pts[k][1] <= vmax + vinc) ||
	      (j == 1 &&
	       umin - sinc <= pts[k][0] && pts[k][0] <= umax + sinc &&
	       vmin - sinc <= pts[k][1] && pts[k][1] <= vmax + sinc) ||
	      j == 2)
	    if (ToLeft(pts[a], pts[b], pts[k]))
	      if (lambda[a] == 2 ||
		  NoIntersects(pts, lambda, sigma, *nNodes, a, b, k))
		kappa[nKappa++] = k;
	}
    
    ab = sqrt((pts[b][0] - pts[a][0])*(pts[b][0] - pts[a][0]) +
	      (pts[b][1] - pts[a][1])*(pts[b][1] - pts[a][1]));

    c = -1;
    kc = kmax = 0;
    for (k=0; k<nKappa; k++) {
      outside = 0;
      for (j=0; j<nKappa; j++)
	if (k != j) {
	  r = Circumradius(pts[a], pts[b], pts[kappa[j]], ab);
	  Circumcenter(pts[a], pts[b], pts[kappa[j]], r, cc);
	  rk = sqrt((pts[kappa[k]][0]-cc[0])*(pts[kappa[k]][0]-cc[0]) +
		    (pts[kappa[k]][1]-cc[1])*(pts[kappa[k]][1]-cc[1]));
	  if (rk > r + EDITOR_TRIM_EXTRA) {
	    outside = 1;
	    break;
	  }
	  else
	    kc++;
	}
      if (!outside) {
	c = kappa[k];
	new_lambda = 2;
	break;
      }
      if (kc > kmax) {
	kmax = kc;
	c = kappa[k];
	new_lambda = 1;
      }
    }
    free_iarray1(kappa);

    if (c > -1 && lambda[c] == 3 || sigma[b] == c || sigma[c] ==a) {
      InsertTriangle(pts, adjs, actual[a], actual[b], actual[c]);
      InsertTriangle(pts, adjs, actual[b], actual[c], actual[a]);
      InsertTriangle(pts, adjs, actual[c], actual[a], actual[b]);
      (*nTriang)++;
    }

    if (c > -1 && lambda[c] == 3) {       /* case 1 */
      sigma[a]  = c;          /* add AC to sigma */
      sigma[c]  = b;          /* add CB to sigma */
      lambda[a] = new_lambda;
      lambda[c] = new_lambda;
      try_again = 0;
    }
    else if (c > -1 && sigma[b] == c) {   /* case 2 */
      sigma[a]  = c;          /* add AC to sigma */
      lambda[a] = new_lambda;
      lambda[b] = 0;          /* b has been triangulated */
      try_again = 0;
    }
    else if (c > -1 && sigma[c] == a) {   /* case 3 */
      sigma[c]  = b;          /* add CB to sigma */
      lambda[c] = new_lambda;
      lambda[a] = 0;          /* a has been triangulated */
      try_again = 0;
    }
    else {                    /* go to next AB */
      if (c > -1 && sigma[b] != c && sigma[c] != a) {
	cross_over = c;
	cross_lambda = new_lambda;
      }

      j = a;
      a = -1;
      for (i=j+1; i < *nNodes; i++)
	if (lambda[i] == 1 || lambda[i] == 2) {
	  a = i;
	  break;
	}
      if (a < 0) {
	if (cross_over) {   /* case 5 */
	  a = j;
	  c = cross_over;
	  new_lambda = cross_lambda;

	  InsertTriangle(pts, adjs, actual[a], actual[b], actual[c]);
	  InsertTriangle(pts, adjs, actual[b], actual[c], actual[a]);
	  InsertTriangle(pts, adjs, actual[c], actual[a], actual[b]);
	  (*nTriang)++;

	  sigma[a] = c;                  /* add AC to sigma, deleting AB */
	  lambda[a] = new_lambda;
	  pts[*nNodes][0] = pts[c][0];   /* duplicate C by C' */
	  pts[*nNodes][1] = pts[c][1];
	  lambda[*nNodes] = new_lambda;
	  actual[*nNodes] = c;
	  for (j=0; j < *nNodes; j++)
	    if (lambda[j] == 1 || lambda[j] == 2)
	      if (sigma[j] ==c) {
		sigma[j] = *nNodes;
		break;
	      }
	  sigma[*nNodes] = b;            /* add C'B to sigma */
	  (*nNodes)++;
	  try_again = 0;
	}
	else
	  return (0);
      }
      else
	try_again = 1;
    }
  }
  return (1);
}

void InsertTriangle(double **pts, AdjType *adjs, int a, int b, int c)
{
  AdjElem *new, *this;
  int quit = 0;

  this = adjs[a].first;
  while (this) {
    if (b == this->adj) {   /* insert c after b in list */
      if (!this->next || c != this->next->adj) {
	new = (AdjElem *)gen_array1(1, sizeof(AdjElem));
	new->adj = c;
	new->prev = this;
	if (new->next = this->next)
	  new->next->prev = new;
	else
	  adjs[a].last = new;
	this->next = new;
	adjs[a].nAdjs += 1;
      }
      quit = 1;
      break;
    }
    this = this->next;
  }

  if (!quit) {
    this = adjs[a].first;
    while (this) {   /* insert b before c in list */
      if (c == this->adj) {
	if (!this->prev || b != this->prev->adj) {
	  new = (AdjElem *)gen_array1(1, sizeof(AdjElem));
	  new->adj = b;
	  new->next = this;
	  if (new->prev = this->prev)
	    new->prev->next = new;
	  else
	    adjs[a].first = new;
	  this->prev = new;
	  adjs[a].nAdjs += 1;
	}
	quit = 1;
	break;
      }
      this = this->next;
    }
  }

  if (!quit) {
    this = adjs[a].last;
    new = (AdjElem *)gen_array1(1, sizeof(AdjElem));
    new->adj = b;
    if (this) {
      new->prev = this;
      this->next = new;
    }
    else
      adjs[a].first = new;
    this = new;
    new = (AdjElem *)gen_array1(1, sizeof(AdjElem));
    new->adj = c;
    new->prev = this;
    this->next = new;
    adjs[a].last = new;
    adjs[a].nAdjs += 2;
  }

  if (!adjs[a].bound)
    if (adjs[a].first->adj == adjs[a].last->adj) {
      this = adjs[a].last;
      adjs[a].last = this->prev;
      adjs[a].last->next = NULL;
      free(this);
      adjs[a].nAdjs -= 1;
    }
}

int ToLeft(double *a, double *b, double *c)
{
  double ca[2], cb[2];

  ca[0] = a[0] - c[0];
  ca[1] = a[1] - c[1];

  cb[0] = b[0] - c[0];
  cb[1] = b[1] - c[1];

  if ((ca[0]*cb[1] - cb[0]*ca[1]) > EDITOR_TRIM_EXTRA)
    return (1);
  else
    return (0);
}

double Circumradius(double *a, double *b, double *c, double ab)
{
  double ca[2], cb[2], sinC, X;

  ca[0] = a[0] - c[0];
  ca[1] = a[1] - c[1];

  cb[0] = b[0] - c[0];
  cb[1] = b[1] - c[1];

  X = ca[0]*cb[1] - cb[0]*ca[1];
  sinC = X/((sqrt(ca[0]*ca[0] + ca[1]*ca[1]) *
	     sqrt(cb[0]*cb[0] + cb[1]*cb[1])));

  return (ab/(2.0*sinC));
}

void Circumcenter(double *a, double *b, double *c, double r, double *cc)
{
  double x1[2], x2[2], y1[2], y2[2], s, t;

  x1[0] = (c[0] + a[0]) / 2.0;
  x1[1] = (c[1] + a[1]) / 2.0;
  x2[0] = x1[0] + c[1] - x1[1];
  x2[1] = x1[1] - c[0] + x1[0];

  y1[0] = (c[0] + b[0]) / 2.0;
  y1[1] = (c[1] + b[1]) / 2.0;
  y2[0] = y1[0] - c[1] + y1[1];
  y2[1] = y1[1] + c[0] - y1[0];
  
  GetIntersect(x1, x2, y1, y2, &s, &t);
  cc[0] = s*x2[0] + (1-s)*x1[0];
  cc[1] = s*x2[1] + (1-s)*x1[1];
}

int NoIntersects(double **pts, int *lambda, int *sigma, int nNodes,
		   int a, int b, int c)
{
  double s1, s2, t1, t2;
  int i, j, no = 1, p1, p2;

  for (i=0; i<nNodes; i++) {
    j = sigma[i];
    if (!lambda[i] && (i != a || j != b)) {
      if ((i == c && j == a) || (i == b && j == c))
	  break;    /* if IJ is CA or BC, then can't intersect */

      if (!Overlap(pts[c], pts[a], pts[i], pts[j]) &&
	  !Overlap(pts[b], pts[c], pts[i], pts[j]))
	continue;   /* if IJ doesn't overlap CA or BC it can't intersect */

      p1 = GetIntersect(pts[c], pts[a], pts[i], pts[j], &s1, &t1);
      p2 = GetIntersect(pts[b], pts[c], pts[i], pts[j], &s2, &t2);
      if (p1 && p2)
	continue;   /* if IJ is parallel to CA or BC it can't intersect */

      if ((s1 <= 0.0 + EDITOR_TRIM_EXTRA || s1 >= 1.0 - EDITOR_TRIM_EXTRA) &&
	  (s2 <= 0.0 + EDITOR_TRIM_EXTRA || s2 >= 1.0 - EDITOR_TRIM_EXTRA))
	continue;   /* if intersects at end points, ok */

      no = 0;      /* intersection */
      break;
    }
  }

  return (no);
}

int Overlap(double *p0, double *p1, double *q0, double *q1)
{
  double pxmin, pxmax, pymin, pymax, qxmin, qxmax, qymin, qymax;

  if (p0[0] < p1[0]) {
    pxmin = p0[0];
    pxmax = p1[0];
  }
  else {
    pxmin = p1[0];
    pxmax = p0[0];
  }

  if (p0[1] < p1[1]) {
    pymin = p0[1];
    pymax = p1[1];
  }
  else {
    pymin = p1[1];
    pymax = p0[1];
  }

  if (q0[0] < q1[0]) {
    qxmin = q0[0];
    qxmax = q1[0];
  }
  else {
    qxmin = q1[0];
    qxmax = q0[0];
  }

  if (q0[1] < q1[1]) {
    qymin = q0[1];
    qymax = q1[1];
  }
  else {
    qymin = q1[1];
    qymax = q0[1];
  }

  if (pxmax <= qxmin || pxmin >= qxmax || pymax <= qymin || pymin >= qymax)
    return (0);

  return (1);
}

int GetIntersect(double *p0, double *p1, double *q0, double *q1, double *s,
		 double *t)
{
  double a, b, c, alpha, beta, denom;
  int parallel = 0;

  a = (p1[0]-p0[0])*(q0[0]-q1[0]) + (p1[1]-p0[1])*(q0[1]-q1[1]);
  b = (q0[0]-q1[0])*(q0[0]-q1[0]) + (q0[1]-q1[1])*(q0[1]-q1[1]);
  c = (p1[0]-p0[0])*(p1[0]-p0[0]) + (p1[1]-p0[1])*(p1[1]-p0[1]);

  alpha = -1.0*((p0[0]-q0[0])*(q0[0]-q1[0]) + (p0[1]-q0[1])*(q0[1]-q1[1]));
  beta  = -1.0*((p0[0]-q0[0])*(p1[0]-p0[0]) + (p0[1]-q0[1])*(p1[1]-p0[1]));
  
  denom = Det(a, b, c, a);
  
  if(denom <= 0.000001) {
    *s = Det(alpha, b, beta, a) / denom;
    *t = Det(a, alpha, c, beta) / denom;
  }
  else
    parallel = 1;

  return (parallel);
}

double Det(double a, double b, double c, double d)
{
  return (a*d - c*b);
}

TrimSurf *InvertTrimSurf(TrimSurf *from, TrimSurf *to)
{
  int i, j, k,m;

  if (!to)
    to = (TrimSurf *)gen_array1(1, sizeof(TrimSurf));

  /* to invert a trimmed surface, invert the boundary surface flag
   * and reverse the orientation of the all of the trimming curves
   */

  to->fgeom = copyfgeom(from->fgeom, NULL);   /* copy untrimmed surface */
  to->bLoop = (from->bLoop ? 0 : 1);          /* invert boundary flag */
  to->nLoops = from->nLoops;
  to->loops = (TrimLoop *)gen_array1(to->nLoops, sizeof(TrimLoop));
 
  /* invert trimming curves */
  /* copy the trimming loops in the reverse order */

  for (i=0,j=to->nLoops-1; j>=0; i++,j--) {
    to->loops[i].nCurves = from->loops[j].nCurves;
    to->loops[i].egeoms = (ParCurv **)gen_array1(to->loops[i].nCurves,
						 sizeof(ParCurv *));

    /* copy each curve of each trimming loop in the reverse order */
    /* remember to reverse the parametrization of the curves */

    for (k=0,m=to->loops[i].nCurves-1; m>=0; k++,m--) {
      to->loops[i].egeoms[k] = copyegeom(from->loops[j].egeoms[m], NULL);
      RevCurvParam(to->loops[i].egeoms[k]);   /* reverse the parametrization */
    }
  }

  return (to);
}

/* Scramble the order of the points sampled in the interior of a trimmed
 * surface, so that the first three points are not collinear.  This is
 * necessary as a precondition for the TRIPACK constrained triangulation.
 */

int ScramblePoints(int nNodes, float *x, float *y, int ncc, int *lcc)
{
  FILE *fp;
  float temp;
  int i, j, k, n1, n2, ok = 0;
  char *tmpfil;

  if (lcc[0] > 3) {   /* at least 3 points in the interior */
    n1 = 0;
    n2 = lcc[0] - 1;
    if (!IsCollinear2(x[n1], y[n1], x[n1+1], y[n1+1], x[n1+2], y[n1+2]))
      ok = 1;
    else
      for (i=n1+3; i<n2; i++)
	if (IsCollinear2(x[n1], y[n1], x[n1+1], y[n1+1], x[i], y[i])) {
	  temp = x[n1+2];
	  x[n1+2] = x[i];
	  x[i] = temp;
	  ok = 1;
	}
  }
  if (!ok)            /* look at the trimming loops */
    for (i=0; i<ncc && !ok; i++) {
      n1 = lcc[i] - 1;
      n2 = ((i < ncc-1) ? lcc[i+1] - 3 : nNodes - 2);
      for (j=n1; j<n2 && !ok; j++) {
	if (!IsCollinear2(x[j], y[j], x[j+1], y[j+1], x[j+2], y[j+2])) {
	  tmpfil = tempnam("/usr/tmp", "trim");
	  fp = fopen(tmpfil, "w");
	  for (k=j; k<n2; k++)
	    fprintf(fp, "%9.7f %9.7f\n", x[k], y[k]);
	  for (k=n1; k<j; k++)
	    fprintf(fp, "%9.7f %f9.7\n", x[k], y[k]);
	  fclose(fp);

	  fp = fopen(tmpfil, "r");
	  for (k=n1; k<n2; k++)
	    fscanf(fp, "%f %f", &x[k], &y[k]);
	  fclose(fp);
	  unlink(tmpfil);

	  ok = 1;
	}
      }
    }

  return ok;
}


/***********************************************
  "free_trim(...)" is modified in v.11.1
**********************************************/
void free_trim(TrimSurf *trim)
{
  int j, k;

  if (trim->fgeom)
    free_fgeom(trim->fgeom);
  if (trim->nLoops) {
    for (j=0; j<trim->nLoops; j++) {
      for (k=0; k<trim->loops[j].nCurves; k++) {
	free_egeom(trim->loops[j].egeoms[k]);
      }
      free(trim->loops[j].egeoms);
    }
    free(trim->loops);

    /********** 11.1 begin **********/
    if(trim->c_loops) { /* C-curve available (explicitly) */
      PrintLog("[Free TrimSurf] \"C-curves are provided explicitly.\"\n");
      for (j=0; j<trim->nLoops; j++) {
	for (k=0; k<trim->c_loops[j].nCurves; k++) {
	  free_egeom(trim->c_loops[j].egeoms[k]);
	}
	free(trim->c_loops[j].egeoms);
      }
      free(trim->c_loops);
    }
    else {
      PrintLog("[Free TrimSurf] \"C-curves are not provided explicitly.\"\n");
    }
    /********** 11.1 end ************/
  }
}

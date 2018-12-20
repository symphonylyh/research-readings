/* Copyright (C) 1993 Massachusetts Institute of Technology, Cambridge, MA
   All rights reserved
*/

/* trim_contour.c - Report regions out of tolerance */

/* ConvertToCos()
   CrossEdge()
   CrossTriangle()
   FindThirdPoint()
   LinearCosNormal()
   NextAdjB()
   TraceContour()
   TraceToBoundary()
   TraceToNext()
   TraceToStart()
   TriFunction()
*/

/* Given a faceted representation of a surface (pts, iAdj, iEnd), find */
/* the intersection curves with a plane of constant z = ct.  The curves */
/* are formed into closed loops.  For each loop, the orientation (clock- */
/* wise or counter-clockwise), a point inside the loop, and the Z value  */
/* at that point (above or below ct) are reported.  The loops and the */
/* NURBS design surface are used to generate trimmed patches. */

#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include "gen.h"
#include "bspl.h"
#include "editor.h"

int TraceContour(FILE *fp, double **pts, int nPts, int *iEnd, int *iAdj,
		 int nTriang, int *count, double ct, int *nBoundary)
{
  double cr[2], u, v;
  int a, b, c, a2, b2, c2;
  int end, i, j, k, k2, nCurves = 0, next, start, *visited;
  int inLoop = 0, quit = -1, skip = 0, stop;

  visited = int_array1(nTriang);

  for (a=0; a<nPts; a++)                       /* find point on boundary */
    if (iAdj[iEnd[a]] == -1)
      break;

  start = a;
  while (1) {
    i = a ? iEnd[a-1]+1 : 0;
    b = iAdj[i];                               /* current edge is AB */
    c = iAdj[i+1];
    k = TriFunction(nPts, iAdj, iEnd, count, a, b, c);
    if (!visited[k-1]) {
      if (CrossEdge(pts, a, b, ct, cr)) {            /* ct intersects AB */
	if (!inLoop) {
	  if (skip)
	    skip = 0;
	  else {                                       /* start new loop */
	    nCurves++;
	    fprintf(fp, "M %f %f %f\n", u = cr[0], v = cr[1], ct);

	    stop  = a;               /* stop loop at A = stop */
	    next = b;                /* start looking for next loop at B */

	    TraceToBoundary(fp, pts, nPts, iEnd, iAdj, count, visited, ct,
			    &a, &b, &c);
	    if (quit < 0)
	      quit = b;                    /* stop looking when A = quit */

	    fprintf(fp, "D %f %f %f\n", pts[c][0], pts[c][1], ct);
	    b = c;

	    inLoop = 1;

	    visited[k-1] = 1;
	  }
	}
	else {                                        /* already in loop */
	  fprintf(fp, "D %f %f %f\n", cr[0], cr[1], ct);
	  TraceToBoundary(fp, pts, nPts, iEnd, iAdj, count, visited, ct,
			  &a, &b, &c);

	  fprintf(fp, "D %f %f %f\n", pts[c][0], pts[c][1], ct);
	  b = c;

	  visited[k-1] = 1;
	}
      }
      else if (pts[a][2] == ct) {           /* contour passes through A */
	if (!inLoop) {
	  if (skip)
	    skip = 0;
	  else {                             /* start new loop */
	    nCurves++;
	    fprintf(fp, "M %f %f %f\n", u = pts[a][0], v = pts[a][1], ct);

	    stop = a;                /* stop loop at A = stop */
	    next = b;                /* start looking for next loop at B */

	    if (CrossEdge(pts, b, c, ct, cr))
	      TraceToBoundary(fp, pts, nPts, iEnd, iAdj, count, visited, ct,
			      &a, &b, &c);
	    else if (pts[b][2] == ct)
	      TraceToBoundary(fp, pts, nPts, iEnd, iAdj, count, visited, ct,
			      &c, &a, &b);
	    else if (pts[c][2] == ct)
	      TraceToBoundary(fp, pts, nPts, iEnd, iAdj, count, visited, ct,
			      &a, &b, &c);
	    else {
	      NextTriangle(pts, nPts, iEnd, iAdj, count, visited, ct,
			   a, &b, &c, 0);
	      TraceToBoundary(fp, pts, nPts, iEnd, iAdj, count, visited, ct,
			      &a, &b, &c);
	    }
	    if (quit < 0)
	      quit = b;                  /* stop looking when A = quit */

	    fprintf(fp, "D %f %f %f\n", pts[c][0], pts[c][1], ct);
	    b = c;

	    inLoop = 1;
	      
	    visited[k-1] = 1;
	  }
	}
	else {                             /* already in loop */
	  fprintf(fp, "D %f %f %f\n", pts[a][0], pts[a][1], ct);

	  if (CrossEdge(pts, b, c, ct, cr))
	    TraceToBoundary(fp, pts, nPts, iEnd, iAdj, count, visited, ct,
			    &a, &b, &c);
	  else if (pts[b][2] == ct)
	    TraceToBoundary(fp, pts, nPts, iEnd, iAdj, count, visited, ct,
			    &c, &a, &b);
	  else if (pts[c][2] == ct)
	    TraceToBoundary(fp, pts, nPts, iEnd, iAdj, count, visited, ct,
			    &a, &b, &c);
	  else {
	    NextTriangle(pts, nPts, iEnd, iAdj, count, visited, ct,
			 a, &b, &c, 0);
	    TraceToBoundary(fp, pts, nPts, iEnd, iAdj, count, visited, ct,
			    &a, &b, &c);
	  }

	  fprintf(fp, "D %f %f %f\n", pts[c][0], pts[c][1], ct);
	  b = c;

	  visited[k-1] = 1;
	}
      }
      else {
	if (inLoop)
	  fprintf(fp, "D %f %f %f\n", pts[b][0], pts[b][1], ct);
      }
    }
    else {
      if (inLoop)
	fprintf(fp, "D %f %f %f\n", pts[b][0], pts[b][1], ct);
    }

    if (b == stop) {                          /* end of current loop */
      fprintf(fp, "D %f %f %f\n", u, v, ct);
      inLoop = 0;
      skip = 1;
      a = next;
    }
    else                                 /* move to next boundary segment */
      a = b;

    if ((nCurves && a == quit) || (!nCurves && a == start))
      break;                     /* have looked at all boundary nodes     */
  }

  *nBoundary = nCurves;

  for (i=0; i<nPts; i++) {
    start = (i ? iEnd[i-1] + 1 : 0);
    end = iEnd[i] + 1;
    for (j=start; j+1 < end && iAdj[j+1] > -1; j++)
      if (iAdj[j] > i && iAdj[j+1] > i) {
	a = i;
	b = iAdj[j];
	c = iAdj[j+1];

	k = TriFunction(nPts, iAdj, iEnd, count, a, b, c);
	if (!visited[k-1]) {
	  if (CrossTriangle(pts, ct, &a, &b, &c)) {
	    if (CrossEdge(pts, a, b, ct, cr)) {
	      nCurves++;
	      fprintf(fp, "M %f %f %f\n", cr[0], cr[1], ct);

	      TraceToStart(fp, pts, nPts, iEnd, iAdj, count, visited, ct, k,
			   a, b, c);
	    }
	    else if (pts[a][2] == ct) {
	      nCurves++;
	      fprintf(fp, "M %f %f %f\n", pts[a][0], pts[a][1], ct);

	      if (CrossEdge(pts, b, c, ct, cr))
		TraceToStart(fp, pts, nPts, iEnd, iAdj, count, visited, ct, k,
			     a, b, c);
	      else if (pts[b][2] == ct)
		TraceToStart(fp, pts, nPts, iEnd, iAdj, count, visited, ct, k,
			     c, a, b);
	      else if (pts[c][2] == ct)
		TraceToStart(fp, pts, nPts, iEnd, iAdj, count, visited, ct, k,
			     c, a, b);
	      else {
		b2 = b;
		k2 = NextTriangle(pts, nPts, iEnd, iAdj, count, visited, ct,
				  a, &b2, &c, 1);
		TraceToStart(fp, pts, nPts, iEnd, iAdj, count, visited, ct, k2,
			     a, b2, c);
	      }
	    }
	  }
	  visited[k-1] = 1;
	}
      }
    if (iAdj[end-1] > -1)
      if (iAdj[end-1] > i && iAdj[start] > i) {
	a = i;
	b = iAdj[end-1];
	c = iAdj[start];

	k = TriFunction(nPts, iAdj, iEnd, count, a, b, c);
	if (!visited[k-1]) {
	  if (CrossTriangle(pts, ct, &a, &b, &c)) {
	    if (CrossEdge(pts, a, b, ct, cr)) {
	      nCurves++;
	      fprintf(fp, "M %f %f %f\n", cr[0], cr[1], ct);
	      
	      TraceToStart(fp, pts, nPts, iEnd, iAdj, count, visited, ct, k,
			   a, b, c);
	    }
	    else if (pts[a][2] == ct) {
	      nCurves++;
	      fprintf(fp, "M %f %f %f\n", pts[a][0], pts[a][1], ct);
	      
	      if (CrossEdge(pts, b, c, ct, cr))
		TraceToStart(fp, pts, nPts, iEnd, iAdj, count, visited, ct, k,
			     a, b, c);
	      else if (pts[b][2] == ct)
		TraceToStart(fp, pts, nPts, iEnd, iAdj, count, visited, ct, k,
			     c, a, b);
	      else if (pts[c][2] == ct)
		TraceToStart(fp, pts, nPts, iEnd, iAdj, count, visited, ct, k,
			     c, a, b);
	      else {
		b2 = b;
		k2 = NextTriangle(pts, nPts, iEnd, iAdj, count, visited, ct,
				  a, &b2, &c, 1);
		TraceToStart(fp, pts, nPts, iEnd, iAdj, count, visited, ct, k2,
			     a, b2, c);
	      }
	    }
	  }
	  visited[k-1] = 1;
	}
      }
  }
  fprintf(fp, "E\n");

  free_iarray1(visited);

  return nCurves;
}

short NextTriangle(double **pts, int nPts, int *iEnd, int *iAdj,
		   int *count, int *visited, double ct, int a, int *b,
		   int *c, int flag)
{
  double cr[2];
  int b2, c2, end;
  int crossed = 0, k, k2 = -1;

  end = *b;
  b2 = *c;
  do {
    c2 = FindThirdPoint(iEnd, iAdj, a, b2);
    k = TriFunction(nPts, iAdj, iEnd, count, a, b2, c2);
    if (!visited[k-1]) {
      if (CrossEdge(pts, b2, c2, ct, cr) || pts[c2][2] == ct) {
	if (!crossed) {
	  *b = b2;
	  *c = c2;
	  k2 = k;
	  crossed = 1;
	}
      }
      else
	visited[k-1] = 1;
    }
    b2 = c2;
  } while ((!flag && iAdj[iEnd[b2]] != -1) || (flag  && b2 != end));

  return k2;
}

void TraceToBoundary(FILE *fp, double **pts, int nPts, int *iEnd,
		     int *iAdj, int *count, int *visited, double ct,
		     int *a, int *b, int *c)
{
  double cr[2];
  int k;

  do {
    *c = FindThirdPoint(iEnd, iAdj, *a, *b);

    k = TriFunction(nPts, iAdj, iEnd, count, *a, *b, *c);

    TraceToNext(fp, pts, nPts, iEnd, iAdj, count, visited, ct, k, a, b, *c);
  } while (iAdj[iEnd[*a]] != -1 || iAdj[iEnd[*b]] != -1 ||
	   iAdj[iEnd[*a]-1] != *b);

  *c = iAdj[*b ? iEnd[*b-1]+1 : 0];
}

void TraceToStart(FILE *fp, double **pts, int nPts, int *iEnd, int *iAdj,
		  int *count, int *visited, double ct, int k, int a,
		  int b, int c)
{
  double cr[2];
  int end, next;

  end = k;
  do {
    next = TraceToNext(fp, pts, nPts, iEnd, iAdj, count, visited, ct, k,
		       &a, &b, c);

    c = FindThirdPoint(iEnd, iAdj, a, b);

    k = TriFunction(nPts, iAdj, iEnd, count, a, b, c);
  } while (k != end && next);
}

int TraceToNext(FILE *fp, double **pts, int nPts, int *iEnd, int *iAdj,
		int *count, int *visited, double ct, int k,
		int *a, int *b, int c)
{
  double cr[2];
  int a2, b2;
  int crossed = 0, k2, next = 0;

  if (!visited[k-1]) {
    next = 1;
    if (CrossEdge(pts, *b, c, ct, cr)) {
      fprintf(fp, "D %f %f %f\n", cr[0], cr[1], ct);
      *a = c;
    }
    else if (CrossEdge(pts, c, *a, ct, cr)) {
      fprintf(fp, "D %f %f %f\n", cr[0], cr[1], ct);
      *b = c;
    }
    else {     /* contour crosses vertex C */
      fprintf(fp, "D %f %f %f\n", pts[c][0], pts[c][1], ct);
      a2 = *a;
      *a = c;
      b2 = NextAdjB(iEnd, iAdj, *a, *b);
      while (b2 != a2) {
	c = FindThirdPoint(iEnd, iAdj, *a, b2);
	k2 = TriFunction(nPts, iAdj, iEnd, count, *a, b2, c);
	if (!visited[k2-1]) {
	  if (!crossed &&
	      (CrossEdge(pts, b2, c, ct, cr) || pts[c][2] == ct)) {
	    *b = b2;
	    crossed = 1;
	  }
	  else
	    visited[k2-1] = 1;
	}
	b2 = c;
      }
      next = crossed;
    }
    visited[k-1] = 1;
  }
  return next;
}

int CrossTriangle(double **pts, double ct, int *a, int *b, int *c)
{
  int k = -1;

  if      ((pts[*a][2] <= ct && ct <  pts[*b][2]) ||
	   (pts[*b][2] <  ct && ct <= pts[*a][2]))
    k = 0;
  else if ((pts[*b][2] <= ct && ct <  pts[*c][2]) ||
	   (pts[*c][2] <  ct && ct <= pts[*b][2])) {
    k  = *a;
    *a = *b;
    *b = *c;
    *c = k;
  }
  else if ((pts[*c][2] <= ct && ct <  pts[*a][2]) ||
	   (pts[*a][2] <  ct && ct <= pts[*c][2])) {
    k  = *a;
    *a = *c;
    *c = *b;
    *b = k;
  }
  return k+1;
}

int CrossEdge(double **pts, int a, int b, double ct, double *cr)
{
  double ba, ca;

  if ((pts[a][2] < ct && ct < pts[b][2]) ||
      (pts[b][2] < ct && ct < pts[a][2])) {
    ca =     ct    - pts[a][2];
    ba = pts[b][2] - pts[a][2];
    cr[0] = pts[a][0] + ca*(pts[b][0] - pts[a][0])/ba;
    cr[1] = pts[a][1] + ca*(pts[b][1] - pts[a][1])/ba;
    return 1;
  }
  return 0;
}

int FindThirdPoint(int *iEnd, int *iAdj, int a, int b)
{
  int i, start;

  start = (a ? iEnd[a-1]+1 : 0);
  for (i=start; i<=iEnd[a]; i++)
    if (iAdj[i] == b)
      break;
  return (i<iEnd[a] ? iAdj[i+1] : iAdj[start]);
}

int NextAdjB(int *iEnd, int *iAdj, int a, int b)
{
  int i, start;

  start = (a ? iEnd[a-1]+1 : 0);
  for (i=start; i<=iEnd[a]; i++)
    if (iAdj[i] == b)
      break;
  return (i < iEnd[a] && iAdj[i+1] == -1 ? iAdj[start] : b);
}

int TriFunction(int nPts, int *iAdj, int *iEnd, int *count,
		int a, int b, int c)
{
  int end, i, j, nTriang, start;

  if (b < a && b < c) {
    i = b;
    b = c;
    c = a;
  }
  else if (c < a && c < b) {
    i = c;
    c = b;
    b = a;
  }
  else
    i = a;

  nTriang = count[i];
  start = (i ? iEnd[i-1] + 1 : 0);
  end = iEnd[i] + 1;
  for (j=start; j+1 < end && iAdj[j+1] > -1; j++)
    if (iAdj[j] > i && iAdj[j+1] > i) {
      nTriang++;
      if (iAdj[j] == b && iAdj[j+1] == c)
	return nTriang;
    }
  if (iAdj[end-1] > -1)
    if (iAdj[end-1] > i && iAdj[start] > i) {
      nTriang++;
      if (iAdj[end-1] == b && iAdj[start] == c)
	return nTriang;
    }
  return -1;
}

void ConvertToCos(FILE *fp, int nCurves, ParCurv **egeoms)
{
  double ct, len, *lengths, u, v, u2, v2;
  char op[8];
  int i, inCurve = 0, j, nCurv = 0, *nPoints, nPts, order = 2, quit = 0;

  lengths = dbl_array1(nCurves);
  nPoints = int_array1(nCurves);

  while (!quit) {
    fscanf(fp, "%s", op);
    if (inCurve) {
      if (op[0] == (unsigned)'D') {
	fscanf(fp, "%lf %lf %lf", &u, &v, &ct);
	lengths[nCurv] += sqrt((u-u2)*(u-u2) + (v-v2)*(v-v2));
	nPts++;
	u2 = u;
	v2 = v;
      }
      else {
	nPoints[nCurv++] = nPts;
	inCurve = 0;
      }
    }
    if (!inCurve) {
      if (op[0] == (unsigned)'M') {
	fscanf(fp, "%lf %lf %lf", &u2, &v2, &ct);
	lengths[nCurv] = 0.0;
	nPts = 1;
	inCurve = 1;
      }
      else
	quit = 1;
    }
  }
  rewind(fp);

  for (i=0; i<nCurves; i++) {
    egeoms[i] = egeomalloc1(order, nPoints[i]);
    egeoms[i]->type = PCurveOpen;
    egeoms[i]->order = order;
    egeoms[i]->ncontpts = nPoints[i];

    egeoms[i]->knots = dbl_array1(order+nPoints[i]);
    egeoms[i]->kmem = order + nPoints[i];

    egeoms[i]->contpts = vec_array1(nPoints[i]);
    egeoms[i]->pmem = nPoints[i];

    egeoms[i]->knots[0] = 0.0;
    for (j=0; j<nPoints[i]; j++) {
      fscanf(fp, "%s %lf %lf %lf", op, &u, &v, &ct);
      if (j)
	len += sqrt((u-u2)*(u-u2) + (v-v2)*(v-v2));
      else
	len = 0.0;

      egeoms[i]->knots[j+1] = len/lengths[i];

      egeoms[i]->contpts[j]->x = u2 = u;
      egeoms[i]->contpts[j]->y = v2 = v;
      egeoms[i]->contpts[j]->z = 0.0;
      egeoms[i]->contpts[j]->w = 1.0;
    }
    egeoms[i]->knots[order+nPoints[i]-1] = 1.0;
  }
  free_darray1(lengths);
  free_iarray1(nPoints);
}

vector *LinearCosNormal(ParCurv *egeom, double t1, double t2)
{
  vector *r1, *r2;

  r1 = evalbsp(egeom, t1);
  r2 = evalbsp(egeom, t2);

  sub_vect1(r2, r1, r1);
  unitvector1(r1, r1);

  r2->x = -r1->y;
  r2->y =  r1->x;
  r2->z =  0.0;
  r2->w =  1.0;

  vectfree(r1);
  
  return r2;
}

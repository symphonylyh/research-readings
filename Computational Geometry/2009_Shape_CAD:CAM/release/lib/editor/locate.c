/* Copyright (C) 1993 Massachusetts Institute of Technology, Cambridge, MA
   All rights reserved
*/

/* locate.c - locate triangle enclosing input x,y
*/

/* GetAdj()
   InterpolateZ()
   LocateTri()
   LocateTri2()
   RightOf()
*/

#include "editor.h"

/***********************************************************************
  LocateTri() - locate triangle enclosing input x,y for concex
                triangulations.  If the triangulation may be
		concave, use LocateTri2()
***********************************************************************/

int LocateTri(double **pts, int *iEnd, int *iAdj, double x, double y,
	      int i, int j, int k, int *a, int *b, int *c)
{
  int n = 0;

  if (RightOf(x, y, pts[i][0], pts[i][1], pts[j][0], pts[j][1])) {
    n = GetAdj(iEnd, iAdj, j, i);
    if (n > -1) {
      k = n;
      n = LocateTri(pts, iEnd, iAdj, x, y, j, i, k, a, b, c);
      return n;       /* return 1 if point is outside, 0 if inside */
    }
  }
  else if (RightOf(x, y, pts[j][0], pts[j][1], pts[k][0], pts[k][1])) {
    n = GetAdj(iEnd, iAdj, k, j);
    if (n > -1) {
      i = n;
      n = LocateTri(pts, iEnd, iAdj, x, y, k, j, i, a, b, c);
      return n;       /* return 1 if point is outside, 0 if inside */
    }
  }
  else if (RightOf(x, y, pts[k][0], pts[k][1], pts[i][0], pts[i][1])) {
    n = GetAdj(iEnd, iAdj, i, k);
    if (n > -1) {
      j = n;
      n = LocateTri(pts, iEnd, iAdj, x, y, i, k, j, a, b, c);
      return n;       /* return 1 if point is outside, 0 if inside */
    }
  }
  *a = i;
  *b = j;
  *c = k;

  return (n == -1);   /* return 1 if point is outside, 0 if inside */
}

int RightOf(double ax, double ay, double bx, double by, double cx, double cy)
{   /* cross product to determine if to right of vector */
  return (((ax-bx)*(cy-by) - (ay-by)*(cx-bx)) > 0.0);
}

int GetAdj(int *iEnd, int *iAdj, int a, int b)
{
  int i, start;

  start = (a ? iEnd[a-1]+1 : 0);
  for (i=start; i<=iEnd[a]; i++)
    if (iAdj[i] == b)
      break;
  return (i < iEnd[a] ? iAdj[i+1] : iAdj[start]);
}

double InterpolateZ(double **pts, double x, double y, int a, int b, int c)
{
  double d, e, f, g;   /* interpolate z value at x,y given corner */
                       /* z values of enclosing triangle          */
  d = (pts[b][1]-pts[a][1])*(pts[c][2]-pts[a][2]) -
      (pts[b][2]-pts[a][2])*(pts[c][1]-pts[a][1]);
  e = (pts[b][2]-pts[a][2])*(pts[c][0]-pts[a][0]) -
      (pts[b][0]-pts[a][0])*(pts[c][2]-pts[a][2]);
  f = (pts[b][0]-pts[a][0])*(pts[c][1]-pts[a][1]) -
      (pts[b][1]-pts[a][1])*(pts[c][0]-pts[a][0]);
  g = d*pts[b][0] + e*pts[b][1] + f*pts[b][2];

  return ((d*x + e*y - g) / -f);
}

/***********************************************************************
  LocateTri2() - locate triangle enclosing input x,y for general concave
                 triangulations.  If the triangulation is known to be
		 convex, use LocateTri()
***********************************************************************/

int LocateTri2(double **pts, int *iEnd, int *iAdj, double x, double y,
	       int i, int j, int k, int *iMark, int *a, int *b, int *c)
{
  int ii, m, n = 0, start;

  start = (i ? iEnd[i-1]+1 : 0);  /* start = index into iAdj and iMark of i */
  for (ii=start; ii<=iEnd[i]; ii++)
    if (iAdj[ii] == j)
      break;                         /* ii = index into iAdj and iMark of j */
  m =  (ii < iEnd[i] ? ii+1 : start); /* m = index into iAdj and iMark of k */
  if (iMark[m]) {   /* have we already visited this triangle */
    *a = i;         /* if so, then point is outside */
    *b = j;
    *c = k;
    return 1;       /* return 1 if point is outside, 0 if inside */
  }
  iMark[m] = 1;     /* mark this triangle as visited */

  if (RightOf(x, y, pts[i][0], pts[i][1], pts[j][0], pts[j][1])) {
    n = GetAdj(iEnd, iAdj, j, i);
    if (n > -1) {
      k = n;
      n = LocateTri2(pts, iEnd, iAdj, x, y, j, i, k, iMark, a, b, c);
      return n;         /* return 1 if point is outside, 0 if inside */
    }
    else {                               /* if triangluation is concave */
      n = (j ? iAdj[iEnd[j-1]+1] : iAdj[0]);   /* jump over concavity */
      if (n != k) {                      /* to rightmost adjacent triangle */
	i = j;
	j = n;   /* make sure that we don't get into an infinite loop */
	k = (i ? iAdj[iEnd[i-1]+2] : iAdj[1]);
	n = LocateTri2(pts, iEnd, iAdj, x, y, i, j, k, iMark, a, b, c);
	return n;   /* return 1 if point is outside, 0 if inside */
      }
      else {
	n = iAdj[iEnd[i]-1];
	if (n != k) {                      /* or leftmost adjacent triangle */
	  j = i;
	  i = n;
	  k = GetAdj(iEnd, iAdj, i, j);
	  n = LocateTri2(pts, iEnd, iAdj, x, y, i, j, k, iMark, a, b, c);
	  return n;     /* return 1 if point is outside, 0 if inside */
	}
      }
    }
  }
  else if (RightOf(x, y, pts[j][0], pts[j][1], pts[k][0], pts[k][1])) {
    n = GetAdj(iEnd, iAdj, k, j);
    if (n > -1) {
      i = n;
      n = LocateTri2(pts, iEnd, iAdj, x, y, k, j, i, iMark, a, b, c);
      return n;
    }
    else {
      n = (k ? iAdj[iEnd[k-1]+1] : iAdj[0]);
      if (n != i) {
	i = k;
	j = n;
	k = (i ? iAdj[iEnd[i-1]+2] : iAdj[1]);
	n = LocateTri2(pts, iEnd, iAdj, x, y, i, j, k, iMark, a, b, c);
	return n;
      }
      else {
	n = iAdj[iEnd[j]-1];
	if (n != i) {
	  i = n;
/*        j = j; */
	  k = GetAdj(iEnd, iAdj, i, j);
	  n = LocateTri2(pts, iEnd, iAdj, x, y, i, j, k, iMark, a, b, c);
	  return n;
	}
      }
    }
  }
  else if (RightOf(x, y, pts[k][0], pts[k][1], pts[i][0], pts[i][1])) {
    n = GetAdj(iEnd, iAdj, i, k);
    if (n > -1) {
      j = n;
      n = LocateTri2(pts, iEnd, iAdj, x, y, i, k, j, iMark, a, b, c);
      return n;
    }
    else {
      n = (i ? iAdj[iEnd[i-1]+1] : iAdj[0]);
      if (n != j) {
/*      i = i; */
	j = n;
	k = (i ? iAdj[iEnd[i-1]+2] : iAdj[1]);
	n = LocateTri2(pts, iEnd, iAdj, x, y, i, j, k, iMark, a, b, c);
	return n;
      }
      else {
	n = iAdj[iEnd[k]-1];
	if (n != j) {
	  i = n;
	  j = k;
	  k = GetAdj(iEnd, iAdj, i, j);
	  n = LocateTri2(pts, iEnd, iAdj, x, y, i, j, k, iMark, a, b, c);
	  return n;
	}
      }
    }
  }
  *a = i;
  *b = j;
  *c = k;

  return 0;   /* return 1 if point is outside, 0 if inside */
}

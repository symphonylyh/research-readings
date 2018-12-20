/* Copyright (C) 1994 Massachusetts Institute of Technology, Cambridge, MA
   All rights reserved
*/

/* alpha_shape.c */

/* AlphaShape()
*/

#include <math.h>
#include "gen.h"
#include "editor.h"

/* Alpha-shape to determine the bounary points of a cloud of 2D points */
/* H. Edelsbrunner et al., "On the Shape of a Set of Points in the Plane,"
   IEEE Transactions on Information Theory IT-29(4):551-559, Julu 1932 */

/* Note: In this program alpha is the postive radius of the discs */
/* In the paper alpha is the inverse of radius and negative */

int AlphaShape(int nPts, double **pts, int *iEnd, int *iAdj,
	       double alpha, int *iAlpha, int *nAlpha)
{
  double ab, d, dmax;
  int **a_extreme, end, i, j, k, cancel = 0, start;

/* find alpha-shape */

/* find a-extreme points  of set S of all points */
/* 1) all points P on convex hull and */
/* 2) interior points P that lie on a disc of radius alpha that does not
   include any other points of S.  The set of centers of these discs is the
   Voronoi region containing P.  To find maximum radius calculate d(P,x)
   for x = vertices of Voronoi region */

  *nAlpha = 0;
  a_extreme = int_array2(nPts, 2);

  for (i=0; i<nPts; i++)
    if (iAdj[iEnd[i]] < 0) {
      a_extreme[i][0] =  1;
      a_extreme[i][1] = -1; /* flag for a-neighbor not set */
    }
    else {
      dmax = 0.0;
      start = (i ? iEnd[i-1]+1 : 0);
      end = iEnd[i];
      for (j=start; j<=end; j++) {
	ab = sqrt((pts[i][0]-pts[iAdj[j]][0])*(pts[i][0]-pts[iAdj[j]][0]) +
		  (pts[i][1]-pts[iAdj[j]][1])*(pts[i][1]-pts[iAdj[j]][1]));
	k = (j < end ? j+1 : start);
	d = Circumradius(pts[i], pts[iAdj[j]], pts[iAdj[k]], ab);
	if (d > dmax) dmax = d;
      }
      if (alpha <= dmax) {
	a_extreme[i][0] =  1;
	a_extreme[i][1] = -1; /* flag for a-neighbor not set */
      }
    }

/* find alpha-neighbors, ie. edges E of alpha-shape */
/* 1) E is convex hull edge and */
/* 2) E is interior edge */

  for (i=0; i<nPts; i++)
    if (a_extreme[i][0]) {
      start = (i ? iEnd[i-1]+1 : 0);
      end = iEnd[i];
      if (iAdj[end] < 0) {
	ab = sqrt((pts[i][0]-pts[iAdj[start]][0])*
		  (pts[i][0]-pts[iAdj[start]][0]) +
		  (pts[i][1]-pts[iAdj[start]][1])*
		  (pts[i][1]-pts[iAdj[start]][1]));
	d = Circumradius(pts[i], pts[iAdj[start]], pts[iAdj[start+1]], ab);
	if (alpha >= d)
	  if (a_extreme[i][1] != iAdj[start]) {
	    a_extreme[i][1] = iAdj[start];
	    (*nAlpha)++;
	  }
      }
      else {
	for (j=start; j<=end; j++)
	  if (a_extreme[iAdj[j]][0]) {
	    k = (j < end ? j+1 : start);
	    ab = sqrt((pts[i][0]-pts[iAdj[j]][0])*
		      (pts[i][0]-pts[iAdj[j]][0]) +
		      (pts[i][1]-pts[iAdj[j]][1])*
		      (pts[i][1]-pts[iAdj[j]][1]));
	    d = Circumradius(pts[i], pts[iAdj[j]], pts[iAdj[k]], ab);

	    k = (j > start ? j-1 : end);
	    ab = sqrt((pts[i][0]-pts[iAdj[k]][0])*
		      (pts[i][0]-pts[iAdj[k]][0]) +
		      (pts[i][1]-pts[iAdj[k]][1])*
		      (pts[i][1]-pts[iAdj[k]][1]));
	    dmax = Circumradius(pts[i], pts[iAdj[k]], pts[iAdj[j]], ab);

	    if (d < dmax && d <= alpha && alpha <= dmax) {
	      if (a_extreme[i][1] != iAdj[j]) {
		a_extreme[i][1] = iAdj[j];
		(*nAlpha)++;
	      }
	    }
	    else if (dmax <= alpha && alpha <= d) {
	      if (a_extreme[iAdj[j]][1] != i) {
		a_extreme[iAdj[j]][1] = i;
		(*nAlpha)++;
	      }
	    }
	  }
      }
    }

  for (i=0; i<nPts; i++)
    if (a_extreme[i][0] && a_extreme[i][1] < 0) {
      cancel = 1;
      break;
    }
  if (!cancel) {
    for (i=0; i<nPts; i++)   /* remove interior loops */
      if (a_extreme[i][0]) { /* find first point on alpha-hull */
	start = i;           /* start is first point */
	j = a_extreme[i][1]; /* j is the next point */
	break;
      }

    for (i=0,k=1; i<(*nAlpha) && j != start; i++,k++)
      j = a_extreme[j][1];   /* count the points on the alpha-hull */
    (*nAlpha) = k;

    for (i=0; i<nPts; i++)
      if (a_extreme[i][0]) { /* find the first point on alpha-hull */
	j = i;               /* j is first point */
	break;
      }

    for (i=0; i<(*nAlpha); i++) {
      iAlpha[i] = j;         /* copy alpha-hull points to iAlpha[i] */
      j = a_extreme[j][1];
    }
  }

  free_iarray2(a_extreme);

  return cancel;
}

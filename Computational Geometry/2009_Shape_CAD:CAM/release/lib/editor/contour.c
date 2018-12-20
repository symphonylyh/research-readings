/* Copyright (C) 1993 Massachusetts Institute of Technology, Cambridge, MA
 *  All rights reserved
 */

/* contour.c */

/* ContourFacets()
 *  IsAbove()
 */

#include <stdio.h>
#include <unistd.h>
#include <malloc.h>
#include "gen.h"
#include "bspl.h"
#include "editor.h"

char *tempnam(const char *, const char *);

int ContourFacets(double **pts, int nPts, int *iEnd, int *iAdj,
		  double ct, ParCurv ***egeoms, int ***inner,
		  int ***outer, int **above, int *nBoundary)
{
  static char *tols[2] = {"out of", "within"};
  static double ray[2] = {1.0, 0.0};

  FILE *fp;
  vector *norm, *r;
  double **point1, **point2, scal, vmin, vmax, z;
  int a, b, c;
  int *count, *orient;
  int end, i, inside, j, k, nCurves, ni, nTriang, start;
  char *tmpfil;

  count = int_array1(nPts);

  nTriang = 0;   /* number of triangles adjacent to each point */
  for (i=start=0; i<nPts; i++) {
    count[i] = nTriang;
    end = iEnd[i] + 1;
    for (j=start; j+1 < end && iAdj[j+1] > -1; j++)
      if (iAdj[j] > i && iAdj[j+1] > i)
	nTriang++;
    if (iAdj[end-1] > -1)
      if (iAdj[end-1] > i && iAdj[start] > i)
	nTriang++;
    start = end;
  }

  tmpfil = tempnam("/usr/tmp", "cont");
  if (fp = fopen(tmpfil, "w")) {
    nCurves = TraceContour(fp, pts, nPts, iEnd, iAdj, nTriang, count, ct,
			   nBoundary);
    fclose(fp);

    if (fp = fopen(tmpfil, "r")) {
      (*egeoms) = (ParCurv **)gen_array1(nCurves, sizeof(ParCurv *));
      ConvertToCos(fp, nCurves, *egeoms);
      fclose(fp);

      (*above) = int_array1(nCurves);
      orient = int_array1(nCurves);
      point1 = dbl_array2(nCurves, 2);
      point2 = dbl_array2(nCurves, 2);

      for (i=0; i<nCurves; i++) {
	vmin =  1.0e30;
	vmax = -1.0e30;
	for (j=0; j<(*egeoms)[i]->ncontpts; j++) {
	  if ((*egeoms)[i]->contpts[j]->y < vmin)
	    vmin = (*egeoms)[i]->contpts[j]->y;
	  if ((*egeoms)[i]->contpts[j]->y > vmax)
	    vmax = (*egeoms)[i]->contpts[j]->y;
	}
	     /* is curve counter-clockwise? */
	orient[i] = RecursiveInt((*egeoms)[i], vmin, vmax, 1);

	     /* find point inside of curve */
	inside = 0;
	for (j=(*egeoms)[i]->order-1;!inside&&j<(*egeoms)[i]->ncontpts-1;j++){
	  r = evalbsp((*egeoms)[i], (*egeoms)[i]->knots[j]);
	  point1[i][0] = r->x/r->w;
	  point1[i][1] = r->y/r->w;
	  vectfree(r);

	  norm = LinearCosNormal((*egeoms)[i], (*egeoms)[i]->knots[j],
				 (*egeoms)[i]->knots[j+1]);

	  for (k=0, scal=0.01; !inside && k<10; k++, scal /= 2) {
	    point2[i][0] = point1[i][0] + scal*norm->x;
	    point2[i][1] = point1[i][1] + scal*norm->y;
	    ni = loop_rayintersect(&(*egeoms)[i], 1, ray, point2[i]);
	    if ((ni/2)*2 != ni) {
	      inside = 1;
	      break;
	    }

	    point2[i][0] = point1[i][0] - scal*norm->x;
	    point2[i][1] = point1[i][1] - scal*norm->y;
	    ni = loop_rayintersect(&(*egeoms)[i], 1, ray, point2[i]);
	    if ((ni/2)*2 != ni) {
	      inside = 1;
	      break;
	    }
	  }
	  vectfree(norm);
	}
	                  /* is inside of curve above or below ct? */
	(*above)[i] = IsAbove((*egeoms)[i], pts, nPts, iEnd, iAdj, ct);
      }

 	    /* is inside of curve interior of trimmed patch? */

      for (i=0; i<nCurves; i++) {
	ni = loop_rayintersect((*egeoms), nCurves, ray, point2[i]);

	if ((ni/2)*2 != ni) {  /* odd number of intersections; inside */
	  if (!orient[i]) {
	    RevCurvParam((*egeoms)[i]);
	    orient[i] = 1;
	  }
	}
	else {
	  if (orient[i]) {     /* even number of intersections; outside */
	    RevCurvParam((*egeoms)[i]);
	    orient[i] = 0;
	  }
	}
      }

      (*inner) = int_array2(nCurves, nCurves);
      (*outer) = int_array2(nCurves, 2);

      for (i=0; i<nCurves; i++)
	if (orient[i])
	  for (j=0; j<nCurves; j++)
	    if (j != i) {
	      ni = loop_rayintersect(&(*egeoms)[i], 1, ray, point1[j]);
	      if ((ni/2)*2 != ni) {
		(*inner)[i][j] = 1;
		(*outer)[i][1]++;
		(*outer)[j][0] = 1;
	      }
	    }
      free_darray2(point1);
      free_darray2(point2);

      *nBoundary = 0;
      for (i=0; i<nCurves; i++)
	if (!(*outer)[i][0])
	  (*nBoundary)++;
	
      free_iarray1(orient);

    }
    unlink(tmpfil);
  }
  free(tmpfil);
  free_iarray1(count);
  
  return nCurves;
}

int IsAbove(ParCurv *egeom, double **pts, int nPts, int *iEnd, int *iAdj,
	    double ct)
{
  double u, umin = 1.0, ut, vmax = -1.0, vmin = 1.0, v, zt;
  int a, b, c;
  int above, haveIt = 0, i, j, *iMark;

  for (i=0; i<egeom->ncontpts; i++) {
    if (egeom->contpts[i]->y < vmin) vmin = egeom->contpts[i]->y;
    if (egeom->contpts[i]->y > vmax) vmax = egeom->contpts[i]->y;
  }

  v = (vmin + vmax)/2.0;

  j = egeom->ncontpts - 1;
  for (i=0; i<egeom->ncontpts; i++) {
    if ((egeom->contpts[j]->y < v && v < egeom->contpts[i]->y) ||
	(egeom->contpts[i]->y < v && v < egeom->contpts[j]->y)) {
      u = (egeom->contpts[j]->x - egeom->contpts[i]->x) * 
	  (          v          - egeom->contpts[i]->y) /
	  (egeom->contpts[j]->y - egeom->contpts[i]->y) +
	  egeom->contpts[i]->x;
      if (u < umin) umin = u;
    }
    j = i;
  }
  j = iEnd[nPts-1] + 1;
  iMark = int_array1(j);
  for (i=0; i<j; i++)
    iMark[i] = 0;     /* mark all triangles as unvisited */
  j = LocateTri2(pts, iEnd, iAdj, umin, v, 0, iAdj[0], iAdj[1], iMark,
		 &a, &b, &c);
  free_iarray1(iMark);

  if ((pts[a][1] < v && v < pts[b][1]) ||
      (pts[b][1] < v && v < pts[a][1])) {
    ut = (pts[a][0] - pts[b][0]) * (v - pts[b][1]) / (pts[a][1] - pts[b][1]) +
         pts[b][0];
    zt = (pts[a][2] - pts[b][2]) * (v - pts[b][1]) / (pts[a][1] - pts[b][1]) +
         pts[b][2];
    if (ut != umin)
      haveIt = 1;
  }

  if (!haveIt) {
    if ((pts[b][1] < v && v < pts[c][1]) ||
	(pts[c][1] < v && v < pts[b][1])) {
      ut = (pts[b][0]-pts[c][0]) * (v - pts[c][1]) / (pts[b][1] - pts[c][1]) +
	   pts[c][0];
      zt = (pts[b][2]-pts[c][2]) * (v - pts[c][1]) / (pts[b][1] - pts[c][1]) +
	   pts[c][2];
      if (ut != umin)
	haveIt = 1;
    }

    if (!haveIt)
      if ((pts[c][1] < v && v < pts[a][1]) ||
	  (pts[a][1] < v && v < pts[c][1])) {
	ut = (pts[c][0]-pts[a][0]) * (v-pts[a][1]) / (pts[c][1] - pts[a][1]) +
  	     pts[a][0];
	zt = (pts[c][2]-pts[a][2]) * (v-pts[a][1]) / (pts[c][1] - pts[a][1]) +
	     pts[a][2];
      }
  }

  if (ct >= 0.0) {
    if (ut < umin)
      above = (zt <= ct);
    else
      above = (zt >= ct);
  }
  else {
    if (ut < umin)
      above = (zt >= ct);
    else
      above = (zt <= ct);
  }

  return above;
}

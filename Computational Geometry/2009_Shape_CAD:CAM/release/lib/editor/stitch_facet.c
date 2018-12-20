/* Copyright (C) 1998 Massachusetts Institute of Technology, Cambridge, MA
 * All rights reserved
 */

/* stitch_facet.c */

/* copyAdj()
 * stitch_facet()
 * updateAdj()
 */

#include "gen.h"
#include "editor.h"

void stitch_facet(int nPts, double **pts, int **iAdj, int *iEnd, int inU)
{
  double xmax1, xmin1, xmax2, xmin2, ymax1, ymin1, ymax2, ymin2;
  int
    atEnd1, atEnd2,
    extra,
    first,            /* first point on boundary */
    i, j, n,
    *jAdj,            /* adjacency matrix */
    *jEnd,            /* adjacency record sizes */
    inside1, inside2,
    l0, l1,           /* first and last points to stitch on left-hand side */
    last,
    ll, lr,           /* lower-left and lower-right points */
    nLeft,            /* number of points between ul and ll */
    nRight,           /* number of points between lr and ur */
    r0, r1,           /* first and last points to stitch on right-hand side */
    that, this,
    ul, ur;           /* upper-left and upper-right points */

  /* find a point on the boundary */

  for (i=0; i<nPts; i++)
    if ((*iAdj)[iEnd[i]] < 0) {
      first = i;                      /* first boundary point */
      break;
    }

  if (inU) {                          /* stitch in U direction */

    /* march along the boundary finding the bottom-most and top-most points */
    /* because the adjacency graph is built with adjacent points given in */
    /* counter-clockwise order, the next point on the boundary is the first */
    /* point adjacent to the current point */

    ymin1 = ymax1 = pts[first][1];
    ll = ul = first;

    this = first;
    while ((this = (*iAdj)[this ? iEnd[this-1]+1 : 0]) != first) {
      if (pts[this][1] < ymin1) {
	ll = this;
	ymin1 = pts[this][1];
      }
      if (pts[this][1] > ymax1) {
	ul = this;
	ymax1 = pts[this][1];
      }
    }

    /* of (the possibly many) points with extreme y values, find the */
    /* left-most and right-most */

    xmin1 = xmax1 = pts[ll][0];
    lr = ll;
    xmin2 = xmax2 = pts[ul][0];
    ur = ul;

    this = first;
    while ((this = (*iAdj)[this ? iEnd[this-1]+1 : 0]) != first) {
      if (pts[this][1] == ymin1) {
	if (pts[this][0] < xmin1) {
	  ll = this;
	  xmin1 = pts[this][0];
	}
	if (pts[this][0] > xmax1) {
	  lr = this;
	  xmax1 = pts[this][0];
	}
      }
      if (pts[this][1] == ymax1) {
	if (pts[this][0] < xmin2) {
	  ul = this;
	  xmin2 = pts[this][0];
	}
	if (pts[this][0] > xmax2) {
	  ur = this;
	  xmax2 = pts[this][0];
	}
      }
    }
  }
  else {                              /* stitch in V direction */

    /* march along the boundary finding the left-most and right-most points */
    /* because the adjacency graph is built with adjacent points given in */
    /* counter-clockwise order, the next point on the boundary is the first */
    /* point adjacent to the current point */

    xmin1 = xmax1 = pts[first][0];
    ll = lr = first;

    this = first;
    while ((this = (*iAdj)[this ? iEnd[this-1]+1 : 0]) != first) {
      if (pts[this][0] < xmin1) {
	ll = this;
	xmin1 = pts[this][0];
      }
      if (pts[this][0] > xmax1) {
	lr = this;
	xmax1 = pts[this][0];
      }
    }

    /* of (the possibly many) points with extreme x values, find the */
    /* top-most and bottom-most */

    ymin1 = ymax1 = pts[ll][1];
    ul = ll;
    ymin2 = ymax2 = pts[lr][1];
    ur = lr;

    this = first;
    while ((this = (*iAdj)[this ? iEnd[this-1]+1 : 0]) != first) {
      if (pts[this][0] == xmin1) {
	if (pts[this][1] < ymin1) {
	  ll = this;
	  ymin1 = pts[this][1];
	}
	if (pts[this][1] > ymax1) {
	  ul = this;
	  ymax1 = pts[this][1];
	}
      }
      if (pts[this][0] == xmax1) {
	if (pts[this][1] < ymin2) {
	  lr = this;
	  ymin2 = pts[this][1];
	}
	if (pts[this][1] > ymax2) {
	  ur = this;
	  ymax2 = pts[this][1];
	}
      }
    }
  }

  /* ll is the lower-left point on the boundary */
  /* ul is the upper-left point on the boundary */
  /* lr is the lower-right point on the boundary */
  /* ur is the upper-right point on the boundary */

  if (inU) {                          /* stitch in U direction */

    /* count the number points (inclusive) between ll and ul, and lr and ur */

    this = ul;
    for (nLeft = 2; (this = (*iAdj)[this ? iEnd[this-1]+1 : 0]) != ll;
	 nLeft++)
      ;
    this = lr;
    for (nRight = 2; (this = (*iAdj)[this ? iEnd[this-1]+1 : 0]) != ur;
	 nRight++)
      ;
  }
  else {                              /* stitch in V direction */

    /* count the number points (inclusive) between ll and lr, and ul and ur */

    this = ll;
    for (nLeft = 2; (this = (*iAdj)[this ? iEnd[this-1]+1 : 0]) != lr;
	 nLeft++)
      ;                               /* nLeft is really counting the bottom */
    this = ur;
    for (nRight = 2; (this = (*iAdj)[this ? iEnd[this-1]+1 : 0]) != ul;
	 nRight++)
      ;                               /* nRight is really counting the top */
  }

  /* stitch the left and right sides of the triangulation together */

  if (nLeft >= nRight) {   /* left side has more points than the right side */

    /* n1 = nLeft - nRight + 1 */
    /* n2 = nRight - 3 */
    /* n1 is the number of new adjacencies at the top of the stitch */
    /* n2 is the number of new adjacencies pairs at the bottom of */
    /* the stitch */

    /* each of n1 adds 2 updates (one for each end point) to the matrix */
    /* each of n2 adds 4 updates (each end point for each of 2 adjacencies) */

    extra = 2*(nLeft - nRight + 1) + 4*(nRight - 3);
  }
  else if (nRight >= nLeft) {

    /* n1 = nRight - nLeft + 1 */
    /* n2 = nLeft - 3 */

    extra = 2*(nRight - nLeft + 1) + 4*(nLeft - 3);
  }

  /* make working copy of the adjacency matrix with room to insert the */
  /* new adjacencies */

  copyAdj(nPts, *iAdj, iEnd, extra, &jAdj, &jEnd);

  if (inU) {

    /* since ll may already be connected to lr (and ul to ur) */
    /* we start and end the stitching one point in from the extreme points */

    r0 = (*iAdj)[lr ? iEnd[lr-1]+1 : 0];    /* first point up from ur */
    r1 = (*iAdj)[iEnd[ur]-1];          /* first point down from upper right */
    l0 = (*iAdj)[iEnd[ll]-1];          /* first point up from lower left */
    l1 = (*iAdj)[ul ? iEnd[ul-1]+1 : 0];    /* first point down from ll */
    last = lr;                         /* pointer to previous "this" value */
    this = r0;                         /* right-side pointer */
    that = ll;                         /* left-side pointer */
  }
  else {

    r0 = (*iAdj)[ur ? iEnd[ur-1]+1 : 0];     /* first point left of ur */
    r1 = (*iAdj)[iEnd[ul]-1];          /* first point right from upper left */
    l0 = (*iAdj)[iEnd[lr]-1];          /* first point left of lower right */
    l1 = (*iAdj)[ll ? iEnd[ll-1]+1 : 0];     /* first point right of ll */
    last = ur;                         /* pointer to previous "this" value */
    this = r0;                         /* right-side pointer */
    that = lr;                         /* left-side pointer */
  }

  /* stitch together the triangulation by alternately incrementing */
  /* the left-side ("that") and right-side ("this") pointers and */
  /* creating an adjacency between "this" and "that" */

  /* whenever "this" or "that" reaches its limit (r1 or l1) */
  /* stop incrementing that pointer */

  /* the "inside?" flags indicate whether the point is on the boundary (0) */
  /* or in the interior (1) AFTER the adjacency is added */

  /* the "atEnd?" flags indicate whether the adjacency for that point */
  /* should be inserted at the end (1) or at the beginning (0) of */
  /* the existing adjacency record */

  inside1 = 0;
  atEnd1 = 1;
  atEnd2 = 1;
  while (last != r1 || that != l1) {
    last = this;

    if (that != l1)
      that = (*iAdj)[iEnd[that]-1];   /* increment left-side pointer */
    if (that == l0 || that == l1)
      inside2 = 0;
    else
      inside2 = 1;
    if (that == l1)
      atEnd2 = 0;

    if (this != r1) {
      this = (*iAdj)[this ? iEnd[this-1]+1 : 0]; /* increment right-side */

      if (that != l1) {
	if (this == r0 || this == r1)
	  inside1 = 0;
	else
	  inside1 = 1;
	atEnd1 = 1;
	updateAdj(this, inside1, atEnd1, that, inside2, atEnd2, nPts, jAdj,
		  jEnd);
      }
    }

    if (last == r0 || last == r1)
      inside1 = 0;
    else
      inside1 = 1;
    if (last == r0)
      atEnd1 = 0;
    else
      atEnd1 = 1;
    updateAdj(last, inside1, atEnd1, that, inside2, atEnd2, nPts, jAdj,
	      jEnd);
  }

  for (i=0; i<nPts; i++)             /* copy new record size array */
    iEnd[i] = jEnd[i];
  free_iarray1(jEnd);                /* reclaim dynamic memory */

  free_iarray1(*iAdj);
  *iAdj = int_array1(iEnd[nPts-1]+1);  /* copy new adjacency array */
  for (i=0; i<=iEnd[nPts-1]; i++)
    (*iAdj)[i] = jAdj[i];
  free_iarray1(jAdj);                /* reclaim dynamic memory */
}

void copyAdj(int nPts, int *iAdj, int *iEnd, int extra, int **jAdj, int **jEnd)
{
  int i;

  *jEnd = int_array1(nPts);
  for (i=0; i<nPts; i++)
    (*jEnd)[i] = iEnd[i];

  *jAdj = int_array1(iEnd[nPts-1]+1 + extra);
  for (i=0; i<=iEnd[nPts-1]; i++)
    (*jAdj)[i] = iAdj[i];
}

void updateAdj(int n1, int inside1, int atEnd1, int n2, int inside2,
	       int atEnd2, int nPts, int *jAdj, int *jEnd)
{
  int i;

  /* if n1 is a boundary point, or if it was a boundary but the boundary */
  /* flag has already been overwritten, we need to slide up by 1 the */
  /* adjacency records for all points > n1 */

  /**************************************************
  printf("# connect %2d(%c,%c) to %2d(%c,%c)\n", n1+1, inside1 ? 'i' : 'b',
	 atEnd1 ? 'E' : 'B', n2+1, inside2 ? 'i' : 'b', atEnd2 ? 'E' : 'B');
  **************************************************/

  if (!inside1 || jAdj[jEnd[n1]] >= 0) {
    for (i=jEnd[nPts-1]; i>jEnd[n1]; i--) /* slide up by 1 the adjacencies */
      jAdj[i+1] = jAdj[i];
    for (i=n1+1; i<nPts; i++)           /* increment the record sizes by 1 */
      jEnd[i] += 1;
  }

  if (!atEnd1) {                        /* insert at beginning of record */
    for (i=jEnd[n1]; i>=(n1 ? jEnd[n1-1]+1 : 0);  i--) /* slide up by 1 */
      jAdj[i+1] = jAdj[i];
    jAdj[n1 ? jEnd[n1-1]+1 : 0] = n2;   /* insert the point */
    jEnd[n1] += 1;
  }
  else {                                /* insert at end of record */
    if (!inside1) {
      jAdj[jEnd[n1]] = n2;
      jEnd[n1] += 1;                    /* remember to move up by 1 the */
      jAdj[jEnd[n1]] = -1;              /* boundary flag */
    }
    else if (jAdj[jEnd[n1]] >= 0) {
      jEnd[n1] += 1;                    /* if we are not overwriting the */
      jAdj[jEnd[n1]] = n2;              /* boundary flag, move up by 1 */
    }
    else
      jAdj[jEnd[n1]] = n2;              /* overwrite the boundary flag */
  }

  if (!inside2 || jAdj[jEnd[n2]] >= 0) {
    for (i=jEnd[nPts-1]; i>jEnd[n2]; i--) /* slide up by 1 the adjacencies */
      jAdj[i+1] = jAdj[i];
    for (i=n2+1; i<nPts; i++)           /* increment the record sizes by 1 */
      jEnd[i] += 1;
  }

  if (!atEnd2) {                        /* insert at beginning of record */
    for (i=jEnd[n2]; i>=(n2 ? jEnd[n2-1]+1 : 0);  i--) /* slide up by 1 */
      jAdj[i+1] = jAdj[i];
    jAdj[n2 ? jEnd[n2-1]+1 : 0] = n1;   /* insert the point */
    jEnd[n2] += 1;
  }
  else {                                /* insert at end of record */
    if (!inside2) {
      jAdj[jEnd[n2]] = n1;
      jEnd[n2] += 1;                    /* remember to move up by 1 the */
      jAdj[jEnd[n2]] = -1;              /* boundary flag */
    }
    else if (jAdj[jEnd[n2]] >= 0) {
      jEnd[n2] += 1;                    /* if we are not overwriting the */
      jAdj[jEnd[n2]] = n1;              /* boundary flag, move up by 1 */
    }
    else
      jAdj[jEnd[n2]] = n1;              /* overwrite the boundary flag */
  }
}

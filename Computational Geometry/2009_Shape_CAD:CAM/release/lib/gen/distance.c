/* Copyright (C) 1997, Massachusetts Institute of Technology
 * All rights reserved
 */

/* linearDistance()
 * pointToLine()
 */

#include <math.h>
#include "gen.h"

/* return the minimum distance from a point (c) to a line segment (ab)
 * this is either the perpindicular distance to the interior of ab
 * or an oblique distance to either a or b
 * the parametric location of the minimum distance projection is 
 * returned in t
 */

double pointToLine(vector *a, vector*b, vector *c, double *t)
{
  double d, dx, dy, dz, rx, ry, rz;

  dx = b->x - a->x;             /* line segment from a to b */
  dy = b->y - a->y;
  dz = b->z - a->z;

  rx = c->x - a->x;             /* vector from point c to a */
  ry = c->y - a->y;
  rz = c->z - a->z;

  *t = (rx*dx + ry*dy + rz*dz) /   /* parametric location of orthogonal */
      (dx*dx + dy*dy + dz*dz);     /* projection of c onto ab */

  if (*t <= 0.0) {
    *t = 0.0;
    d = sqrt(rx*rx + ry*ry + rz*rz);   /* vector from point to end of line */
  }
  else {
    if (*t >= 1.0) {
      *t = 1.0;
      dx = c->x - b->x;         /* vector from point to end of line */
      dy = c->y - b->y;
      dz = c->z - c->z;
    }
    else {
      rx = a->x + (*t)*dx;      /* cartesian location of orthogonal */
      ry = a->y + (*t)*dy;      /* projection of c onto ab */
      rz = a->z + (*t)*dz;
      dx = c->x - rx;           /* vector from point to line */
      dy = c->y - ry;
      dz = c->z - rz;
    }
    d = sqrt(dx*dx + dy*dy + dz*dz);
  }
  return d;
}

/* return the minimum distance from a point (c) to a linear B-spline
 * curve (egeom)
 * this should be the perpindicular distance to the interior of 
 * the linear span between a pair of control points, or the 
 * perpindicular distance to a control point
 * the parametric location of the minimum distance projection is
 * return in t
 */

double linearDistance(ParCurv *egeom, vector *c, double *t)
{
  double d, dmin, u;
  int i;

  dmin = 1.0e+30;
  for (i=1; i<egeom->ncontpts; i++) {

    /* find distance from point c to the line segment between the
     * the (i-1)th and ith control points */

    d = pointToLine(egeom->contpts[i-1], egeom->contpts[i], c, &u);

    if (d < dmin) {
      /* if this distance is smaller than the current distance, save it */
      dmin = d;

      /* the parametric value of the closest point is the knot value of
       * the first control point, plus the proportion of the knot span 
       */
      *t = egeom->knots[i] + u*(egeom->knots[i+1] - egeom->knots[i]);
    }
  }
  return dmin;
}

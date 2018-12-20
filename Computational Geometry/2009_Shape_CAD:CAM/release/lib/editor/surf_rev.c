/************************************************************************
 *									*
			    Copyright (C) 1992
     Massachusetts Institute of Technology, Cambridge, Massachusetts.
			   All rights reserved.

     This portion of source code was prepared in the Ocean Engineering
     Design Laboratory by Bradley A. Moran.

     bug fix 1-4-1991; changed 0 to nrb->vcontpts-1 on lines 432-437
     (change made by E. C. Sherbrooke)

 *									*
 ************************************************************************/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <malloc.h>
#include "gen.h"
#include "bspl.h"
#include "editor.h"

/* only use single precision to give a small tolerance */
#define M_2_PI          6.2831853
#define M_PI		3.14159265358979323846
#define M_PI_2		1.57079632679489661923
#define M_SQRT2		1.41421356237309504880
#define M_SQRT1_2	0.70710678118654752440

void qsort(void *, size_t, size_t, int (*)(double *, double *));

static double mat[4][4];

#define VORDER   3
#define VCONTPTS 12

ParSurf *SurfRev_to_ParSurf(SurfRev *rev, ParSurf *fgeom)
{
  double ang, v_split;
  ParCurv *prof;
  ParSurf *keeper, *trash;
  unsigned i, j, k;

  prof = copyegeom(rev->de, 0);
/*
 * Allocate the first nurbs surface, which represents the full sweep of the
 * profile geometry.  The v direction is the rational parameter, in the
 * circular sweep.
*/
  if (fgeom)
    keeper = fgeom;
  else
    keeper = fgeomalloc1(prof->order, 3, prof->ncontpts, 9);
  keeper->type = PSurfaceOpen;

  keeper->vknots[0] = keeper->vknots[1] = keeper->vknots[2] = 0.0;
  keeper->vknots[3] = keeper->vknots[4] = 0.25;
  keeper->vknots[5] = keeper->vknots[6] = 0.50;
  keeper->vknots[7] = keeper->vknots[8] = 0.75;
  keeper->vknots[9] = keeper->vknots[10] = keeper->vknots[11] = 1.0;

  memcpy(keeper->uknots, prof->knots,
	 (prof->order + prof->ncontpts)*sizeof(double));
/*
 * Compute the homogeneous transformation matrix that rotates the profile
 * curve about the srf_rev axis through the start angle.  Perform the
 * actual transformation and store the resulting curve in the isoparametric
 * lines v = 0 and v = 1 of the full sweep surface.
 */
  init_tmat(rev->axis, rev->start, 1.0);
  for (i = 0; i < prof->ncontpts; ++i) {
    post4x4(prof->contpts[i], keeper->contpts[i][0]);
    copyvector(keeper->contpts[i][0], keeper->contpts[i][8]);
  }
/*
 * Now do the same type of transformation for the remaining principal
 * orientations, i.e., sa + pi/2, sa + pi and sa + 3pi/2;
 */
  init_tmat(rev->axis, rev->start + M_PI/2, 1.0);
  for (i = 0; i < prof->ncontpts; ++i)
    post4x4(prof->contpts[i], keeper->contpts[i][2]);
  init_tmat(rev->axis, rev->start + M_PI, 1.0);
  for (i = 0; i < prof->ncontpts; ++i)
    post4x4(prof->contpts[i], keeper->contpts[i][4]);
  init_tmat(rev->axis, rev->start + 3*M_PI/2, 1.0);
  for (i = 0; i < prof->ncontpts; ++i)
    post4x4(prof->contpts[i], keeper->contpts[i][6]);
/*
 * The intermediate sections get a similar transformation, but the
 * homogeneous coordinates must be appropriately scaled for circular
 * geometry.
 */
  init_tmat(rev->axis, rev->start + M_PI/4, M_SQRT2);
  for (i = 0; i < prof->ncontpts; ++i) {
    post4x4(prof->contpts[i], keeper->contpts[i][1]);
    sethomogeq(keeper->contpts[i][1], M_SQRT1_2);
  }
  init_tmat(rev->axis, rev->start + 3*M_PI/4, M_SQRT2);
  for (i = 0; i < prof->ncontpts; ++i) {
    post4x4(prof->contpts[i], keeper->contpts[i][3]);
    sethomogeq(keeper->contpts[i][3], M_SQRT1_2);
  }
  init_tmat(rev->axis, rev->start + 5*M_PI/4, M_SQRT2);
  for (i = 0; i < prof->ncontpts; ++i) {
    post4x4(prof->contpts[i], keeper->contpts[i][5]);
    sethomogeq(keeper->contpts[i][5], M_SQRT1_2);
  }
  init_tmat(rev->axis, rev->start + 7*M_PI/4, M_SQRT2);
  for (i = 0; i < prof->ncontpts; ++i) {
    post4x4(prof->contpts[i], keeper->contpts[i][7]);
    sethomogeq(keeper->contpts[i][7], M_SQRT1_2);
  }
  trash = fgeomalloc1(keeper->uorder, VORDER, keeper->ucontpts, VCONTPTS);

/*
 * Compute the value of the v split and create the knot vectors for the new
 * surface.
 */
  v_split = split_val((ang = rev->term - rev->start) > M_2_PI ? M_2_PI : ang);
  memcpy(trash->uknots, keeper->uknots,
	 (keeper->uorder+keeper->ucontpts)*sizeof(double));
  memcpy(trash->vknots, keeper->vknots,	 VCONTPTS*sizeof(double));
  for (i = 0; i < VORDER; ++i)
    trash->vknots[VCONTPTS+i] = v_split;
  qsort(trash->vknots, VORDER+VCONTPTS, sizeof(double), surfrev_compare);
  
/*
 * Subdivide the surface.
 */
  surfoslo3(keeper, trash, 4);

  free_fgeom(keeper);
/*
 * Count the number of knots that are not greater than v_split.
 */
  k = 0;
  while (k < VORDER+VCONTPTS && trash->vknots[k] <= v_split)
    ++k;

  keeper = fgeomalloc1(trash->uorder, VORDER, trash->ucontpts, k - VORDER);
  keeper->type = PSurfaceOpen;
  memcpy(keeper->uknots, trash->uknots,
	 (trash->uorder+trash->ucontpts)*sizeof(double));
  memcpy(keeper->vknots, trash->vknots, k*sizeof(double));
  for (i = 0; i < trash->ucontpts; ++i)
    for (j = 0; j < k - VORDER; ++j)
      copyvector(trash->contpts[i][j], keeper->contpts[i][j]);
  free_fgeom(trash);
  
  knot_normalize(keeper->uknots, keeper->uorder + keeper->ucontpts);
  knot_normalize(keeper->vknots, keeper->vorder + keeper->vcontpts);
  free_egeom(prof);

  return (keeper);
}

SurfRev *ParSurf_to_SurfRev(ParSurf *nrb, SurfRev *sgeom)
{
  SurfRev *rev;
  ParCurv *gen;
  register unsigned i;

  if (sgeom)
    rev = sgeom;
  else
    rev = (SurfRev *)gen_array1(1, sizeof(SurfRev));
/*
 * If the surface of revolution is successfully allocated, build a
 * generatrix curve as a nurbs by extracting from the nurbs surface the
 * isoparametric curve v = 0.
 */
  gen = egeomalloc1(nrb->uorder, nrb->ucontpts);
  gen->type = PCurveOpen;

  memcpy(gen->knots, nrb->uknots, nrb->ukmem*sizeof(double));
  for (i = 0; i < nrb->ucontpts; ++i)
    memcpy(gen->contpts[i], nrb->contpts[i][0], sizeof(vector));
/*
 * If type specifies that the generatrix be a power basis spline, IGES type
 * 112, convert the nurbs and free memory.  Otherwise, gen is the
 * generatrix.
 */
  rev->de = gen;
/*
 * Compute the axis, start and terminate angles.  The start angle will
 * arbitrarily assigned at sa = 0.  The terminate angle must be computed
 * from the geometry
 */
  rev->start = 0.0;
  rev->term = axis_and_terminate(nrb, rev->axis);

  return (rev);
}

void init_tmat(double *axis, double start, double scale)
{
  double dx, dy, dz, mg, th, bt;
  double arr[23], *v = arr-1;

  dx = axis[3] - axis[0];
  dy = axis[4] - axis[1];
  dz = axis[5] - axis[2];
  mg = sqrt(dx*dx + dy*dy + dz*dz);

  th = atan2(dy, dx);

  if (!(th == th))	/* check for NaN, set to 0.0 if NaN */
    th = 0.0;

  bt = asin(dz / mg);

  v[1] = sin(th);
  v[2] = sin(bt);
  v[3] = sin(start);
  v[4] = cos(th);
  v[5] = cos(start);
  v[6] = -scale*(v[2]*v[3]*v[4] + v[5]*v[1]);
  v[7] = cos(bt);
  v[8] = v[7]*v[7];
  v[9] = scale*(v[3]*v[1] - v[2]*v[5]*v[4]);
  v[10] = v[8]*v[4] - v[2]*v[9];
  v[11] = v[4]*v[10] - v[1]*v[6];
  v[12] = scale*(v[5]*v[4] - v[2]*v[3]*v[1]);
  v[13] = -scale*(v[2]*v[5]*v[1] + v[3]*v[4]);
  v[14] = v[8]*v[1] - v[2]*v[13];
  v[15] = v[4]*v[14] - v[1]*v[12];
  v[16] = v[7]*v[2]*(1 - v[5]*scale);
  v[17] = v[16]*v[4] - v[7]*v[3]*scale*v[1];
  v[18] = v[1]*v[10] + v[4]*v[6];
  v[19] = v[1]*v[14] + v[4]*v[12];
  v[20] = v[16]*v[1] + v[7]*v[3]*scale*v[4];
  v[21] = v[7]*(v[9] + v[2]*v[4]);
  v[22] = v[7]*(v[13] + v[2]*v[1]);
  v[23] = v[8]*v[5]*scale + v[2]*v[2];

  mat[0][0] = v[11];
  mat[0][1] = v[15];
  mat[0][2] = v[17];
  mat[0][3] = 0.0;
  
  mat[1][0] = v[18];
  mat[1][1] = v[19];
  mat[1][2] = v[20];
  mat[1][3] = 0.0;
  
  mat[2][0] = v[21];
  mat[2][1] = v[22];
  mat[2][2] = v[23];
  mat[2][3] = 0.0;
  
  mat[3][0] = -v[21]*axis[2] - v[18]*axis[1] - v[11]*axis[0] + axis[0];
  mat[3][1] = -v[22]*axis[2] - v[19]*axis[1] + axis[1] - v[15]*axis[0];
  mat[3][2] = -v[23]*axis[2] + axis[2] - v[20]*axis[1] - v[17]*axis[0];
  mat[3][3] = 1.0;
}

void post4x4(vector *in, vector *out)
{
  vector tmp;
  
  tmp.x = mat[0][0]*in->x + mat[1][0]*in->y + mat[2][0]*in->z +
    mat[3][0]*in->w;
  tmp.y = mat[0][1]*in->x + mat[1][1]*in->y + mat[2][1]*in->z +
    mat[3][1]*in->w;
  tmp.z = mat[0][2]*in->x + mat[1][2]*in->y + mat[2][2]*in->z +
    mat[3][2]*in->w;
  tmp.w = mat[0][3]*in->x + mat[1][3]*in->y + mat[2][3]*in->z +
    mat[3][3]*in->w;
  
  memcpy(out, &tmp, sizeof(vector));
}

double split_val(double delta)
{
  double s, val;
  int k = 0;
  
  while (delta > M_PI_2) {
    delta -= M_PI_2;
    ++k;
  }
  s = tan(delta / 2);
  val = s*M_SQRT2 / (1 + s*(M_SQRT2 - 1));
  val *= 0.25;
  val += 0.25*k;
  return (val);
}

int surfrev_compare(double *a, double *b)
{
  if (*a < *b)
    return -1;
  else if (*a > *b)
    return 1;
  else
    return 0;
}

double axis_and_terminate(ParSurf *nrb, double *axis)
{
  double alpha, denom;
  int internal_knots;
  vector a0, b0, da, db, *eval, v0, v1, v2, a, b, n;

/*
 * Start by selecting three points on the circle.  The number of internal
 * knots is used to determine how many quadrants the arc spans.  Ideally,
 * the three points selected will form a well-behaved triangle.
 */
  memcpy(&v0, nrb->contpts[0][0], sizeof(vector));

  internal_knots = nrb->vcontpts - nrb->vorder;
  switch (internal_knots) {

  case 0:
  case 2:
    evalrsurf(nrb, 0.0, 0.5, 0, &eval);
    memcpy(&v1, eval, sizeof(vector));
    memcpy(&v2, nrb->contpts[0][nrb->vcontpts-1], sizeof(vector));
    break;
    
  case 4:
    memcpy(&v1, nrb->contpts[0][2], sizeof(vector));
    memcpy(&v2, nrb->contpts[0][nrb->vcontpts-1], sizeof(vector));
    break;
    
  case 6:
    memcpy(&v1, nrb->contpts[0][2], sizeof(vector));
    memcpy(&v2, nrb->contpts[0][4], sizeof(vector));
    break;
    
  default:
    fputs("geometry error in compute_axis\n", stderr);
  }
/*
 * Using the three points, determine a normal vector to the plane formed.
 * The orientation of the normal is chosen to point upward.
 */
  sub_vect1(&v2, &v1, &db);
  sub_vect1(&v0, &v1, &da);
  cross1(&db, &da, &n);
  unitvector1(&n, &n);

/*
 * Compute the midpoints of two of the sides of the triangle.
 */
  add_vect1(&v0, &v1, &a);
  scale4(0.5, &a);
  add_vect1(&v1, &v2, &b);
  scale4(0.5, &b);

/*
 * The vectors a0 and b0 are perpendicular bisectors of the sides of the
 * triangle.  Use the normal to get their directions.
 */
  cross1(&da, &n, &a0);
  cross1(&n, &db, &b0);
/*
 * Compute the point O by intersecting the lines given by a+alpha*a0 and
 * b+beta*b0.
 */
  denom = a0.y*b0.x - a0.x*b0.y;
  if (fabs(denom) > 1.0e-06)
    alpha = (b0.x*(b.y - a.y) - b0.y*b.x + a.x*b0.y) / denom;
  else
    alpha = (b0.x*(b.z - a.z) - b0.z*b.x + a.x*b0.z) /
      (a0.z*b0.x - a0.x*b0.z);

  axis[0] = a.x + alpha*a0.x;
  axis[1] = a.y + alpha*a0.y;
  axis[2] = a.z + alpha*a0.z;
  axis[3] = axis[0] + n.x;
  axis[4] = axis[1] + n.y;
  axis[5] = axis[2] + n.z;
/*
 * Determine the terminate angle by taking the inverse cosine of the dot
 * product of vectors which correspond to R(0,0) - O and R(0,1) - O.  Check
 * the quadrant of the solution.
 */
  da.x = nrb->contpts[0][0]->x / nrb->contpts[0][0]->w - axis[0];
  da.y = nrb->contpts[0][0]->y / nrb->contpts[0][0]->w - axis[1];
  da.z = nrb->contpts[0][0]->z / nrb->contpts[0][0]->w - axis[2];
  unitvector1(&da, &da);
  db.x = nrb->contpts[0][nrb->vcontpts-1]->x /
    nrb->contpts[0][nrb->vcontpts-1]->w - axis[0];
  db.y = nrb->contpts[0][nrb->vcontpts-1]->y /
    nrb->contpts[0][nrb->vcontpts-1]->w - axis[1];
  db.z = nrb->contpts[0][nrb->vcontpts-1]->z /
    nrb->contpts[0][nrb->vcontpts-1]->w - axis[2];
  unitvector1(&db, &db);
  
  alpha = acos(dot(&da, &db));

  if (internal_knots == 4 || internal_knots == 6)
    alpha = 2*M_PI - alpha;

  return (alpha);
}

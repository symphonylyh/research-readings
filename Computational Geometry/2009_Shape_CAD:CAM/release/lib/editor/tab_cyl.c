/************************************************************************
 *									*
			    Copyright (C) 1992
     Massachusetts Institute of Technology, Cambridge, Massachusetts.
			   All rights reserved.

     This portion of source code was prepared in the Ocean Engineering
     Design Laboratory by Bradley A. Moran.
 *									*
 ************************************************************************/

#include <stdio.h>
#include <string.h>
#include <malloc.h>
#include "gen.h"
#include "editor.h"

ParSurf *TabCyl_to_ParSurf(TabCyl *cyl, ParSurf *fgeom)
{
  double dx, dy, dz, w;
  unsigned i;
  ParCurv *dir;
  ParSurf *nrb;

  dir = cyl->de;
/*
 * Allocate the nurbs surface, linear in v.  Set the type to aperiodic.
 */
  if (fgeom)
    nrb = fgeom;
  else
    nrb = fgeomalloc1(dir->order, 2, dir->ncontpts, 2);
  nrb->type = PSurfaceOpen;

/*
 * The u knot vector of the nurbs will be the same as the knot vector of
 * the directrix.  The v knot vector will be 0,0,1,1.
 */

  memcpy(nrb->uknots, dir->knots, dir->kmem*sizeof(double));
  nrb->vknots[0] = nrb->vknots[1] = 0.0;
  nrb->vknots[2] = nrb->vknots[3] = 1.0;

/*
 * These three delta values represent the translation required to shift the
 * directrix to the other end of the line.
 */

  dx = cyl->lx - dir->contpts[0]->x / dir->contpts[0]->w;
  dy = cyl->ly - dir->contpts[0]->y / dir->contpts[0]->w;
  dz = cyl->lz - dir->contpts[0]->z / dir->contpts[0]->w;
     
/*
 * Compute the control points:  the points corresponding to v = 0 come
 * directly from the directrix; the ones for v = 1 are simply a translation
 * of the former.
 */

  for (i = 0; i < nrb->ucontpts; ++i) {
    memcpy(nrb->contpts[i][0], dir->contpts[i], sizeof(vector));
    memcpy(nrb->contpts[i][1], dir->contpts[i], sizeof(vector));
    nrb->contpts[i][1]->x += dx*nrb->contpts[i][1]->w;
    nrb->contpts[i][1]->y += dy*nrb->contpts[i][1]->w;
    nrb->contpts[i][1]->z += dz*nrb->contpts[i][1]->w;
  }

/*
 * The final step is to determine whether or not to free the directrix.
 * It will be freed only if its pointer value is different from the pointer
 * value of the original directrix, which signifies that a conversion took
 * place.
 */

  if (dir != (ParCurv *) cyl->de)
    free_egeom (dir);

  return nrb;
}

TabCyl *ParSurf_to_TabCyl(ParSurf *nrb, TabCyl *tgeom)
{
  TabCyl *cyl;
  ParCurv *dir;
  register unsigned i;

  if (tgeom)
    cyl = tgeom;
  else
    cyl = (TabCyl *)gen_array1(1, sizeof(TabCyl));

/*
 * If cyl is successfully allocated, build a nurbs directrix curve from the
 * isoparametric line v = 0.
 */
  dir = egeomalloc1(nrb->uorder, nrb->ucontpts);
  dir->type = PCurveOpen;

  memcpy(dir->knots, nrb->uknots, nrb->ukmem*sizeof(double));
  for (i = 0; i < nrb->ucontpts; ++i)
    memcpy(dir->contpts[i], nrb->contpts[i][0], sizeof(vector));

/*
 * If type specifies that the directrix be a power basis spline, type 112,
 * convert the nurbs and free memory.  Otherwise, dir is the directrix.
 */
  
  cyl->de = dir;

/*
 * The terminate point of the line segment is the control point of the
 * original surface corresponding to u = 0, v = 1.
 */
  cyl->lx = nrb->contpts[0][1]->x / nrb->contpts[0][1]->w;
  cyl->ly = nrb->contpts[0][1]->y / nrb->contpts[0][1]->w;
  cyl->lz = nrb->contpts[0][1]->z / nrb->contpts[0][1]->w;

  return cyl;
}

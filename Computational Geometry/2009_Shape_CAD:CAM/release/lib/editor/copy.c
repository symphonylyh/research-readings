/* Copyright (C) 1992 Massachusetts Institute of Technology, Cambridge, MA
   All rights reserved
*/

/* copy.c */

/* CopyFacetSurf()
   CopyGridSurf()
   CopyListCurv()
   CopyParUv()
   FreeFacetSurf()
*/

#include <string.h>
#include <malloc.h>
#include "gen.h"
#include "editor.h"

ListCurv *CopyListCurv(ListCurv *from, ListCurv *to)
{
  int i, j;

  if (!to)
    to = (ListCurv *)gen_array1(1, sizeof(ListCurv));
  memcpy(to, from, sizeof(ListCurv));

  to->pts = dbl_array2(to->npoints, 3);
  for (i=0; i<to->npoints; i++)
    for (j=0; j<3; j++)
      to->pts[i][j] = from->pts[i][j];

  return (to);
}

GridSurf *CopyGridSurf(GridSurf *from, GridSurf *to)
{
  int i, j, k;

  if (!to)
    to = (GridSurf *)gen_array1(1, sizeof(GridSurf));
  memcpy(to, from, sizeof(GridSurf));

  to->pts = dbl_array3(to->ncol, to->nrow, 3);
  for (i=0; i<to->ncol; i++)
    for (j=0; j<to->nrow; j++)
      for (k=0; k<3; k++)
	to->pts[i][j][k] = from->pts[i][j][k];

  return (to);
}

ParUv *CopyParUv(ParUv *from, ParUv *to)
{
  int i, j;

  if (!to)
    to = (ParUv *)gen_array1(1, sizeof(ParUv));
  memcpy(to, from, sizeof(ParUv));

  to->pts = dbl_array2(to->npoints, 2);
  for (i=0; i<to->npoints; i++)
    for (j=0; j<2; j++)
      to->pts[i][j] = from->pts[i][j];

  return (to);
}

FacetSurf *CopyFacetSurf(FacetSurf *from, FacetSurf *to)
{
  int i, j;

  if (!to)
    to = (FacetSurf *)gen_array1(1, sizeof(FacetSurf));
  memcpy(to, from, sizeof(FacetSurf));

  to->pts = dbl_array2(to->nPts, 3);
  for (i=0; i<to->nPts; i++)
    for (j=0; j<3; j++)
      to->pts[i][j] = from->pts[i][j];
  to->iEnd = int_array1(to->nPts);
  for (i=0; i<to->nPts; i++)
    to->iEnd[i] = from->iEnd[i];
  to->iAdj = int_array1(to->nAdj);
  for (i=0; i<to->nAdj; i++)
    to->iAdj[i] = from->iAdj[i];

  return (to);
}

void FreeFacetSurf(FacetSurf *tgeom)
{
  free_darray2(tgeom->pts);
  free_iarray1(tgeom->iEnd);
  free_iarray1(tgeom->iAdj);
  free_garray1((char *)tgeom);
}

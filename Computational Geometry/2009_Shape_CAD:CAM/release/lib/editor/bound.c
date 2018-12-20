/* Copyright (C) 1992 Massachusetts Institute of Technology, Cambridge, MA
   All rights reserved
*/

/* bound.c */

/* CompareBoundBox()
   GetCurvBoundBox()
   GetGridBoundBox()
   GetListBoundBox()
   GetMaxRange()
   GetSurfBoundBox()
   GetUvBoundBox()
   GetVectBoundBox()
*/

#include "gen.h"
#include "editor.h"

void GetListBoundBox(ListCurv *list, float *box)
{
  int i;
  
  box[0] = box[2] = box[4] = 1.0e30;
  box[1] = box[3] = box[5] = -1.0e30;
  for (i=0; i<list->npoints; i++) {
    if (list->pts[i][0] < box[0]) box[0] = list->pts[i][0];
    if (list->pts[i][0] > box[1]) box[1] = list->pts[i][0];
    if (list->pts[i][1] < box[2]) box[2] = list->pts[i][1];
    if (list->pts[i][1] > box[3]) box[3] = list->pts[i][1];
    if (list->pts[i][2] < box[4]) box[4] = list->pts[i][2];
    if (list->pts[i][2] > box[5]) box[5] = list->pts[i][2];
  }
}

void GetGridBoundBox(GridSurf *grid, float *box)
{
  int i, j;
  
  box[0] = box[2] = box[4] = 1.0e30;
  box[1] = box[3] = box[5] = -1.0e30;
  for (i=0; i<grid->ncol; i++)
    for (j=0; j<grid->nrow; j++) {
      if (grid->pts[i][j][0] < box[0]) box[0] = grid->pts[i][j][0];
      if (grid->pts[i][j][0] > box[1]) box[1] = grid->pts[i][j][0];
      if (grid->pts[i][j][1] < box[2]) box[2] = grid->pts[i][j][1];
      if (grid->pts[i][j][1] > box[3]) box[3] = grid->pts[i][j][1];
      if (grid->pts[i][j][2] < box[4]) box[4] = grid->pts[i][j][2];
      if (grid->pts[i][j][2] > box[5]) box[5] = grid->pts[i][j][2];
    }
}

void GetCurvBoundBox(ParCurv *egeom, float *box)
{
  int i;
  
  box[0] = box[2] = box[4] = 1.0e30;
  box[1] = box[3] = box[5] = -1.0e30;
  for (i=0; i<egeom->ncontpts; i++) {
    if (egeom->contpts[i]->x < box[0]) box[0] = egeom->contpts[i]->x;
    if (egeom->contpts[i]->x > box[1]) box[1] = egeom->contpts[i]->x;
    if (egeom->contpts[i]->y < box[2]) box[2] = egeom->contpts[i]->y;
    if (egeom->contpts[i]->y > box[3]) box[3] = egeom->contpts[i]->y;
    if (egeom->contpts[i]->z < box[4]) box[4] = egeom->contpts[i]->z;
    if (egeom->contpts[i]->z > box[5]) box[5] = egeom->contpts[i]->z;
  }
}

void GetSurfBoundBox(ParSurf *fgeom, float *box)
{
  int i, j;
  
  box[0] = box[2] = box[4] = 1.0e30;
  box[1] = box[3] = box[5] = -1.0e30;
  for (i=0; i<fgeom->ucontpts; i++)
    for (j=0; j<fgeom->vcontpts; j++) {
      if (fgeom->contpts[i][j]->x < box[0]) box[0] = fgeom->contpts[i][j]->x;
      if (fgeom->contpts[i][j]->x > box[1]) box[1] = fgeom->contpts[i][j]->x;
      if (fgeom->contpts[i][j]->y < box[2]) box[2] = fgeom->contpts[i][j]->y;
      if (fgeom->contpts[i][j]->y > box[3]) box[3] = fgeom->contpts[i][j]->y;
      if (fgeom->contpts[i][j]->z < box[4]) box[4] = fgeom->contpts[i][j]->z;
      if (fgeom->contpts[i][j]->z > box[5]) box[5] = fgeom->contpts[i][j]->z;
    }
}

void GetUvBoundBox(ParUv *ugeom, float *box)
{
  int i;
  
  box[0] = box[2] = 1.0e30;
  box[1] = box[3] = -1.0e30;
  box[4] = box[5] = 0.0;
  for (i=0; i<ugeom->npoints; i++) {
    if (ugeom->pts[i][0] < box[0]) box[0] = ugeom->pts[i][0];
    if (ugeom->pts[i][0] > box[1]) box[1] = ugeom->pts[i][0];
    if (ugeom->pts[i][1] < box[2]) box[2] = ugeom->pts[i][1];
    if (ugeom->pts[i][1] > box[3]) box[3] = ugeom->pts[i][1];
  }
}

void GetVectBoundBox(double **pts, float *box, int nPts)
{
  int i;
  
  box[0] = box[2] = box[4] = 1.0e30;
  box[1] = box[3] = box[5] = -1.0e30;
  for (i=0; i<nPts; i++) {
    if (pts[i][0] < box[0]) box[0] = pts[i][0];
    if (pts[i][0] > box[1]) box[1] = pts[i][0];
    if (pts[i][1] < box[2]) box[2] = pts[i][1];
    if (pts[i][1] > box[3]) box[3] = pts[i][1];
    if (pts[i][2] < box[4]) box[4] = pts[i][2];
    if (pts[i][2] > box[5]) box[5] = pts[i][2];
  }
}

void CompareBoundingBox(float *update, float *box)
{
  int i;

  for (i=0; i<6; i += 2)
    if (box[i] < update[i]) update[i] = box[i];
  for (i=1; i<6; i += 2)
    if (box[i] > update[i]) update[i] = box[i];
}

float GetMaxRange(float *box)
{
  float m = 0.0, r;

  if ((r = box[1] - box[0]) > m) m = r;
  if ((r = box[3] - box[2]) > m) m = r;
  if ((r = box[5] - box[4]) > m) m = r;

  return (m);
}

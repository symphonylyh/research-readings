/* Transformright (C) 1992 Massachusetts Institute of Technology, Cambridge, MA
   All rights reserved
*/

/* transform.c */

/* transformegeom()
   transformfgeom()
   TransformGridSurf()
   TransformListCurv()
   TransformPts()
   TransformTrimSurf()
*/

#include <string.h>
#include <malloc.h>
#include "gen.h"
#include "editor.h"

ListCurv *TransformListCurv(ListCurv *from, ListCurv *to, double *tr)
{
  vector *frv, *tov;
  double **m;
  int i;

  if (!to)
    to = (ListCurv *)gen_array1(1, sizeof(ListCurv));
  memcpy(to, from, sizeof(ListCurv));

  m = compute_transfm(tr);
  frv = vectalloc();
  tov = vectalloc();
  frv->w = 1.0;

  to->pts = dbl_array2(to->npoints, 3);
  for (i=0; i<to->npoints; i++) {
    frv->x = tr[6]*from->pts[i][0];   /* scaling */
    frv->y = tr[7]*from->pts[i][1];
    frv->z = tr[8]*from->pts[i][2];
    mult4x4(m, frv, tov);             /* rotations and translations */
    to->pts[i][0] = tov->x/tov->w;
    to->pts[i][1] = tov->y/tov->w;
    to->pts[i][2] = tov->z/tov->w;
  }
  vectfree(frv);
  vectfree(tov);
  free_darray2(m);

  return (to);
}

GridSurf *TransformGridSurf(GridSurf *from, GridSurf *to, double *tr)
{
  vector *frv, *tov;
  double **m;
  int i, j, k;

  if (!to)
    to = (GridSurf *)gen_array1(1, sizeof(GridSurf));
  memcpy(to, from, sizeof(GridSurf));

  m = compute_transfm(tr);
  frv = vectalloc();
  tov = vectalloc();
  frv->w = 1.0;

  to->pts = dbl_array3(to->ncol, to->nrow, 3);
  for (i=0; i<to->ncol; i++)
    for (j=0; j<to->nrow; j++) {
      frv->x = tr[6]*from->pts[i][j][0];   /* scaling */
      frv->y = tr[7]*from->pts[i][j][1];
      frv->z = tr[8]*from->pts[i][j][2];
      mult4x4(m, frv, tov);                /* rotations and translations */
      to->pts[i][j][0] = tov->x/tov->w;
      to->pts[i][j][1] = tov->y/tov->w;
      to->pts[i][j][2] = tov->z/tov->w;
    }
  vectfree(frv);
  vectfree(tov);
  free_darray2(m);

  return (to);
}

ParCurv *transformegeom(ParCurv *old, ParCurv *new, double *tr)
{
  vector *frv;
  double **m;
  int i, nbytes;

  if (!new) {
    new = (ParCurv *)gen_array1(1, sizeof(ParCurv));

    if (new) {
      memcpy(new, old, sizeof(ParCurv));
      new->knots = dbl_array1(old->kmem);
      new->contpts = vec_array1(old->pmem);
    }
  }
  else {
    if (new->kmem < old->kmem) {
      free_darray1(new->knots);
      new->knots = dbl_array1(old->kmem);
      new->kmem = old->kmem;
    }
    if (new->pmem < old->pmem) {
      free_varray1(new->contpts, new->ncontpts);
      new->contpts = vec_array1(old->pmem);
      new->pmem = old->pmem;
    }
    new->type = old->type;
    new->ncontpts = old->ncontpts;
    new->order = old->order;
  }
  nbytes = sizeof(double) * old->kmem;
  memcpy(new->knots, old->knots, nbytes);
  
  m = compute_transfm(tr);
  frv = vectalloc();
  for (i=0; i<old->ncontpts; i++) {
    frv->x = tr[6]*old->contpts[i]->x;   /* scaling */
    frv->y = tr[7]*old->contpts[i]->y;
    frv->z = tr[8]*old->contpts[i]->z;
    frv->w =       old->contpts[i]->w;
    mult4x4(m, frv, new->contpts[i]);    /* rotations and translations */
  }
  vectfree(frv);
  free_darray2(m);

  return new;
}

ParSurf *transformfgeom(ParSurf *old, ParSurf *new, double *tr)
{
  vector *frv;
  double **m;
  int i, j;

  if (!new) {
    new = (ParSurf *)gen_array1(1, sizeof(ParSurf));

    if (new) {
      memcpy(new, old, sizeof(ParSurf));
      new->uknots = dbl_array1(old->ukmem);
      new->vknots = dbl_array1(old->vkmem);
      new->contpts = vec_array2(old->upmem, old->vpmem);
    }
  }
  else {
    if (new->ukmem < old->ukmem) {
      free_darray1(new->uknots);
      new->uknots = dbl_array1(old->ukmem);
      new->ukmem = old->ukmem;
    }
    if (new->vkmem < old->vkmem) {
      free_darray1(new->vknots);
      new->vknots = dbl_array1(old->vkmem);
      new->vkmem = old->vkmem;
    }
    if (new->vpmem < old->vpmem || new->upmem < old->upmem) {
      free_varray2(new->contpts, new->ucontpts, new->vcontpts);
      new->contpts = (vector ***)ptr_array2(old->upmem, old->vpmem);
      new->upmem = old->upmem, new->vpmem = old->vpmem;
    }
    new->type = old->type;
    new->uorder = old->uorder, new->ucontpts = old->ucontpts;
    new->vorder = old->vorder, new->vcontpts = old->vcontpts;
  }
  memcpy(new->uknots, old->uknots, old->ukmem*sizeof(double));
  memcpy(new->vknots, old->vknots, old->vkmem*sizeof(double));

  m = compute_transfm(tr);
  frv = vectalloc();
  for (i=0; i<old->ucontpts; i++)
    for (j=0; j<old->vcontpts; j++) {
      frv->x = tr[6]*old->contpts[i][j]->x;   /* scaling */
      frv->y = tr[7]*old->contpts[i][j]->y;
      frv->z = tr[8]*old->contpts[i][j]->z;
      frv->w =       old->contpts[i][j]->w;
      mult4x4(m, frv, new->contpts[i][j]);    /* rotations and translations */ 
    }
  vectfree(frv);
  free_darray2(m);

  return new;
}

TrimSurf *TransformTrimSurf(TrimSurf *from, TrimSurf *to, double *tr)
{
  int i, j;

  if (!to)
    to = (TrimSurf *)gen_array1(1, sizeof(TrimSurf));

  to->fgeom = transformfgeom(from->fgeom, NULL, tr);
  to->bLoop = from->bLoop;
  to->nLoops = from->nLoops;
  to->loops = (TrimLoop *)gen_array1(to->nLoops, sizeof(TrimLoop));
  for (i=0; i<to->nLoops; i++) {
    to->loops[i].nCurves = from->loops[i].nCurves;
    to->loops[i].egeoms = (ParCurv **)gen_array1(to->loops[i].nCurves,
						 sizeof(ParCurv *));
    for (j=0; j<to->loops[i].nCurves; j++)
      to->loops[i].egeoms[j] = copyegeom(from->loops[i].egeoms[j], NULL);
  }

  /******************* 11.1 begin *********************/
  if(!from->c_loops) { /* No explicit C-curve info */
    to->c_loops = NULL;
    PrintLog("--> \"Original trimmed surface does not have C-curve info (explicitly).\"\n");
    PrintLog("    \"Its transform won't have it either.\"\n");
  }
  else { /* "from" has explicit C-curve info. */
    PrintLog("--> \"Original trimmed surface has C-curve info (explicitly).\"\n");
    PrintLog("    \"C-curve(s) should be transformed too.\"\n");

    to->c_loops = (TrimLoop *)gen_array1(to->nLoops, sizeof(TrimLoop));
    for (i=0; i<to->nLoops; i++) {
      to->c_loops[i].nCurves = from->c_loops[i].nCurves;
      to->c_loops[i].egeoms = (ParCurv **)gen_array1(to->c_loops[i].nCurves,
						     sizeof(ParCurv *));
      for (j=0; j<to->c_loops[i].nCurves; j++)
	to->c_loops[i].egeoms[j] = transformegeom(from->c_loops[i].egeoms[j], NULL, tr);
    }
  }
  /******************* 11.1 end *********************/

  return (to);
}

double **TransformPts(double **old, double **new, double *tr, int nPts)
{
  vector n, o;
  double **m;
  int i;

  if (!new)
    new = dbl_array2(nPts, 3);
  
  m = compute_transfm(tr);
  for (i=0; i<nPts; i++) {
    o.x = tr[6]*old[i][0];   /* scaling */
    o.y = tr[7]*old[i][1];
    o.z = tr[8]*old[i][2];
    o.w = 1.0;               /* rotations and translations */
    mult4x4(m, &o, &n);
    new[i][0] = n.x/n.w;
    new[i][1] = n.y/n.w;
    new[i][2] = n.z/n.w;
  }
  free_darray2(m);

  return new;
}

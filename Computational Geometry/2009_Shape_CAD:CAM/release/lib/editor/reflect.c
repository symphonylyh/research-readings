/* Transformright (C) 1995 Massachusetts Institute of Technology, Cambridge, MA
   All rights reserved
*/

/* reflect.c */

/* reflectegeom()
   reflectfgeom()
   ReflectGridSurf()
   ReflectListCurv()
   ReflectPts()
   ReflectTrimSurf()
*/

#include <string.h>
#include <malloc.h>
#include "gen.h"
#include "editor.h"

ListCurv *ReflectListCurv(ListCurv *from, ListCurv *to, int plane)
{
  short i;

  if (!to)
    to = (ListCurv *)gen_array1(1, sizeof(ListCurv));
  memcpy(to, from, sizeof(ListCurv));

  for (i=0; i<to->npoints; i++) {
    if (plane == 1)   /* reflect about yz-plane */
      to->pts[i][0] = -from->pts[i][0];

    if (plane == 2)   /* reflect about xz-plane */
      to->pts[i][1] = -from->pts[i][1];

    if (plane == 0)   /* reflect about xy-plane */
      to->pts[i][2] = -from->pts[i][2];
  }

  return (to);
}

GridSurf *ReflectGridSurf(GridSurf *from, GridSurf *to, int plane)
{
  short i, j;

  if (!to)
    to = (GridSurf *)gen_array1(1, sizeof(GridSurf));
  memcpy(to, from, sizeof(GridSurf));

  for (i=0; i<to->ncol; i++)
    for (j=0; j<to->nrow; j++) {
      if (plane == 1) /* reflect about yz-plane */
	to->pts[i][j][0] = -from->pts[i][j][0];
      
      if (plane == 2) /* reflect about xz-plane */
	to->pts[i][j][1] = -from->pts[i][j][1];

      if (plane == 0) /* reflect about xy-plane */
	to->pts[i][j][2] = -from->pts[i][j][2];
    }

  return (to);
}

ParCurv *reflectegeom(ParCurv *old, ParCurv *new, int plane)
{
  short i, nbytes;

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
      new->contpts = (vector **)ptr_array1(old->pmem);
      new->pmem = old->pmem;
    }
    new->type = old->type;
    new->ncontpts = old->ncontpts;
    new->order = old->order;
  }
  nbytes = sizeof(double) * old->kmem;
  memcpy(new->knots, old->knots, nbytes);
  
  for (i=0; i<old->ncontpts; i++) {
    if (plane == 1) {   /* reflect about yz-plane */
      new->contpts[i]->x = -old->contpts[i]->x;
      new->contpts[i]->y =  old->contpts[i]->y;
      new->contpts[i]->z =  old->contpts[i]->z;
    }

    if (plane == 2) {   /* reflect about xz-plane */
      new->contpts[i]->x =  old->contpts[i]->x;
      new->contpts[i]->y = -old->contpts[i]->y;
      new->contpts[i]->z =  old->contpts[i]->z;
    }

    if (plane == 0) {   /* reflect about xy-plane */
      new->contpts[i]->x =  old->contpts[i]->x;
      new->contpts[i]->y =  old->contpts[i]->y;
      new->contpts[i]->z = -old->contpts[i]->z;
    }
  }

  return new;
}

ParSurf *reflectfgeom(ParSurf *old, ParSurf *new, int plane)
{
  short i, j, nbytes;

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

  for (i=0; i<old->ucontpts; i++)
    for (j=0; j<old->vcontpts; j++) {
      if (plane == 1) {   /* reflect about yz-plane */
	new->contpts[i][j]->x = -old->contpts[i][j]->x;
	new->contpts[i][j]->y =  old->contpts[i][j]->y;
	new->contpts[i][j]->z =  old->contpts[i][j]->z;
      }

      if (plane == 2) {   /* reflect about yz-plane */
	new->contpts[i][j]->x =  old->contpts[i][j]->x;
	new->contpts[i][j]->y = -old->contpts[i][j]->y;
	new->contpts[i][j]->z =  old->contpts[i][j]->z;
      }

      if (plane == 0) {  /* reflect about xy-plane */
	new->contpts[i][j]->x =  old->contpts[i][j]->x;
	new->contpts[i][j]->y =  old->contpts[i][j]->y;
	new->contpts[i][j]->z = -old->contpts[i][j]->z; /* plane */
      }
    }

  return new;
}

TrimSurf *ReflectTrimSurf(TrimSurf *from, TrimSurf *to, int plane)
{
  short i, j;

  if (!to)
    to = (TrimSurf *)gen_array1(1, sizeof(TrimSurf));

  to->fgeom = reflectfgeom(from->fgeom, NULL, plane);
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
    PrintLog("    \"Its reflection won't have it either.\"\n");
  }
  else { /* "from" has explicit C-curve info. */
    PrintLog("--> \"Original trimmed surface has C-curve info (explicitly).\"\n");
    PrintLog("    \"C-curve(s) should be reflected too.\"\n");

    to->c_loops = (TrimLoop *)gen_array1(to->nLoops, sizeof(TrimLoop));
    for (i=0; i<to->nLoops; i++) {
      to->c_loops[i].nCurves = from->c_loops[i].nCurves;
      to->c_loops[i].egeoms = (ParCurv **)gen_array1(to->c_loops[i].nCurves,
						     sizeof(ParCurv *));
      for (j=0; j<to->c_loops[i].nCurves; j++)
	to->c_loops[i].egeoms[j] = reflectegeom(from->c_loops[i].egeoms[j], NULL, plane);
    }
  }
  /******************* 11.1 end *********************/

  return (to);
}

double **ReflectPts(double **old, double **new, int plane, int nPts)
{
  vector n, o;
  short i;

  if (!new)
    new = dbl_array2(nPts, 3);
  memcpy(new, old, nPts*3*sizeof(double));

  for (i=0; i<nPts; i++) {
    if (plane == 1)   /* reflect about yz-plane */
      new[i][0] = -old[i][0];

    if (plane == 2)   /* reflect about xz-plane */
      new[i][1] = -old[i][1];

    if (plane == 0)   /* reflect about xy-plane */
      new[i][2] = -old[i][2];
  }

  return new;
}

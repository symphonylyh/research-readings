/* Copyright (C) 1993 Massachusetts Institute of Technology, Cambridge, MA
   All rights reserved
*/

/* knot_rem.c */

/* extract_knots()
   extract_knots_cv()
   get_cv_removal()
   get_iso_removal()
   knot_rem_and_pert()
   knot_rem_cv()
   rem_1st_knot()
*/

#include "gen.h"
#include "editor.h"

void knot_rem_cv(ParCurv **egeoms, int np, ParCurv *bspls)
{
  ParCurv *curv;
  int *uc;
  short i;
  
  uc = int_array1(np);

  for (i=0; i<np; i++)
    uc[i] = egeoms[i]->ncontpts;

  curv = egeomalloc1(bspls->order, bspls->ncontpts);

  extract_knots_cv(bspls, np, curv, uc);
  get_cv_removal(bspls, np, curv, uc);
  bspls = copyegeom(curv, bspls);

  free_egeom(curv);
  free_iarray1(uc);
}

void get_cv_removal(ParCurv *bspls, int np, ParCurv *curv, int *uc)
{
  short j, k, n, m;

  for (j=m=n=0; j< np; j++,m++)
    for (k=0; k<uc[j]-1; k++,m++,n++)
      curv->contpts[n] = copyvector(bspls->contpts[m], curv->contpts[n]);

  curv->contpts[n] = copyvector(bspls->contpts[m-1], curv->contpts[n]);
}

void extract_knots_cv(ParCurv *bspls, int np, ParCurv *curv, int *uc)
{
  short i, j, k, m, n;

  n = uc[0]+bspls->order-1;
  for (i=k=m=0; i<np; i++) {
    for (j=m; j<n; j++,k++)
      curv->knots[k] = bspls->knots[j];
    m = n+1; 
    n = n+uc[i+1];
  }
  curv->knots[k] = 1.0;
  curv->ncontpts = k - curv->order + 1;
}

double knot_rem_and_pert(ParSurf ***fgeom, int nps, int npt, ParSurf *bspls,
			 double PERT, int join_key)
{
  ParSurf *surf;
  double err;
  int *uc, *vc;
  short i, j, k;

  uc = int_array1(nps);
  vc = int_array1(npt);

  for(i=0; i<nps; i++)
    uc[i] = fgeom[0][i]->ucontpts;
  for(i=0; i<npt; i++)
    vc[i] = fgeom[i][0]->vcontpts;

  surf = fgeomalloc1(bspls->uorder, bspls->vorder,
		     bspls->ucontpts, bspls->vcontpts);

  extract_knots(surf, nps, npt, uc, vc, bspls);
  rem_1st_knot(surf, nps, npt, uc, vc, bspls);
  bspls = copyfgeom(surf, bspls);

  err = get_error(fgeom, nps, npt, bspls, 0, 1);

  if (join_key > 1)
    perturb_knots(bspls, nps, npt, uc, vc, PERT);

  free_iarray1(uc);
  free_iarray1(vc);

  return err;
}

void extract_knots(ParSurf *surf, int nps, int npt, int *uc, int *vc,
		   ParSurf *bspls)
{
  short i, j, k, m, n;

  n = uc[0]+bspls->uorder-1;
  for (i=k=m=0; i<nps; i++) {
    for (j=m; j<n; j++,k++)
      surf->uknots[k] = bspls->uknots[j];
    m = n+1; 
    n = n+uc[i+1];
  }
  surf->uknots[k] = 1.0;
  surf->ucontpts = k - surf->uorder + 1;

  n = vc[0]+bspls->vorder-1;
  for(i=k=m=0; i<npt; i++) {
    for(j=m; j<n; j++,k++)
      surf->vknots[k] = bspls->vknots[j];
    m = n+1; 
    n = n+vc[i+1];
  }
  surf->vknots[k] = 1.0;
  surf->vcontpts = k - surf->vorder + 1;
}

void rem_1st_knot(ParSurf *surf, int nps, int npt, int *uc, int *vc,
		  ParSurf *bspls)
{
  short i, j, k, m, n;
 
  for (j=m=n=0; j< npt; j++,m++)
    for (k=0; k<vc[j]-1; k++,m++,n++)
      get_iso_removal(n, m, surf, bspls, nps, uc);

  get_iso_removal(n, m-1, surf, bspls, nps, uc);
}

void get_iso_removal(int i1, int i2, ParSurf *surf, ParSurf *bspls, 
		     int nps, int *uc)
{
  short j, k, m, n;

  for (j=m=n=0; j< nps; j++,m++)
    for (k=0; k<uc[j]-1; k++,m++,n++)
      surf->contpts[n][i1] = 
	copyvector(bspls->contpts[m][i2], surf->contpts[n][i1]);
  surf->contpts[n][i1] =
    copyvector(bspls->contpts[m-1][i2], surf->contpts[n][i1]);
}

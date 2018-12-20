/* Copyright (C) 1993 Massachusetts Institute of Technology, Cambridge, MA
   All rights reserved
*/

/* get_error.c */

/* extract_cv()
   extract_patch()
   get_error()
   get_error_cv()
   recover_cv()
   recover_patch()
*/

#include "gen.h"
#include "bspl.h"
#include "editor.h"

double get_error_cv(ParCurv **egeom, int np, ParCurv *bspl, int PERTUR,
		    int NMK)
{
  ParCurv *negeom, *neg;
  double errors, FINAL_ERROR=0.0;
  int nmk, *uc;
  short i, j;

  uc = int_array1(np+1);

  uc[0] = 0;
  uc[1] = egeom[0]->ncontpts + 1;

  if(PERTUR == 1)
    nmk = NMK - 2;
  else
    nmk = NMK;

  for(i=2; i<=np; i++)
    uc[i] = uc[i-1] + egeom[i-1]->ncontpts - nmk;

  uc[np] = bspl->ncontpts;

  negeom = recover_cv(bspl, np, uc, NMK);
  
  if (PERTUR == 0)
    for (i=1; i<=np; i++)
      uc[i] = uc[i-1] + egeom[i-1]->ncontpts;
  else {
    for (i=2; i<np; i++)
      uc[i] = uc[i-1] + egeom[i-1]->ncontpts + 2;
    uc[np] = uc[np-1] + egeom[np-1]->ncontpts + 1;
  }

  for (i=0; i<np; i++) {
    neg = extract_cv(negeom, i, uc);
    errors = ConvexHullTest(neg, egeom[i]);
    if(errors > FINAL_ERROR)
      FINAL_ERROR = errors;
    free_egeom(neg);
  }
  free_egeom(negeom);
  free_iarray1(uc);

  return FINAL_ERROR;
}

ParCurv *extract_cv(ParCurv *negeom, int m, int *uc) 
{
  ParCurv *neg;
  double a, b, c;
  short i, j;

  neg = egeomalloc1(negeom->order, uc[m+1]-uc[m]);
  
  for (i=uc[m],j=0; i<uc[m+1]+negeom->order; i++,j++)
    neg->knots[j] = negeom->knots[i];

  for (i=0; i<neg->ncontpts; i++)
    neg->contpts[i] = copyvector(negeom->contpts[i+uc[m]], neg->contpts[i]);

  a = neg->knots[neg->ncontpts] - neg->knots[0];
  c = neg->knots[0];  
  for (i=0; i<neg->ncontpts+neg->order; i++)
    neg->knots[i] = (neg->knots[i] - c)/a;

  return neg;
}

ParCurv *recover_cv(ParCurv *bspl, int np, int *uc, int NMK)
{
  ParCurv *negeom;
  double **inpoly, **outpoly;
  short i, j, k;

  negeom = egeomalloc1(bspl->order, bspl->ncontpts+(np-1)*NMK);

  for (i=k=0; i< np; i++) {
    for (j=uc[i]; j<uc[i+1]; j++,k++)
      negeom->knots[k] = bspl->knots[j];
    for (j=0; j<NMK; j++,k++)
      negeom->knots[k] = bspl->knots[uc[i+1]];
  }
  for (i=1; i<bspl->order+1; i++)
    negeom->knots[negeom->ncontpts+negeom->order-i] = 1.0;

  inpoly = dbl_array2(bspl->ncontpts,4);
  outpoly = dbl_array2(negeom->ncontpts, 4);
  for (i=0; i<bspl->ncontpts; i++) {
    inpoly[i][0] = bspl->contpts[i]->x;
    inpoly[i][1] = bspl->contpts[i]->y;
    inpoly[i][2] = bspl->contpts[i]->z;
    inpoly[i][3] = bspl->contpts[i]->w;
  }

  curve_oslo1(bspl->order, bspl->ncontpts, negeom->order+negeom->ncontpts,
	      bspl->knots, negeom->knots, inpoly, outpoly);

  for (i=0; i<negeom->ncontpts; i++) {
    negeom->contpts[i]->x = outpoly[i][0];
    negeom->contpts[i]->y = outpoly[i][1];
    negeom->contpts[i]->z = outpoly[i][2];
    negeom->contpts[i]->w = outpoly[i][3];
  }

  free_darray2(inpoly);
  free_darray2(outpoly);

  return negeom;
}

double get_error(ParSurf ***fgeom, int nps, int npt, ParSurf *bspl_surf,
		 int PERTUR, int NMK)
{
  ParSurf *nfgeom, *nfg;
  double errors, FINAL_ERROR=0.0;
  int nmk, *uc, *vc;
  short i, j;

  uc = int_array1(nps+1);
  vc = int_array1(npt+1);

  uc[0] = 0;
  vc[0] = 0;  
  uc[1] = fgeom[0][0]->ucontpts + 1;
  vc[1] = fgeom[0][0]->vcontpts + 1;

  if (PERTUR == 1)
    nmk = NMK - 2;
  else
    nmk = NMK;

  for (i=2; i<=nps; i++)
    uc[i] = uc[i-1] + fgeom[0][i-1]->ucontpts - nmk;

  for (i=2; i<=npt; i++)
    vc[i] = vc[i-1] + fgeom[i-1][0]->vcontpts - nmk;

  uc[nps] = bspl_surf->ucontpts;
  vc[npt] = bspl_surf->vcontpts;

  nfgeom = recover_patch(bspl_surf, nps, npt, uc, vc, NMK);

  if (PERTUR == 0) {
    for (i=1; i<=nps; i++)
      uc[i] = uc[i-1] + fgeom[0][i-1]->ucontpts;
    for (i=1; i<=npt; i++)
      vc[i] = vc[i-1] + fgeom[i-1][0]->vcontpts;
  }
  else {
    for (i=2; i<nps; i++)
      uc[i] = uc[i-1] + fgeom[0][i-1]->ucontpts + 2;
/***** only add extra 1 if there are multiple surfaces in u direction */
    uc[nps] = uc[nps-1] + fgeom[0][nps-1]->ucontpts + (nps > 1 ? 1 : 0);
/***** the following two lines were wrong */
/** for (i=2; i<=npt; i++) ****************************/
/**** vc[i] = vc[i-1] + fgeom[i-1][0]->vcontpts + 1; **/
/***** this is what they should really look like */
    for (i=2; i<npt; i++)
      vc[i] = vc[i-1] + fgeom[i-1][0]->vcontpts + 2;
/***** only add extra 1 if there are multiple surfaces in v direction */
    vc[npt] = vc[npt-1] + fgeom[npt-1][0]->vcontpts + (npt > 1 ? 1 : 0);
  }

  for (i=0; i<nps; i++)
    for (j=0; j<npt; j++) {
      nfg = extract_patch(nfgeom, i, j, uc, vc);
      errors = ConvexHullSurf(nfg, fgeom[j][i], 500, 0);

      if (errors > FINAL_ERROR)
	FINAL_ERROR = errors;
      free_fgeom(nfg);
    }

  free_fgeom(nfgeom);
  free_iarray1(uc);
  free_iarray1(vc);

  return FINAL_ERROR;
}

ParSurf *extract_patch(ParSurf *nfgeom, int m, int n, int *uc, int *vc) 
{
  ParSurf *nfg;
  int i, j;
  double a, b, c;

  nfg = fgeomalloc1(nfgeom->uorder, nfgeom->vorder, uc[m+1]-uc[m], 
		    vc[n+1]-vc[n]);
  
  for (i=uc[m],j=0; i<uc[m+1]+nfgeom->uorder; i++,j++)
    nfg->uknots[j] = nfgeom->uknots[i];

  j=0;

  for (i=vc[n],j=0; i<vc[n+1]+nfgeom->vorder; i++,j++)
    nfg->vknots[j] = nfgeom->vknots[i];

  for (i=0; i<nfg->ucontpts; i++)
    for (j=0; j<nfg->vcontpts; j++)
      nfg->contpts[i][j] = copyvector(nfgeom->contpts[i+uc[m]][j+vc[n]],
				      nfg->contpts[i][j]);

  a = nfg->uknots[nfg->ucontpts] - nfg->uknots[nfg->uorder-1];

  c = nfg->uknots[nfg->uorder-1];  
  for (i=0; i<nfg->ucontpts+nfg->uorder; i++) {
    nfg->uknots[i] = (nfg->uknots[i] - c)/a;

  }
  c = nfg->vknots[nfg->vorder-1];
  b = nfg->vknots[nfg->vcontpts] - nfg->vknots[nfg->vorder-1];
  for (i=0; i<nfg->vcontpts+nfg->vorder; i++)
    nfg->vknots[i] = (nfg->vknots[i] - c)/b;

  return nfg;
}

ParSurf *recover_patch(ParSurf *bspls, int nps, int npt, int *uc, int *vc,
		       int NMK)
{
  ParSurf *nfgeom;
  short i, j, k;

  nfgeom = fgeomalloc1(bspls->uorder, bspls->vorder,
		       bspls->ucontpts+(nps-1)*NMK,
		       bspls->vcontpts+(npt-1)*NMK);

  for (i=k=0; i< nps; i++) {
    for(j=uc[i]; j<uc[i+1]; j++,k++)
      nfgeom->uknots[k] = bspls->uknots[j];

    for (j=0; j<NMK; j++,k++)
      nfgeom->uknots[k] = bspls->uknots[uc[i+1]];
  }
  for (i=1; i<bspls->uorder+1; i++)
    nfgeom->uknots[nfgeom->ucontpts+nfgeom->uorder-i] = 1.0;

  for (i=k=0; i< npt; i++) {
    for (j=vc[i]; j<vc[i+1]; j++,k++)
      nfgeom->vknots[k] = bspls->vknots[j];
    for (j=0; j<NMK; j++,k++)
      nfgeom->vknots[k] = bspls->vknots[vc[i+1]];
  }
  for (i=1; i<bspls->vorder+1; i++)
    nfgeom->vknots[nfgeom->vcontpts+nfgeom->vorder-i] = 1.0;

  surfoslo3(bspls, nfgeom, 4);

  return nfgeom;
}

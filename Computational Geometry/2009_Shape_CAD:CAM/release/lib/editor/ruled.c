/* Copyright (C) 1992 Massachusetts Institute of Technology, Cambridge, MA
   All rights reserved
*/

/* ruled.c */

/* MergeKnots()
   RulSurf_to_ParSurf()
*/

/* base on: /usr/deslab/epraxiteles/praxlib/editor/rs_conv.c */

#include <stdio.h>
#include "gen.h"
#include "bspl.h"
#include "editor.h"

ParSurf *RulSurf_to_ParSurf(RulSurf *rgeom, ParSurf *fgeom)
{
  ParCurv *egeom[2];
  double eps = 1.0e-12;
  int order, nknots;
  short i;

  MergeKnots(rgeom->de1, rgeom->de2, rgeom->dirflg, egeom, &nknots, &order);

  if (!fgeom)
    fgeom = fgeomalloc1(order, 2, egeom[0]->ncontpts, 2);
  
  for (i=0; i<2; i++) {
    fgeom->vknots[i] = 0.0;
    fgeom->vknots[i+2] = 1.0;
  }

  for (i=0; i<nknots; i++)
    fgeom->uknots[i] = egeom[0]->knots[i];
  for (i=0; i < nknots - order; i++)
    fgeom->contpts[i][0] = copyvector(egeom[0]->contpts[i],
				      fgeom->contpts[i][0]);

  for (i=0; i < nknots - order; i++)
    fgeom->contpts[i][1] = copyvector(egeom[1]->contpts[i],
				      fgeom->contpts[i][1]);

  free_egeom(egeom[0]);
  free_egeom(egeom[1]);

  return (fgeom);
}

void MergeKnots(ParCurv *de1, ParCurv *de2, short dirflg, ParCurv **egeom,
		int *nknots, int *order)
{
  ParCurv *bgeom[2];
  double eps = 1.0e-12, *knot, **lim;
  int indx[2], ncontpts;
  short i, k;

  bgeom[0] = de1;
  bgeom[1] = de2;

  *order = MAX(bgeom[0]->order, bgeom[1]->order);
  ncontpts = MAX(bgeom[0]->ncontpts, bgeom[1]->ncontpts);
  indx[0] = MAX(1, MAX((*order)-bgeom[0]->order, (*order)-bgeom[1]->order));
  knot = dbl_array1(((*order)+ncontpts)*2*indx[0]);
  
  egeom[0] = egeomalloc1((*order), 2*indx[0]*ncontpts);
  egeom[1] = egeomalloc1((*order), 2*indx[0]*ncontpts);
  
  copyegeom(bgeom[0], egeom[0]);
  copyegeom(bgeom[1], egeom[1]);

  if (dirflg == 1)
    for (i=0; i<bgeom[1]->ncontpts; i++)
      egeom[1]->contpts[ncontpts-1-i] = 
	copyvector(bgeom[1]->contpts[i], egeom[1]->contpts[ncontpts-1-i]);

  (*nknots) = (*order) + ncontpts;

  if ((*order) > egeom[0]->order) {
    lim = dbl_array2(((*order)+ncontpts)*2*indx[0], ncontpts*2*indx[0]);
    for (i=0; i < (*order) - bgeom[0]->order; i++) {
      k = compute_lim(egeom[0], lim, knot, eps);
      raise_byone(egeom[0], lim, knot, k);
    }
    free_darray2(lim);
  }
  else if ((*order) > egeom[1]->order) {
    lim = dbl_array2(((*order)+ncontpts)*2*indx[0], ncontpts*2*indx[0]);      
    for (i=0; i < (*order) - bgeom[1]->order; i++) {
      k = compute_lim(egeom[1], lim, knot, eps);
      raise_byone(egeom[1], lim, knot, k);
    }
    free_darray2(lim);
  }
  (*nknots) = merge_knotv(egeom, 0, 2, knot, indx, eps);
  addpoints(egeom, 0, 2, knot, (*nknots));

  free_darray1(knot);
}

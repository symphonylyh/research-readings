/* Copyright (C) 1993 Massachusetts Institute of Technology, Cambridge, MA
   All rights reserved
*/

/* perturb.c */

/* max_pert()
   max_pert_cv()
   perturb_knots()
   perturb_knots_cv()
*/

#include "gen.h"
#include "editor.h"

double max_pert_cv(ParCurv **egeoms, int np, ParCurv *bspl, double err,
		   double factor)
{
  double maxPer=0.0, pert;
  short i, j;

  for (i=j=0; i<np; i++) {
    j += egeoms[i]->ncontpts;
    pert = (bspl->knots[j] - bspl->knots[j-1])/2.0;
    if (pert > maxPer)
      maxPer = pert;
  }
  maxPer = MIN(maxPer, err*factor);

  return maxPer;
}

void perturb_knots_cv(ParCurv **egeoms, int np, ParCurv *bspls, double PERT)
{
  int *uc;
  short i, j;

  uc = int_array1(np);

  uc[0] = egeoms[0]->ncontpts+1;
  for (i=1; i<np-1; i++)
    uc[i] = uc[i-1] + egeoms[i]->ncontpts - 1;
  for (i=0; i<np-1; i++) {
    j = uc[i];
    bspls->knots[j-1] -= PERT;
    bspls->knots[j+1] += PERT; 
  }
  free_iarray1(uc);
}

double max_pert(ParSurf ***fgeom, int nps, int npt, ParSurf *bspls, double err,
		double factor)
{
  double maxPer=0.0, pert;
  short i, j;
  
  for (i=j=0; i<nps; i++) {
    j += fgeom[0][i]->ucontpts;
    pert = (bspls->uknots[j] - bspls->uknots[j-1])/2.0;
    if(pert > maxPer)
      maxPer = pert;
  }
    
  for (i=j=0; i<npt; i++) {
    j += fgeom[i][0]->vcontpts;
    pert = (bspls->vknots[j] - bspls->vknots[j-1])/2.0;
    if(pert > maxPer)
      maxPer = pert;
  }
  maxPer = MIN(maxPer, err*factor);

  return maxPer;
}

void perturb_knots(ParSurf *surf, int nps, int npt, int *uc, int *vc,
		   double PERT)
{
  short i, j;
  
  uc[0]++;
  for(i=1; i<nps-1; i++)
    uc[i] = uc[i-1] + uc[i] - 1;

  for (i=0; i<nps-1; i++) {
    j = uc[i];
    surf->uknots[j-1] -= PERT;
    surf->uknots[j+1] += PERT; 
  }
  vc[0]++;
  for(i=1; i<npt-1; i++)
    vc[i] = vc[i-1] + vc[i] - 1;

  for (i=0; i<npt-1; i++) {
    j = vc[i] ;
    surf->vknots[j-1] -= PERT;
    surf->vknots[j+1] += PERT;
  }
}

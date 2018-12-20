/* Copyright (C) 1993 Massachusetts Institute of Technology, Cambridge, MA
   All rights reserved
*/

/* split_surf.c */

/* SplitSurf()
   SubdivSurfNxM()
*/

#include "gen.h"
#include "bspl.h"
#include "editor.h"

void SubdivSurfNxM(ParSurf *fgeom, int nu, int nv, double *ui, double *vi,
		   ParSurf ***fgeoms)
{
  ParSurf *ugeom, *vgeom, *utemp1, *utemp2, *vtemp1, *vtemp2;
  double u, v;
  short i, j;

  ugeom = copyfgeom(fgeom, NULL);
  for (i=0; i<=nu; i++) {
    if (i < nu) {
      u = i ? (ui[i] - ui[i-1]) / (1.0 - ui[i-1]) : ui[i];
      SplitSurf(ugeom, 1, u, &utemp1, &utemp2);
      free_fgeom(ugeom);
      ugeom = copyfgeom(utemp1, NULL);
      free_fgeom(utemp1);
    }
    vgeom = copyfgeom(ugeom, NULL);
    for (j=0; j<=nv; j++) {
      if (j < nv) {
	v = j ? (vi[j] - vi[j-1]) / (1.0 - vi[j-1]) : vi[j];
	SplitSurf(vgeom, 0, v, &vtemp1, &vtemp2);
	free_fgeom(vgeom);
	vgeom = copyfgeom(vtemp1, NULL);
	free_fgeom(vtemp1);
      }
      fgeoms[i][j] = vgeom;
      if (j < nv) {
	vgeom = copyfgeom(vtemp2, NULL);
	free_fgeom(vtemp2);
      }
    }
    if (i < nu) {
      free_fgeom(ugeom);
      ugeom = copyfgeom(utemp2, NULL);
      free_fgeom(utemp2);
    }
  }
  free_fgeom(ugeom);
}

void SplitSurf(ParSurf *surf, int uFlag, double param, ParSurf **fgeom1,
	       ParSurf **fgeom2)
{
  ParSurf *fgeom, *temp;
  short i, j, k, n, indx, uKnots, vKnots;

  if (uFlag) {
    uKnots = surf->uorder;
    vKnots = 0;
  }
  else {
    uKnots = 0;
    vKnots = surf->vorder;
  }

  fgeom = copyfgeom(surf, 0);
  temp = copyfgeom(fgeom, 0);

  free_fgeom(fgeom);
  fgeom = fgeomalloc1(temp->uorder, temp->vorder, temp->ucontpts + uKnots,
		      temp->vcontpts + vKnots);
  copyfgeom(temp, fgeom);
  fgeom->ucontpts = temp->ucontpts + uKnots;
  fgeom->vcontpts = temp->vcontpts + vKnots;

  if (uKnots) {
    fgeom->uknots[0] = temp->uknots[0];
    for (i=1,j=0; i < temp->ucontpts + fgeom->uorder; i++) {
      if(temp->uknots[i] >= param && temp->uknots[i-1] < param)
	for(j=0; j<uKnots; j++)
	  fgeom->uknots[i+j] = param;
      fgeom->uknots[i+j] = temp->uknots[i];
    }
  }

  if (vKnots) {
    fgeom->vknots[0] = temp->vknots[0];
    for (i=1,j=0; i < temp->vcontpts + fgeom->vorder; i++) {
      if(temp->vknots[i] >= param && temp->vknots[i-1] < param)
	for(j=0; j<vKnots; j++)
	  fgeom->vknots[i+j] = param;
      fgeom->vknots[i+j] = temp->vknots[i];
    }
  }
  surfoslo3(temp, fgeom, 4);
  free_fgeom(temp);

  if (uKnots) {
    for(i=fgeom->uorder; i< fgeom->ucontpts; i++) {
/* allow the surface to have multiple internal knots, other than at param */
      if(fgeom->uknots[i] == param &&
	 fgeom->uknots[i] == fgeom->uknots[i+1]) {
	indx = 2;
	for (j=i+2; j<i+fgeom->uorder && j<fgeom->ucontpts; j++)
	  if (fgeom->uknots[i] == fgeom->uknots[j])
	    indx++;
      }
      if(indx == fgeom->uorder) {
	indx = i + fgeom->uorder;

	(*fgeom1) = fgeomalloc1(fgeom->uorder, fgeom->vorder,
				   indx - fgeom->uorder, fgeom->vcontpts);
	(*fgeom2) = fgeomalloc1(fgeom->uorder, fgeom->vorder,
				   fgeom->ucontpts - indx + fgeom->uorder,
				   fgeom->vcontpts);
	for(j=0; j<indx; j++)
	  (*fgeom1)->uknots[j] = fgeom->uknots[j];
	k=0;
	for(j=indx-fgeom->uorder; j<fgeom->uorder+fgeom->ucontpts; j++) {
	  (*fgeom2)->uknots[k] = fgeom->uknots[j];
	  k++;
	}
	for(n=0; n<fgeom->vcontpts; n++)
	  for(j=0;j<indx - fgeom->uorder;j++)
	    copyvector(fgeom->contpts[j][n],(*fgeom1)->contpts[j][n]);
	    
	for(n=0;n<fgeom->vcontpts;n++) {
	  k = fgeom->ucontpts - indx + fgeom->uorder;
	  for(j = fgeom->ucontpts - 1;
	      j >= indx - fgeom->uorder;
	      j--) {
	    k--;
	    copyvector(fgeom->contpts[j][n], (*fgeom2)->contpts[k][n]);
	  }
	}
	  
	for(n=0;n<fgeom->vcontpts + fgeom->vorder;n++){
	  (*fgeom1)->vknots[n] = fgeom->vknots[n];
	  (*fgeom2)->vknots[n] = fgeom->vknots[n];
	}
	knot_normalize((*fgeom1)->uknots, (*fgeom1)->ukmem);
	knot_normalize((*fgeom2)->uknots, (*fgeom2)->ukmem);
	knot_normalize((*fgeom1)->vknots, (*fgeom1)->vkmem);
	knot_normalize((*fgeom2)->vknots, (*fgeom2)->vkmem);
	
	break;
      }
    }
  }
  else if (vKnots) {
    for (i=fgeom->vorder; i< fgeom->vcontpts; i++) {
/* allow the surface to have multiple internal knots, other than at param */
      if (fgeom->vknots[i] == param &&
	  fgeom->vknots[i] == fgeom->vknots[i+1]) {
	indx = 2;
	for (j=i+2; j<i+fgeom->vorder && j<fgeom->vcontpts; j++)
	  if (fgeom->vknots[i] == fgeom->vknots[j])
	    indx++;
      }
      if(indx == fgeom->vorder) {
	indx = i + fgeom->vorder;

	(*fgeom1) = fgeomalloc1(fgeom->uorder, fgeom->vorder,
				   fgeom->ucontpts, indx - fgeom->vorder);
	(*fgeom2) = fgeomalloc1(fgeom->uorder,fgeom->vorder,
				   fgeom->ucontpts, fgeom->vcontpts-indx +
				   fgeom->vorder);

	for (j=0; j<indx; j++)
	  (*fgeom1)->vknots[j] = fgeom->vknots[j];
	k=0;
	for (j=indx-fgeom->vorder; j<fgeom->vorder+fgeom->vcontpts; j++) {
	  (*fgeom2)->vknots[k] = fgeom->vknots[j];
	  k++;
	}
	for (n=0; n<fgeom->ucontpts; n++)
	  for (j=0; j<indx - fgeom->vorder; j++)
	    copyvector(fgeom->contpts[n][j],(*fgeom1)->contpts[n][j]);
	    
	for (n=0;n<fgeom->ucontpts;n++) {
	  k = fgeom->vcontpts - indx + fgeom->vorder;
	  for(j = fgeom->vcontpts - 1;
	      j >= indx - fgeom->vorder;
	      j--) {
	    k--;
	    copyvector(fgeom->contpts[n][j], (*fgeom2)->contpts[n][k]);
	  }
	}
	
	for(n=0; n<fgeom->ucontpts + fgeom->uorder; n++){
	  (*fgeom1)->uknots[n] = fgeom->uknots[n];
	  (*fgeom2)->uknots[n] = fgeom->uknots[n];
	}
	knot_normalize((*fgeom1)->vknots, (*fgeom1)->vkmem);
	knot_normalize((*fgeom2)->vknots, (*fgeom2)->vkmem);
	knot_normalize((*fgeom1)->uknots, (*fgeom1)->ukmem);
	knot_normalize((*fgeom2)->uknots, (*fgeom2)->ukmem);
	
	break;
      }
    }
  }
  free_fgeom(fgeom);
}

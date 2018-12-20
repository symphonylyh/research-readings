/* Copyright (C) 1993 Massachusetts Institute of Technology, Cambridge, MA
   All rights reserved
*/

/* subdivde.c */

/* subdiv_bez()
   subdiv_cv()
*/

#include "gen.h"
#include "bspl.h"
#include "editor.h"

/* iu = 0     do not subdiv */
/* iu = 1     subdiv at the first span of left hand side */
/* iu = 2     subdiv at the last  span of left hand side */

void subdiv_cv(ParCurv *oegeom, ParCurv *negeom, int iu)
{
  double **inpoly, **outpoly;
  short i;

  if (iu==0) return;

  if (iu == 1) {
    for (i=0; i<oegeom->order; i++)
      negeom->knots[i] = oegeom->knots[i];

    for (i=oegeom->order; i<oegeom->order*2; i++)
      negeom->knots[i] = oegeom->knots[oegeom->order];

    for (i=oegeom->order*2; i<oegeom->order*2+oegeom->ncontpts-1; i++)
      negeom->knots[i] = oegeom->knots[i-oegeom->order+1];
  }

  if (iu == 2) {
    for (i=0; i<oegeom->ncontpts; i++)
      negeom->knots[i] = oegeom->knots[i];

    for (i=oegeom->ncontpts; i<oegeom->order+oegeom->ncontpts-1; i++)
      negeom->knots[i] = oegeom->knots[oegeom->ncontpts-1];

    for (i=oegeom->ncontpts; i<oegeom->order+oegeom->ncontpts; i++)
      negeom->knots[i+oegeom->order-1] = oegeom->knots[i];
  }

  inpoly = dbl_array2(oegeom->ncontpts,4);
  outpoly = dbl_array2(negeom->ncontpts, 4);
  for (i=0; i<oegeom->ncontpts; i++) {
    inpoly[i][0]=oegeom->contpts[i]->x;
    inpoly[i][1]=oegeom->contpts[i]->y;
    inpoly[i][2]=oegeom->contpts[i]->z;
    inpoly[i][3]=oegeom->contpts[i]->w;
  }

  curve_oslo1(oegeom->order, oegeom->ncontpts, negeom->order+negeom->ncontpts,
	      oegeom->knots, negeom->knots, inpoly, outpoly);

  for(i=0; i<negeom->ncontpts; i++) {
    negeom->contpts[i]->x=outpoly[i][0];
    negeom->contpts[i]->y=outpoly[i][1];
    negeom->contpts[i]->z=outpoly[i][2];
    negeom->contpts[i]->w=outpoly[i][3];
  }

  free_darray2(inpoly);
  free_darray2(outpoly);
}

/* iu = 0     do not subdiv along u-dir */
/* iu = 1     subdiv at the first span of left hand side */
/* iu = 2     subdiv at the last  span of left hand side */

void subdiv_bez(ParSurf *ofgeom, ParSurf *nfgeom, int iu, int iv)
{
  short i;

  if (iu == 1) {
    for(i=0; i<ofgeom->uorder; i++)
      nfgeom->uknots[i] = ofgeom->uknots[i];

    for(i=ofgeom->uorder; i<ofgeom->uorder*2; i++)
      nfgeom->uknots[i] = ofgeom->uknots[ofgeom->uorder];

    for(i=ofgeom->uorder*2; i<ofgeom->uorder*2+ofgeom->ucontpts-1; i++)
      nfgeom->uknots[i] = ofgeom->uknots[i-ofgeom->uorder+1];

    for(i=0; i<ofgeom->vorder+ofgeom->vcontpts; i++)
      nfgeom->vknots[i] = ofgeom->vknots[i];
  }

  if (iu == 2) {
    for(i=0; i<ofgeom->ucontpts; i++)
      nfgeom->uknots[i] = ofgeom->uknots[i];

    for(i=ofgeom->ucontpts; i<ofgeom->uorder+ofgeom->ucontpts-1; i++)
      nfgeom->uknots[i] = ofgeom->uknots[ofgeom->ucontpts-1];

    for(i=ofgeom->ucontpts; i<ofgeom->uorder+ofgeom->ucontpts; i++)
      nfgeom->uknots[i+ofgeom->uorder-1] = ofgeom->uknots[i];

    for(i=0; i<ofgeom->vorder+ofgeom->vcontpts; i++)
      nfgeom->vknots[i] = ofgeom->vknots[i];
  }

  if (iv == 1) {
    for(i=0; i<ofgeom->vorder; i++)
      nfgeom->vknots[i] = ofgeom->vknots[i];

    for(i=ofgeom->vorder; i<ofgeom->vorder*2; i++)
      nfgeom->vknots[i] = ofgeom->vknots[ofgeom->vorder];

    for(i=ofgeom->vorder*2; i<ofgeom->vorder*2+ofgeom->vcontpts-1; i++)
      nfgeom->vknots[i] = ofgeom->vknots[i-ofgeom->vorder+1];

    for(i=0; i<ofgeom->uorder+ofgeom->ucontpts; i++)
      nfgeom->uknots[i] = ofgeom->uknots[i];
  }


  if (iv == 2) {
    for(i=0; i<ofgeom->vcontpts; i++)
      nfgeom->vknots[i] = ofgeom->vknots[i];

    for(i=ofgeom->vcontpts; i<ofgeom->vorder+ofgeom->vcontpts-1; i++)
      nfgeom->vknots[i] = ofgeom->vknots[ofgeom->vcontpts-1];

    for(i=ofgeom->vcontpts; i<ofgeom->vorder+ofgeom->vcontpts; i++)
      nfgeom->vknots[i+ofgeom->vorder-1] = ofgeom->vknots[i];

    for(i=0; i<ofgeom->uorder+ofgeom->ucontpts; i++)
      nfgeom->uknots[i] = ofgeom->uknots[i];
  }

  surfoslo3(ofgeom, nfgeom, 4);
}

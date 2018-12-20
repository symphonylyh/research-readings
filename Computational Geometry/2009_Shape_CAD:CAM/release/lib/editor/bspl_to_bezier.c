/* Copyright (C) 1995 Massachusetts Institute of Technology, Cambridge, MA
   All rights reserved
*/

/* bspl_to_bezier.c */

/* BezierContpts()
   BezierKnots()
   BsplineToBezierKnots()
   BsplineToBezierSurf()
   copy_contpts()
   copy_knots()
   CopyUniqueKnots()
   FindUniqueKnots()
   FreeBezierContpts()
*/

#include <malloc.h>
#include <math.h>
#include "gen.h"
#include "bspl.h"
#include "editor.h"

ParSurf *** BsplineToBezierSurf(ParSurf *fgeom, int *nsegu, int *nsegv)
{
  ParSurf *fgm;
  ParSurf ***bezier;
  vector *****contpts;
  double *vknots, *uknots, **uknots_s, **vknots_s;
  int i, j, u2, v2;
  short n_ucontpts, n_uknots, n_vcontpts, n_vknots;

  if (fgeom->uorder == fgeom->ucontpts && fgeom->vorder == fgeom->vcontpts) {
    bezier = (ParSurf ***)gen_array2(1, 1, sizeof(ParSurf *));
    bezier[0][0] = copyfgeom(fgeom, NULL);

    *nsegu = *nsegv = 1;
  }
  else {
    uknots = BsplineToBezierKnots(fgeom->uknots, fgeom->uorder, 
				  fgeom->ukmem, &n_ucontpts);
      
    vknots = BsplineToBezierKnots(fgeom->vknots, fgeom->vorder, 
				  fgeom->vkmem, &n_vcontpts);
      
    fgm = fgeomalloc1(fgeom->uorder, fgeom->vorder, n_ucontpts, n_vcontpts);
      
    n_uknots = n_ucontpts + fgeom->uorder;
    n_vknots = n_vcontpts + fgeom->vorder;
      
    *nsegu = n_ucontpts / fgeom->uorder;
    *nsegv = n_vcontpts / fgeom->vorder;
      
    copy_knots(fgm->uknots, uknots, n_uknots);
    copy_knots(fgm->vknots, vknots, n_vknots);
      
    surfoslo3(fgeom, fgm, 4);  
      
    uknots_s = BezierKnots(uknots, fgeom->uorder, n_ucontpts);
    vknots_s = BezierKnots(vknots, fgeom->vorder, n_vcontpts);
    free_darray1(uknots);
    free_darray1(vknots);

    contpts = BezierContpts(fgm->contpts, *nsegu, *nsegv, fgm->uorder,
			    fgm->vorder);

    bezier = (ParSurf ***)gen_array2(*nsegu, *nsegv, sizeof(ParSurf *));
    for (i=0; i< *nsegu; i++)
      for (j=0; j < *nsegv; j++)
	bezier[i][j] = fgeomalloc1(fgeom->uorder, fgeom->vorder,
				   fgeom->uorder, fgeom->vorder);
    u2 = 2*(fgeom->uorder);
    v2 = 2*(fgeom->vorder);
      
    for (i=0; i < *nsegu; i++)
      for (j=0; j < *nsegv; j++) {
	copy_knots(bezier[i][j]->uknots, uknots_s[i], u2);
	copy_knots(bezier[i][j]->vknots, vknots_s[j], v2);
      }
    free_darray2(uknots_s);
    free_darray2(vknots_s);

    for (i=0; i < *nsegu; i++)
      for (j=0; j < *nsegv; j++)
	copy_contpts(bezier[i][j]->contpts, contpts[i][j], fgeom->uorder,
		     fgeom->vorder);
      
    FreeBezierContpts(contpts, *nsegu, *nsegv, fgm->uorder, fgm->vorder);
    free_fgeom(fgm);
  }
  return bezier;
}

vector *****BezierContpts(vector ***contpts, short nsegu, short nsegv,
			  short uorder, short vorder)
{
  vector *****contpts_n;
  short i, j, k, l;

  contpts_n = (vector *****)calloc(nsegu, sizeof(vector ****));
  for (i=0; i<nsegu; i++)
    contpts_n[i]=(vector ****)calloc(nsegv, sizeof(vector ***));

  for (i=0; i<nsegu; i++)
    for (j=0; j<nsegv; j++)
      contpts_n[i][j] = vec_array2(uorder, vorder);

  for (i=0; i<nsegu; i++)
    for (j=0; j<nsegv; j++)
      for (k=0; k<uorder; k++)
	for (l=0; l<vorder; l++)
	  copyvector(contpts[i*uorder + k][j*vorder + l],
		     contpts_n[i][j][k][l]);

  return contpts_n;
}

void FreeBezierContpts(vector *****contpts, short nsegu, short nsegv,
		       short uorder, short vorder)
{
  short i, j;

  for (i=0; i<nsegu; i++)
    for (j=0; j<nsegv; j++)
      free_varray2(contpts[i][j], uorder, vorder);

  for (i=0; i<nsegu; i++)
    free(contpts[i]);

  free(contpts);
}

double **BezierKnots(double *knots, short order, short ncontpts)
{
  double **knots_s;
  short i, j, nseg1, nseg2;

  nseg1 = ncontpts/order;
  nseg2 = 2*order;
  knots_s = dbl_array2(nseg1, nseg2);

  for (i=0; i<nseg1; i++)
    for (j=0; j<nseg2; j++)
      knots_s[i][j] = knots[i*order + j];

  return knots_s;
}

double *BsplineToBezierKnots(double *knots, short order, short nknots,
			     short *ncontpts)
{
  double *knots_b, *knots_n;
  short nUnique;

  knots_n = FindUniqueKnots(knots, nknots, &nUnique);
  knots_b = CopyUniqueKnots(knots_n, order, ncontpts, nUnique);
  free_darray1(knots_n);

  return knots_b;
}

double *FindUniqueKnots(double *knots, short nknots, short *nUnique)
{
  double *knots_n;
  int i, j, *nk;

  nk = int_array1(nknots);
  nk[0] = 0;

  for (i=j=1; i<nknots; i++)
/** if (fabs(knots[i] - knots[i-1]) > ZERO) { **/
/** make this consistent with curve routine checknot() in bezierop.c **/
    if (fabs(knots[i] - knots[i-1]) >= EDITOR_TOLKNOT) {
      nk[j] = i;
      j++;
    }

  *nUnique = j;
  knots_n = dbl_array1(*nUnique);
  for (i=0; i < *nUnique; i++)
    knots_n[i] = knots[nk[i]];
  free_iarray1(nk);

  return knots_n;
}

double *CopyUniqueKnots(double *knots_n, short order, short *new_ncontpts,
			short nUnique)
{
  double *knots_b;
  int i, j, nknots;

  nknots = order*nUnique;

  *new_ncontpts = nknots - order;
  knots_b = dbl_array1(nknots);

  for (i=0; i<nUnique; i++)
    for (j=0; j<order; j++)
      knots_b[i*order + j] = knots_n[i];

  return knots_b;
}

void copy_contpts(vector ***result, vector ***original_contpts,
		  short ucontpts, short vcontpts)
{
  short i , j;

  for (i=0; i<ucontpts; i++)
    for (j=0; j<vcontpts; j++)
      copyvector(original_contpts[i][j], result[i][j]);
}

void copy_knots(double *knots, double *new_knots, short nknots)
{
  short i;

  for (i=0; i<nknots; i++)
    knots[i] = new_knots[i];
}

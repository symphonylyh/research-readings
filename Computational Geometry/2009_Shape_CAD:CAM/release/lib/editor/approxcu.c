/************************************************************************
 *									*
			Copyright (C) 1992 by
	Massachusetts Institute of Technology, Cambridge, MA
			 All rights reserved
 *									*
 ************************************************************************/

#include <stdio.h>
#include <math.h>
#include "gen.h"
#include "bspl.h"
#include "appr.h"
#include "editor.h"

void e02baf_(int *, int *, double *, double *,double *, double *,
	     double *, double *, double *, double *, int *);

void approx_cubic_curv1(ParCurv *egeom, double **xyz, double *w, double *u,
			int npoints, int ncap, short ip, double *maxerr)
{
  vector *vect;
  double2 *wk1, *x, *c;
  double2 ss, **cxyz, **wk2, xx[1];
  int ifail = 0;
  short i, j;

  c = dbl_array1(npoints+4);
  x = dbl_array1(npoints);
  wk1 = dbl_array1(npoints);
  cxyz = dbl_array2(npoints, 3);
  wk2 = dbl_array2(npoints+4, 4);

  calc_par_chord(xyz, u, npoints);

  egeom->ncontpts = ncap-4;

  for (i=0; i<4; i++)
    egeom->knots[i] = 0.0;
  for (i=4; i<ncap-4; i++)
    egeom->knots[i] = ((double2)i-3)/((double2)ncap-7);
  for (i=ncap-4; i<ncap; i++)
    egeom->knots[i] = 1.0;

  vect = vectalloc();

  for (j=0; j<3; j++) {
    for (i=0; i<npoints; i++) 
      x[i] = xyz[i][j];
    e02baf_(&npoints,&ncap,u,x,w,egeom->knots,wk1,&wk2[0][0],c,&ss,&ifail);
    for (i=0; i<egeom->ncontpts; i++) 
      cxyz[i][j] = c[i];
  }

  for (j=0; j<egeom->ncontpts; j++) {
    egeom->contpts[j]->x = cxyz[j][0];
    egeom->contpts[j]->y = cxyz[j][1];
    egeom->contpts[j]->z = cxyz[j][2];
    egeom->contpts[j]->w = 1.0;
  }  

  free_darray1(c);
  free_darray1(x);
  free_darray1(wk1);
  free_darray2(cxyz);
  free_darray2(wk2);
}

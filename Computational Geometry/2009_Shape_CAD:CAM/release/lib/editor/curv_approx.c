/* Copyright (C) 1992 Massachusetts Institute of Technology, Cambridge, MA
   All rights reserved
*/

/* approxcurve.c */

/* approx_user_knots
*/

#include <malloc.h>
#include "gen.h"
#include "appr.h"
#include "editor.h"

void e02baf_(int*, int*, double*, double*, double*, double*, double*,
	     double*, double*, double*, int*);

short approx_user_knots(ParCurv *egeom, double **xyz, double *w, int npoints,
			int ncap, int ip, double *knots, double *maxerr)
{
  vector *vect;
  double2 *c, **cxyz, ss, *u, *wk1, **wk2, *x;
  int ifail;
  short i, j;
  
  vect = vectalloc();

  c = dbl_array1(ncap);
  x = dbl_array1(npoints);
  u = dbl_array1(npoints);
  wk1 = dbl_array1(npoints);
  cxyz = dbl_array2(ncap-4, 3);
  wk2 = dbl_array2(ncap, 4);

  if (ip == 1)
    calc_par_chord(xyz, u, npoints);
  else if (ip == 2)
    calc_par_uniform(xyz, u, npoints);
  else if (ip == 3)
    calc_par_foley(xyz, u, npoints);
  else if (ip == 4)
    calc_par_hartley(xyz, u, npoints);

  for (j=0; j<3; j++){
    ifail = 1;
    for (i=0; i<npoints; i++)
      x[i] = xyz[i][j];
    e02baf_(&npoints, &ncap, u, x, w, knots, wk1, &wk2[0][0], c, &ss, &ifail);
    if(ifail == 5) { /* this checks for singular matrices */
      egeom->type = -1;
      free_darray1(u);
      free_darray1(c);  
      free_darray1(x);
      free_darray1(wk1);
      free_darray2(wk2);
      free_darray2(cxyz);
      vectfree(vect);
      return (-1);
    }
    for (i=0; i<ncap-4; i++)
      cxyz[i][j] = c[i];
  }

  for (j=0; j<ncap; j++)
    egeom->knots[j] = knots[j];

  for (j=0; j<ncap-4; j++){
    vect->x = cxyz[j][0];
    vect->y = cxyz[j][1];
    vect->z = cxyz[j][2];
    vect->w = 1.0;
    copyvector(vect, egeom->contpts[j]);
  }

  find_error_curv1(egeom, xyz, u, npoints, maxerr);

  egeom->type = PCurveOpen;

  free_darray1(u);
  free_darray1(c);
  free_darray1(x);
  free_darray1(wk1);
  free_darray2(wk2);
  free_darray2(cxyz);
  vectfree(vect);

  return (0);
}

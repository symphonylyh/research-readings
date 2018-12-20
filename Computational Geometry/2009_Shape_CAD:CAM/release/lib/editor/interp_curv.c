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
#include "editor.h"

ParCurv *int_curve(double **cq, short ordera, short n, short cntr)
{
  ParCurv *egeom;
  double *xx,*yy,*wk1,**wk2,*cc,*tau,*ww;
  double **nij,**rhs,**sol,**aa,**bb,*wk;
  short i;

  egeom = egeomalloc1(ordera, n);
  tau = dbl_array1(n);
  if (cntr) {
    ww = dbl_array1(n);
    cc = dbl_array1(n+7);
    xx = dbl_array1(n+7);
    yy = dbl_array1(n);
    wk1 = dbl_array1(n); 
    wk2 = dbl_array2(n+7, 4);

    for (i=0; i<n; i++)
      ww[i] = 1.0;  
  
    for (i=0; i<ordera; i++)
      egeom->knots[i] = 0.0;
    opt_param(egeom, cq, tau, n);
    interp_points(egeom, cq, tau, xx, yy, ww, wk1, wk2, cc, n);

    free_darray1(xx);
    free_darray1(ww);
    free_darray1(cc);   
    free_darray1(yy);
    free_darray1(wk1);
    free_darray2(wk2);
  }
  else {
    nij = dbl_array2(n, n); 
    rhs = dbl_array2(3, n);
    sol = dbl_array2(3, n); 
    aa = dbl_array2(n, n);
    bb = dbl_array2(3, n);
    wk = dbl_array1(n+1);

    egeom->knots[0] = 0.0;
    opt_param_per(egeom, cq, tau, n);
    interp_points_per(egeom, cq, tau, nij, rhs, sol, aa, bb, wk, n);

    free_darray2(nij);
    free_darray2(rhs);
    free_darray1(wk);
    free_darray2(sol);
    free_darray2(aa);
    free_darray2(bb);
  }
  free_darray1(tau);

  return (egeom);
}

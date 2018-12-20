/************************************************************************
 *									*
			Copyright (C) 1992 by
	Massachusetts Institute of Technology, Cambridge, MA
			 All rights reserved
 *									*
 ************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include "gen.h"
#include "bspl.h"
#include "editor.h"

void f04aef_(double *, int *, double *, int *, int *, int *, double *,
	     int *, double *, double *, int *, double *, int *, int *);

void opt_param_per(ParCurv *egeom, double **cq, double *tau, short n)
{
  struct vector *v1,*v2;
  double s,x,y,z;
  short i,j,n1,k;

  egeom->knots[n] = 1.0;
  n1 = n-1;
  v1 = vectalloc();
  v1->x = cq[0][0]/cq[0][3];
  v1->y = cq[0][1]/cq[0][3];
  v1->z = cq[0][2]/cq[0][3];
  v1->w = 1.0;
  v2 = vectalloc();
  for (i=1; i<=n; i++) { 
    j = i%n; 
    v2->x = cq[j][0]/cq[j][3];
    v2->y = cq[j][1]/cq[j][3]; 
    v2->z = cq[j][2]/cq[j][3];
    v2->w = 1.0;
    x = v2->x-v1->x;
    y = v2->y-v1->y;
    z = v2->z-v1->z;
    tau[i-1] = sqrt(x*x+y*y+z*z);
    copyvector(v2,v1); 
  }
  vectfree(v1);
  vectfree(v2);
  s = 0.0; 
  for (i=0; i<n; i++) {
    for (j=0; j<3; j++) {
      k = (i+j)%n;
      s += tau[k];
    }
    if (i<n1)
      egeom->knots[i+1] = s;
  }
  for (i=1; i<n; i++)
    egeom->knots[i] /= s;
  nodes_bsp_per(egeom->knots, n, 4, tau);
}

void interp_points_per(ParCurv *egeom, double **cq, double *tau,
		       double **nij, double **rhs, double **sol,
		       double **aa,  double **bb, double *wk, int n)
{
  int i,j,k,ifail;
  
  n_matrix_per(nij, egeom->knots, tau, 4, n);
  
  for (i=0; i<n; i++) {
    j = (i+2)%n;
    rhs[0][i] = cq[j][0]; 
    rhs[1][i] = cq[j][1];
    rhs[2][i] = cq[j][2];
  }
  j = n;
  k = 3;
  ifail = 0;
  f04aef_(&nij[0][0], &j, &rhs[0][0], &j, &n, &k, &sol[0][0], &j, wk,
	  &aa[0][0], &j, &bb[0][0], &j, &ifail);
  if (ifail) {
    printf("interp_points_per: ifail = %d from f04aef!\n", ifail);
    exit(1);
  }
  for (i=0; i<n; i++) {
    egeom->contpts[i]->x = sol[0][i];
    egeom->contpts[i]->y = sol[1][i];
    egeom->contpts[i]->z = sol[2][i];
    egeom->contpts[i]->w = 1.0;
  }
}

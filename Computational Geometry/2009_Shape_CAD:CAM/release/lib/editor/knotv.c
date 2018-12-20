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

void e02baf_(int *, int *, double *, double *, double *, double *, double *,
	     double*, double *, double *, int *);

void opt_param(ParCurv *egeom, double **cq, double *tau, short n)
{
  vector *v1,*v2;
  double s,x,y,z;
  short i,j;

  for(i=0; i<4; i++)
    egeom->knots[i+n] = 1.0 ;
  
  v1 = vectalloc();
  v1->x = cq[0][0]/cq[0][3];
  v1->y = cq[0][1]/cq[0][3];
  v1->z = cq[0][2]/cq[0][3];
  v1->w = 1.0;
  v2 = vectalloc();
  for (i=1; i<n; i++) {
    v2->x = cq[i][0]/cq[i][3];
    v2->y = cq[i][1]/cq[i][3]; 
    v2->z = cq[i][2]/cq[i][3];
    v2->w = 1.0;
    x = v2->x-v1->x;
    y = v2->y-v1->y;
    z = v2->z-v1->z;
    tau[i-1] = sqrt(x*x+y*y+z*z);
    copyvector(v2, v1);
  }
  vectfree(v1);
  vectfree(v2);
  s = 0.0 ;
  for (i=0;i<n-3;i++) {
    for (j=0; j<3; j++)
      s += tau[i+j] ;
    if (i<n-4)
      egeom->knots[i+4] = s;
  }
  for (i=4; i<n; i++)
    egeom->knots[i] /= s;
  nodes_bsp(egeom->knots, n, 4, tau);
}

void interp_points(ParCurv *egeom, double **cq, double *tau, double *xx,
		   double *yy, double *ww, double *wk1, double **wk2,
		   double *cc, int n)
{
  double s;
  int ncap,i,j,ifail;
  
  ncap = n+4;
  for (i=0; i<ncap; i++)
    xx[i] = egeom->knots[i];
  for (j=0; j<3; j++) {
    for (i=0;  i<n; i++)
      yy[i] = cq[i][j]; 
    ifail = 0;
    e02baf_(&n, &ncap, tau, yy, ww, xx, wk1, &wk2[0][0], cc, &s, &ifail);
    if (ifail) {
      printf("interp_points: ifail = %d from e02baf!\n", ifail);
      exit(1);
    }
    for (i=0; i<n; i++) {
      if (j == 0)
	egeom->contpts[i]->x = cc[i];
      else if(j == 1)
	egeom->contpts[i]->y = cc[i];
      else
	egeom->contpts[i]->z = cc[i]; 
      egeom->contpts[i]->w = 1.0;
    }
  }
}

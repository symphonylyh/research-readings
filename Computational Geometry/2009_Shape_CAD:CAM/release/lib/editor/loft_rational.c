/************************************************************************
 *									*
			Copyright (C) 1991 by
	Massachusetts Institute of Technology, Cambridge, MA
			 All rights reserved
 *									*
 ************************************************************************/

#include <stdio.h>
#include <math.h>
#include "gen.h"
#include "bspl.h"
#include "editor.h"

#define MACHPREC 1.1102230246251567E-13	  /* changed on 25 JUN 90 */
#define MACHPR 1.1*MACHPREC

/*
void f01lbf_(int *, int *, int *, double *, int *, double *, int *, int *,
	     int *, int *);
void f04ldf_(int *, int *, int *, int *, double *, int *, double *, int *,
	     int *, double *, int *, int *);
*/
void f07bdf_(int *, int *, int *, int *, double *, int *, int *, int *);
void f07bef_(char *, int *, int *, int *, int *, double *, int *, int *,
             double *, int *, int *);


ParSurf *tol_loft_rational(ParCurv **egeom, double *U, short uorder, int nsec)

/* Form a "lofted" rational B-spline surface patch out of the nsec
   sections egeom[k], 0<=k<nsec, where each egeom[k] is a rational B-spline
   curve. All egeom's have the same order and knot vector. The parametric
   lines of the "lofting" patch at u(sub)k, 0<=k<nsec are egeom[k]. The
   parameter values u(sub)k are contained in the double array U of
   dimension nsec. uorder is the order of the "lofting surface in the
   u-direction */

{
  ParSurf *fgeom;
  double t,**al,**rhs,**aij,**N;
  int i,j,k,j1,j2,na,no,ia,il,m1,m2,tl,vorder,vpts,*in;
  double **aijn;
  int *ipiv,ldab,kl,ku;
  char trans;


/* Meaning of the most important variables (in alphabetical order)

aij      =  2-d array containing the linear system matrix (see NAG
            Manual for details).
al       =  2-d auxiliary array required by the NAG routines.
fgeom    =  lofting rational B-spline patch (returned by this routine).
ia       =  see NAG Manual (f01lbf).
il       =  see NAG Manual (f01lbf).
in       =  integer array required by f01lbf, f04ldf (see NAG Manual).
m1, m2   =  lower and upper bandwidths of the matrix of the linear
            system.
N        =  2-d array containing the B-spline basis functions evaluated
            at the nodes of eg3d's knot vector.
rhs      =  the right hand side of the linear system.
vorder   =  v-order of fgeom.
vpts     =  number of fgeom's control points in the v-direction.
*/

/* Assign the v-knot vector of fgeom. Calculate the u-knot vector. */

  vpts = egeom[0]->ncontpts;
  vorder = egeom[0]->order;
  fgeom = fgeomalloc1(uorder, vorder, nsec, vpts);
  for (i=0; i<vorder+vpts; i++)
    fgeom->vknots[i] = egeom[0]->knots[i];
  for (i=0; i<uorder; i++) {
    fgeom->uknots[i] = U[0];
    fgeom->uknots[nsec+i] = U[nsec-1];
  }
  t = (double)(uorder-1);
  for (i=uorder; i<nsec; i++) {
    fgeom->uknots[i] = 0.0;
    for (j=i-uorder+1; j<i; j++)
      fgeom->uknots[i] += U[j];
    fgeom->uknots[i] /= t;
  }

/* Compute the nodes of eg3d and the bandwidth of the Gram matrix.
   Allocate storage */

  tol_calc_bandwidth(fgeom->uknots, U, nsec, uorder, &m1, &m2) ;
  ia = MIN(nsec, m1+m2+1);
  il = MAX(1, m1);
  in = int_array1(nsec);
  aij = dbl_array2(nsec, ia);
  al = dbl_array2(nsec, il);
  
/* Build the matrix of the linear system (see NAG Manual for details) */

  N = dbl_array2(uorder+1, uorder+1);
  for (i=0; i<nsec; i++) {
    if (fabs(U[i]-1.0) < MACHPREC)
      t = U[i]-MACHPR;
    else
      t = U[i];
    tl = find(uorder+nsec, fgeom->uknots, t);
    nbasisd(tl, uorder, t, fgeom->uknots, N);
    j1 = MAX(i-m1, 0);
    j2 = MIN(i+m2,nsec);
    k = 0 ; na = uorder-tl;
    for (j=j1; j<=j2; j++) {
      no = na+j;
      if (no<0 || no>uorder)
	aij[i][k++] = 0.0; 
      else
	aij[i][k++] = N[no][uorder];
    }
  }
  free_darray2(N);
  j = 0;

  kl = m1;
  ku = m2;

  ldab = 2*kl+ku+1;
  aijn = dbl_array2((unsigned)nsec,(unsigned)ldab);
  ipiv = int_array1(nsec);
  trans = 'N';

  for (i=0;i<nsec;i++) {
    for (j=0;j<nsec;j++) {
      if (i<kl) {
	aijn[j][kl+ku+i-j]=aij[i][j];
      }
      else {
	if ((kl+j-i)>=0) {
	  aijn[j][kl+ku+i-j]=aij[i][kl+j-i];
	}
      }
    }
  }
  /*
  f01lbf_(&nsec, &m1, &m2, &aij[0][0], &ia, &al[0][0], &il, in, &i, &j);
  */

  f07bdf_(&nsec,&nsec,&kl,&ku,&aijn[0][0],&ldab,ipiv,&j);  
  if (j) {
    printf("Failure in f07bdf - called from tol_loft_rational!\n");
    return (NULL);
  }
 
/* Calculate the right hand side, consisting of 4*vpts columns */

  rhs = dbl_array2(4*vpts, nsec);
  for (i=0; i<nsec; i++) {
    k = 0;
    for (j=0; j<vpts; j++) {
      rhs[k++][i] = egeom[i]->contpts[j]->x;
      rhs[k++][i] = egeom[i]->contpts[j]->y;
      rhs[k++][i] = egeom[i]->contpts[j]->z;
      rhs[k++][i] = egeom[i]->contpts[j]->w;
    }
  }

/* Solve the linear system */

  k = 4*vpts;
  j = 0;
  /*
  f04ldf_(&nsec, &m1, &m2, &k, &aij[0][0], &ia, &al[0][0], &il, in,
	  &rhs[0][0], &nsec, &j);
  */

  f07bef_(&trans,&nsec,&kl,&ku,&k,&aijn[0][0],&ldab,ipiv,&rhs[0][0],&nsec,&j);  

  if (j) {
    printf("Failure in f07bef - called from tol_<loft_rational!\n");
    return (NULL);
  }

/* Build fgeom */

  for (i=0; i<nsec; i++) {
    k = 0;
    for (j=0; j<vpts; j++) {
      fgeom->contpts[i][j]->x = rhs[k++][i];
      fgeom->contpts[i][j]->y = rhs[k++][i];
      fgeom->contpts[i][j]->z = rhs[k++][i];
      fgeom->contpts[i][j]->w = rhs[k++][i];
    }
  }
  free_darray2(rhs);
  free_darray2(al);
  free_darray2(aij);
  free_iarray1(in);
  free_darray2(aijn);
  free_iarray1(ipiv);

  return (fgeom);
}

void tol_calc_bandwidth(double *knotv, double *chsi, int ktpts, int order,
			int *m1, int *m2)

/* Calculate the lower and upper bandwidth of the Gram matrix.
   knotv = knot vector.
   chsi  = array containing the points, where the B-spline is evaluated.
   ktpts = number of vertices ofthe control polygon.
   order = order of the B-spline.
   m1    = lower bandwidth.
   m2    = upper bandwidth.
*/
{
  int i,j,b,l,r;

  *m1 = 0;
  *m2 = 0;
  for (i=0; i<ktpts; i++) {
    l = 0;
    r = 0;
    for (j=0; j<ktpts; j++) {
      b = chsi[j]>knotv[i]+MACHPREC && chsi[j]<knotv[i+order]-MACHPREC;
      l += b && j<i;
      r += b && j>i;
    }
    if (*m1<l) *m1 = l;
    if (*m2<r) *m2 = r;
  }
}

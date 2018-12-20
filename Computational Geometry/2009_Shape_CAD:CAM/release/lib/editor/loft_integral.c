/************************************************************************
 *									*
			Copyright (C) 1992 by
	Massachusetts Institute of Technology, Cambridge, MA
			 All rights reserved
 *									*
 ************************************************************************/

#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include "gen.h"
#include "bspl.h"
#include "appr.h"
#include "editor.h"

#define MACHPREC 1.1102230246251567E-13	  /* changed on 25 JUN 90 */
#define MACHPR 1.1*MACHPREC
#define MACPREC 1.5*MACHPREC
#define MAXSECT 500
#define MAXEQ 800
#define MAXBW 150
#define MAXHB 75

/*
void f01lbf_(int *, int *, int *, double *, int *, double *, int *, int *,
	     int *, int *);
void f04ldf_(int *, int *, int *, int *, double *, int *, double *, int *,
	     int *, double *, int *, int *);
*/

void f07bdf_(int *, int *, int *, int *, double *, int *, int *, int *);
void f07bef_(char *, int *, int *, int *, int *, double *, int *, int *,
             double *, int *, int *);


ParSurf *loft_integral(ParCurv **eg, double *U, double *V, short uorder,
		       short nsec, short vsec)

/* Construct an integral B-spline surface patch by interpolating nsec
   planar curves, eg[k], 1<=k<=nsec. egeom[k] is expressed as a
   rational B-spline curve. All eg's have the same order and knot
   vector.
   The output geometry  has the same v-order and knot vector as egeoms
   and interpolates eg's at the nodes of their common knot vector.
   Its u-order is uorder and the u-knot vector is
   computed in the routine by averaging the elements of U (see below).
   The u-parameter values at the section locations are contained in the
   array U, which contains nsec elements. The parameter values of the
   data points along the circumference are contained in the array V. Their
   number is vsec.
   The routine determines the vertices of the output geometry's control
   polyhedron by solving a linear system. An optimum numbering scheme for
   the Gram matrix has been devised, which minimizes the total bandwidth.
   aij is a 2-d array containing the linear system matrix (see NAG
   Manual for details).
   al is a 2-d auxiliary array required by the NAG routines.
   rhs is the the right hand side of the linear system.
   The arrays aij,al and rhs are allocated in the FORTRAN main program and
   passed to loft_integral. If space is allocated by dbl_array2 within
   loft_integral in the usual manner the program craches!! */

{
  ParSurf *canal ;
  vector *v;
  double **c,**z,*chsi,*zeta,*knotu,*knotv,d,**aij,**al,**rhs ;
  int kspan,lspan,i,j,k,l,m,n,e1,e2,*kb,*lb,n0,nu,nv,no,mo,ia,il,*in,ij ;
  int n3,n6,k3,m6,m3,tk,tl,i0,j0,ii,jj,e20,nno,mmo,jin,e,b1,b2,m1,m2,neq ;
  double **aijn;
  int *ipiv,ldab,kl,ku;
  char trans;
  

/* Meaning of the most important local variables (in alphabetical order)

b1       =  total bandwidth of the Gram matrix.
c        =  2-D array containing the B-spline basis functions evaluated
            at the u- or v-parameter values.
canal    =  approximating canal surface (the output geometry).
chsi     =  u-parameter values at the data points.
ia       =  see NAG Manual (f01lbf).
il       =  see NAG Manual (f01lbf).
in       =  integer array required by f01lbf, f04ldf (see NAG Manual).
kb       =  integer array of dimension 2. 
            kb[0] = maximum number of data points belonging to the
                    support of the basis u-B-spline function corresponding
                    to one data point, d, which have a u-parameter value
                    smaller than d.
            kb[1] = maximum number of data points belonging to the
                    support of the basis u-B-spline function corresponding
                    to one data point, d, which have a u-parameter value
                    greater than d.
knotu    =  knot vector along the "faster increasing" direction.
knotv    =  knot vector along the "slower increasing" direction.
lb       =  integer array of dimension 2. 
            lb[0] = maximum number of data points belonging to the
                    support of the basis v-B-spline function corresponding
                    to one data point, d, which have a v-parameter value
                    smaller than d.
            lb[1] = maximum number of data points belonging to the
                    support of the basis v-B-spline function corresponding
                    to one data point, d, which have a v-parameter value
                    greater than d.
m        =  number of vertices in the "slower increasing" direction.
m1, m2   =  lower and upper bandwidths of the matrix of the linear
            system.
mo       =  order in the "slower" increasing direction.
mmo      =  number of knot elements in the "slower increasing" direction.
n        =  number of vertices in the "faster increasing" direction.
n0       =  flag indicating vertex and data points storage sequence.
            0 : first index (u-vertices) increases faster.
            1 : second index (v-vertices) increases faster.
neq      =  rank of the Gram matrix.
nno      =  number of knot elements in the "faster increasing" direction.
no       =  order in the "faster" increasing direction.
nu       =  number of u-vertices of the control polyhedron.
nv       =  number of v-vertices of the control polyhedron.
z        =  2-D array containing the B-spline basis functions evaluated
            at the v- or u-parameter values.
zeta     =  v-parameter values at the data points (nodes of eg).

*/

  kb = int_array1(2);
  lb = int_array1(2) ;
  chsi = dbl_array1(MAXSECT);
  zeta = dbl_array1(MAXSECT);

/* Calculate canal's knot vectors */

  nu = nsec;
  nv = vsec;
  mo = uorder;
  no = eg[0]->order;
  canal = fgeomalloc1(mo, no, nu, nv);
  for (i=0; i<mo; i++) {
    canal->uknots[i] = U[0];
    canal->uknots[nu+i] = U[nu-1];
  }
  d = (double)(mo-1);
  for (i=mo; i<nu; i++) {
    canal->uknots[i] = 0.0;
    for (j=i-mo+1; j<i; j++)
      canal->uknots[i] += U[j];
    canal->uknots[i] /= d;
  }
  for (i=0; i<no; i++) {
    canal->vknots[i] = V[0];
    canal->vknots[nv+i] = V[nv-1] ;
  }
  d = (double)(no-1);
  for (i=no; i<nv; i++) {
    canal->vknots[i] = 0.0;
    for (j=i-no+1; j<i; j++)
      canal->vknots[i] += V[j];
    canal->vknots[i] /= d;
  }
  for (i=0; i<nu; i++)
    chsi[i] = U[i];
  d = canal->uknots[mo+nu-1];
  for (i=0; i<nu; i++)
    if (fabs(d-chsi[i]) < MACHPREC)
      chsi[i] -= MACHPR;
  for (i=0; i<nv; i++)
    zeta[i] = V[i];
  e = eg[0]->ncontpts;
  d = eg[0]->knots[no+e-1];
  for (i=0; i<nv; i++)
    if (fabs(d-zeta[i]) < MACHPREC)
      zeta[i] -= MACHPR ;

/* Establish optimum numbering sequence to obtain minimum bandwidth */

  intcyl_calc_bounds_offs(canal, chsi, zeta, kb, lb);
  kspan = kb[0]+kb[1];
  lspan = lb[0]+lb[1];
  b1 = nv*kspan+lspan+1;
  b2 = nu*lspan+kspan+1;
  n0 = b1<b2;
  n = (n0)? nv : nu;
  m = (n0)? nu : nv;
  i = (nu<nv)? nv : nu;
  if (n0-1) {
    swap_vectors(chsi, zeta, i);
    swap_ivectors(kb, lb, 2);
    b1 = b2;
  }
  if (b1 > MAXBW) {
    printf("loft_integral: maximum permissible bandwidth exceeded!\n");
    printf("               computed bandwidth = %d\n",b1);
    return (NULL);
  }
  no = (n0) ? canal->vorder : canal->uorder;
  nno = n+no;
  mo = (n0) ? canal->uorder : canal->vorder;
  mmo = m+mo;
  knotu = dbl_array1(nno);
  knotv = dbl_array1(mmo);
  for (i=0; i<nno; i++)
    knotu[i]=(n0) ? canal->vknots[i] : canal->uknots[i];
  for (i=0; i<mmo; i++)
    knotv[i]=(n0) ? canal->uknots[i] : canal->vknots[i];

/* Compute the bandwidth. Allocate storage. */

  n6 = n;
  n3 = n;
  m6 = m;
  m1 = kb[0]*n6+lb[0];
  m2 = kb[1]*n6+lb[1];
  neq = n6*m6;
  m3 = m;
  free_iarray1(lb);
  free_iarray1(kb);
  if (neq > MAXEQ) {
    printf("loft_integral: maximum permissible linear system rank exceeded!\n");
    printf("               computed value = %d\n\n",neq) ;
    return (NULL);
  }
  if (m1 > MAXHB) {
    printf("loft_integral: too many subdiagonals..computed value = %d\n",m1);
    return (NULL);
  }
  i0 = mo+1;
  c = dbl_array2(i0, i0);
  i0 = no+1;
  z = dbl_array2(i0, i0);
  
  ia = b1;
  il = m1;
  aij = dbl_array2(neq, ia);
  al  = dbl_array2(neq, il);
  rhs = dbl_array2(3, neq);

/* Compute the Gram matrix and the rhs of the linear system*/

  for (i=0; i<neq; i++)
    for (j=0; j<b1; j++)
      aij[i][j] = 0.0;
  for (k=0; k<m3; k++) {
    k3 = n6*k;
    tk = find(mmo, knotv, chsi[k]);
    i0 = tk-mo;
    nbasisd(tk, mo, chsi[k], knotv,c);
    e20 = n6*i0;
    ij = -1;
    while (c[++ij][mo] < MACHPREC);
    for (l=0; l<n3; l++) {
      tl = find(nno, knotu, zeta[l]);
      nbasisd(tl, no, zeta[l], knotu, z);
      j0 = tl-no;
      e1 = k3+l;
      jin = -1;
      while (z[++jin][no] < MACHPREC)
	;

/* Compute data points */

      v = (n0) ? rbspeval(eg[k], zeta[l],0) : rbspeval(eg[l], chsi[k],0);

/* Compute the non-zero elements of the Gram matrix */

      i = ij;
      while (i<=mo && c[i][mo] > MACHPREC) {
	ii = i+i0;
	e2 = e20+n6*i;
	j = jin ;
	if (ii>=0 && ii<m6)
	  while (j<=no && z[j][no] > MACHPREC) {
	    jj = j0+j;
	    if (jj>=0 && jj<n6) {
	      e = (e1>m1)? e2+jj+m1-e1 : e2+jj;
	      aij[e1][e] = c[i][mo]*z[j][no];
	    }
	    j++;
	  }
	i++;
      }
      d = v->w;
      rhs[0][e1]=v->x/d;
      rhs[1][e1]=v->y/d;
      rhs[2][e1]=v->z/d;
      vectfree(v);
    }
  }

/* Solve the linear system */

  free_darray2(z);
  free_darray2(c);
  free_darray1(zeta);
  free_darray1(chsi);
  free_darray1(knotu);
  free_darray1(knotv);
  in = int_array1(neq);
  k = 3;
  e = l = 0;
  e1 = neq;

  /*
  f01lbf_(&neq, &m1, &m2, &aij[0][0], &ia, &al[0][0], &il, in, &l, &e);
  if (e) {
    printf("failure in f01lbf - called from loft_integral!\n") ;
    return (NULL);
  }
  f04ldf_(&neq, &m1, &m2, &k, &aij[0][0], &ia, &al[0][0], &il, in,
	  &rhs[0][0], &e1, &e);
  if (e) {
    printf("failure in f04ldf - called from loft_integral!\n") ;
    return (NULL);
  }
  */

  kl = m1;
  ku = m2;
  ldab = 2*kl+ku+1;
  aijn = dbl_array2((unsigned)neq,(unsigned)ldab);
  ipiv = int_array1(neq);
  trans = 'N';

   for (i=0;i<neq;i++) {
    for (j=0;j<neq;j++) {
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
  
  f07bdf_(&neq,&neq,&kl,&ku,&aijn[0][0],&ldab,ipiv,&e);
  if (e) {
    printf(" Failure in f07bdf - called from loft_integral!\n");
    return (NULL);
  }  
  f07bef_(&trans,&neq,&kl,&ku,&k,&aijn[0][0],&ldab,ipiv,&rhs[0][0],&neq,&e);  
  if (e) {
    printf(" Failure in f07bef - called from loft_integral!\n");
    return (NULL);
  }   
 


  for (k=0; k<m3; k++) {
    k3 = k*n6;
    for (l=0; l<n3; l++) {
      e = k3+l;
      e1 = (n0)? k : l;
      e2 = (n0)? l : k;
      j = 0; 
      canal->contpts[e1][e2]->x = rhs[j++][e];
      canal->contpts[e1][e2]->y = rhs[j++][e];
      canal->contpts[e1][e2]->z = rhs[j++][e];
      canal->contpts[e1][e2]->w = 1.0;
    }
  }
  free_darray2(aij);
  free_darray2(al);
  free_darray2(rhs);
  free_iarray1(in);
  free_iarray1(ipiv);
  free_darray2(aijn);
  return(canal);
}

   
void intcyl_calc_bounds_offs(ParSurf *fgeom, double *chsi, double *zeta,
			     int *kb, int *lb)

/* Calculate the maximum number of nodes affected by one u-vertex and by
   one v-vertex. The routine computes the vectors kb,lb. For the meaning
   of the input variables see above (canal_to_int). */

{
  int i,j,k,l,m1,m2;
  
  kb[0] = kb[1] = lb[0] = lb[1] = 0;
  l = fgeom->ucontpts;
  for (i=0; i<l; i++) {
    j = 0;
    m1 = m2 = 0;
    while (chsi[j]<fgeom->uknots[i]+MACPREC &&
	   (chsi[j]>MACHPREC || i>=fgeom->uorder))
      j++;
    k = i+fgeom->uorder;
    while (chsi[j]<chsi[i]-MACPREC) {
      if (chsi[j]>MACHPREC && i!=j)
	m1++;
      j++;
    }
    while (chsi[j]<fgeom->uknots[k]-MACPREC && j<l) { 
      if (chsi[j]>MACHPREC && i!=j)
	m2++;
      j++;
    }
    if (kb[0]<m1)
      kb[0] = m1;
    if (kb[1]<m2)
      kb[1] = m2;
  }
  
  l = fgeom->vcontpts ;
  for (i=0; i<l; i++) {
    j = 0;
    m1 = m2 = 0; 
    while (zeta[j]<fgeom->vknots[i]+MACPREC &&
	   (zeta[j]>MACHPREC || i>=fgeom->vorder))
      j++ ;
    k = i+fgeom->vorder; 
    while (zeta[j]<zeta[i]-MACPREC) {
      if (zeta[j]>MACHPREC && i!=j)
	m1++;
      j++;
    }
    while (zeta[j]<fgeom->vknots[k]-MACPREC && j<l) {
      if (zeta[j]>MACHPREC && i!=j)
	m2++;
      j++;
    }
    if (lb[0]<m1)
      lb[0] = m1;
    if (lb[1]<m2)
      lb[1] = m2;
  }
}

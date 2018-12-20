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

void alloc_fgeompts(ParSurf *);
/*
void f01lbf_(int *, int *, int *, double *, int *, double *, int *,
	     int *, int *, int *);
void f04ldf_(int *, int *, int *, int *, double *, int *, double *,
	     int *, int *, double *, int *, int *);
*/

void f07bdf_(int *, int *, int *, int *, double *, int *, int *, int *);
void f07bef_(char *, int *, int *, int *, int *, double *, int *, int *,
             double *, int *, int *);


#define MACHPREC 1.0e-13
#define MACPREC 1.5*MACHPREC
#define MACHPR 1.1*MACHPREC

#define MAXSECT 500
#define MAXBW 150
#define MAXEQ 800
#define MAXHB 75

ParSurf *integral_offset_surf(ParSurf *fgeom, double dist, double **knv,
			      int *nk, double eps, int *noff)

/* Build the offset to an integral B-spline patch. The offset is an
   integral B-spline with the same orders as the progenitor.
   dist  is the value of the offset distance (positive or negative).
   knv[][0]  is the initial u-knot vector, knv[][1]  the initial v-knot
   vector.
   nk[0]  is the number of u-knots, nk[1]  is the number of v-knots.
   eps  is the maximum allowable position difference between
   approximating and true surfaces.
   noff[0], noff[1]  are the number of sampling points on the
   approximating and the true offset surfaces of fgeom in the u and v
   direction, respectively.
*/
{
  ParSurf *offset;
  double maxError;
  short i, ku, kv, again, uo, vo, up, vp;

/* Meaning of the most important local variables (in alphabetical order)

ku       =  number of u-knots of fgeom - 1.
kv       =  number of v-knots of fgeom - 1.
offset   =  the offset surface, returned by this routine.
again    =  logical variable. If true, iterate to build approximating
            offset surface with more data points.
uo       =  u-order of fgeom
up       =  number of u-control polyhedron vertices of offset.
vo       =  v-order of fgeom
vp       =  number of v-control polyhedron vertices of offset.
*/

  uo = fgeom->uorder;
  vo = fgeom->vorder;
  up = nk[0] - uo;
  vp = nk[1] - vo;
  ku = uo + fgeom->ucontpts - 1;
  kv = vo + fgeom->vcontpts - 1;
  offset = fgeomalloc2(fgeom->uorder, fgeom->vorder, MAXSECT, MAXSECT);
  offset->uorder = uo;
  offset->vorder = vo;
  offset->ucontpts = up;
  offset->vcontpts = vp;
  for (i=0; i<nk[0]; i++)
    offset->uknots[i] = knv[i][0];
  for (i=0; i<nk[1]; i++)
    offset->vknots[i] = knv[i][1];

  again = 1;
  while (again) {
    alloc_fgeompts(offset);

    if (integral_build_offset(offset, fgeom, dist))
      return (NULL);

    if ((again = integral_sample_offset(fgeom, offset, noff[0], noff[1],
					eps, dist, &maxError)) < 0)
      return (NULL);

    printf("Error is %g; %2d points added\n", maxError, again);
    fflush(stdout);
  }
  return(offset) ;
}

int integral_build_offset(ParSurf *offset, ParSurf *fgeom, double dist)

/* Calculate the control polygon vertices of offset.
   fgeom  is the progenitor (integral B-spline patch).
   dist  is the offset distance (positive or negative).
   An optimum numbering scheme for the Gram matrix has been
   devised, which minimizes the total bandwidth.
   aij  is a 2-d array containing the linear system matrix (see NAG
   Manual for details).
   al  is a 2-d auxiliary array required by the NAG routines.
   rhs  is the the right hand side of the linear system.
*/
{
  vector *nrm, *v, *vup;
  double **c,**z,*chsi,*zeta,*knotu,*knotv,d,dn,r2,**aij,**al,**rhs;
  int kspan,lspan,i,j,k,l,m,n,e1,e2,*kb,*lb,n0,nu,nv,no,mo,ia,il,*in,ij;
  int n3,n6,k3,m6,m3,tk,tl,i0,j0,ii,jj,e20,nno,mmo,jin,e,b1,b2,m1,m2,neq;
  double **aijn;
  int *ipiv,ldab,kl,ku;
  char trans;

/* Meaning of the most important local variables (in alphabetical order)

b1       =  total bandwidth of the Gram matrix.
c        =  2-D array containing the B-spline basis functions evaluated
            at the u- or v-parameter values.
chsi     =  u-parameter values at the data points. chsi and zeta (see
            below) can be swapped to achieve a smaller bandwidth.
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
nrm      =  normal vector to fgeom at the data points (scaled by dist).
nno      =  number of knot elements in the "faster increasing" direction.
no       =  order in the "faster" increasing direction.
nrm      =  normal vector to fgeom at a data point.
nu       =  number of u-vertices of the control polyhedron.
nv       =  number of v-vertices of the control polyhedron.
vup      =  position vector of one data point lying on the offset.
z        =  2-D array containing the B-spline basis functions evaluated
            at the v- or u-parameter values.
zeta     =  v-parameter values at the data points. chsi (see above) and zeta
            can be swapped to achieve a smaller bandwidth.
*/

  kb = int_array1(2);
  lb = int_array1(2);
  chsi = dbl_array1(MAXSECT);
  zeta = dbl_array1(MAXSECT);

/* Establish optimum numbering sequence to obtain minimum bandwidth */

  nu = offset->ucontpts;
  nv = offset->vcontpts;
  mo = offset->uorder;
  nodes_bsp(offset->uknots, nu, mo, chsi);
  d = offset->uknots[mo+nu-1];
  for (i=0; i<nu; i++)
    if (fabs(d - chsi[i]) < MACHPREC)
      chsi[i] -= MACHPR;
  no = offset->vorder;
  nodes_bsp(offset->vknots, nv ,no, zeta);
  d = offset->vknots[no+nv-1];
  for (i=0; i<nv; i++)
    if (fabs(d - zeta[i]) < MACHPREC)
      zeta[i] -= MACHPR;
  integral_calc_bounds_offs(offset, chsi, zeta, kb, lb);
  kspan = kb[0] + kb[1];
  lspan = lb[0] + lb[1];
  b1 = nv*kspan+lspan+1;
  b2 = nu*lspan+kspan+1;
  n0 = b1<b2;
  n = (n0) ? nv : nu;
  m = (n0) ? nu : nv;
  i = (nu<nv) ? nv : nu;
  if (n0-1) {
    swap_vectors(chsi, zeta, i);
    swap_ivectors(kb, lb, 2);
    b1 = b2;
  }
  if (b1 > MAXBW)
    return (1);

  no = (n0)? offset->vorder : offset->uorder;
  nno = n + no;
  mo = (n0)? offset->uorder : offset->vorder;
  mmo = m + mo;
  knotu = dbl_array1(nno);
  knotv = dbl_array1(mmo);
  for (i=0; i<nno; i++)
    knotu[i] = (n0) ? offset->vknots[i] : offset->uknots[i];
  for (i=0; i<mmo; i++)
    knotv[i] = (n0) ? offset->uknots[i] : offset->vknots[i];
  
  /* Compute the bandwidth. Allocate storage. */

  n6 = n;
  n3 = n;
  m6 = m;
  m1 = kb[0]*n6+lb[0];
  m2 = kb[1]*n6+lb[1];
  neq = n6*m6;
  m3 = m;
  free_iarray1(lb);
  free_iarray1(kb) ;
  if (neq > MAXEQ)
    return (2);

  if (m1 > MAXHB) {
    printf("integral_build_offset: too many subdiagonals!\n");
    printf("                       computed value = %d\n", m1) ;

    return (1);
  }
  i0 = mo+1;
  c = dbl_array2(i0, i0);
  i0 = no+1;
  z = dbl_array2(i0, i0);
  
/* allocate memory */

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
    if (tk > m-1)
      tk = m - 1;     /* takes care of bug in find() */
    i0 = tk-mo;
    nbasisd(tk, mo, chsi[k], knotv, c);
    e20 = n6*i0 ;
    ij = -1;
    while (c[++ij][mo] < MACHPREC)
      ;
    for (l=0; l<n3; l++) {
      tl = find(nno, knotu, zeta[l]);
      if (tl > n-1)
	tl = n - 1;   /* takes care of bug in find() */
      nbasisd(tl, no, zeta[l], knotu, z);
      j0 = tl-no;
      e1 = k3+l;
      jin = -1;
      while (z[++jin][no] < MACHPREC)
	;

/* Compute data points */

      nrm = (n0) ? normalsurf(fgeom, chsi[k], zeta[l]) :
	normalsurf(fgeom, zeta[l], chsi[k]);
      scale_vect1(dist, nrm, nrm) ;
      v = (n0) ? revalderivsurf(fgeom, chsi[k], zeta[l], 0, 0) :
	revalderivsurf(fgeom, zeta[l], chsi[k], 0, 0);
      vup = add_vect(v, nrm);
      vectfree(v);
      vectfree(nrm);
      
/* Compute the non-zero elements of the Gram matrix */

      i = ij;
      while (i<=mo && c[i][mo]>MACHPREC) {
	ii = i+i0;
	e2 = e20+n6*i;
	j = jin;
	if (ii>=0 && ii<m6)
	  while (j<=no && z[j][no]>MACHPREC) {
	  jj = j0+j ; if (jj>=0 && jj<n6) {
	    e = (e1>m1)? e2+jj+m1-e1 : e2+jj;
	    aij[e1][e] = c[i][mo]*z[j][no];
	  }
	  j++;
        }
	i++;
      }
      rhs[0][e1] = vup->x;
      rhs[1][e1] = vup->y;
      rhs[2][e1] = vup->z;
      vectfree(vup);
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
  e = -1;
  l = 0;
  e1 = neq;

  kl = m1;
  ku = m2;
  ldab = 2*kl+ku+1;
  aijn = dbl_array2((unsigned)neq,(unsigned)ldab);
  ipiv = int_array1(neq);
  trans = 'N';
  /*
  f01lbf_(&neq, &m1, &m2, &aij[0][0], &ia, &al[0][0], &il, in, &l, &e);
  */

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
  if (e)
    return (3);

  f07bef_(&trans,&neq,&kl,&ku,&k,&aijn[0][0],&ldab,ipiv,&rhs[0][0],&neq,&e);  
  /*
  f04ldf_(&neq, &m1, &m2, &k, &aij[0][0], &ia, &al[0][0], &il, in, &rhs[0][0],
	  &e1, &e);
  */
  if (e)
    return (4);
  
  for (k=0; k<m3; k++) {
    k3 = k*n6;
    for (l=0; l<n3; l++) {
      e = k3+l;
      e1 = (n0) ? k : l;
      e2 = (n0) ? l : k;
      j = 0 ; 
      offset->contpts[e1][e2]->x = rhs[j++][e];
      offset->contpts[e1][e2]->y = rhs[j++][e];
      offset->contpts[e1][e2]->z = rhs[j++][e];
      offset->contpts[e1][e2]->w = 1.0;
    }
  }
  free_darray2(aij);
  free_darray2(al);
  free_darray2(rhs);
  free_iarray1(in) ;
  free_iarray1(ipiv) ; 
  free_darray2(aijn);   
  
  return (0);
}

int integral_sample_offset(ParSurf *fgeom, ParSurf *offset, short nu, short nv,
			   double eps, double dist, double *maxError)

/* This routine samples uniformly the approximating offset surface.
   The position difference between isoparametric points of the true and
   the approximating offsets is computed and compared against eps.
   fgeom  is the progenitor surface and dist the offset distance.
   nu,nv  are the number of sampled points in the u,v direction,
   respectively.
   eps  is the maximum allowable position difference between
   isoparametric points of the progenitor and the offset.
   The routine updates the knot vectors of offset.
   The total number of added points is returned. */
 
{
  vector *v1,*nrm,*vq;
  double2 d,du,dv,u,v,*kn;
  int i,j,k,l,b,tu,tv,uo,vo,up,vp,un,vn,b0,**ll;

/* Meaning of the most important local variables (in alphabetical order)

b0       =  number of added data points, value returned by this routine.
du       =  sampling step along the u-direction.
dv       =  sampling step along the v-direction.
kn       =  the refined u- or v-knot vector of offset.
ll       =  2-dimensional array.
            ll[i][j] = 0 : no failure.
            ll[i][j] = 1 : one added point in the rectangular interval
                           (u[i],u[i+1])x(v[j],v[j+1]) ,
            where u,v are the u- and v-knot vectors.
tu       =  index indicating the relative position of u (see below) with
            respect to par[0][.]
tv       =  index indicating the relative position of v (see below) with
            respect to par[1][.]
u        =  u-parameter value of a sampled point.
v        =  v-parameter value of a sampled point.
vq       =  position difference between isoparametric points of
            the true and the approximating offsets.

*/

  du = (fgeom->uknots[fgeom->ucontpts] - fgeom->uknots[fgeom->uorder-1]) /
       (double2)(nu-1);
  dv = (fgeom->vknots[fgeom->vcontpts] - fgeom->vknots[fgeom->vorder-1]) /
       (double2)(nv-1);
  b0 = 0;
  uo = fgeom->uorder;
  vo = fgeom->vorder;
  up = offset->ucontpts;
  vp = offset->vcontpts;
  un = up+uo;
  vn = vp+vo;
  ll = int_array2(MAXSECT, MAXSECT);
  *maxError = 0.0;
  
/* Sample uniformly */

  for (i=0; i<up; i++)
    for (j=0; j<vp; j++)
      ll[i][j] = 0;

  u = fgeom->uknots[fgeom->uorder-1];
  for (i=0; i<nu; i++) {
    if (fabs(1.0-u) < MACHPREC)
      u -= MACHPR ;
    tu = find(un, offset->uknots, u);
    v = fgeom->vknots[fgeom->vorder-1];
    for (j=0; j<nv; j++) {
      if (fabs(1.0-v) < MACHPREC) v -= MACHPR;
      tv = find(vn, offset->vknots, v);

/* Compute position vectors to the offsets at (u,v) */

      if (!ll[tu][tv]) {
	nrm = normalsurf(fgeom, u, v);
	scale_vect1(dist, nrm, nrm);
	v1 = revalderivsurf(fgeom, u, v, 0, 0);
	vq = add_vect(v1, nrm);
	vectfree(v1);
	vectfree(nrm) ;
	
/* Compute position difference */

	if (!ll[tu][tv]) {
	  v1 = revalderivsurf(offset, u, v, 0, 0);
	  d = distance(v1, vq);
	  ll[tu][tv] = d>eps;
	  b0 += ll[tu][tv];
	  if (*maxError < d)
	    *maxError = d;
	  vectfree(v1) ;
	}
	vectfree(vq);
      }
      v += dv;
    }
    u += du;
  }
  
/* Compute new knot vectors of offset. */

  if (b0) {
    kn = dbl_array1(MAXSECT);

    k = uo-1;
    for (i=0; i<uo-1; i++)
      kn[i] = offset->uknots[i];
    for (i=uo-1; i<up; i++) {
      kn[k++] = offset->uknots[i];
      if (k >= MAXSECT)
	return (-1);

      l = 1;
      for (j=vo-1; j<vp; j++)
	if (l && ll[i][j]) {
	  kn[k++] = (offset->uknots[i]+offset->uknots[i+1])/2.0;
	  l = 0 ;
	  if (k >= MAXSECT)
	    return (-2);
	}
    }
    offset->ucontpts = k;
    for (i=0; i<k; i++)
      offset->uknots[i] = kn[i];
    for (i=k; i<k+uo; i++)
      offset->uknots[i] = fgeom->uknots[fgeom->ucontpts];

    k = vo-1;
    for (j=0; j<vo-1; j++)
      kn[j] = offset->vknots[j];
    for (j=vo-1; j<vp; j++) {
      kn[k++] = offset->vknots[j];
      if (k >= MAXSECT)
	return (-3);

      l = 1;
      for (i=uo-1; i<up; i++)
	if (l && ll[i][j]) {
	  kn[k++] = (offset->vknots[j]+offset->vknots[j+1])/2.0;
	  l = 0;
	  if (k >= MAXSECT)
	    return (-4);
	}
    }
    offset->vcontpts = k;
    for (j=0; j<k; j++)
      offset->vknots[j] = kn[j];
    for (j=k; j<k+vo; j++)
      offset->vknots[j] = fgeom->vknots[fgeom->vcontpts];

    free_darray1(kn); 
  }
  free_iarray2(ll);
  
  return (b0);
}

void integral_calc_bounds_offs(ParSurf *fgeom, double *chsi, double *zeta,
			       int *kb, int *lb)

/* Calculate the maximum number of nodes affected by one u-vertex and by
   one v-vertex. The routine computes the vectors kb,lb. For the meaning
   of the input variables see above (build_offset_...). */

{
  int i, j, k, l, m1, m2;

  kb[0] = kb[1] = lb[0] = lb[1] = 0;
  l = fgeom->ucontpts;
  for (i=0; i<l; i++) {
    j = 0;
    m1 = m2 = 0;
    while ((chsi[j] < fgeom->uknots[i] + MACPREC) &&
	   (chsi[j] > MACHPREC || i >= fgeom->uorder))
      j++;
    k = i + fgeom->uorder;
    while (chsi[j] < chsi[i]-MACPREC) {
      if (chsi[j] > MACHPREC && i!=j)
	m1++;
      j++;
    }
    while (chsi[j] < fgeom->uknots[k] - MACPREC && j < l) { 
      if (chsi[j] > MACHPREC && i != j)
	m2++;
      j++ ;
    }
    if (kb[0]<m1) kb[0] = m1;
    if (kb[1]<m2) kb[1] = m2;
  }
  
  l = fgeom->vcontpts;
  for (i=0; i<l; i++) {
    j = 0;
    m1 = m2 = 0; 
    while ((zeta[j] < fgeom->vknots[i] + MACPREC) &&
	   (zeta[j] > MACHPREC || i >= fgeom->vorder))
      j++;
    k = i+fgeom->vorder; 
    while (zeta[j] < zeta[i]-MACPREC) {
      if (zeta[j]>MACHPREC && i!=j)
	m1++;
      j++;
    }
    while (zeta[j] < fgeom->vknots[k]-MACPREC && j<l) {
      if (zeta[j] > MACHPREC && i != j)
	m2++;
      j++;
    }
    if (lb[0] < m1) lb[0] = m1;
    if (lb[1] < m2) lb[1] = m2;
  }
}

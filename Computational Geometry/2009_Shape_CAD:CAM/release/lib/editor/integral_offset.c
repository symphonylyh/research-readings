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

#define MAXEQ 800
#define MAXSECT 500
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

short integral_offset(ParSurf *plate, ParSurf **offset, double tz,
		      int *nsec, double eps, int *noff)

/* Build the tolerance region around a deformed plate, represented by
   the rational B-spline surface patch plate.
   tz is the value of the tolerance.
   offset are the approximating B-splines to plate's offsets.
   offset[0] = "upper" offset (in the direction to plate's normal).
   offset[1] = "lower" offset (in a direction opposite to that
   of plate's normal).
   nsec[0] is the initial number of uniformly distributed internal knots
   in the u direction, nsec(1) the corresponding number for the v direction.
   eps is the maximum allowable position difference between
   approximating and true surfaces.
   noff[0], noff[1] are the number of sampling points on the
   approximating and the true offset surfaces of plate in the u and v
   direction, respectively. */

{
  double du,dv;
  int i,l,s,uo,vo,up,vp,ku,kv;

/* Meaning of the most important local variables (in alphabetical order)
du       =  uniform u-knots' spacing.
dv       =  uniform v-knots' spacing.
l        =  logical variable. If true, iterate to build approximating
            offset surface with more data points.
s        =  loop control variable.
            s = 0 characterizes a v = ct. edge.
            s = 1 characterizes a u = ct. edge.
uo       =  u-order of plate
up       =  number of u-control polyhedron vertices of offset.
vo       =  v-order of plate
vp       =  number of v-control polyhedron vertices of offset.
*/

/* Allocate memory for offsets. */

  uo = plate->uorder; 
  vo = plate->vorder;
  up = uo+nsec[0];
  vp = vo+nsec[1];
  du = (plate->uknots[plate->ucontpts]-plate->uknots[plate->uorder-1])/
       (double)(nsec[0]+1);
  dv = (plate->vknots[plate->vcontpts]-plate->vknots[plate->vorder-1])/
       (double)(nsec[1]+1);
  for (s=0; s<2; s++) {
    offset[s] = fgeomalloc2(plate->uorder, plate->vorder, MAXSECT, MAXSECT);
    offset[s]->uorder = uo;
    offset[s]->vorder = vo;
    offset[s]->ucontpts = up;
    offset[s]->vcontpts = vp;
    for (i=0; i<uo; i++) {
      offset[s]->uknots[i] = plate->uknots[plate->uorder-1];
      offset[s]->uknots[i+up] = plate->uknots[plate->ucontpts];
    }
    for (i=uo; i<up; i++)
      offset[s]->uknots[i] = offset[s]->uknots[i-1]+du;
    for (i=0; i<vo; i++) {
      offset[s]->vknots[i] = plate->vknots[plate->vorder-1];
      offset[s]->vknots[i+vp] = plate->vknots[plate->vcontpts];
    }
    for (i=vo; i<vp; i++)
      offset[s]->vknots[i] = offset[s]->vknots[i-1]+dv;
  }
  
/* Compute approximating offset surfaces. Iterate until convergence. */

  printf(" Compute and sample approximating offset surfaces\n");
  l = 1;
  while (l) {
    for (s=0; s<2; s++)
      alloc_fgeompts(offset[s]);
    if (build_offsets_int(offset, plate, tz))
      return (1);

    if ((l = tol_sample_offsets(plate,offset,noff[0],noff[1],eps,tz)) == -1)
      return (1);
  }

  return (0);
}

short build_offsets_int(ParSurf **offset, ParSurf *plate, double tz)

/* Builds the upper and lower offset surfaces to plate
   (offset[0],offset[1], respectively). Upper means in the direction of the
   vector normal to plate, lower in a direction opposite to this.
   plate is a rational B-spline surface patch representing the deformed
   plate.
   pars is a two-dimensional array of dimensions 2 x MAXSECT containing
   the parametric values at the data points.
   pars[0][.] = u-values at the data points.
   pars[1][.] = v-values at the data points.
   tz is the value of the tolerance.
   It is assumed that the knot vectors and the first perimetrical row of
   vertices of offset's control polyhedron are known. The routine
   determines the remaining vertices of the control polyhedron by solving a
   linear system. An optimum numbering scheme for the Gram matrix has been
   devised, which minimizes the total bandwidth. */

{
  vector *ru,*rv,*v,*nrm,*vup,*vdw;
  double **c,**z,*chsi,*zeta,*knotu,*knotv,d,dn,r2;
  double **aij, **al, **rhs;
  double **aijn;
  int kspan,lspan,i,j,k,l,m,n,e1,e2,*kb,*lb,n0,nu,nv,no,mo,ia,il,*in,ij;
  int n3,n6,k3,m6,m3,tk,tl,i0,j0,ii,jj,e20,nno,mmo,jin,e,b1,b2,m1,m2,neq;
  int *ipiv,ldab,kl,ku;
  char trans;

/* aij is a 2-d array containing the linear system matrix (see NAG
   Manual for details).
   al is a 2-d auxiliary array required by the NAG routines.
   rhs is the the right hand side of the linear system. */

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
nrm      =  normal vector to plate at the data points (scaled by tz).
nno      =  number of knot elements in the "faster increasing" direction.
no       =  order in the "faster" increasing direction.
nrm      =  normal vector to plate at a data point.
nu       =  number of u-vertices of the control polyhedron.
nv       =  number of v-vertices of the control polyhedron.
vdw      =  position vector of one data point lying on the lower offset.
vup      =  position vector of one data point lying on the upper offset.
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

  nu = offset[0]->ucontpts;
  nv = offset[0]->vcontpts;
  mo = offset[0]->uorder;
  nodes_bsp(offset[0]->uknots, nu, mo, chsi);
  d = offset[0]->uknots[mo+nu-1];
  for (i=0; i<nu; i++)
    if (fabs(d-chsi[i])<MACHPREC)
      chsi[i] -= MACHPR;
  no = offset[0]->vorder;
  nodes_bsp(offset[0]->vknots, nv, no, zeta);
  d = offset[0]->vknots[no+nv-1];
  for (i=0; i<nv; i++)
    if (fabs(d-zeta[i])<MACHPREC)
      zeta[i] -= MACHPR;
  calc_bounds_offs(offset[0], chsi, zeta, kb, lb);
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
  no = (n0) ? offset[0]->vorder : offset[0]->uorder;
  nno = n+no;
  mo = (n0) ? offset[0]->uorder : offset[0]->vorder;
  mmo = m+mo;
  knotu = dbl_array1((unsigned)nno);
  knotv = dbl_array1((unsigned)mmo);
  for (i=0; i<nno; i++)
    knotu[i]=(n0) ? offset[0]->vknots[i] : offset[0]->uknots[i];
  for (i=0; i<mmo; i++)
    knotv[i] = (n0) ? offset[0]->uknots[i] : offset[0]->vknots[i];

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
    printf(" build_offsets_int: maximum permissible linear system rank exceeded!");
    printf("                    computed value = %d\n", neq);
    return (1);
  }
  
  i0 = mo+1;
  c = dbl_array2((unsigned)i0, (unsigned)i0);
  i0 = no+1;
  z = dbl_array2((unsigned)i0, (unsigned)i0);

/* allocate memory for aij, al, and rhs */

  ia = b1;
  il = m1;

  aij = dbl_array2((unsigned)neq, (unsigned)ia);
  al  = dbl_array2((unsigned)neq, (unsigned)il);
  rhs = dbl_array2((unsigned)6, (unsigned)neq);

/* Compute the Gram matrix and the rhs of the linear system*/
  
  for (i=0; i<neq; i++)
    for (j=0; j<b1; j++)
      aij[i][j] = 0.0;
  for (k=0; k<m3; k++) {
    k3 = n6*k;
    tk = find(mmo, knotv, chsi[k]);
    i0 = tk-mo;
    nbasisd(tk, mo, chsi[k], knotv, c);
    e20 = n6*i0;
    ij = -1;
    while (c[++ij][mo]<MACHPREC);
    for (l=0; l<n3; l++) {
      tl = find(nno,knotu, zeta[l]);
      nbasisd(tl, no, zeta[l], knotu, z);
      j0 = tl-no;
      e1 = k3+l;
      jin = -1;
      while (z[++jin][no]<MACHPREC);

/* Compute data points */

      nrm = (n0) ? normalsurf(plate, chsi[k], zeta[l]) :
	normalsurf(plate, zeta[l], chsi[k]);
      scale_vect1(tz, nrm, nrm);
      v = (n0) ? revalderivsurf(plate, chsi[k], zeta[l], 0, 0) :
	revalderivsurf(plate, zeta[l], chsi[k], 0, 0);
      vup = add_vect(v,nrm);
      vdw = sub_vect(v,nrm);
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
	    jj = j0+j;
	    if (jj>=0 && jj<n6)
	      e = (e1>m1)? e2+jj+m1-e1 : e2+jj; aij[e1][e] = c[i][mo]*z[j][no];
	    j++;
	  }
	i++;
      }
      rhs[0][e1] = vup->x; rhs[1][e1] = vup->y; rhs[2][e1] = vup->z;
      rhs[3][e1] = vdw->x; rhs[4][e1] = vdw->y; rhs[5][e1] = vdw->z;
      vectfree(vup);
      vectfree(vdw);
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
  k = 6;
  e = l = 0;

  e1 = neq;

  kl = m1;
  ku = m2;
  ldab = 2*kl+ku+1;
  aijn = dbl_array2((unsigned)neq,(unsigned)ldab);
  ipiv = int_array1(neq);
  trans = 'N';

  /*
  f01lbf_(&neq,&m1,&m2,&aij[0][0],&ia,&al[0][0],&il,in,&l,&e);
  if (e) {
    printf(" Failure in f01lbf - called from build_offsets!\n");
    return (1);
  }
  f04ldf_(&neq,&m1,&m2,&k,&aij[0][0],&ia,&al[0][0],&il,in,&rhs[0][0],&e1,&e);
  if (e) {
    printf(" Failure in f04ldf - called from build_offsets!\n");
    return (1);
  }
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
  if (e) {
    printf(" Failure in f07bdf - called from build_offsets!\n");
    return (1);
  }  
  f07bef_(&trans,&neq,&kl,&ku,&k,&aijn[0][0],&ldab,ipiv,&rhs[0][0],&neq,&e);  
  if (e) {
    printf(" Failure in f07bef - called from build_offsets!\n");
    return (1);
  }   

  for (k=0; k<m3; k++) {
    k3 = k*n6;
    for (l=0; l<n3; l++) {
      e = k3+l;
      e1 = (n0) ? k : l;
      e2 = (n0) ? l : k;
      j = 0; 
      for (i=0; i<2; i++) {
	offset[i]->contpts[e1][e2]->x = rhs[j++][e];
	offset[i]->contpts[e1][e2]->y = rhs[j++][e];
	offset[i]->contpts[e1][e2]->z = rhs[j++][e];
	offset[i]->contpts[e1][e2]->w = 1.0;
      }
    }
  }
  free_iarray1(in);
  free_iarray1(ipiv);
  free_darray2(aij);
  free_darray2(aijn);
  free_darray2(al);
  free_darray2(rhs);

  return (0);
}

short tol_sample_offsets(ParSurf *plate, ParSurf **offset, int nu, int nv,
			 double eps, double tz)

/* This routine samples uniformly the approximating offset surfaces to
   plate, offset[0] and offset[1]. offset[0] is the "upper" offset, i.e. in
   the direction of plate's normal, offset[1] is the "lower" offset, i.e.
   in the opposite direction.
   The position difference between isoparametric points of the true and
   the approximating offsets is computed and compared against eps.
   plate is the progenitor surface (deformed plate geometry) and tz the
   offset distance (tolerance).
   nu,nv are the number of sampled points in the u,v direction,
   respectively.
   eps is the maximum allowable position difference between
   isoparametric points of the progenitor and the offsets.
   The routine updates the knot vectors of offset.
   The total number of added points is returned. */
 
{
 vector *v1,*nrm,**vq;
 double du,dv,u,v,*kn;
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
vdw      =  position vector to a point lying on the "lower" offset.
vq       =  vector array.
            vq[0] = position difference between isoparametric points of
                    the true and the approximating "upper" offsets.
            vq[1] = position difference between isoparametric points of
                    the true and the approximating "lower" offsets.
vup      =  position vector to a point lying on the "upper" offset.
*/

 du = (plate->uknots[plate->ucontpts]-plate->uknots[plate->uorder-1])/
      (double)(nu-1);
 dv = (plate->vknots[plate->vcontpts]-plate->vknots[plate->vorder-1])/
      (double)(nv-1);
 u = plate->uknots[plate->uorder-1];
 b0 = 0;
 uo = plate->uorder;
 vo = plate->vorder;
 up = offset[0]->ucontpts;
 vp = offset[0]->vcontpts;
 un = up+uo; vn = vp+vo;
 ll = int_array2(MAXSECT, MAXSECT);
 vq = vec_array1(2);

/* Sample uniformly */

 for (i=0; i<up; i++)
   for (j=0; j<vp; j++)
     ll[i][j] = 0;
 for (i=0; i<nu; i++) {
   if (fabs(1.0-u)<MACHPREC)
     u -= MACHPR;
   tu = find(un, offset[0]->uknots, u);
   v = plate->vknots[plate->vorder-1];
   for (j=0; j<nv; j++) {
     if (fabs(1.0-v)<MACHPREC)
       v -= MACHPR;
     tv = find(vn, offset[0]->vknots, v);

/* Compute position vectors to the offsets at (u,v) */

     if (!ll[tu][tv]) {
       nrm = normalsurf(plate, u, v);
       scale_vect1(tz, nrm, nrm);
       v1 = revalderivsurf(plate, u, v, 0, 0);
       vq[0] = add_vect(v1, nrm);
       vq[1] = sub_vect(v1, nrm);
       vectfree(v1);
       vectfree(nrm);

/* Compute position difference */

       for (b=0; b<2; b++)
	 if (!ll[tu][tv]) {
	   v1 = revalderivsurf(offset[b], u, v, 0, 0);
	   ll[tu][tv] = distance(v1, vq[b])>eps;
	   b0 += ll[tu][tv];
	   vectfree(v1);
	   vectfree(vq[b]);
	 }
     }
     v += dv;
   }
   u += du;
 }

/* Compute new knot vectors of offset. */

 free_varray1(vq, 2);
 if (b0) {
   kn = dbl_array1(MAXSECT);
   k = uo-1;
   for (i=0; i<uo-1; i++)
     kn[i] = offset[0]->uknots[i];
   for (i=uo-1; i<up; i++) {
     kn[k++] = offset[0]->uknots[i];
     if (k>=MAXSECT) {
       printf(" Max no of u ctpts exceeded (1) in sample_offsets!\n");
       return (-1);
     }
     l = 1;
     for (j=vo-1; j<vp; j++)
       if (l && ll[i][j]) {
	 kn[k++] = (offset[0]->uknots[i]+offset[0]->uknots[i+1])/2.0;
	 l = 0;
	 if (k>=MAXSECT) {
	   printf(" Max no of u ctpts exceeded (2) in sample_offsets!\n");
	   return (-1);
	 }
       }
   }
   for (b=0; b<2; b++) {
     offset[b]->ucontpts = k;
     for (i=0; i<k; i++)
       offset[b]->uknots[i] = kn[i];
     for (i=k; i<k+uo; i++)
       offset[b]->uknots[i] = plate->uknots[plate->ucontpts];
   }
   k = vo-1;
   for (j=0; j<vo-1; j++)
     kn[j] = offset[0]->vknots[j];
   for (j=vo-1; j<vp; j++) {
     kn[k++] = offset[0]->vknots[j];
     if (k>=MAXSECT) {
       printf(" Max no of v ctpts exceeded (1) in sample_offsets!\n");
       return (-1);
     }
     l = 1;
     for (i=uo-1; i<up; i++)
       if (l && ll[i][j]) {
	 kn[k++] = (offset[0]->vknots[j]+offset[0]->vknots[j+1])/2.0;
	 l = 0;
	 if (k>=MAXSECT) {
	   printf(" Max no of v ctpts exceeded (2) in sample_offsets!\n");
	   return (-1);
	 }
       }
   }
   for (b=0; b<2; b++) {
     offset[b]->vcontpts = k;
     for (j=0; j<k; j++)
       offset[b]->vknots[j] = kn[j];
     for (j=k; j<k+vo; j++)
       offset[b]->vknots[j] = plate->vknots[plate->vcontpts];
   }
   free_darray1(kn); 
 }
 free_iarray2(ll); 

 return (b0);
}

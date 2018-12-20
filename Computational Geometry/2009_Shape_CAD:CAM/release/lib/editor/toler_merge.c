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
#include "editor.h"

#define MAXSECT 500
#define MACHPREC 1.1102230246251567E-13	  /* changed on 25 JUN 90 */
#define MACHPR 1.1*MACHPREC
#define MACPREC 1.5*MACHPREC

short merge_can_off(ParSurf **canals, ParSurf **offset)

/* Merge the knot vectors of the lateral canal surfaces and the offset
   surfaces to the plate geometry.
   canals is an array of ParSurf containing the canal surfaces.
   canals[0] is the canal surface with spine the line v=0 of the
             deformed plate geometry.
   canals[1] is the canal surface with spine the line v=1 of the
             deformed plate geometry.
   canals[2] is the canal surface with spine the line u=0 of the
             deformed plate geometry.
   canals[3] is the canal surface with spine the line u=1 of the
             deformed plate geometry.
   offset[0] is the "upper" offset, offset[1] the lower offset to the
   deformed plate geometry. */

{
  ParSurf *fg,*ofg;
  ParCurv **eg;
  double *kn,ep;
  int i,k,s,j,p,q,r,*index;

/* Meaning of the most important local variables (in alphabetical order)

eg       =  array comprising two adjacent edges of offset and canals.
fg       =  auxiliary face accepting the refined u-knot vector of
            canals.
index    =  auxiliary space required by merge_knotv.
k        =  total number of knots resulted after merging of knot
            vectors.
kn       =  knot vector resulting after merging operations.
ofg      =  auxiliary face accepting the refined u- and v-knot vectors of
            offset.
q        =  2*s (see below).
s        =  loop variable indicating either a canal surface having spine
            the edge v=ct. or u=ct. of the original plate geometry or the
            upper, lower offsets.

*/

/* Allocate working space */

  eg = (ParCurv **)gen_array1(2,sizeof(ParCurv *));
  ofg = fgeomalloc2(offset[0]->uorder,offset[0]->vorder, MAXSECT, MAXSECT);
  kn = dbl_array1(2*MAXSECT);
  index = int_array1(2);
  ofg->uorder = offset[0]->uorder;
  ofg->vorder = offset[0]->vorder;
  j = MAX(offset[0]->uorder, offset[0]->vorder);
  for (i=0; i<2; i++)
    eg[i] = egeomalloc2((short)j, MAXSECT);
  ep = (double)MACHPREC;
  
/* Merge offset and canal geometries */

  for (s=0; s<2; s++) {
    p = 1-s;
    q = 2*s;

/* Merge knot vectors */

    extract_edge(eg[0], canals[q], 1, 0);
    extract_edge(eg[1], offset[0], p, 0);
    k = merge_knotv(eg, 0, 2, kn, index, ep);
    if (k>MAXSECT) {
      printf(" Max no of sections exceeded in merge_can_off!\n");
      return (1);
    }

/* Apply Oslo algorithm to canal surfaces */

    for (r=0; r<2; r++) {
      j = q+r;
      p = k-canals[j]->uorder;
      fg=fgeomalloc1(canals[j]->uorder, canals[j]->vorder, (short)p,
		     canals[j]->vcontpts);
      for (i=0; i<k; i++)
	fg->uknots[i] = kn[i];
      p = canals[j]->vorder+canals[j]->vcontpts;
      for (i=0; i<p; i++)
	fg->vknots[i] = canals[j]->vknots[i];
      surfoslo3(canals[j], fg, 4);
      free_fgeom(canals[j]);
      canals[j] = copyfgeom(fg, 0);
      free_fgeom(fg);
    }

/* Assign new knot vector to auxiliary face ofg */
  
    if (s) {
      ofg->vcontpts = k-offset[0]->vorder;
      for (i=0; i<k; i++)
	ofg->vknots[i] = kn[i];
    }
    else   {
      ofg->ucontpts = k-offset[0]->uorder;
      for (i=0; i<k; i++) 
	ofg->uknots[i] = kn[i];
    }
  }

/* Apply Oslo algorithm to offsets */

  free_earray1(eg, 2);
  free_darray1(kn);
  free_iarray1(index);
  alloc_fgeompts(ofg);
  for (s=0; s<2; s++) {
    surfoslo3(offset[s], ofg, 3);
    free_fgeom(offset[s]);
    for (i=0; i<ofg->ucontpts; i++)
      for (j=0; j<ofg->vcontpts; j++) 
	ofg->contpts[i][j]->w = 1.0;
    offset[s] = copyfgeom(ofg, 0);
  }
  free_fgeom(ofg);

  return (0);
}

void merge_tol_edges(ParSurf *plate, ParSurf **canals, ParSurf **offset,
		     double tz)

/* Merge the long edges of the tolerance region. For each edge a new
   control polygon having as vertices the midpoints of the segments
   connecting corresponding vertices of canals, offset is built.
   plate is the deformed plate geometry.
   canals is an array of ParSurf containing the canal surfaces.
   canals[0] is the canal surface with spine the line v=0 of the
             deformed plate geometry.
   canals[1] is the canal surface with spine the line v=1 of the
             deformed plate geometry.
   canals[2] is the canal surface with spine the line u=0 of the
             deformed plate geometry.
   canals[3] is the canal surface with spine the line u=1 of the
             deformed plate geometry.
   offset[0] is the "upper" offset, offset[1] the lower offset to the
   deformed plate geometry.
   tz is the value of the tolerance. */

{
  vector *v0,*v1,*v2,*v3;
  double h;
  int i,cm,ud,*p,pu,pv;

/* Meaning of the most important local variables (in alphabetical order)

cm       =  2*s + uv (see below).
            cm = 0 corresponds to the edge v = 0 .
            cm = 1 corresponds to the edge v = 1 .
            cm = 2 corresponds to the edge u = 0 .
            cm = 3 corresponds to the edge u = 1 .
p        =  integer array.
            p[0] = 0.
            p[1] = number of v-control points of canals (-1).
pu       =  number of u-control points of offset (-1).
pv       =  number of v-control points of offset (-1).
ud       =  loop index.
            ud=0 corresponds to the upper offset.
            ud=1 corresponds to the upper offset.
v0       =  control polygon vertex of the edge resulting after the
            merging operation.
v1,v2,v3 =  edge control polygon vertices.

*/

  p = int_array1(2);
  v0 = vectalloc();
  p[0] = 0;
  p[1] = canals[0]->vcontpts-1;
  h = (double)0.5;
  pu = offset[0]->ucontpts-1;
  pv = offset[0]->vcontpts-1;
  
/* Loop over the long edges of the tolerance region */

  for (cm=0; cm<4; cm++)
    for (ud=0; ud<2; ud++)
      for (i=1; i<canals[cm]->ucontpts-1; i++) {
	v1 = canals[cm]->contpts[i][p[ud]];
	switch (cm) {
	case 0 :
	  v2 = offset[ud]->contpts[i][0];
	  break;
	case 1 :
	  v2 = offset[ud]->contpts[i][pv];
	  break;
	case 2 :
	  v2 = offset[ud]->contpts[0][i];
	  break;
	case 3 :
	  v2 = offset[ud]->contpts[pu][i];
	  break;
	}
	add_vect1(v1, v2, v0);
	scale_vect1(h, v0, v0);
	v1 = copyvector(v0, v1);
	v2 = copyvector(v0, v2);
      }
  
/* Compute averages of the corner vertices */
  
  h = (double)(1.0/3.0);
  for (ud=0; ud<2; ud++) {
    v1 = canals[0]->contpts[0][p[ud]];
    v2 = canals[2]->contpts[0][p[ud]];
    v3 = offset[ud]->contpts[0][0];
    add_vect1(v1, v2, v0);
    add_vect1(v3, v0, v0);
    scale_vect1(h, v0, v0);
    v1 = copyvector(v0, v1);
    v2 = copyvector(v0, v2);
    v3 = copyvector(v0, v3);
    v1 = canals[0]->contpts[pu][p[ud]];
    v2 = canals[3]->contpts[0][p[ud]];
    v3 = offset[ud]->contpts[pu][0];
    add_vect1(v1, v2, v0);
    add_vect1(v3, v0, v0);
    scale_vect1(h, v0, v0);
    v1 = copyvector(v0, v1);
    v2 = copyvector(v0, v2);
    v3 = copyvector(v0, v3);
    v1 = canals[1]->contpts[0][p[ud]];
    v2 = canals[2]->contpts[pv][p[ud]];
    v3 = offset[ud]->contpts[0][pv];
    add_vect1(v1, v2, v0);
    add_vect1(v3, v0, v0);
    scale_vect1(h, v0, v0);
    v1 = copyvector(v0, v1);
    v2 = copyvector(v0, v2);
    v3 = copyvector(v0, v3);
    v1 = canals[1]->contpts[pu][p[ud]];
    v2 = canals[3]->contpts[pv][p[ud]];
    v3 = offset[ud]->contpts[pu][pv];
    add_vect1(v1, v2, v0);
    add_vect1(v3, v0, v0);
    scale_vect1(h, v0, v0);
    v1 = copyvector(v0, v1);
    v2 = copyvector(v0, v2);
    v3 = copyvector(v0, v3);
  }
  vectfree(v0);
  free_iarray1(p);
}

short sample_mod_surf(ParSurf *plate, ParSurf **canals, ParSurf **offset,
		      int *nm, double eps, double tz)

/* Sample uniformly exact and approximating surfaces along the first
   nontrivial span of the knot vector in a direction normal to the long
   edges of the tolerance region. This routine is called after merging of
   the above edges by merge_tol_edges. The routine returns flase (0) if eps
   is exceeded at any sampled point.
   plate is the deformed plate geometry.
   canals is an array of ParSurf containing the canal surfaces.
   canals[0] is the canal surface with spine the line v=0 of the
             deformed plate geometry.
   canals[1] is the canal surface with spine the line v=1 of the
             deformed plate geometry.
   canals[2] is the canal surface with spine the line u=0 of the
             deformed plate geometry.
   canals[3] is the canal surface with spine the line u=1 of the
             deformed plate geometry.
   offset[0] is the "upper" offset, offset[1] the lower offset to the
   deformed plate geometry.
   nm is an array containing the number of sampling points.
   nm[0] is the number of sampling points along the edge of the
         tolerance region.
   nm[1] is the number of sampling points in a direction normal to the
         edge of the tolerance region.
   eps is the maximum allowable positional error between exact and
   approximating surfaces.
   tz is the value of the tolerance. */

{
  ParCurv *eg;
  vector *v1,*v2,*nrm;
  double u,v,du,dv,ue,ve,uoe,voe,u0,v0,t0,*d,dm0,dm1,f,uo0,vo0;
  int s,cm,ud,uv,i,j,p,q,l,uo,vo,up,vp,ou,ov,nm0,nm1;

/* Meaning of the most important local variables (in alphabetical order)

cm       =  2*s + uv (see below).
            cm = 0 corresponds to the edge v = 0 .
            cm = 1 corresponds to the edge v = 1 .
            cm = 2 corresponds to the edge u = 0 .
            cm = 3 corresponds to the edge u = 1 .
d        =  auxiliary array containing the sampling steps in the u and
            the v-direction of the canal surfaces.
dm0      =  nm[0]-1 converted to a double.
dm1      =  nm[1]-1 converted to a double.
du       =  sampling step in the u-direction of offset.
dv       =  sampling step in the v-direction of offset.
eg       =  the generatrix geometry.
l        =  value returned by this routine. True if eps is not exceeded
            at any sampled point, false otherwise.
nrm      =  unit normal vector to plate.
ou       =  last index of the u-knot vector of offset.
ov       =  last index of the v-knot vector of offset.
p        =  number of v-knots of canals (-1).
q        =  v-order of canals. Ocasionally used as a loop index.
s        =  loop index.
            s=0 corresponds to a v=ct. edge of plate, offset.
            s=1 corresponds to a u=ct. edge of plate, offset.
u        =  u-parameter value of a sampling point.
u0       =  u-parameter value of the first sampling point of canals.
ud       =  loop index.
            ud=0 corresponds to the upper offset.
            ud=1 corresponds to the lower offset.
ue       =  u-parameter value of the last sampling point of canals.
uo       =  u-order of plate, offset.
uoe      =  u-parameter value of the last sampling point of offset.
up       =  number of u-control points (-1) of offset.
uv       =  loop index.
            uv=0 corresponds to a u=0 or v=0 edge.
            uv=1 corresponds to a u=1 or v=1 edge.
v        =  v-parameter value of a sampling point.
v0       =  v-parameter value of the first sampling point of canals.
ve       =  v-parameter value of the last sampling point of canals.
voe      =  v-parameter value of the last sampling point of offset.
vo       =  v-order of plate, offset.
vp       =  number of v-control points (-1) of offset.

*/

  q = canals[0]->vorder;
  p = canals[0]->vcontpts-1;
  l = 1;
  uo = plate->uorder;
  vo = plate->vorder;
  up = offset[0]->ucontpts-1;
  vp = offset[0]->vcontpts-1;
  ou = uo+up; ov = vo+vp;
  u0 = canals[0]->uknots[canals[0]->uorder-1];
  v0 = canals[0]->vknots[canals[0]->vorder-1];
  ve = canals[0]->vknots[p+q];
  uo0 = offset[0]->uknots[offset[0]->uorder-1];
  vo0 = offset[0]->vknots[offset[0]->vorder-1];
  dm0 = (double)(nm[0]-1);
  dm1 = (double)(nm[1]-1);
  d = dbl_array1(2);
  d[0] = (canals[0]->vknots[q]-v0)/dm1;
  d[1] = (ve-canals[0]->vknots[p])/dm1;

/* Sample canal surfaces */

  for (s=0; s<2; s++)
    if (l)
      for (uv=0; uv<2; uv++)
	if (l) {
	  cm = 2*s+uv; 
	  ue = canals[cm]->uknots[canals[cm]->uorder+canals[cm]->ucontpts-1];
	  du = (ue-u0)/dm0;
	  u = u0;
	  for (i=0; i<nm[0]; i++) {
	    for (q=0; q<2; q++) if (l) {
	      dv = d[q]; v = (q) ? ve : v0;
	      for (j=0; j<nm[1]; j++) if (l) {
		if (fabs(ue-u)<MACHPREC) u -= MACHPR;
		if (fabs(ve-v)<MACHPREC) v -= MACHPR;
		if (fabs(u-u0)<MACPREC) u = u0;
		if (fabs(v-v0)<MACPREC) v = v0;
	  
/* Compute position vector on the generalized cylinder */

		t0 = (double)uv;
		eg = generatrix(plate, u, t0, s, tz);
		v1 = rbspeval(eg, v, 0);
		free_egeom(eg);
		v2 = revalderivsurf(canals[cm], u, v, 0,0);
		l = distance(v1, v2)<eps;
		vectfree(v2);
		vectfree(v1);
		if (q)
		  v -= dv;
		else
		  v += dv;
	      }
	    }
	    u += du;
	  }
    
/* Sample offset surfaces */

	  switch (cm) {
	  case 0 :
	    uoe = offset[0]->uknots[ou];
	    voe = offset[0]->vknots[vo];
	    du = (uoe-offset[0]->uknots[offset[0]->uorder-1])/dm0;
	    dv = (voe-offset[0]->vknots[offset[0]->vorder-1])/dm1;
	    u = offset[0]->uknots[offset[0]->uorder-1];
	    nm0 = nm[0];
	    nm1 = nm[1];
	    break;
	  case 1 :
	    uoe = offset[0]->uknots[ou];
	    voe = offset[0]->vknots[ov];
	    du = (uoe-offset[0]->uknots[offset[0]->uorder-1])/dm0;
	    dv = (voe-offset[0]->vknots[vp])/dm1;
	    u = offset[0]->uknots[offset[0]->uorder-1];
	    nm0 = nm[0];
	    nm1 = nm[1];
	    break;
	  case 2 :
	    uoe = offset[0]->uknots[uo];
	    voe = offset[0]->vknots[ov];
	    du = (uoe-offset[0]->uknots[offset[0]->uorder-1])/dm1;
	    dv = (voe-offset[0]->vknots[offset[0]->vorder-1])/dm0;
	    u = offset[0]->uknots[offset[0]->uorder-1];
	    nm0 = nm[1];
	    nm1 = nm[0];
	    break;
	  case 3 :
	    uoe = offset[0]->uknots[ou];
	    voe = offset[0]->vknots[ov];
	    du = (uoe-offset[0]->uknots[up])/dm1;
	    dv = (voe-offset[0]->vknots[offset[0]->vorder-1])/dm0;
	    u = uoe; nm0 = nm[1]; 
	    nm1 = nm[0];
	    break;
	  }
	  for (i=0;i<nm0;i++) if (l) {
	    switch (cm) {
	    case 0 :
	      v = offset[0]->vknots[offset[0]->vorder-1]; 
	      break;
	    case 1 :
	      v = voe; 
	      break;
	    case 2 :
	      v = offset[0]->vknots[offset[0]->vorder-1];
	      break;
	    case 3 :
	      v = offset[0]->vknots[offset[0]->vorder-1];
	      break;
	    }
	    for (j=0;j<nm1;j++) if (l) {
	      if (fabs(uoe-u)<MACHPREC) u -= MACHPR;
	      if (fabs(voe-v)<MACHPREC) v -= MACHPR;
	      if (fabs(u-uo0)<MACPREC) u = uo0;
	      if (fabs(v-vo0)<MACPREC) v = vo0;
	      for (ud=0; ud<2; ud++)
		if (l) {
		  nrm = normalsurf(plate, u, v);
		  f = tz*(double)(1-2*ud);
		  scale_vect1(f, nrm, nrm);
		  v1 = revalderivsurf(plate, u, v, 0, 0);
		  add_vect1(nrm, v1, v1);
		  v2 = revalderivsurf(offset[ud], u, v, 0, 0);
		  l = distance(v1,v2)<eps;
		  vectfree(v2);
		  vectfree(v1);
		}
	      if (cm==1)
		v -= dv;
	      else
		v += dv;
	    }
	    if (cm==3)
	      u -= du;
	    else
	      u += du;
	  }
	} 
  free_darray1(d);
  
  return(l);
}

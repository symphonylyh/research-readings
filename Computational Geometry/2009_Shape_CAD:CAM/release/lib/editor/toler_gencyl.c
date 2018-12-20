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

ParSurf *canal_to_rat(ParSurf *plate, double t0, int ed, int ns, int nc,
		      double eps, double *pars, ParCurv **gen, int *m,
		      int m0, double tz)

/* Approximate a generalized cylinder by a rational B-spline surface
   patch.
   plate is a rational B-spline surface patch representing the deformed
   plate.
   The spine is a parametric line of plate with u=t0 if
   ed=1, v=t0 if ed=0. The output geometry interpolates the generalized
   cylinder at m sections, positioned at parameter values contained in the
   array pars (m and pars are computed in this routine).
   gen is a 2-d array containing the interpolated sections (generatrices).
   The iterative aproximation procedure starts with m0 equidistant
   sections between consecutive interpolated sections.
   ns is the number of sampled sections between consecutive interpolated
   sections.
   nc is the number of sampling points on each section (generatrix).
   eps is the maximum allowable positional error between the
   approximating B-spline and the generalized cylinder.
   tz is the value of the tolerance. */

{
  ParSurf *canal;
  double t,d;
  int i,order;

/* Meaning of the most important local variables (in alphabetical order)

canal    =  The output geometry (rational B-spline surface patch).

*/

  *m = m0;
  d = (ed) ? plate->vknots[plate->vcontpts]-plate->vknots[plate->vorder-1] :
             plate->uknots[plate->ucontpts]-plate->uknots[plate->uorder-1];
  order = (ed) ? plate->vorder : plate->uorder;
  d /= (double)(*m-1);
  if (m0<order) {
    printf(" Number of sampling points smaller than spine order!\n");
    return (NULL);
  }
  pars[0] = (ed) ? plate->vknots[plate->vorder-1] :
                   plate->uknots[plate->uorder-1];
  for (i=0; i<*m; i++) {
    if (i)
      pars[i] = pars[i-1]+d;
    gen[i] = generatrix(plate, pars[i], t0, ed, tz);
  }
  
/* Iterate as required for convergence */

  printf("\n Iterate until convergence");
  fflush(stdout);
  do {
    if (*m>MAXSECT) {
      printf(" Maximum number of sections exceeded in canal_to_rat!\n");
      return (NULL);
    }
    if (*m != m0) free_fgeom(canal);
    if ((canal = tol_loft_rational(gen, pars, order, *m)) == NULL)
      return (NULL);
  } while (sample_gencyl(plate, t0, ed, ns, nc, eps, canal, gen, pars, m, tz));
  
  return(canal);
}

short sample_gencyl(ParSurf *plate, double t0, int ed, int ns, int nc,
		    double eps, ParSurf *canal, ParCurv **gen, double *pars,
		    int *m, double tz)

/* Sample uniformly the generalized cylinder. canal is the approximating
   rational B-spline patch, gen are the plane sections (generatrices) at m
   parameter values of the spine, which are contained in pars. The routine
   returns the number of new sections required for better approximation
   (one new section midway between existing, consecutive sections for each
   failed interval). pars, m and gen are updated.
   plate,t0,ed,ns,nc,eps,tz have the same meaning as in canal_to_rat. */

{
  ParCurv *eg;
  vector *v1,*v2;
  double ds,d,f,f0,df,mg;
  int i,j,k,uo,up,vo,vp,*l,b,b0;

/* Meaning of the most important local variables (in alphabetical order)

b0       =  total number of failed intervals (=number of new sections).
eg       =  the generatrix.
l        =  integer array. For each interval :
            1 if some point failed between gen[i] and gen[i-1].
            0 otherwise.
vo       =  v-order of plate.
vp       =  number of v-vertices of plate's control polyhedron.
uo       =  u-order of plate.
up       =  number of u-vertices of plate's control polyhedron.

*/

  uo = plate->uorder;
  vo = plate->vorder;
  l = int_array1(MAXSECT);
  up = plate->ucontpts;
  vp = plate->vcontpts;
  b0 = 0;
  f0 = (double)(nc-1);

/* Do for each interval between consecutive interpolated sections */

  for (i=1;i<*m;i++) {
    ds = (pars[i]-pars[i-1])/(double)(ns+1);
    l[i] = 0; d = pars[i-1];

/* Generate new sections */
  
    for (j=0; j<ns; j++)
      if (!l[i]) {
	d += ds;
	if (i>1 || j)
	  free_egeom(eg);
	eg = generatrix(plate, d, t0, ed, tz);
	df = (eg->knots[eg->order+eg->ncontpts-1]-eg->knots[0])/f0;
	f = eg->knots[0];
	if (fabs(1.0-d)<MACHPREC)
	  d -= MACHPR;

/* Sample uniformly for each new section */

      for (k=0;k<nc; k++)
	if (!l[i]) {
	  if (fabs(1.0-f)<MACHPREC)
	    f -= MACHPR;
	  v1 = rbspeval(eg,f,0);
	  v2 = revalderivsurf(canal, d, f, 0, 0);
	  sub_vect1(v1, v2, v1);
	  vectfree(v2);
	  mg = mag(v1);
	  vectfree(v1);
	  l[i] = l[i] || (mg>eps);
	  f += df;
	}
      }  
    b0 += l[i];
  }
  
/* Generate one section for each failed interval. Rearrange indices. */

  b = b0;
  free_egeom(eg);
  if (b)
    for (i = *m-1; i>0; i--)
      if (b) {
	j = i+b; gen[j] = copyegeom(gen[i],gen[j]);
	if (l[i]) {
	  b--;
	  j--;
	  if (j<*m)
	    free_egeom(gen[j]);
	  d = (pars[i]+pars[i-1])/2.0;
	  gen[j] = generatrix(plate, d, t0, ed, tz);
	}
      }

/* Compute new pars */

  b = b0;
  if (b)
    for (i = *m-1; i>0;i--)
      if (b) {
	j = i+b;
	pars[j] = pars[i];
	if (l[i]) {
	  b--;
	  j--;
	  pars[j] = (pars[i]+pars[i-1])/2.0;
	}
      }
  free_iarray1(l);
  *m += b0;

  return(b0);
}

short merge_canals(ParSurf *plate, ParSurf **canals, ParCurv ***gen,
		   double **pars, double **par, int *m, int *mp, double tz)

/* Add generatrices to opposite pairs of canal surfaces, such that each
   pair is composed of an equal number of interpolated sections
   (generatrices), which are located at the same parameter values along the
   respective spines.
   plate is the deformed plate geometry.
   canals[0] = canal surface along the edge v=0.
   canals[1] = canal surface along the edge v=1.
   canals[2] = canal surface along the edge u=0.
   canals[3] = canal surface along the edge u=1.
   gen[i][.] is a two-dimensional array containing the sections
   (generatrices) of canals[i].
   pars[i][.] is a two-dimensional array containing the parameter values
   of canals[i]'s spine corresponding to the location of interpolated
   sections.
   m[i] contains the number of elements of pars[i][.].
   par[0][.] contains the parameter values resulting after merging
   pars[0][.] and pars[1][.] .
   par[1][.] contains the parameter values resulting after merging
   pars[2][.] and pars[3][.] .
   mp[i] is the number of elements of par[i].
   tz is the value of the tolerance (radius of the sweeped circular
   sections).
   All input variables are updated by this routine (except plate,tz). */

{
  double d,t0;
  int i,j,k,n,i1,i2,p;

/* Meaning of the most important local variables (in alphabetical order)
p        =  0 corresponds to the edges u = ct.
            1 corresponds to the edges v = ct.
t0       =  the constant u or v value along one edge of plate.
*/

/* Merge opposite pairs of parameter values (pars's) */

  for (i=0; i<4; i += 2) {
    p = i/2;
    n = i+1;
    i1 = i2 = k = 0;
    while (i1<m[i] && i2<m[n]) {
      if (fabs(pars[i][i1]-pars[n][i2])<MACHPREC) {
	par[p][k++] = pars[i][i1];
	i1++;
	i2++;
      }
      else
	par[p][k++]=(pars[i][i1]<pars[n][i2]) ? pars[i][i1++] : pars[n][i2++];
      if (k>=MAXSECT) {
	printf(" Max no of sections exceeded in merge_canals!\n");
	return (1);
      }
    }
    if (i1<m[i])
      for (j=i1; j<m[i]; j++) {
	par[p][k++] = pars[i][i1++];
	if (k>=MAXSECT) {
	  printf(" Max no of sections exceeded in merge_canals!\n");
	  return (1);
	}
      }
    else if (i2<m[n])
      for (j=i2; j<m[n]; j++) {
      par[p][k++] = pars[n][i2++];
      if (k>=MAXSECT) {
	printf(" Max no of sections exceeded in merge_canals!\n");
	return (1);
      }
    }
    mp[p] = k;
  }

/* Calculate the approximating B-splines to the canal surfaces */
  
  for (n=0; n<4; n++) {
    t0 = (double)(n%2);
    p = n/2;
    if (addsections_canal(plate, canals[n], gen[n], pars[n], &m[n], par[p],
			  mp[p], t0, p, tz))
      return (1);
  }
  
  return (0);
}

short addsections_canal(ParSurf *plate, ParSurf *canal, ParCurv **gen,
			double *pars, int *m, double *par, int m0,
			double t0, int ed, double tz)

/* Add new sections (generatrices) to achieve a closer approximation of a
   generalized cylinder.
   plate is the underlying deformed plate geometry and tz the tolerance
   value (radius of the generatrix).
   pars contains the locations of the existing sections along the spine.
   The sections (generatrices) are described by gen. 
   par contains the new locations, which are to be added to pars.
   m, m0 are the number of elements of pars, par, respectively.
   canal is the rational B-spline approximating the generalized
   cylinder. It is updated so that the added sections are also
   interpolated.
   pars and m are also updated. */

{
  double d;
  int i,j,k,b,b0,*l,**ll,order;

/* Meaning of the most important local variables (in alphabetical order)
b0       =  total number of added sections.
l        =  integer array indicating the number of new sections in each
            parameter interval (pars[j],pars[j+1])
ll       =  pointer to par used to position the new section within pars
            properly.
*/

/* Compute pointers to new section locations */

  l = int_array1(MAXSECT);
  ll = int_array2(MAXSECT, MAXSECT);
  b0 = 0;
  for (i=1; i<*m; i++) {
    for (j=0; j<m0; j++)
      ll[j][i] = 0;
    l[i] = 0;
    for (j=0; j<m0; j++) {
      l[i] += par[j]>pars[i-1]+MACHPREC && par[j]<pars[i]-MACHPREC;
      k = l[i]-1;
      if (l[i] && !ll[k][i])
	ll[k][i] = j+1;
    }
    b0 += l[i];
  }
  
/* Compute new sections (generatrices) */

  b = b0;
  if (b) {
    for (i = *m-1;i>0;i--) if (b) {
      j = i+b;
      gen[j] = copyegeom(gen[i], gen[j]); 
      pars[j] = pars[i];
      if (l[i])
	for (k=0; k<l[i]; k++) {
	  b--;
	  j--;
	  if (j<*m)
	    free_egeom(gen[j]);
	  d=par[ll[k][i]-1];
	  pars[j]=d;
	  gen[j]=generatrix(plate, d, t0, ed, tz);
	}
    }
    *m += b0;
    if(*m>MAXSECT) {
      printf(" Max no of sections exceeded in addsections_canal!\n");
      return (1);
    }

/* Interpolate sections to build a rational B-spline */

    free_fgeom(canal);
    order = canal->uorder;
    canal = tol_loft_rational(gen, pars, order, *m);
  }
  free_iarray2(ll);
  free_iarray1(l);
  
  return (0);
}

ParCurv *generatrix(ParSurf *plate, double t, double t0, int ed, double tz)

/* This routine returns the generatrix at the parameter value t of the
   spine. The spine is the line u=t0 of plate if ed=1, the line v=t0 if
   ed=0. The generatrix is trimmed in position using the Frenet's
   trihedron at t, so that the local y axis points along the normal and the
   local x axis along the binormal. The routine generatrix_lk is called,
   which returns the generatrix at point t, defined as a rational B-spline
   in the local coordinate system. */

{
  ParCurv *eg;
  vector *v1,**v0;
  double sg,u,v,**tt; 
  int i;

/* Meaning of the most important local variables (in alphabetical order)
eg       =  the generatrix expressed as a rational B-spline in its local
            system. eg is transformed to the global system and returned by this
            routine.
sg       =  control variable passed as an argument to generatrix_lk
            guaranteeing that the canal surface is built on the proper side with
            respect to the plate.
tt       =  transformation matrix (4 x 4). The control polygon
            coordinates of the generatrix at its final position are obtained by
            multiplying [t] with [x y 0 w]T, where (x,y,w) are the homogeneous
            coordinates in the local system.
v0       =  the Frenet trihedron.
*/

  if (ed) {
    u = t0;
    v = t;
  }
  else {
    u = t;
    v = t0;
  }
  v1 = revalderivsurf(plate, u, v, 0, 0);
  v0 = frenet_tr(plate, t, t0, ed);
  sg = fabs(t0+(double)ed);
  if (fabs(sg-1.0)<MACHPREC)
    sg = tz;
  else
    sg = -tz;
  eg = generatrix_lk(t, sg);
  tt = dbl_array2(4, 4);
  for (i=0; i<3; i++) {
    tt[0][i] = v0[i]->x;
    tt[1][i] = v0[i]->y;
    tt[2][i] = v0[i]->z;
  }
  u = v1->w;
  tt[0][3] = v1->x/u;
  tt[1][3] = v1->y/u;
  tt[2][3] = v1->z/u;
  tt[3][0] = tt[3][1] = tt[3][2] = 0.0;
  tt[3][3] = 1.0;
  free_varray1(v0, 3);
  vectfree(v1);
  for (i=0; i<eg->ncontpts; i++)
    mult4x4(tt,eg->contpts[i], eg->contpts[i]);
  free_darray2(tt);

  return(eg);
}

ParCurv *generatrix_lk(double t, double sg)

/* this routine returns the generatrix at the parameter value t of the
   spine, expressed as a rational B-spline in its local system. In its
   current form the routine returns a hemicircle.
   sg is the radius of the hemicircle. */

{
  ParCurv *eg;
  double tz,f1,d,f;
  int i;

  eg = egeomalloc1(3, 5);
  for (i=0; i<3; i++) {
    eg->knots[i] = 0.0;
    eg->knots[i+5] = 1.0;
  }
  eg->knots[3] = eg->knots[4] = 0.5;
  d = 1.0/sqrt((double)2.0);
  tz = fabs(sg);
  f = d*tz; f1 = d*sg;
  for (i=0; i<5; i++)
    eg->contpts[i]->z = 0.0;
  eg->contpts[0]->x = 0.0;
  eg->contpts[0]->y = tz;
  eg->contpts[0]->w = 1.0;
  eg->contpts[1]->y = f;
  eg->contpts[1]->x = f1;
  eg->contpts[1]->w = d;
  eg->contpts[2]->y = 0.0;
  eg->contpts[2]->x = sg;
  eg->contpts[2]->w = 1.0;
  eg->contpts[3]->y = -f;
  eg->contpts[3]->x = f1;
  eg->contpts[3]->w = d;
  eg->contpts[4]->x = 0.0;
  eg->contpts[4]->y = -tz;
  eg->contpts[4]->w = 1.0;

  return(eg);
}

vector **frenet_tr(ParSurf *plate, double t, double t0, int ed)

/* Compute the components of the Frenet trihedron at the parameter t
   (u=t or v=t according to the value of ed). t0 is the constant u or v
   value at the edge of the plate. Output is th Frenet trihedron in the
   order (X(s),Y(s),Z(s)). Y(s) is taken along the normal to plate, 
   Z(s) along the tangent. */

{
  vector **v0,*v1;
  double u,v;
  int ef;
  
  ef = 1-ed;
  if (ed) {
    u = t0;
    v = t;
  }
  else {
    u = t;
    v = t0;
  }
  v0 = vec_array1(3);
  v1=revalderivsurf(plate, u, v, ef, ed);
  v0[2]=unitvector(v1);
  vectfree(v1);
  v0[1] = normalsurf(plate, u, v);
  v0[0] = cross(v0[1], v0[2]);

  return(v0);
}

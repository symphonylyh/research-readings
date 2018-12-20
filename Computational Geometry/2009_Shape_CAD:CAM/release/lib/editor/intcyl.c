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

#define MACHPREC 1.1102230246251567E-13	  /* changed on 25 JUN 90 */
#define MACHPR 1.1*MACHPREC
#define MAXSECT 500

ParSurf *canal_to_int(ParCurv *directrix, short ns, short nc, double eps,
		      double *pars, double *circ, ParCurv **gen, int *m,
		      int *n, short m0, short n0, double tz, double theta)

/* Approximate a generalized cylinder by an integral B-spline surface
   patch.
   plate is a rational B-spline surface patch representing the deformed
   plate.
   The spine is a parametric line directrix.
   The output geometry interpolates the generalized
   cylinder at m sections, positioned at parameter values contained in the
   array pars and at n points along the generatric circumference, whose
   parameter values are contained in circ (m, n, pars and circ are
   updated by this routine).
   gen is an array containing the interpolated sections (generatrices).
   The iterative aproximation procedure starts with m0 equidistant
   sections between consecutive interpolated sections.
   ns is the number of sampled sections between consecutive interpolated
   sections.
   nc is the number of sampling points on each section (generatrix).
   eps is the maximum allowable positional error between the
   approximating B-spline and the generalized cylinder.
   tz is the value of the tolerance.
*/

{
  ParSurf *canal;
  double t,d ;
  int i,order,vo,vp,np;

/* Meaning of the most important local variables (in alphabetical order)

canal    =  The output geometry (rational B-spline surface patch).
uo       =  u-order of plate.
up       =  number of u-vertices of plate's control polyhedron.
vo       =  v-order of plate or gen.
vp       =  number of v-vertices of plate's or gen's control polyhedron.

*/

  *m = m0 ;
  order = directrix->order;
  np = directrix->ncontpts;
  d = directrix->knots[np+order-1] - directrix->knots[0];
  d /= (double)(*m-1) ;
  if (m0<order) {
    printf("canal_to_int: initial number of sections smaller than spine order!\n");
    return (NULL);
  }
  pars[0] = directrix->knots[0];
  for (i=0;i<*m;i++) {
    if (i) pars[i] = pars[i-1]+d;
    gen[i] = cyl_generatrix(directrix, pars[i], tz, theta) ;
  } 
  vo = gen[0]->order ; vp = gen[0]->ncontpts ; *n = n0 ;
  if (n0<vo) {
    printf("canal_to_int: initial number of sections smaller than spine order!\n");
    return (NULL);
  }
  d = (gen[0]->knots[vo+vp-1]-gen[0]->knots[0])/(double)(n0-1) ;
  circ[0] = gen[0]->knots[0];
  for (i=1; i<n0; i++)
    circ[i] = circ[i-1]+d ;

/* Iterate as required for convergence */

  do {
    if (*m != m0)
      free_fgeom(canal);
    if ((canal = loft_integral(gen,pars,circ,order,*m,*n)) == NULL)
      return (NULL);
  } while (int_sample_gencyl(directrix, ns, nc, eps, canal, gen, pars, circ,
			     m, n, tz, theta));
  return(canal);
}

int int_sample_gencyl(ParCurv *directrix, short ns, short nc, double eps,
		      ParSurf *canal, ParCurv **gen, double *pars,
		      double *circ, int *m, int *n, double tz, double theta)

/* Sample uniformly the generalized cylinder. canal is the approximating
   rational B-spline patch, gen are the plane sections (generatrices) at m
   parameter values of the spine, which are contained in pars. n is the
   number of data points along the circumference of the generatrix and circ
   contains the respective parameter values. The routine returns the number
   of new sections (one new section midway between existing, consecutive
   sections for each failed interval) and the number and parameter values
   of data points along the circumference required for better approximation.
   pars, circ, m, n and gen are updated.
   directrix,ns,nc,eps,tz have the same meaning as in canal_to_rat. */

{
  ParCurv *eg;
  vector *v1,*v2;
  double ds,d,f,f0,df,mg;
  int i,j,k,*l,**ll,b,b0,b1,tv;

/* Meaning of the most important local variables (in alphabetical order)

b0       =  total number of failed intervals along the spine
            (=number of new sections).
b1       =  total number of failed intervals along the circumference of
            the generatrix.
eg       =  the generatrix.
l        =  integer array. For each interval :
            1 if some point failed between gen[i] and gen[i-1].
            0 otherwise.
ll       =  integer array. For each interval between consecutive
            parameter values along the circumference :
            1 if some point failed between circ[j-1] and circ[j] and
              between gen[i] and gen[i-1].
            0 otherwise.

*/

  l = int_array1(MAXSECT);
  ll = int_array2(MAXSECT, MAXSECT);
  f0 = (double)(nc-1);
  b0 = b1 = 0;

/* Do for each interval between consecutive interpolated sections */

  for (i=1; i<*m; i++) {
    ds = (pars[i]-pars[i-1])/(double)(ns+1);
    l[i] = 0;
    for (k=0; k<*n; k++)
      ll[i][k] = 0;
    d = pars[i-1] ;

/* Generate new sections */
  
    for (j=0;j<ns;j++) {
      d += ds;
      if (i>1 || j)
	free_egeom(eg);
      eg = cyl_generatrix(directrix, d, tz, theta);
      df = (eg->knots[eg->order+eg->ncontpts-1]-eg->knots[0])/f0;
      f = eg->knots[0] ;
      if (fabs(1.0-d) < MACHPREC) d -= MACHPR;

/* Sample uniformly for each new section */

      for (k=0; k<nc; k++) {
	if (fabs(1.0-f) < MACHPREC)
	  f -= MACHPR ;
	tv = find(*n,circ,f);
	if (!ll[i][tv]) {
	  v1 = rbspeval(eg, f, 0);
	  v2 = revalderivsurf(canal, d, f, 0, 0);
	  sub_vect1(v1, v2, v1);
	  vectfree(v2);
	  mg = mag(v1);
	  vectfree(v1);
	  b = mg>eps;
	  l[i] = l[i] || b;
	  ll[i][tv] = ll[i][tv] || b;
        }
	f += df;
      }
    }  
    b0 += l[i];
  }
  if (*m+b0 > MAXSECT)
    errormsg(6, "max no of sections exceeded in canal_to_int");

/* Determine failed intervals along the generatrix circumference */

  for (j=0; j<*n-1; j++) {
    ll[0][j] = 0;
    for (i=1; i<*m; i++)
      ll[0][j] = ll[0][j] || ll[i][j];
    b1 += ll[0][j];
  }
  if (*n+b1 > MAXSECT)
    errormsg(6,"max no of gen pnts exceeded in canal_to_int") ;

/* Generate one section for each failed interval. Rearrange indices. */

  b = b0;
  free_egeom(eg);
  if (b)
    for (i = *m-1; i>0; i--)
      if (b) {
	j = i+b;
	gen[j] = copyegeom(gen[i], gen[j]);
	if (l[i]) {
	  b--;
	  j--;
	  if (j<*m)
	    free_egeom(gen[j]);
	  d = (pars[i]+pars[i-1])/2.0;
	  gen[j] = cyl_generatrix(directrix, d, tz, theta);
	}
      }
  
/* Compute new pars */

  b = b0;
  if (b)
    for (i = *m-1; i>0; i--)
      if (b) {
	j = i+b;
	pars[j] = pars[i];
	if (l[i]) {
	  b--;
	  j--;
	  pars[j] = (pars[i]+pars[i-1])/2.0;
	}
      }

/* Compute new circ */

  b = b1;
  if (b)
    for (i = *n-1; i>0; i--)
      if (b) {
	j = i+b;
	circ[j] = circ[i];
	if (ll[0][i-1]) {
	  b--;
	  j--;
	  circ[j] = (circ[i]+circ[i-1])/2.0;
	}
      }
  free_iarray2(ll);
  free_iarray1(l);
  *m += b0;
  *n += b1;

  return(b0);
}

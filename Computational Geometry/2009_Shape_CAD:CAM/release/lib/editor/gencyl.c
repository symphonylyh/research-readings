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

ParSurf *cyl_canal_to_rat(ParCurv *directrix, short ns, short nc, double eps,
			  double *pars, ParCurv **gen, int *m, short m0,
			  double tz, double theta)

/* Approximate a generalized cylinder by a rational B-spline surface
   patch.
   plate is a rational B-spline surface patch representing the deformed
   plate.
   The output geometry interpolates the generalized
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
  int i,order,np;
  double t,d ;

/* Meaning of the most important local variables (in alphabetical order)

canal    =  The output geometry (rational B-spline surface patch).

*/

  *m = m0 ;
  order = directrix->order;
  np = directrix->ncontpts;

  d = directrix->knots[np+order-1] - directrix->knots[0];
  d /= (double)(*m-1) ;

  if (m0<order) {
    printf("Number of sampling points smaller than spine order!\n");
    return (NULL);
  }
  pars[0] = directrix->knots[0];

  for (i=0; i<*m; i++) {
    if (i)
      pars[i] = pars[i-1] + d;
    gen[i] = cyl_generatrix(directrix, pars[i], tz, theta);
  }

/* Iterate as required for convergence */

  do {
    if (*m>MAXSECT) {
      printf("Maximum number of sections exceeded in canal_to_rat!\n") ;
      return (NULL);
    }
    if (*m != m0)
      free_fgeom(canal) ;
    if ((canal = tol_loft_rational(gen, pars, order, *m)) == NULL) {
      return (NULL);
    }
  } while (cyl_sample_gencyl(directrix, ns, nc, eps, canal, gen, pars, m, tz,
			     theta));

  return(canal);
}

short cyl_sample_gencyl(ParCurv *directrix, short ns, short nc, double eps,
			ParSurf *canal, ParCurv **gen, double *pars, int *m,
			double tz, double theta)

/* Sample uniformly the generalized cylinder. canal is the approximating
   rational B-spline patch, gen are the plane sections (generatrices) at m
   parameter values of the spine, which are contained in pars. The routine
   returns the number of new sections required for better approximation
   (one new section midway between existing, consecutive sections for each
   failed interval). pars, m and gen are updated.
   directrix,ns,nc,eps,tz,theta have the same meaning as in canal_to_rat. */

{
  ParCurv *eg;
  vector *v1,*v2;
  double ds,d,f,f0,df,mg;
  int i,j,k,*l,b,b0;

/* Meaning of the most important local variables (in alphabetical order)

b0       =  total number of failed intervals (=number of new sections).
eg       =  the generatrix.
l        =  integer array. For each interval :
            1 if some point failed between gen[i] and gen[i-1].
            0 otherwise.
*/

  l = int_array1(MAXSECT);
  b0 = 0;
  f0 = (double)(nc-1);

/* Do for each interval between consecutive interpolated sections */

  for (i=1; i<*m; i++) {
    ds = (pars[i]-pars[i-1])/(double)(ns+1);
    l[i] = 0;
    d = pars[i-1];

/* Generate new sections */
  
    for (j=0; j<ns; j++)
      if (!l[i]) {
	d += ds;
	if (i>1 || j)
	  free_egeom(eg);
	eg = cyl_generatrix(directrix, d, tz, theta);
	df = (eg->knots[eg->order+eg->ncontpts-1]-eg->knots[0])/f0;
	f = eg->knots[0];
	if (fabs(1.0-d) < MACHPREC)
	  d -= MACHPR;

/* Sample uniformly for each new section */

	for (k=0;k<nc;k++)
	  if (!l[i]) {
	    if (fabs(1.0-f) < MACHPREC)
	      f -= MACHPR ;
	    v1 = rbspeval(eg, f, 0);
	    v2 = revalderivsurf(canal, d, f, 0, 0);
	    sub_vect1(v1,v2,v1);
	    vectfree(v2);
	    mg = mag(v1);
	    vectfree(v1) ;
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
	j = i+b;
	gen[j] = copyegeom(gen[i],gen[j]);
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
  free_iarray1(l);
  *m += b0;
  return (b0);
}

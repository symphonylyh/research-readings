/* Copyright (C) 1992 Massachusetts Institute of Technology, Cambridge, MA
   All rights reserved
*/

/* powsurf.c */
/* based on: /usr/deslab/epraxiteles/praxlib/iges/nrb_to_pow.c */
/*           /usr/deslab/epraxiteles/praxlib/iges/pow_to_nrb.c */

/* berntomono2()
 * calc_points_surf()
 * monotobern2()
 * ParSurf_to_PowSurf()
 * ParSurf_to_PowSurf2()
 * PowParSurf_compare()
 * PowParSurf_compare2()
 * PowSurf_error()
 * PowSurf_iso()
 * PowSurf_to_ParSurf()
 * PowSurf_to_ParSurf_int()
 * PowSurf_to_ParSurf_loft()
 * RemoveSurfKnotU()
 * RemoveSurfKnotV()
*/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <malloc.h>
#include "iges.h"   /* must be before "gen.h" is referenced */
#include "appr.h"
#include "bspl.h"
#include "editor.h"

/********* calc_points_surf() *********
* 1     Purpose
* 
*       Evaluate points on a parametric spline surface.
* 
* 2     Specification
* 
*       #include "iges.h"
*       void calc_points_surf(double *unodes, double *vnodes, PowSurf *sgeom,
*                             double **xyz, double *u, double *v, int nc);
* 
* 3     Description
* 
*       This function evaluates a parametric spline surface at its nodal values.
* 
* 4     References
* 
*       [1]   Digital Representation for Communication of Product Definition
*             Data, US PRO/IPO-100, Initial Graphics Exchange Specification
* 	    (IGES) 5.2, IGES/PDES Organization, U.S. Product Data Association,
* 	    Fairfax, VA, November 1993.
* 
* 5     Parameters
* 
*           1.  double * unodes
*               On entry: the address of an array containing the nodes of the
* 	      surface in the u direction.
* 
*           2.  double * vnodes
*               On entry: the address of an array containing the nodes of the
* 	      surface in the v direction.
* 
*           3.  PowSurf * sgeom
*               On entry:  the pointer to a structure containing the parametric
* 	      spline surface.
* 
*           4.  double *** xyz
*               On entry:  the address of a triply-dimensioned array that will
* 	      contain the evaluated points.
* 
*               On exit:  the address of a triply-dimensioned array containing
* 	      the evaluated points.
* 
*           5.  double * u
*               On exit:  the address of an array containing the nodes of the
* 	      surface in the u direction.
* 
*           6.  double * v
*               On exit:  the address of an array containing the nodes of the
* 	      surface in the v direction.
* 
*           7.  int nc
*               On entry:  flag specifying the continuity of the surface.
* 
*                 nc   Continuity
*                 --   ----------
*                  0   Curvature
*                  1   Tangent
*                  2   Position
* 
* Functions referenced by calc_points_surf() are:
*  PowSurf_eval()
*  vectfree()
*
* Functions that reference calc_points_surf() are:
*  PowSurf_to_ParSurf_int()
*/

void calc_points_surf(double *unodes, double *vnodes, PowSurf *sgeom,
		      double **xyz, double *u, double *v, int nc)
{
  vector *vect;
  int i, j, m, n;

  m = (sgeom->vsegmts-1)*nc + 4;
  n = (sgeom->usegmts-1)*nc + 4;
  for (i=0; i<m; i++)
    for (j=0; j<n; j++) {
      vect = PowSurf_eval(sgeom, unodes[j], vnodes[i],0,0);
      xyz[i*n + j][0] = vect->x;
      xyz[i*n + j][1] = vect->y;
      xyz[i*n + j][2] = vect->z;
      xyz[i*n + j][3] = 1.0;
      u[i*n + j] = unodes[j];
      v[i*n + j] = vnodes[i];
      vectfree(vect);
    }
}

/********* ParSurf_to_PowSurf() *********
* 1     Purpose
* 
*       Convert NURBS surface to parametric spline surface.
* 
* 2     Specification
* 
*       #include "iges.h"
*       PowSurf *ParSurf_to_PowSurf(ParSurf *fgeom, PowSurf *sgeom,
*                                   double *maxerr);
* 
* 3     Description
* 
*       This function converts a NURBS surface into a parametric spline surface.
* 
* 4     References
* 
*       [1]   Digital Representation for Communication of Product Definition
*             Data, US PRO/IPO-100, Initial Graphics Exchange Specification
* 	    (IGES) 5.2, IGES/PDES Organization, U.S. Product Data Association,
* 	    Fairfax, VA, November 1993.
* 
* 5     Parameters
* 
*           1.  ParSurf * fgeom
*               On entry:  the address of a structure containing the NURBS
* 	      surface.
* 
*           2.  PowSurf * sgeom
*               On entry:  the address of a structure that will contain the
* 	      parametric spline surface.  If specified as NULL then a new
* 	      structure will be allocated by this function.
* 
*               On exit:  the address of a structure that contains the
* 	      parametric surface.  This same address is returned by the
* 	      function.
* 
*           3.  double * maxerr
*               On exit:  the addrees of a variable that will contain the
* 	      global conversion error calculated as the maximum distance 
* 	      between points on the parametric spline and the NURBS surface
* 	      evaluated at the breakpoints and the midpoints of the
* 	      breakpoints.
* 
* 6     Return Values, Error Indicators and Warnings
* 
*       The address of the structure containing the parametric spline surface
*       is returned.  If sgeom is not NULL on entry, then its address is
*       returned; otherwise, the address of the a newly allocated structure
*       is returned.
* 
* 8     Further Comments
* 
*       The returned structure should be deallocated with free_sgeom().
* 
* Functions referenced by ParSurf_to_PowSurf() are:
*  copyvector()
*  evalderivsurf()
*  PowParSurf_compare()
*  scale_vect1()
*  sgeomalloc()
*  vectfree()
*
* Functions that reference ParSurf_to_PowSurf() are:
*  SaveIgesSurf()
*/

PowSurf *ParSurf_to_PowSurf(ParSurf *fgeom, PowSurf *sgeom, double *maxerr)
{
  static double fact03[] = {1.0, 1.0, 2.0, 6.0};
  vector *v;
  double f;
  int ubkpts = 0, vbkpts = 0;
  int i, j, m, n;
     
  if (!sgeom)
    sgeom = sgeomalloc(fgeom->ucontpts - fgeom->uorder + 1,
		       fgeom->vcontpts - fgeom->vorder + 1);

  sgeom->uknots[0] = fgeom->uknots[0];
  sgeom->vknots[0] = fgeom->vknots[0];

  for (i=1; i<fgeom->uorder + fgeom->ucontpts; i++) 
    if (fgeom->uknots[i] != fgeom->uknots[i-1])
      sgeom->uknots[++ubkpts] = fgeom->uknots[i];
     
  for (i=1; i<fgeom->vorder + fgeom->vcontpts; i++) 
    if (fgeom->vknots[i] != fgeom->vknots[i-1])
      sgeom->vknots[++vbkpts] = fgeom->vknots[i];
     
  for (i=0; i<ubkpts; i++)
    for (j=0; j<vbkpts; j++)
      for (m=0; m<4; m++) 
	for (n=0; n<4; n++) {
	  v = evalderivsurf(fgeom, sgeom->uknots[i], sgeom->vknots[j], n, m);
	  v->w = 1.0;
	  f = 1 / (fact03[m]*fact03[n]);
	  scale_vect1(f, v, v);
	  copyvector(v, sgeom->contpts[4*m + n][i][j]);
	  vectfree(v);
	}
  PowParSurf_compare (fgeom, sgeom, maxerr, 1);

  return (sgeom);
}

/********* ParSurf_to_PowSurf2() *********
* 1     Purpose
* 
*       Convert NURBS surface to parametric spline surface.
* 
* 2     Specification
* 
*       #include "iges.h"
*       PowSurf *ParSurf_to_PowSurf2(ParSurf *fgeom, PowSurf *sgeom,
*                                   double *maxerr);
* 
* 3     Description
* 
*       This function converts a NURBS surface into a parametric spline surface.
* 
* 4     References
* 
*       [1]   Digital Representation for Communication of Product Definition
*             Data, US PRO/IPO-100, Initial Graphics Exchange Specification
* 	    (IGES) 5.2, IGES/PDES Organization, U.S. Product Data Association,
* 	    Fairfax, VA, November 1993.
* 
* 5     Parameters
* 
*           1.  ParSurf * fgeom
*               On entry:  the address of a structure containing the NURBS
* 	      surface.
* 
*           2.  PowSurf * sgeom
*               On entry:  the address of a structure that will contain the
* 	      parametric spline surface.  If specified as NULL then a new
* 	      structure will be allocated by this function.
* 
*               On exit:  the address of a structure that contains the
* 	      parametric surface.  This same address is returned by the
* 	      function.
* 
*           3.  double * maxerr
*               On exit:  the addrees of a variable that will contain the
* 	      global conversion error calculated as the maximum distance 
* 	      between points on the parametric spline and the NURBS surface
* 	      evaluated at the breakpoints and the midpoints of the
* 	      breakpoints.
* 
* 6     Return Values, Error Indicators and Warnings
* 
*       The address of the structure containing the parametric spline surface
*       is returned.  If sgeom is not NULL on entry, then its address is
*       returned; otherwise, the address of the a newly allocated structure
*       is returned.
* 
* 8     Further Comments
* 
*       The returned structure should be deallocated with free_sgeom().
* 
* Functions referenced by ParSurf_to_PowSurf2() are:
*  copyfgeom()
*  dbl_array2()
*  fgeomalloc1()
*  free_darray2()
*  free_fgeom()
*  PowParSurf_compare()
*  sgeomalloc()
*  surfoslo3()
*
* Functions that reference ParSurf_to_PowSurf2() are:
*  SaveIgesSurf()
*/

PowSurf *ParSurf_to_PowSurf2(ParSurf *fgeom, PowSurf *sgeom, double *maxerr)
{
  ParSurf
    *bernstein;
  int
    i, ii, j, jj, k, kk, m, n,
    nInternU, nInternV,
    nKnotU, nKnotV;
  double
    **c,
    coeff;

  /* find number of unique internal knot values */

  nInternU = 0;
  for (i=fgeom->uorder+1; i<=fgeom->ucontpts; i++)
    if (i == fgeom->ucontpts || fgeom->uknots[i] > fgeom->uknots[i-1])
      ++nInternU;
  nInternV = 0;
  for (i=fgeom->vorder+1; i<=fgeom->vcontpts; i++)
    if (i == fgeom->vcontpts || fgeom->vknots[i] > fgeom->vknots[i-1])
      ++nInternV;

  /* allocate memory if necessary */

  if (!sgeom)
    sgeom = sgeomalloc(nInternU + 1, nInternV + 1);

  sgeom->uknots[0] = fgeom->uknots[0];
  sgeom->vknots[0] = fgeom->vknots[0];
  sgeom->uorder = fgeom->uorder;
  sgeom->vorder = fgeom->vorder;

  /* find number of unique internal knot values */

  nInternU = 0;
  for (i=fgeom->uorder+1; i<=fgeom->ucontpts; i++)
    if (i == fgeom->ucontpts || fgeom->uknots[i] > fgeom->uknots[i-1])
      sgeom->uknots[++nInternU] = fgeom->uknots[i-1];
  sgeom->uknots[nInternU+1] = fgeom->uknots[fgeom->ucontpts];

  nInternV = 0;
  for (i=fgeom->vorder+1; i<=fgeom->vcontpts; i++)
    if (i == fgeom->vcontpts || fgeom->vknots[i] > fgeom->vknots[i-1])
      sgeom->vknots[++nInternV] = fgeom->vknots[i-1];
  sgeom->vknots[nInternV+1] = fgeom->vknots[fgeom->vcontpts];

  /* allocate space to allow for "order" knot multiplicity of the
   * internal knots
   */
  bernstein = fgeomalloc1(fgeom->uorder, fgeom->vorder,
			  fgeom->uorder + fgeom->uorder*nInternU,
			  fgeom->vorder + fgeom->vorder*nInternV);

  if (nInternU || nInternV) {

    /* new knot vectors with "order" knot multiplicity */

    /* "uorder" knots (u=0) at beginning */

    for (m=0; m<fgeom->uorder; m++)
      bernstein->uknots[m] = fgeom->uknots[m];

    /* U internal knots */

    for (i=fgeom->uorder+1, k=1; i<=fgeom->ucontpts; i++) {
      bernstein->uknots[m++] = fgeom->uknots[i-1];
      if (i == fgeom->ucontpts || fgeom->uknots[i] > fgeom->uknots[i-1]) {
	for (j=0; j<fgeom->uorder-k;j++) {
	  bernstein->uknots[m++] = fgeom->uknots[i-1];
	}
	k = 1;
      }
      else
	k++;         /* count U knot multiplicity */
    }

    /* "uorder" knots (u=1) at end */

    for (i=fgeom->ucontpts; i<fgeom->uorder+fgeom->ucontpts; i++)
      bernstein->uknots[m++] = fgeom->uknots[i];

    /* "vorder" knots (v=0) at beginning */

    for (m=0; m<fgeom->vorder; m++)
      bernstein->vknots[m] = fgeom->vknots[m];

    /* V internal knots */

    for (i=fgeom->vorder+1, k=1; i<=fgeom->vcontpts; i++) {
      bernstein->vknots[m++] = fgeom->vknots[i-1];
      if (i == fgeom->vcontpts || fgeom->vknots[i] > fgeom->vknots[i-1]) {
	for (j=0; j<fgeom->vorder-k;j++) {
	  bernstein->vknots[m++] = fgeom->vknots[i-1];
	}
	k = 1;
      }
      else
	k++;         /* count V knot multiplicity */
    }

    /* "vorder" knots (v=1) at end */

    for (i=fgeom->vcontpts; i<fgeom->vorder+fgeom->vcontpts; i++)
      bernstein->vknots[m++] = fgeom->vknots[i];

    nKnotU = fgeom->uorder*2 + fgeom->uorder*nInternU;
    nKnotV = fgeom->vorder*2 + fgeom->vorder*nInternV;

    /* use the Oslo algorithm to insert knots so that every internal
     * knot has "order" multiplicity;
     * this converts the original B-spline surface  to Bernstein basis
     */

    surfoslo3(fgeom, bernstein, 4);
  }
  else {

    /* if 1 internal knot in both U and V, then Bezier surface which is
     * already in the Bernstein basis
     */

    copyfgeom(fgeom, bernstein);
  }

  /* convert from the Bernstein to the monomial basis */

  c = dbl_array2(bernstein->uorder, bernstein->vorder);
  for (k=0; k<=nInternU; k++) {
    m = k*4;

    for (kk=0; kk<=nInternV; kk++) {
      n = kk*4;

      for (j=0; j<3; j++) {     /* convert each dimension (X,Y,Z) separately */
	for (i=0; i<bernstein->uorder; i++)
	  for (ii=0; ii<bernstein->vorder; ii++) {
	    if (j == 0)
	      c[i][ii] = bernstein->contpts[m+i][n+ii]->x;
	    else if (j == 1)
	      c[i][ii] = bernstein->contpts[m+i][n+ii]->y;
	    else
	      c[i][ii] = bernstein->contpts[m+i][n+ii]->z;
	    c[i][ii] /= bernstein->contpts[m+i][n+ii]->w;
	  }

	/* multiply the coefficient by (nInternU+1)^i
	 * to scale the parameter range from [0,1] to the
	 * portion of the parameter space between the bounding
	 * knot values
	 */

	for (i=0; i<bernstein->uorder; i++) {
	  for (ii=0; ii<bernstein->vorder; ii++) {
	    coeff = berntomono2(bernstein->uorder-1, bernstein->vorder-1, i,
				ii, c)*pow(nInternU+1, i)*pow(nInternV+1, ii);
	    jj = bernstein->vorder*ii + i;
	    if (j == 0)
	      sgeom->contpts[jj][k][kk]->x = coeff;
	    else if (j == 1)
	      sgeom->contpts[jj][k][kk]->y = coeff;
	    else {
	      sgeom->contpts[jj][k][kk]->z = coeff;
	      sgeom->contpts[jj][k][kk]->w = 1.0;
	    }
	  }
	}
      }
    }
  }
  free_darray2(c);
  free_fgeom(bernstein);

  PowParSurf_compare(fgeom, sgeom, maxerr, 1);

  return sgeom;
}

/* Converts i,j - Bernstein coefficient to monomial basis.
 * m,n	  : degrees in u,v
 * c(i,j) : coeff of Bi,m(u)Bj,n(v)
 */

/*****************************************************************************
*                               berntomono2()
******************************************************************************
* 
* 1   Purpose
*     This function converts a polynomial in Bernstein form to one 
*     in the monomial form.  This function is called once for each i =
*     0,1,...,m-1, and j=0,1,...,n-1, for each spatial coordinate.
* 
* 2   Specification
*     #include "iges.h"
*     double berntomono2(int m, int n, int i, int j, double **c)
* 
* 3   Description
*     This routine converts Bernstein basis representation of an algebraic 
*     surface the monomial basis.
* 
* 4   References
*     Not applicable
* 
* 5   Parameters
*       1. int m
*         On entry: the maximum degree in u of the input Bernstein polynomial
*       2. int n
*         On entry: the maximum degree in v of the input Bernstein polynomial
*       3. int i
*         On entry: degree in u for which the monomial coefficient is required.
*       4. int j
*         On entry: degree in v for which the monomial coefficient is required.
*       5. double ** c
*         On entry: array of length m+1 by n+1 holding the control polygon
* 	   of the Bernstein polynomial
* 
* 6   Return Values, Error Indicators and Warnings
*     Return the monomial coefficient of the monomial corresponding to
*     the input Bernstein polynomial.
* 
* 7   Accuracy
*     Not applicable.
* 
* 8   Further Comments
*     Not applicable.
* 
* 9   Functions referenced by berntomono2() are:
*     factorial()
* 
* 10  Functions that reference berntomono2() are: ParSurf_to_PowSurf2()
* 
******************************************************************************/

double berntomono2(int m, int n, int i, int j, double **c)
{
  int k,l;
  double C, f1, f2, fmi, fnj, fn, fm, minus1, exp;

  fm = factorial(m);	
  fn = factorial(n);
  fmi = factorial(m-i);	
  fnj = factorial(n-j);

  C = 0;
  for (k=0; k<=i; k++) {
    f1 = 1.0 / factorial(i-k) / factorial(k);
    for (l=0; l<=j; l++) {
      f2 = 1.0 / factorial(j-l) / factorial(l);
      exp = i + j - k - l;
      minus1 = pow(-1.0, exp); 
      C += (minus1*f1*f2*c[k][l]);
    }
  }
     
  C *= (fm*fn/fmi/fnj);
     
  return C;
}

/********* PowParSurf_compare() *********
* 1     Purpose
* 
*       Position error between parametric spline surface and NURBS surface.
* 
* 2     Specification
* 
*       #include "iges.h"
*       void PowParSurf_compare(ParSurf *fgeom, PowSurf *sgeom, double *maxerr,
*                               int nc);
* 
* 3     Description
* 
*       This function calculates the position error between a parametric
*       spline surface and a NURBS surface.
* 
* 4     References
* 
*       [1]   Digital Representation for Communication of Product Definition
*             Data, US PRO/IPO-100, Initial Graphics Exchange Specification
* 	    (IGES) 5.2, IGES/PDES Organization, U.S. Product Data Association,
* 	    Fairfax, VA, November 1993.
* 
* 5     Parameters
* 
*           1.  ParSurf * fgeom
*               On entry:  the address of a structure containing the NURBS
* 	      surface.
* 
*           2.  PowSurf * sgeom
*               On entry:  the address of a structure containing the parametric
* 	      spline surface.
* 
*           3.  double * maxerr
*               On exit:  the addrees of a variable that will contain the
* 	      global conversion error calculated as the maximum distance
* 	      between points on the parametric spline and the NURBS surface
* 	      breakpoints evaluated at the breakpoints and the midpoints of
*             the breakpoints.
* 
*           4.  int nc
*               On entry: flag specifying the degree of continuity at the
* 	      breakpoints of the parametric polynomials.
* 
*                 nc   Continuity
*                 --   ----------
*                  0   Surfature
*                  1   Tangent
*                  2   Position
* 
* Functions referenced by PowParSurf_compare() are:
*  distance()
*  PowSurf_eval()
*  revalderivsurf()
*  vectfree()
*
* Functions that reference PowParSurf_compare() are:
*  ParSurf_to_PowSurf()
*  PowSurf_to_ParSurf_int()
*  PowSurf_to_ParSurf_loft()
*/

void PowParSurf_compare(ParSurf *fgeom, PowSurf *sgeom, double *maxerr,
			int nc)
{
  vector *vect1, *vect2, *vect3, *vect4;
  double2 u1, u2, v1, v2, u3, v3, u4, v4, dist1, dist2, max;
  int k, l, m, n;
     
  *maxerr = 0.0;
  
  for (k = 0; k < sgeom->usegmts; k++) {
    m = k*nc;
    for (l = 0; l < sgeom->vsegmts; l++) {
      n = l*nc;
      u1 = fgeom->uknots[m+3] + 0.5*(fgeom->uknots[m+4] - fgeom->uknots[m+3]);
      u2 = fgeom->uknots[m+4] - ZERO;
      v1 = fgeom->vknots[n+3] + 0.5*(fgeom->vknots[n+4] - fgeom->vknots[n+3]);
      v2 = fgeom->vknots[n+4] - ZERO;
      u3 = sgeom->uknots[k] + 0.5*(sgeom->uknots[k+1] - sgeom->uknots[k]);
      v3 = sgeom->vknots[l] + 0.5*(sgeom->vknots[l+1] - sgeom->vknots[l]);
      u4 = sgeom->uknots[k+1] - ZERO;
      v4 = sgeom->vknots[l+1] - ZERO;
      vect1 = revalderivsurf(fgeom, u1, v1, 0, 0);
      vect2 = revalderivsurf(fgeom, u2, v2, 0, 0);
      vect3 = PowSurf_eval(sgeom, u3, v3, 0, 0);
      vect4 = PowSurf_eval(sgeom, u4, v4, 0, 0);
      dist1 = distance(vect1, vect3);
      dist2 = distance(vect2, vect4);
      max = MAX(dist1, dist2);
      *maxerr = MAX(max, *maxerr);
      vectfree(vect1);
      vectfree(vect2);
      vectfree(vect3);
      vectfree(vect4);
    }
  }
}

/********* PowParSurf_compare2() *********
* 1     Purpose
* 
*       Position error between parametric spline surface and NURBS surface.
* 
* 2     Specification
* 
*       #include "iges.h"
*       void PowParSurf_compare2(ParSurf *fgeom, PowSurf *sgeom,
*                                double *maxerr, int nc);
* 
* 3     Description
* 
*       This  function calculates the position error between a parametric
*       spline  surface  and  a NURBS surface.
* 
* 4     References
* 
*       [1]   Digital Representation for Communication of Product Definition
*             Data, US PRO/IPO-100, Initial Graphics Exchange Specification
* 	    (IGES) 5.2, IGES/PDES Organization, U.S. Product Data Association,
* 	    Fairfax, VA, November 1993.
* 
* 5     Parameters
* 
*           1.  ParSurf * fgeom
*               On entry:  the address of a structure containing the NURBS
* 	      surface.
* 
*           2.  PowSurf * sgeom
*               On entry:  the address of a structure containing the parametric
* 	      spline surface.
* 
*           3.  double * maxerr
*               On exit:  the address of a variable that will contain the
* 	      global conversion error calculated as the maximum distance
* 	      between points on the parametric spline and the NURBS surface
* 	      evaluated at the breakpoints and the midpoints of the
*             breakpoints.
* 
*           4.  int nc
*               On entry: flag specifying the degree of continuity at the
* 	      breakpoints of the parametric polynomials.
* 
*                 nc   Continuity
*                 --   ----------
*                  0   Surfature
*                  1   Tangent
*                  2   Position
* 
* Functions referenced by PowParSurf_compare2() are:
*  distance()
*  PowSurf_eval()
*  revalderivsurf()
*  vectfree()
*
* Functions that reference PowParSurf_compare2() are:
*  ParSurf_to_PowSurf2()
*  PowSurf_to_ParSurf_int()
*  PowSurf_to_ParSurf_loft()
*/

void PowParSurf_compare2(ParSurf *fgeom, PowSurf *sgeom, double *maxerr,
			int nc)
{
  vector *vect1, *vect2, *vect3, *vect4;
  double2 u1, u2, v1, v2, u3, v3, u4, v4, dist1, dist2, max;
  int i, j, k, l, m, n;
     
  *maxerr = 0.0;
  
  for (k = 0; k < sgeom->usegmts; k++) {
    m = k*nc;
    i = fgeom->uorder;
    for (l = 0; l < sgeom->vsegmts; l++) {
      n = l*nc;
      j = fgeom->vorder;
      u1 = fgeom->uknots[i-1] + 0.5*(fgeom->uknots[i] - fgeom->uknots[i-1]);
      u2 = fgeom->uknots[i] - ZERO;
      v1 = fgeom->vknots[j-1] + 0.5*(fgeom->vknots[j] - fgeom->vknots[j-1]);
      v2 = fgeom->vknots[j] - ZERO;
      u3 = sgeom->uknots[k] + 0.5*(sgeom->uknots[k+1] - sgeom->uknots[k]);
      v3 = sgeom->vknots[l] + 0.5*(sgeom->vknots[l+1] - sgeom->vknots[l]);
      u4 = sgeom->uknots[k+1] - ZERO;
      v4 = sgeom->vknots[l+1] - ZERO;
      vect1 = revalderivsurf(fgeom, u1, v1, 0, 0);
      vect2 = revalderivsurf(fgeom, u2, v2, 0, 0);
      vect3 = PowSurf_eval(sgeom, u3, v3, 0, 0);
      vect4 = PowSurf_eval(sgeom, u4, v4, 0, 0);
      dist1 = distance(vect1, vect3);
      dist2 = distance(vect2, vect4);
      max = MAX(dist1, dist2);
      *maxerr = MAX(max, *maxerr);
      vectfree(vect1);
      vectfree(vect2);
      vectfree(vect3);
      vectfree(vect4);
    }
  }
}

/********* PowSurf_error() *********
* 1     Purpose
* 
*       Discontinuity errors of parametric spline surface.
* 
* 2     Specification
* 
*       #include "iges.h"
*       void PowSurf_error(PowSurf *sgeom, double errors[2]);
* 
* 3     Description
* 
*       This  function  calculates  the  positional  and  tangent
*       discontinuity  errors  of  a  parametric spline surface.
* 
* 4     References
* 
*       [1]   Digital Representation for Communication of Product Definition
*             Data, US PRO/IPO-100, Initial Graphics Exchange Specification
* 	    (IGES) 5.2, IGES/PDES Organization, U.S. Product Data Association,
* 	    Fairfax, VA, November 1993.
* 
* 5     Parameters
* 
*           1.  PowSurf * sgeom
*               On entry:  the address of a structure containing the parametric
* 	      spline surface.
* 
*           2.  double errors[2]
*               On exit:  the address of an array that will contain the
* 	      discontinuity errors.  errors[0] contains the positional
* 	      discontinuity.  errors[1] contains the tangent discontinuity,
* 	      in radians.
* 
* Functions referenced by PowSurf_error() are:
*  add_vect1()
*  cross1()
*  distance()
*  dot()
*  scale_vect1()
*  unitvector1()
*
* Functions that reference PowSurf_error() are:
*  ReadIgesPowerSurf()
*/

void PowSurf_error(PowSurf *pgm, double *errors)
{
  vector vt0, vt1, vt2, vt3, val, vs, vt, vn;
  double dist, u, v;
  int j, k;
  
  errors[0] = errors[1] = 0.0;

  for (j=0; j<pgm->usegmts-1; j++) { 
    for (k=0; k<pgm->vsegmts-1; k++) { 
      u = pgm->uknots[j+1] - pgm->uknots[j];   /* compute far end of */
      v = pgm->vknots[k+1] - pgm->vknots[k];   /* j-th,k-th segment */

      scale_vect1(u, pgm->contpts[15][j][k], &vt3);
      add_vect1(&vt3, pgm->contpts[14][j][k], &vt3);
      scale_vect1(u, &vt3, &vt3);
      add_vect1(&vt3, pgm->contpts[13][j][k], &vt3);
      scale_vect1(u, &vt3, &vt3);
      add_vect1(&vt3, pgm->contpts[12][j][k], &vt3);
      
      scale_vect1(u, pgm->contpts[11][j][k], &vt2);
      add_vect1(&vt2, pgm->contpts[10][j][k], &vt2);
      scale_vect1(u, &vt2, &vt2);
      add_vect1(&vt2, pgm->contpts[ 9][j][k], &vt2);
      scale_vect1(u, &vt2, &vt2);
      add_vect1(&vt2, pgm->contpts[ 8][j][k], &vt2);
      
      scale_vect1(u, pgm->contpts[ 7][j][k], &vt1);
      add_vect1(&vt1, pgm->contpts[ 6][j][k], &vt1);
      scale_vect1(u, &vt1, &vt1);
      add_vect1(&vt1, pgm->contpts[ 5][j][k], &vt1);
      scale_vect1(u, &vt1, &vt1);
      add_vect1(&vt1, pgm->contpts[ 4][j][k], &vt1);
      
      scale_vect1(u, pgm->contpts[ 3][j][k], &vt0);
      add_vect1(&vt0, pgm->contpts[ 2][j][k], &vt0);
      scale_vect1(u, &vt0, &vt0);
      add_vect1(&vt0, pgm->contpts[ 1][j][k], &vt0);
      scale_vect1(u, &vt0, &vt0);
      add_vect1(&vt0, pgm->contpts[ 0][j][k], &vt0);
      
      scale_vect1(v, &vt3, &val);
      add_vect1(&vt2, &val, &val);
      scale_vect1(v, &val, &val);
      add_vect1(&vt1, &val, &val);
      scale_vect1(v, &val, &val);    /* and compare to beginning of */
      add_vect1(&vt0, &val, &val);   /* (j+1)th,(k+1)th segment */
					
      dist = distance(&val, pgm->contpts[0][j+1][k+1]);
      errors[0] = MAX(dist, errors[0]);

      scale_vect1(1.5*v, &vt3, &vt);   /* find 1st deriv wrt t */
      add_vect1(&vt2, &vt, &vt);
      scale_vect1 (2.0*v, &vt, &vt);
      add_vect1(&vt1, &vt, &vt);

      
      scale_vect1(1.5*u, pgm->contpts[15][j][k], &vt3);   /* 1st deriv wrt s */
      add_vect1(&vt3, pgm->contpts[14][j][k], &vt3);
      scale_vect1(2.0*u, &vt3, &vt3);
      add_vect1(&vt3, pgm->contpts[13][j][k], &vt3);
      
      scale_vect1(1.5*u, pgm->contpts[11][j][k], &vt2);
      add_vect1(&vt2, pgm->contpts[10][j][k], &vt2);
      scale_vect1(2.0*u, &vt2, &vt2);
      add_vect1(&vt2, pgm->contpts[ 9][j][k], &vt2);
      
      scale_vect1(1.5*u,pgm->contpts[ 7][j][k], &vt1);
      add_vect1(&vt1, pgm->contpts[ 6][j][k], &vt1);
      scale_vect1(2.0*u, &vt1, &vt1);
      add_vect1(&vt1,   pgm->contpts[ 5][j][k], &vt1);
      
      scale_vect1(1.5*u,pgm->contpts[ 3][j][k], &vt0);
      add_vect1(&vt0, pgm->contpts[ 2][j][k], &vt0);
      scale_vect1(2.0*u, &vt0, &vt0);
      add_vect1(&vt0, pgm->contpts[ 1][j][k], &vt0);
      
      scale_vect1(v, &vt3, &vs);
      add_vect1 (&vt2, &vs, &vs);
      scale_vect1(v, &vs, &vs);
      add_vect1 (&vt1, &vs, &vs);
      scale_vect1(v, &vs, &vs);
      add_vect1 (&vt0, &vs, &vs);

      cross1(&vs,&vt, &val);
      unitvector1(&val,&val);
      cross1(pgm->contpts[1][j+1][k+1], pgm->contpts[4][j+1][k+1], &vn);
      unitvector1(&vn, &vn);
      dist = acos(dot(&val, &vn));

      if (dist < 1.0/ZERO)
	errors[1] = MAX(dist, errors[1]);
    }
  }
}

/********* PowSurf_iso() *********
* 1     Purpose
* 
*       Extract isoparameter curve from parametric spline surface.
* 
* 2     Specification
* 
*       #include "iges.h"
*       PowCurv *PowSurf_iso(PowSurf *sgeom, int index, double param,
*                            PowCurv *pc);
* 
* 3     Description
* 
*       This  function  extracts  an  isoparameter Curve  (represented  as  a
*       parametric  spline  curve) from a parametric spline surface.
* 
* 4     References
* 
*       [1]   Digital Representation for Communication of Product Definition
*             Data, US PRO/IPO-100, Initial Graphics Exchange Specification
* 	    (IGES) 5.2, IGES/PDES Organization, U.S. Product Data Association,
* 	    Fairfax, VA, November 1993.
* 
* 5     Parameters
* 
*           1.  PowSurf * sgeom
*               On entry:  the address of a structure containing the parametric
* 	      spline surface.
* 
*           2.  int index
*               On entry:  flag specifying the direction of the isoparameter
* 	      curve.  If set equal to 0, then the curve is constant in u;
* 	      otherwise, it is constant in v.
* 
*           3.  double param
*               On entry:  the parametric value of the curve.  If index is
* 	      equal to 0 then param is a constant in u; otherwise it is a
* 	      constant in v.
* 
*           4.  PowCurv * pc
*               On  entry:  the  address  of  a  structure  that  will  contain
* 	      the  isoparameter Curve.   If specified as NULL then a new
* 	      structure will be allocated by this function.
* 
*               On exit:  the address of a structure that contains the
* 	      isoparameter curve.  This same address is returned by the
* 	      function.
* 
* 6     Return Values, Error Indicators and Warnings
* 
*       The  address  of  a  structure  containing  the  isoparameter Curve  is
*       returned.   If  pc  is  not NULL on entry, then its address is returned;
*       otherwise, the address of the a newly allocated structure is returned.
* 
* 8     Further Comments
* 
*       The returned structure should be deallocated with free_pgeom().
* 
* Functions referenced by PowSurf_iso() are:
*  add_vect1()
*  find()
*  pgeomalloc()
*  scale_vect1()
*
* Functions that reference PowSurf_iso() are:
*  PowSurf_to_ParSurf_loft()
*/

PowCurv *PowSurf_iso(PowSurf *sgeom, int index, double param, PowCurv *pc)
{
  vector vt0, vt1, vt2, vt3, *val;
  double s, t;
  int j, k;
 
  if (!pc)
    pc = pgeomalloc(sgeom->vsegmts);

  if (index == 0) {
    j = find(sgeom->usegmts + 1, sgeom->uknots, param);
    if (j >= sgeom->usegmts)
      j = sgeom->usegmts - 1;
  
    s = param - sgeom->uknots[j];

    for (k=0; k<sgeom->vsegmts; k++) {
      pc->knots[k] = sgeom->vknots[k];

      scale_vect1(s, sgeom->contpts[15][j][k], pc->contpts[k][3]);
      add_vect1(pc->contpts[k][3], sgeom->contpts[14][j][k],
		pc->contpts[k][3]);
      scale_vect1(s, pc->contpts[k][3], pc->contpts[k][3]);
      add_vect1(pc->contpts[k][3], sgeom->contpts[13][j][k],
	    	pc->contpts[k][3]);
      scale_vect1(s, pc->contpts[k][3], pc->contpts[k][3]);
      add_vect1(pc->contpts[k][3], sgeom->contpts[12][j][k],
		pc->contpts[k][3]);
      
      scale_vect1(s, sgeom->contpts[11][j][k], pc->contpts[k][2]);
      add_vect1(pc->contpts[k][2], sgeom->contpts[10][j][k],
		pc->contpts[k][2]);
      scale_vect1(s, pc->contpts[k][2], pc->contpts[k][2]);
      add_vect1(pc->contpts[k][2], sgeom->contpts[9][j][k], pc->contpts[k][2]);
      scale_vect1(s, pc->contpts[k][2], pc->contpts[k][2]);
      add_vect1(pc->contpts[k][2], sgeom->contpts[8][j][k], pc->contpts[k][2]);
      
      scale_vect1(s, sgeom->contpts[7][j][k], pc->contpts[k][1]);
      add_vect1(pc->contpts[k][1], sgeom->contpts[6][j][k], pc->contpts[k][1]);
      scale_vect1(s, pc->contpts[k][1], pc->contpts[k][1]);
      add_vect1(pc->contpts[k][1], sgeom->contpts[5][j][k], pc->contpts[k][1]);
      scale_vect1(s, pc->contpts[k][1], pc->contpts[k][1]);
      add_vect1(pc->contpts[k][1], sgeom->contpts[4][j][k], pc->contpts[k][1]);
      
      scale_vect1(s, sgeom->contpts[3][j][k], pc->contpts[k][0]);
      add_vect1(pc->contpts[k][0], sgeom->contpts[2][j][k], pc->contpts[k][0]);
      scale_vect1(s, pc->contpts[k][0], pc->contpts[k][0]);
      add_vect1(pc->contpts[k][0], sgeom->contpts[1][j][k], pc->contpts[k][0]);
      scale_vect1(s, pc->contpts[k][0], pc->contpts[k][0]);
      add_vect1(pc->contpts[k][0], sgeom->contpts[0][j][k], pc->contpts[k][0]);
    }
    pc->knots[sgeom->vsegmts] = sgeom->vknots[sgeom->vsegmts];
  }
  else {
    j = find(sgeom->vsegmts + 1, sgeom->vknots, param);
    if (j >= sgeom->vsegmts)
      j = sgeom->vsegmts - 1;
  
    s = param - sgeom->vknots[j];

    for(k=0; k<sgeom->usegmts; k++) {
      pc->knots[k] = sgeom->uknots[k];

      scale_vect1(s, sgeom->contpts[15][k][j], pc->contpts[k][3]);
      add_vect1(pc->contpts[k][3], sgeom->contpts[11][k][j],
		pc->contpts[k][3]);
      scale_vect1(s, pc->contpts[k][3], pc->contpts[k][3]);
      add_vect1(pc->contpts[k][3], sgeom->contpts[7][k][j], pc->contpts[k][3]);
      scale_vect1(s, pc->contpts[k][3], pc->contpts[k][3]);
      add_vect1(pc->contpts[k][3], sgeom->contpts[3][k][j], pc->contpts[k][3]);
      
      scale_vect1(s, sgeom->contpts[14][k][j], pc->contpts[k][2]);
      add_vect1(pc->contpts[k][2], sgeom->contpts[10][k][j],
		pc->contpts[k][2]);
      scale_vect1(s, pc->contpts[k][2], pc->contpts[k][2]);
      add_vect1(pc->contpts[k][2], sgeom->contpts[6][k][j], pc->contpts[k][2]);
      scale_vect1(s, pc->contpts[k][2], pc->contpts[k][2]);
      add_vect1(pc->contpts[k][2], sgeom->contpts[2][k][j], pc->contpts[k][2]);
      
      scale_vect1(s, sgeom->contpts[13][k][j], pc->contpts[k][1]);
      add_vect1(pc->contpts[k][1], sgeom->contpts[9][k][j], pc->contpts[k][1]);
      scale_vect1(s, pc->contpts[k][1], pc->contpts[k][1]);
      add_vect1(pc->contpts[k][1], sgeom->contpts[5][k][j], pc->contpts[k][1]);
      scale_vect1(s, pc->contpts[k][1], pc->contpts[k][1]);
      add_vect1(pc->contpts[k][1], sgeom->contpts[1][k][j], pc->contpts[k][1]);
      
      scale_vect1(s, sgeom->contpts[12][k][j], pc->contpts[k][0]);
      add_vect1(pc->contpts[k][0], sgeom->contpts[8][k][j], pc->contpts[k][0]);
      scale_vect1(s, pc->contpts[k][0], pc->contpts[k][0]);
      add_vect1(pc->contpts[k][0], sgeom->contpts[4][k][j], pc->contpts[k][0]);
      scale_vect1 (s, pc->contpts[k][0], pc->contpts[k][0]);
      add_vect1(pc->contpts[k][0], sgeom->contpts[0][k][j], pc->contpts[k][0]);
    }
    pc->knots[sgeom->usegmts] = sgeom->vknots[sgeom->usegmts];
  }
  return (pc);
}

/********* PowSurf_to_ParSurf() *********
* 1     Purpose
* 
*       Convert parametric spline surface to NURBS surface.
* 
* 2     Specification
* 
*       #include "iges.h"
*       ParSurf *PowSurf_to_ParSurf(PowSurf *sgeom, ParSurf *fgeom,
*                                   double *maxerr, int byLofting, int nc);
* 
* 3     Description
* 
*       This function converts a parametric spline surface into a NURBS surface
*       with a specified degree of continuity at the breakpoints of the
*       parametric polynomials.
* 
* 4     References
* 
*       [1]   Digital Representation for Communication of Product Definition
*             Data, US PRO/IPO-100, Initial Graphics Exchange Specification
* 	    (IGES) 5.2, IGES/PDES Organization, U.S. Product Data Association,
* 	    Fairfax, VA, November 1993.
* 
* 5     Parameters
* 
*           1.  PowSurf * sgeom
*               On entry:  the address of a structure containing the parametric
* 	      spline surface.
* 
*           2.  ParSurf * fgeom
*               On entry: the address of a structure that will contain the NURBS
* 	      surface.  If specified as NULL then a new structure will be
* 	      allocated by this function.
* 
*               On  exit:  the  address  of  a  structure  that  contains  the
* 	      NURBS surface.   This same address is returned by the function.
* 
*           3.  double * maxerr
*               On exit:  the addrees of a variable that will contain the global
* 	      conversion error calculated as the maximum distance between
* 	      points on the parametric spline and the NURBS surface evaluated
* 	      at the breakpoints and the midpoints of the breakpoints.
* 
*           4.  int byLofting
* 	      On entry:  flag specifying the conversion method.  If set equal
* 	      to 1 then the surface is converted by lofting.  Each of a set of
* 	      isoparameter curves of the parametric spline surface is converted
* 	      to a NURBS curve and then the curves are lofted to form
* 	      a NURBS surface.  (See PowSurf_to_ParSurf_loft().)
* 
*               If  set  equal  to  0  then  the  surface  is  converted
* 	      directly  by  interpolating  a  surface through the node points
* 	      of the parametric spline surface. (See PowSurf_to_ParSurf_int().)
* 
*           5.  int nc
*               On entry: flag specifying the degree of continuity at the
* 	      breakpoints of the parametric polynomials.
* 
*                 nc   Continuity
*                 --   ----------
*                  0   Surfature
*                  1   Tangent
*                  2   Position
* 
* 6     Return Values, Error Indicators and Warnings
* 
*       The address of the structure containing the NURBS surface is returned.
*       If fgeom is not NULL on entry, then its address is returned; otherwise,
*       the address of the a newly allocated structure is returned.
* 
* 8     Further Comments
* 
*       The returned structure should be deallocated with free_fgeom().
* 
* Functions referenced by PowSurf_to_ParSurf() are:
*  fgeomalloc1()
*  PowSurf_to_ParSurf_int()
*  PowSurf_to_ParSurf_loft()
*
* Functions that reference PowSurf_to_ParSurf() are:
*  ReadIgesPowerSurf()
*/

ParSurf *PowSurf_to_ParSurf(PowSurf *sgeom, ParSurf *fgeom, double *maxerr,
			    int byLofting, int nc)
{
  if (!fgeom)
    fgeom = fgeomalloc1(4, 4, (sgeom->usegmts-1)*nc + 4,
			(sgeom->vsegmts-1)*nc + 4);

  if (byLofting)
    PowSurf_to_ParSurf_loft(sgeom, fgeom, maxerr, nc);
  else 
    PowSurf_to_ParSurf_int(sgeom, fgeom, maxerr, nc);

  return(fgeom);
}  

/********* PowSurf_to_ParSurf_int() *********
* 1     Purpose
* 
*       Convert parametric spline surface to NURBS surface by interpolation.
* 
* 2     Specification
* 
*       #include "iges.h"
*       ParSurf  *PowSurf_to_ParSurf_int(PowSurf  *sgeom,  ParSurf  *fgeom,
*                                        double  *maxerr,  int nc);
* 
* 3     Description
* 
*       This function converts a parametric spline surface into a NURBS with a
*       specified degree of continuity at the breakpoints of the parametric
*       polynomials.  The conversion method uses interpolation.
* 
* 4     References
* 
*       [1]   Digital Representation for Communication of Product Definition
*             Data, US PRO/IPO-100, Initial Graphics Exchange Specification
* 	    (IGES) 5.2, IGES/PDES Organization, U.S. Product Data Association,
* 	    Fairfax, VA, November 1993.
* 
* 5     Parameters
* 
*           1.  PowSurf * sgeom
*               On entry:  the address of a structure containing the parametric
* 	      spline surface.
* 
*           2.  ParSurf * fgeom
*               On entry: the address of a structure that will contain the NURBS
* 	      surface.  If specified as NULL then a new structure will be
* 	      allocated by this function.
* 
*               On  exit:  the  address  of  a  structure  that  contains  the
* 	      NURBS surface.  This same address is returned by the function.
* 
*           3.  double * maxerr
*               On exit:  the addrees of a variable that will contain the
* 	      global conversion error calculated as the maximum distance
* 	      between points on the parametric spline and the NURBS surface
* 	      evaluated at the breakpoints and the midpoints of the
* 	      breakpoints.
* 
*           4.  int nc
*               On entry: flag specifying the degree of continuity at the
* 	      breakpoints of the parametric polynomials.
* 
*                 nc   Continuity
*                 --   ----------
*                  0   Surfature
*                  1   Tangent
*                  2   Position
* 
* 6     Return Values, Error Indicators and Warnings
* 
*       The address of the structure containing the NURBS surface is returned.
*       If fgeom is not NULL on entry, then its address is returned; otherwise,
*       the address of the a newly allocated structure is returned.
* 
* 8     Further Comments
* 
*       The returned structure should be deallocated with free_fgeom().
* 
* Functions referenced by PowSurf_to_ParSurf_int() are:
*  calc_gram_surf()
*  calc_points_surf()
*  copyvector()
*  dbl_array1()
*  dbl_array2()
*  free_darray1()
*  free_darray2()
*  knot_normalize()
*  make_knots()
*  nodes_bsp()
*  PowParSurf_compare()
*  solve_gram()
*
* Functions that reference PowSurf_to_ParSurf_int() are:
*  PowSurf_to_ParSurf()
*/

void PowSurf_to_ParSurf_int(PowSurf *sgeom, ParSurf *fgeom, double *maxerr,
			    int nc)
{
  double2 *u, *v, *unodes, *vnodes;
  double2 **xyz, **basis;
  int m1, m2, mpoints, npoints, tpoints;
  int i, j;
  
  fgeom->type = PSurfaceOpen;

  npoints = (sgeom->usegmts-1)*nc + 4;
  mpoints = (sgeom->vsegmts-1)*nc + 4;
  tpoints = npoints*mpoints;
  
  u = dbl_array1(tpoints);
  v = dbl_array1(tpoints);
  unodes = dbl_array1(npoints);
  vnodes = dbl_array1(mpoints);
  xyz = dbl_array2(tpoints,4);
  basis = dbl_array2(tpoints, tpoints);
  
  make_knots(fgeom->uknots, sgeom->uknots, sgeom->usegmts, nc);
  nodes_bsp(fgeom->uknots, fgeom->ucontpts, 4, unodes);
  make_knots(fgeom->vknots, sgeom->vknots, sgeom->vsegmts, nc);
  nodes_bsp(fgeom->vknots, fgeom->vcontpts, 4, vnodes);

  calc_points_surf(unodes, vnodes, sgeom, xyz, u, v, nc);

  knot_normalize(fgeom->uknots, fgeom->uorder + fgeom->ucontpts);
  knot_normalize(u, tpoints);
  knot_normalize(fgeom->vknots, fgeom->vorder + fgeom->vcontpts);
  knot_normalize(v, tpoints);

  calc_gram_surf(u, v, fgeom->uknots, fgeom->vknots, 4, 4, mpoints, npoints,
		 basis, &m1, &m2);
  solve_gram(basis, xyz, tpoints, m1, m2);

  for (i=0; i<npoints; i++)
    for (j=0; j<mpoints; j++)
      copyvector((vector *)xyz[i*mpoints+j], fgeom->contpts[i][j]);

  PowParSurf_compare(fgeom, sgeom, maxerr, nc);

  free_darray1(unodes);
  free_darray1(vnodes);
  free_darray2(xyz);
  free_darray2(basis);
  free_darray1(u);
  free_darray1(v);
}

/********* PowSurf_to_ParSurf_loft() *********
* 1     Purpose
* 
*       Convert parametric spline surface to NURBS surface by lofting.
* 
* 2     Specification
* 
*       #include "iges.h"
*       ParSurf  *PowSurf_to_ParSurf_loft(PowSurf  *sgeom,  ParSurf  *fgeom,
*                                         double  *maxerr,  int nc);
* 
* 3     Description
* 
*       This function converts a parametric spline surface into a NURBS with a
*       specified degree of continuity at the breakpoints of the parametric
*       polynomials.  The conversion method first extracts and converts a series
*       of isoparameter curves and then lofts a surface through the curves.
* 
* 4     References
* 
*       [1]   Digital Representation for Communication of Product Definition
*             Data, US PRO/IPO-100, Initial Graphics Exchange Specification
* 	    (IGES) 5.2, IGES/PDES Organization, U.S. Product Data Association,
* 	    Fairfax, VA, November 1993.
* 
* 5     Parameters
* 
*           1.  PowSurf * sgeom
*               On entry:  the address of a structure containing the parametric
* 	      spline surface.
* 
*           2.  ParSurf * fgeom
*               On entry: the address of a structure that will contain the NURBS
* 	      surface.  If specified as NULL then a new structure will be
* 	      allocated by this function.
* 
*               On  exit:  the  address  of  a  structure  that  contains  the
* 	      NURBS surface.  This same address is returned by the function.
* 
*           3.  double * maxerr
*               On exit:  the addrees of a variable that will contain the
* 	      global conversion error calculated as the maximum distance
* 	      between points on the parametric spline and the NURBS surface
* 	      evaluated at the breakpoints and the midpoints of the
* 	      breakpoints.
* 
*           4.  int nc
*               On entry: flag specifying the degree of continuity at the
* 	      breakpoints of the parametric polynomials.
* 
*                 nc   Continuity
*                 --   ----------
*                  0   Surfature
*                  1   Tangent
*                  2   Position
* 
* 6     Return Values, Error Indicators and Warnings
* 
*       The address of the structure containing the NURBS surface is returned.
*       If fgeom is not NULL on entry, then its address is returned; otherwise,
*       the address of the a newly allocated structure is returned.
* 
* 8     Further Comments
* 
*       The returned structure should be deallocated with free_fgeom().
* 
* Functions referenced by PowSurf_to_ParSurf_loft() are:
*  copyfgeom()
*  dbl_array1()
*  free_darray1()
*  free_egeom()
*  free_fgeom()
*  free_pgeom()
*  gen_array1()
*  knot_normalize()
*  loft_rational()
*  make_knots()
*  nodes_bsp()
*  pgeomalloc()
*  PowCurv_to_ParCurv1()
*  PowParSurf_compare()
*  PowSurf_iso()
*
* Functions that reference PowSurf_to_ParSurf_loft() are:
*  PowSurf_to_ParSurf()
*/

void PowSurf_to_ParSurf_loft(PowSurf *sgeom, ParSurf *fgeom, double *maxerr,
			     int nc)
{
  ParSurf *lgeom;
  ParCurv **egeom;
  PowCurv *pc;
  double *nodes;
  int index, nseg;
  int i, j;

  nodes = dbl_array1((sgeom->usegmts-1)*nc + 4);
  make_knots(fgeom->uknots, sgeom->uknots, sgeom->usegmts, nc);
  nodes_bsp(fgeom->uknots, fgeom->ucontpts, 4, nodes);
  nseg = (sgeom->usegmts-1)*nc+ 4;

  egeom = (ParCurv **)gen_array1(nseg, sizeof(ParCurv *));

  pc = pgeomalloc(sgeom->vsegmts);
  for (i=0; i<nseg; i++) {
    PowSurf_iso(sgeom, 0, nodes[i], pc);
    printf(" Isoparameter curve %d extracted;", i);
    egeom[i] = PowCurv_to_ParCurv1(pc, 0, maxerr, nc);
    printf(" conversion error = %g\n", *maxerr);    
  }
  free_pgeom(pc);

  knot_normalize (fgeom->uknots, 4 + nseg);
  knot_normalize (nodes, nseg);
  lgeom = loft_rational(egeom, nodes, fgeom->uknots, 4, nseg);
  copyfgeom(lgeom, fgeom);

  free_fgeom(lgeom);
  fgeom->type = PSurfaceOpen;
  for (i=0; i<nseg; i++)
    free_egeom(egeom[i]);

  PowParSurf_compare(fgeom, sgeom, maxerr, nc);

  free_darray1(nodes);
}

/********* PowSurf_to_ParSurf2() *********
* 1     Purpose
* 
*       Convert parametric spline surface to NURBS surface.
* 
* 2     Specification
* 
*       #include "iges.h"
*       ParSurf *PowSurf_to_ParSurf2(PowSurf *sgeom, ParSurf *fgeom,
*                                    double *maxerr, int nc);
* 
* 3     Description
* 
*       This function converts a parametric spline surface into a NURBS surface
*       with a specified degree of continuity at the breakpoints of the
*       parametric polynomials.
* 
* 4     References
* 
*       [1]   Digital Representation for Communication of Product Definition
*             Data, US PRO/IPO-100, Initial Graphics Exchange Specification
* 	    (IGES) 5.2, IGES/PDES Organization, U.S. Product Data Association,
* 	    Fairfax, VA, November 1993.
* 
* 5     Parameters
* 
*           1.  PowSurf * sgeom
*               On entry:  the address of a structure containing the parametric
* 	      spline surface.
* 
*           2.  ParSurf * fgeom
*               On entry: the address of a structure that will contain the NURBS
* 	      surface.  If specified as NULL then a new structure will be
* 	      allocated by this function.
* 
*               On  exit:  the  address  of  a  structure  that  contains  the
* 	      NURBS surface.   This same address is returned by the function.
* 
*           3.  double * maxerr
*               On exit:  the addrees of a variable that will contain the global
* 	      conversion error calculated as the maximum distance between
* 	      points on the parametric spline and the NURBS surface evaluated
* 	      at the breakpoints and the midpoints of the breakpoints.
* 
*           4.  int nc
*               On entry: flag specifying the degree of continuity at the
* 	      breakpoints of the parametric polynomials.
* 
*                 nc   Continuity
*                 --   ----------
*                  0   Surfature
*                  1   Tangent
*                  2   Position
* 
* 6     Return Values, Error Indicators and Warnings
* 
*       The address of the structure containing the NURBS surface is returned.
*       If fgeom is not NULL on entry, then its address is returned; otherwise,
*       the address of the a newly allocated structure is returned.
* 
* 8     Further Comments
* 
*       The returned structure should be deallocated with free_fgeom().
* 
* Functions referenced by PowSurf_to_ParSurf2() are:
*  monotobern2()
*
* Functions that reference PowSurf_to_ParSurf2() are:
*  ReadIgesPowerSurf()
*/


ParSurf *PowSurf_to_ParSurf2 (PowSurf *sgeom, ParSurf *fgeom, double *maxerr,
			      int nc)
{
  ParSurf *egm;
  double **c, coeff;
  int i, ii, j, jj, k, kk, m, mm, n, nn;

  if (!fgeom)			/* allocate nurbs if nil */
    fgeom = fgeomalloc1(sgeom->uorder, sgeom->vorder,
			sgeom->uorder*sgeom->usegmts,
			sgeom->vorder*sgeom->vsegmts);
  fgeom->type = PSurfaceOpen;

  /* convert monomial break points to Bernstein knot vector */

  for (k=0; k<=sgeom->usegmts; k++)
    for (i=0; i<sgeom->uorder; i++)
      fgeom->uknots[k*sgeom->uorder + i] = sgeom->uknots[k];
  for (k=0; k<=sgeom->vsegmts; k++)
    for (i=0; i<sgeom->vorder; i++)
      fgeom->vknots[k*sgeom->vorder + i] = sgeom->vknots[k];

  c = dbl_array2(sgeom->uorder, sgeom->vorder);
  m = sgeom->uorder - 1;
  mm= sgeom->vorder - 1;

  /* convert the monomial coefficients to Bernstein control polygon */

  for (k=0; k<sgeom->usegmts; k++) {
    n = k*fgeom->uorder;

    for (kk=0; kk<sgeom->vsegmts; kk++) {
      nn = kk*fgeom->vorder;

      for (j=0; j<3; j++) {
	for (i=0; i<sgeom->uorder; i++)
	  for (ii=0; ii<sgeom->vorder; ii++) {
	    jj = sgeom->vorder*ii + i;

	    if (j == 0)
	      c[i][ii] = sgeom->contpts[jj][k][kk]->x;
	    else if (j == 1)
	      c[i][ii] = sgeom->contpts[jj][k][kk]->y;
	    else
	      c[i][ii] = sgeom->contpts[jj][k][kk]->z;

	    /* rescale the parameter space of to [0,1] for each segment */
	    
	    c[i][ii] /= pow(sgeom->usegmts, i);
	    c[i][ii] /= pow(sgeom->vsegmts, ii);
	  }

	for (i=0; i<sgeom->uorder; i++)
	  for (ii=0; ii<sgeom->vorder; ii++) {
	    coeff = monotobern2(m, mm, i, ii, c);

	    if (j == 0)
	      fgeom->contpts[n+i][nn+ii]->x = coeff;
	    else if (j == 1)
	      fgeom->contpts[n+i][nn+ii]->y = coeff;
	    else {
	      fgeom->contpts[n+i][nn+ii]->z = coeff;
	      fgeom->contpts[n+i][nn+ii]->w = 1.0;
	    }
	  }
      }
    }
  }

  free_darray2(c);

  /* remove knots so that each unique internal knot value has multiplicity
   * of 1
   */

  RemoveSurfKnotU(fgeom);
  RemoveSurfKnotV(fgeom);
  PowParSurf_compare2(fgeom, sgeom, maxerr, nc);

  return fgeom;
}

/*****************************************************************************
*                               monotobern2()
******************************************************************************
* 
* 1   Purpose
*     This function converts a monomial to the Bernstein form.
*     This function is called once for each i = 0,1,...,m-1, and j=0,1,...,n-1,
*     for each spatial coordinate.
* 
* 2   Specification
*     #include "iges.h"
*     double monotobern2(int m, int n, int i, int j, double **c)
* 
* 3   Description
*     This routine converts a monomial basis representation of an algebraic 
*     surface to the Bernstein basis.
* 
* 4   References
*     Not applicable
* 
* 5   Parameters
*       1. int m
*         On entry: the maximum degree in u of the input monomial
*       2. int n
*         On entry: the maximum degree in v of the input monmial
*       3. int i
*         On entry: degree in u for which the Bernstein coefficient is
*                   required.
*       4. int j
*         On entry: degree in v for which the Bernstein coefficient is
*                   required.
*       5. double ** c
*         On entry: array of length m+1 by n+1 holding the control polygon
* 	   of the monomial
* 
* 6   Return Values, Error Indicators and Warnings
*     Return the Bernstein coefficient corresponding to
*     the input monomial.
* 
* 7   Accuracy
*     Not applicable.
* 
* 8   Further Comments
*     Not applicable.
* 
* 9   Functions referenced by monotobern2() are:
*     factorial()
* 
* 10  Functions that reference monotobern2() are: ParSurf_to_PowSurf2()
* 
******************************************************************************/

double monotobern2(int m, int n, int i, int j, double **c)
/* Converts i,j - monomial coefficient to bernstein basis.
 * n,m : degrees in u,v
 * i,j : coeff. of pow(u,i),  pow(v,j).
*/
{
  int k,l;
  double2 C,f1,f2,fi,fj,fn,fm;

  fm = factorial(m); fn = factorial(n);
  fi = factorial(i); fj = factorial(j);

  C=0;
  for (k=0; k<=i; k++) {
    f1 = factorial(m-k)/factorial(i-k);
    for (l=0; l<=j; l++) {
      f2 = factorial(n-l)/factorial(j-l);
      C += (f1*f2*c[k][l]);
    }
  }
     
  C *= (fi*fj/(fm*fn));
     
  return C;
}

/*****************************************************************************
*                               RemoveSurfKnotU()
******************************************************************************
* 
* 1   Purpose
*     This function removes multiple copies of U knots from
*     a B-spline surface.
* 
* 2   Specification
*     #include "iges.h"
*     int RemoveSurfKnotU(ParSurf *fgeom);
* 
* 3   Description
*     This function removes multiple copies U knots from
*     a B-spline surface.
* 
* 4   References
*     [1] W. Tiller, "Knot-Removal Algorithms for NURBS Curves and Surfaces,"
*         Computer Aided Design, 24(8):445-453, August 1992.
* 
* 5   Parameters
*       1. ParSurf * fgeom
*         On entry: the initial B-spline surface.
*         On exit: the modified B-spline surface.
* 
* 6   Return Values, Error Indicators and Warnings
*     The actual number of U knots removed is returned.
* 
* 7   Accuracy
*     Not applicable.
* 
* 8   Further Comments
*     Not applicable.
* 
* 9   Functions referenced by RemoveSurfKnotU() are:
*     addv()
*     copyvector()
*     dbl_array1()
*     distance()
*     multv()
*     subv()
*     vect_array2()
* 
* 10  Functions that reference RemoveSurfKnotU() are:
*     PowSurf_to_ParSurf2()
* 
******************************************************************************/

int RemoveSurfKnotU(ParSurf *fgeom)
{
  int ord, hispan, gap, s, mult, fout, last, first, bgap, agap, remflag;
  double hiu, u, *alfas, *uKt, alfi, alfj;
  double TOL = 1.0e-6;
  int t, i, j, k, lf2, m, ki, kj, n, p, r;
  vector ***temp;

  m = fgeom->vcontpts-1;
  n = fgeom->ucontpts-1;
  p = fgeom->uorder-1;
  r = fgeom->uorder+fgeom->ucontpts-1;

  alfas = dbl_array1(2*p+1);
  temp = vec_array2(m+1,2*p+1);
  uKt = dbl_array1(r+1);

  uKt = fgeom->uknots;

  ord=p+1;
  hispan = r-ord;
  if(hispan < ord) return 0; 
  hiu = uKt[hispan];
  gap = 0;
  u = uKt[ord];
  s = ord;
  while (u == uKt[s+1]) s++;
  mult = s-p;
  fout = (2*s-mult-p)/2;
  last = s-mult;
  first = mult;
  bgap = s;
  agap = bgap + 1;
  while (1) {
    for(t=0;t<mult;first--,last++,t++) {
      remflag = 1;
      lf2 = last-first+2;
      i = first;
      j=last;
      while(j-i>t) {
	alfas[i-first+1] = (u-uKt[i])/(uKt[i+ord+gap+t]-uKt[i]);
	alfas[j-first+1] = (u-uKt[j-t])/(uKt[j+ord+gap]-uKt[j-t]);
	i++;
	j--;
      }
      if (j-i == t) {
	alfas[i-first+1] = (u-uKt[i])/(uKt[i+ord+gap+t]-uKt[i]);
      }	
      for(k=0;k<=m;k++) {
	ki = first-1;
	kj = last+1;
	/*
	printf("kj = %d, k = %d, lf2 = %d\n",kj,k,lf2);
	printf("fgeom->contpts[%d][%d] = %lf, %lf, %lf\n",kj,k,
	fgeom->contpts[kj][k]->x,fgeom->contpts[kj][k]->y,
	fgeom->contpts[kj][k]->z);
	*/
	copyvector(fgeom->contpts[ki++][k],temp[k][0]);
	copyvector(fgeom->contpts[kj--][k],temp[k][lf2]);
	i=1;
	j=lf2-1;
	while(j-i>t) {
	  alfi = alfas[i];
	  alfj = alfas[j];
	  copyvector(multv(1/alfi,subv(fgeom->contpts[ki++][k],
			      multv((1.-alfi),temp[k][i-1]))),temp[k][i]);
	  copyvector(multv(1/(1.-alfj),subv(fgeom->contpts[kj--][k],
			      multv(alfj,temp[k][j+1]))),temp[k][j]);
	  i++;
	  j--;
	}
	if(j-i<t) {
	  if(distance(temp[k][i-1], temp[k][j+1])>TOL) {
	    remflag = 0;
	    break;
	  }
	}
	else {
	  alfi = alfas[i];
	  if(distance(fgeom->contpts[ki][k],addv(multv(alfi,temp[k][i+t+1]),
				     multv((1.-alfi),temp[k][i-1])))>TOL) {
	    remflag = 0;
	    break;
	  }
	}
      }
      if(remflag == 0) break;
      for(k=0;k<=m;k++) {
	ki = first;
	kj = last;
	i = 1;
	j = lf2-1;
	while(j-i>t){
	  copyvector(temp[k][i++],fgeom->contpts[ki++][k]);
	  copyvector(temp[k][j--],fgeom->contpts[kj--][k]);
	}
      }
    }
    if(t>0) {
      j = fout;
      i = j;
      for(k=1;k<t;k++)
	if(k%2 == 1) i++; else j--;
      for(i=i+1,k=0;k<=m;k++)
	for(kj=j,ki=i;ki<=bgap;ki++)
	  copyvector(fgeom->contpts[ki][k],fgeom->contpts[kj++][k]);
      bgap -= t;
    }
    if(u==hiu) {
      gap +=t;
      break;
    }
    else{
      j = i = s-t+1;
      k = s+gap+1;
      u = uKt[k];
      while (u == uKt[k]) uKt[i++] = uKt[k++];
      mult = i-j;
      s = i-1;
      gap += t;
      for(ki = bgap+1,k=0;k<=m;k++)
	for(i=ki,j=0;j<mult;j++)
	  copyvector(fgeom->contpts[agap+j][k],fgeom->contpts[i++][k]);
      bgap += mult; agap += mult;
      fout = (2*s-p-mult)/2;
      last = s-mult;
      first = s-p;
    }
  }
  if(gap>0)
    for(i=hispan+1,k=i-gap,j=1;j<=ord;j++)
      uKt[k++] = uKt[i++];

  fgeom->ucontpts -= gap;

  return gap;
}

/*****************************************************************************
*                               RemoveSurfKnotV()
******************************************************************************
* 
* 1   Purpose
*     This function removes multiple copies of V knots from
*     a B-spline surface.
* 
* 2   Specification
*     #include "iges.h"
*     int RemoveSurfKnotV(ParSurf *fgeom);
* 
* 3   Description
*     This function removes multiple copies of V knots from
*     a B-spline surface.
* 
* 4   References
*     [1] W. Tiller, "Knot-Removal Algorithms for NURBS Curves and Surfaces,"
*         Computer Aided Design, 24(8):445-453, August 1992.
* 
* 5   Parameters
*       1. ParSurf * fgeom
*         On entry: the initial B-spline surface.
*         On exit: the modified B-spline surface.
* 
* 6   Return Values, Error Indicators and Warnings
*     The actual number of V knots removed is returned.
* 
* 7   Accuracy
*     Not applicable.
* 
* 8   Further Comments
*     Not applicable.
* 
* 9   Functions referenced by RemoveSurfKnotV() are:
*     addv()
*     copyvector()
*     dbl_array1()
*     distance()
*     multv()
*     subv()
*     vect_array2()
* 
* 10  Functions that reference RemoveSurfKnotV() are:
*     PowSurf_to_ParSurf2()
* 
******************************************************************************/

int RemoveSurfKnotV(ParSurf *fgeom)
{
  int ord, hispan, gap, s, mult, fout, last, first, bgap, agap, remflag;
  double hiv, v, *alfas, *vKt, alfi, alfj;
  double TOL = 1.0e-6;
  int t, i, j, k, lf2, m, ki, kj, n, p, r;
  vector ***temp;

  m = fgeom->ucontpts-1;
  n = fgeom->vcontpts-1;
  p = fgeom->vorder-1;
  r = fgeom->vorder+fgeom->vcontpts-1;

  alfas = dbl_array1(2*p+1);
  temp = vec_array2(m+1,2*p+1);
  vKt = dbl_array1(r+1);

  vKt = fgeom->vknots;

  ord=p+1;
  hispan = r-ord;
  if(hispan < ord) return 0; 
  hiv = vKt[hispan];
  gap = 0;
  v = vKt[ord];
  s = ord;
  while (v == vKt[s+1]) s++;
  mult = s-p;
  fout = (2*s-mult-p)/2;
  last = s-mult;
  first = mult;
  bgap = s;
  agap = bgap + 1;
  while (1) {
    for(t=0;t<mult;first--,last++,t++) {
      remflag = 1;
      lf2 = last-first+2;
      i = first;
      j=last;
      while(j-i>t) {
	alfas[i-first+1] = (v-vKt[i])/(vKt[i+ord+gap+t]-vKt[i]);
	alfas[j-first+1] = (v-vKt[j-t])/(vKt[j+ord+gap]-vKt[j-t]);
	i++;
	j--;
      }
      if (j-i == t) {
	alfas[i-first+1] = (v-vKt[i])/(vKt[i+ord+gap+t]-vKt[i]);
      }	
      for(k=0;k<=m;k++) {
	ki = first-1;
	kj = last+1;
	/*
	printf("kj = %d, k = %d, lf2 = %d\n",kj,k,lf2);
	printf("fgeom->contpts[%d][%d] = %lf, %lf, %lf\n",k,kj,
	fgeom->contpts[k][kj]->x,fgeom->contpts[k][kj]->y,
	fgeom->contpts[k][kj]->z);
	*/
	copyvector(fgeom->contpts[k][ki++],temp[k][0]);
	copyvector(fgeom->contpts[k][kj--],temp[k][lf2]);
	i=1;
	j=lf2-1;
	while(j-i>t) {
	  alfi = alfas[i];
	  alfj = alfas[j];
	  copyvector(multv(1/alfi,subv(fgeom->contpts[k][ki++],
			      multv((1.-alfi),temp[k][i-1]))),temp[k][i]);
	  copyvector(multv(1/(1.-alfj),subv(fgeom->contpts[k][kj--],
			      multv(alfj,temp[k][j+1]))),temp[k][j]);
	  i++;
	  j--;
	}
	if(j-i<t) {
	  if(distance(temp[k][i-1], temp[k][j+1])>TOL) {
	    remflag = 0;
	    break;
	  }
	}
	else {
	  alfi = alfas[i];
	  if(distance(fgeom->contpts[k][ki],addv(multv(alfi,temp[k][i+t+1]),
				     multv((1.-alfi),temp[k][i-1])))>TOL) {
	    remflag = 0;
	    break;
	  }
	}
      }
      if(remflag == 0) break;
      for(k=0;k<=m;k++) {
	ki = first;
	kj = last;
	i = 1;
	j = lf2-1;
	while(j-i>t){
	  copyvector(temp[k][i++],fgeom->contpts[k][ki++]);
	  copyvector(temp[k][j--],fgeom->contpts[k][kj--]);
	}
      }
    }
    if(t>0) {
      j = fout;
      i = j;
      for(k=1;k<t;k++)
	if(k%2 == 1) i++; else j--;
      for(i=i+1,k=0;k<=m;k++)
	for(kj=j,ki=i;ki<=bgap;ki++)
	  copyvector(fgeom->contpts[k][ki],fgeom->contpts[k][kj++]);
      bgap -= t;
    }
    if(v==hiv) {
      gap +=t;
      break;
    }
    else{
      j = i = s-t+1;
      k = s+gap+1;
      v = vKt[k];
      while (v == vKt[k]) vKt[i++] = vKt[k++];
      mult = i-j;
      s = i-1;
      gap += t;
      for(ki = bgap+1,k=0;k<=m;k++)
	for(i=ki,j=0;j<mult;j++)
	  copyvector(fgeom->contpts[k][agap+j],fgeom->contpts[k][i++]);
      bgap += mult; agap += mult;
      fout = (2*s-p-mult)/2;
      last = s-mult;
      first = s-p;
    }
  }
  if(gap>0)
    for(i=hispan+1,k=i-gap,j=1;j<=ord;j++)
      vKt[k++] = vKt[i++];

  fgeom->vcontpts -= gap;

  return gap;
}

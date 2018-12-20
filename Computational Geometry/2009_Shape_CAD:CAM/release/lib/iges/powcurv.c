/* Copyright (C) 1992 Massachusetts Institute of Technology, Cambridge, MA
   All rights reserved
*/

/* powcurv.c */
/* based on: /usr/deslab/epraxiteles/praxlib/iges/nrb_to_pow.c */
/*           /usr/deslab/epraxiteles/praxlib/iges/pow_to_nrb.c */
/*           /usr/deslab/epraxiteles/praxlib/iges/geom112.c */

/* berntomono1()
 * calc_points_curv()
 * comb()
 * make_knots()
 * monotobern1()
 * ParCurv_to_PowCurv()
 * ParCurv_to_PowCurv1()
 * PowCurv_error()
 * PowCurv_to_ParCurv()
 * PowCurv_to_ParCurv1()
 * PowParCurv_compare()
 * RemoveCurveKnot()
*/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <malloc.h>
#include "iges.h"   /* must be before "gen.h" is referenced */
#include "appr.h"
#include "bspl.h"
#include "editor.h"
#include "gen.h"

/********* calc_points_curv() *********
* 1     Purpose
* 
*       Evaluate points on a parametric spline curve.
* 
* 2     Specification
* 
*       #include "iges.h"
*       void calc_points_curv(double *nodes, PowCurv *pgeom, double **xyz,
*                             int nc);
* 
* 3     Description
* 
*       This function evaluates a parametric spline curve at its nodal values.
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
*           1.  double * nodes
*               On entry:  the address of an array containing the nodes of the
* 	      curve.
* 
*           2.  PowCurv * pgeom
*               On entry:  the pointer to a structure containing the parametric
* 	      spline curve.
* 
*           3.  double *** xyz
*               On entry:  the address of a triply-dimensioned array that will
* 	      contain the evaluated points.
* 
*               On exit:  the address of a triply-dimensioned array containing
* 	      the evaluated points.
* 
*           4.  int nc
*               On entry:  flag specifying the continuity of the curve.
* 
*                 nc   Continuity
*                 --   ----------
*                  0   Curvature
*                  1   Tangent
*                  2   Position
* 
* Functions referenced by calc_points_curv() are:
*  copyvector()
*  PowCurv_eval()
*  vectfree()
*
* Functions that reference calc_points_curv() are:
*  PowCurv_to_ParCurv()
*/

void calc_points_curv(double *nodes, PowCurv *pgeom, double **xyz, int nc)
{
  vector *v;
  int i;
    
  for (i=0; i<(pgeom->nsegmts-1)*nc + 4; i++) {
    v = PowCurv_eval (pgeom, nodes[i], 0);
    copyvector (v, (vector *)xyz[i]);
    vectfree(v);
  }
}

/********* make_knots() *********
* 1     Purpose
* 
*       Copy knots from parametric spline to NURBS curve.
* 
* 2     Specification
* 
*       #include "iges.h"
*       void make_knots(double *bknots, double *pknots, int nseg, int nc);
* 
* 3     Description
* 
*       This function assigns the knot vector of a NURBS curve based on the
*       break points of a parametric spline curve.
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
*           1.  double * bknots
*               On entry:  the address of the knot vector of the NURBS curve.
* 
*               On exit:  the address of the knot vector of the NURBS curve,
* 	      to which values have been assigned.
* 
* 
*           2.  double * pknots
*               On entry:  the address of the knot vector of the parametric
* 	      spline curve.
* 
*           3.  int nseg
*               On entry:  the number of break points of the parametric spline
* 	      curve.
* 
*           4.  int nc
*               On entry:  the degree of continuity at the knot values of the
* 	      NURBS curve.
* 
*                 nc   Continuity
*                 --   ----------
*                  0   Curvature
*                  1   Tangent
*                  2   Position
* 
* 8     Further Comments
* 
*       The knot vector values are assigned as follows:
* 
*       bknots[i]               = pknots[0]    for 0 <= i <= 3
*       bknots[i + nc(i-4) + j] = pknots[i-3]  for 4 <= i <= nseg+2 and
*                                                  0 <= j <= nc
*       bknots[i + nc(nseg-1)]  = pknots[nseg] for nseg+3 <= i <= nseg+6
* 
* Functions that reference make_knots() are:
*  PowCurv_to_ParCurv()
*  PowSurf_to_ParSurf_int()
*  PowSurf_to_ParSurf_loft()
*/

void make_knots(double *bknots, double *pknots, int nseg, int nc)
{
  int i, j, k;
  
  for (i=k=0; i<4; i++,k++)
    bknots[k] = pknots[0];
  for (i=4; i<nseg+3; i++)
    for (j=0; j<nc; j++,k++)
      bknots[k] = pknots[i-3];
  for (i=nseg+3; i<nseg+7; i++,k++)
    bknots[k] = pknots[nseg];
}

/********* ParCurv_to_PowCurv() *********
* 1     Purpose
* 
*       Convert NURBS curve to parametric spline curve.
* 
* 2     Specification
* 
*       #include "iges.h"
*       PowCurv *ParCurv_to_PowCurv(ParCurv *egeom, PowCurv *pgeom,
*                                   double *maxerr);
* 
* 3     Description
* 
*       This function converts a NURBS curve into a parametric spline curve.
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
*           1.  ParCurv * egeom
*               On entry:  the address of a structure containing the NURBS
* 	      curve.
* 
*           2.  PowCurv * pgeom
*               On entry:  the address of a structure that will contain the
* 	      parametric spline curve.  If specified as NULL then a new
* 	      structure will be allocated by this function.
* 
*               On  exit:  the address of a structure that contains the
* 	      parametric curve.  This same address is returned by the
* 	      function.
* 
*           3.  double * maxerr
*               On exit:  the address of a variable that will contain the
* 	      global conversion error calculated as the maximum distance
* 	      between points on the parametric spline and the NURBS curve
* 	      evaluated at the breakpoints and the midpoints of the
* 	      breakpoints.
* 
* 6     Return Values, Error Indicators and Warnings
* 
*       The address of the structure containing the parametric spline curve
*       is returned.  If pgeom is not NULL on entry, then its address is
*       returned; otherwise, the address of the a newly allocated structure
*       is returned.
* 
* 8     Further Comments
* 
*       The returned structure should be deallocated with free_pgeom().
* 
* Functions referenced by ParCurv_to_PowCurv() are:
*  copyvector()
*  dbl_array1()
*  free_darray1()
*  free_varray2()
*  pgeomalloc()
*  PowParCurv_compare()
*  rbspeval()
*  scale_vect1()
*  vectfree()
*  vec_array2()
*
* Functions that reference ParCurv_to_PowCurv() are:
*  SaveIgesCos()
*  SaveIgesCurv()
*  SaveIgesSurf()
*/

PowCurv *ParCurv_to_PowCurv(ParCurv *egeom, PowCurv *pgeom, double *maxerr)
{
  static double fact03[] = {1.0, 1.0, 2.0, 6.0};
  vector *v;
  double f;
  int i, j, nbkpts = 0;

  if (!pgeom)
    pgeom = pgeomalloc(egeom->ncontpts - egeom->order + 1);
  else if (pgeom->kmem < egeom->ncontpts - egeom->order + 2 || 
	   pgeom->pmem < egeom->ncontpts - egeom->order + 1) {
    if (pgeom->kmem)
      free_darray1(pgeom->knots);
    pgeom->knots = dbl_array1(egeom->ncontpts - egeom->order + 2);
    pgeom->kmem = egeom->ncontpts - egeom->order + 2;
    if (pgeom->pmem)
      free_varray2(pgeom->contpts, egeom->ncontpts - egeom->order + 1, 4);
    pgeom->contpts = vec_array2(egeom->ncontpts - egeom->order + 1, 4);
    pgeom->pmem = egeom->ncontpts - egeom->order + 1;
  }

  pgeom->knots[0] = egeom->knots[0];
  pgeom->order = egeom->order;

  for (i=1; i<egeom->order + egeom->ncontpts; i++) 
    if (egeom->knots[i] != egeom->knots[i-1]) {
      pgeom->knots[++nbkpts] = egeom->knots[i];

      for (j=0; j<4; j++) {
	v = rbspeval(egeom, egeom->knots[i-1], j);
	scale_vect1(1.0/fact03[j], v, v);
	copyvector(v, pgeom->contpts[nbkpts-1][j]);
	vectfree(v);
      }
    }
  PowParCurv_compare(egeom, pgeom, maxerr, 1);

  return (pgeom);
}

/********* ParCurv_to_PowCurv1() *********
* 1     Purpose
* 
*       Convert NURBS curve to parametric spline curve.
* 
* 2     Specification
* 
*       #include "iges.h"
*       PowCurv *ParCurv_to_PowCurv1(ParCurv *egeom, PowCurv *pgeom,
*                                    double *maxerr);
* 
* 3     Description
* 
*       This function converts a NURBS curve into a parametric spline curve.
*       First, convert the B-spline into the Bernstein basis by adding
*       extra knots at each internal knot value so that each internal
*       knot has "order" multiplicity.  Each Bernstein polynomial segment is
*       thus equivalent to a Bezier polynomial.  Then, convert each
*       Bernstein polynomial to the equivalent monomial basis.  Finally,
*       scale the parameter range of each monomial segment to put
*       the PowCurv into the IGES representation.
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
*           1.  ParCurv * egeom
*               On entry:  the address of a structure containing the NURBS
* 	      curve.
* 
*           2.  PowCurv * pgeom
*               On entry:  the address of a structure that will contain the
* 	      parametric spline curve.  If specified as NULL then a new
* 	      structure will be allocated by this function.
* 
*               On  exit:  the address of a structure that contains the
* 	      parametric curve.  This same address is returned by the
* 	      function.
* 
*           3.  double * maxerr
*               On exit:  the address of a variable that will contain the
* 	      global conversion error calculated as the maximum distance
* 	      between points on the parametric spline and the NURBS curve
* 	      evaluated at the breakpoints and the midpoints of the
* 	      breakpoints.
* 
* 6     Return Values, Error Indicators and Warnings
* 
*       The address of the structure containing the parametric spline curve
*       is returned.  If pgeom is not NULL on entry, then its address is
*       returned; otherwise, the address of the a newly allocated structure
*       is returned.
* 
* 8     Further Comments
* 
*       The returned structure should be deallocated with free_pgeom().
* 
* Functions referenced by ParCurv_to_PowCurv1() are:
*  berntomono1()
*  copyegeom()
*  curve_oslo1()
*  dbl_array1()
*  dbl_array2()
*  egeomalloc1()
*  free_darray1()
*  free_egeom()
*  free_varray2()
*  pgeomalloc()
*  PowParCurv_compare1()
*  vec_array2()
*
* Functions that reference ParCurv_to_PowCurv1() are:
*  SaveIgesCos()
*  SaveIgesCurv()
*  SaveIgesSurf()
*/

PowCurv *ParCurv_to_PowCurv1(ParCurv *egeom, PowCurv *pgeom, double *maxerr)
{
  ParCurv
    *bernstein;
  int
    i, j, k, m,
    nInternal,
    nKnots;
  double
    *c,
    coeff,
    **inpoly,
    **outpoly;

  /* find number of unique internal knot values */

  nInternal=0;
  for (i=egeom->order+1; i<=egeom->ncontpts; i++)
    if (i == egeom->ncontpts || egeom->knots[i] > egeom->knots[i-1])
      ++nInternal;

  /* allocate memory if necessary */

  if (!pgeom)
    pgeom = pgeomalloc(nInternal + 1);
  else if (pgeom->kmem < nInternal + 2 || 
	   pgeom->pmem < nInternal + 1) {
    if (pgeom->kmem)
      free_darray1(pgeom->knots);
    pgeom->knots = dbl_array1(nInternal + 2);
    pgeom->kmem = nInternal + 2;
    if (pgeom->pmem)
      free_varray2(pgeom->contpts, nInternal + 1, 4);
    pgeom->contpts = vec_array2(nInternal + 1, 4);
    pgeom->pmem = nInternal + 1;
  }

  pgeom->knots[0] = egeom->knots[0];
  pgeom->order = egeom->order;

  /* find number of unique internal knot values */

  nInternal=0;
  for (i=egeom->order+1; i<=egeom->ncontpts; i++)
    if (i == egeom->ncontpts || egeom->knots[i] > egeom->knots[i-1])
      pgeom->knots[++nInternal] = egeom->knots[i-1];
  pgeom->knots[nInternal+1] = egeom->knots[egeom->ncontpts];

  /* allocate space to allow for "order" knot multiplicity of the
   * internal knots
   */
  bernstein = egeomalloc1(egeom->order, egeom->order + egeom->order*nInternal);

  if (nInternal) {

    /* new knot vector with "order" knot multiplicity */

    /* "order" knots (t=0) at beginning */

    for (m=0; m<egeom->order; m++)
      bernstein->knots[m] = egeom->knots[m];

    /* internal knots */

    for (i=egeom->order+1, k=1; i<=egeom->ncontpts; i++) {
      bernstein->knots[m++] = egeom->knots[i-1];
      if (i == egeom->ncontpts || egeom->knots[i] > egeom->knots[i-1]) {
	for (j=0; j<egeom->order-k;j++) {
	  bernstein->knots[m++] = egeom->knots[i-1];
	}
	k = 1;
      }
      else
	k++;         /* count knot multiplicity */
    }

    /* "order" knots (t=1) at end */

    for (i=egeom->ncontpts; i<egeom->order+egeom->ncontpts; i++)
      bernstein->knots[m++] = egeom->knots[i];

    nKnots = egeom->order*2 + egeom->order*nInternal;

    /* use the Oslo algorithm to insert knots so that every internal
     * knot has "order" multiplicity;
     * this converts the original B-spline curve to Bernstein basis
     */

    inpoly = dbl_array2(egeom->ncontpts, 4);
    for (i=0; i<egeom->ncontpts; i++) {
      inpoly[i][0] = egeom->contpts[i]->x;
      inpoly[i][1] = egeom->contpts[i]->y;
      inpoly[i][2] = egeom->contpts[i]->z;
      inpoly[i][3] = egeom->contpts[i]->w;
    }

    outpoly = dbl_array2(nKnots-egeom->order, 4);

    curve_oslo1(egeom->order, egeom->ncontpts, nKnots, egeom->knots,
		bernstein->knots, inpoly, outpoly);

    for (i=0; i<nKnots - egeom->order; i++) {
      bernstein->contpts[i]->x = outpoly[i][0];
      bernstein->contpts[i]->y = outpoly[i][1];
      bernstein->contpts[i]->z = outpoly[i][2];
      bernstein->contpts[i]->w = outpoly[i][3];
    }
  }
  else {

    /* if 1 internal knot, then Bezier curve which is already in
     * the Bernstein basis
     */

    copyegeom(egeom, bernstein);
  }

  /* convert from the Bernstein to the monomial basis */

  c = dbl_array1(bernstein->order);
  for (k=0; k<=nInternal; k++) {
    m = k*4;

    for (j=0; j<3; j++) {     /* convert each dimension (X,Y,Z) separately */
      for (i=0; i<bernstein->order; i++) {
	if (j == 0)
	  c[i] = bernstein->contpts[m+i]->x;
	else if (j == 1)
	  c[i] = bernstein->contpts[m+i]->y;
	else
	  c[i] = bernstein->contpts[m+i]->z;
	c[i] /= bernstein->contpts[m+i]->w;
      }

      /* multiply the coefficient by (nInternal+1)^i
       * to scale the parameter range from [0,1] to the
       * portion of the parameter space between the bounding
       * knot values
       */

      for (i=0; i<bernstein->order; i++) {
	coeff = berntomono1(bernstein->order-1, i, c)*pow(nInternal+1, i);
	if (j == 0)
	  pgeom->contpts[k][i]->x = coeff;
	else if (j == 1)
	  pgeom->contpts[k][i]->y = coeff;
	else {
	  pgeom->contpts[k][i]->z = coeff;
	  pgeom->contpts[k][i]->w = 1.0;
	}
      }
    }
  }
  free_darray1(c);
  free_egeom(bernstein);

  PowParCurv_compare1(egeom, pgeom, maxerr, 1);

  return (pgeom);
}

/*****************************************************************************
*                               berntomono1()
******************************************************************************
* 
* 1   Purpose
*     This function converts a polynomial in Bernstein form to one 
*     in the monomial form.  This function is called once for each i =
*     0,1,...,m-1, for each spatial coordinate.
* 
* 2   Specification
*     #include "iges.h"
*     double berntomono1(int m, int i, double *c)
* 
* 3   Description
*     This routine converts Bernstein basis representation of an algebraic 
*     curve the monomial basis.
* 
* 4   References
*     Not applicable
* 
* 5   Parameters
*       1. int m
*         On entry: the maximum degree of the input Bernstein polynomial
*       2. int i
*         On entry: degree for which the monomial coefficient is required.
*       3. double * c
*         On entry: array of length m+1 holding the control polygon
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
* 9   Functions referenced by berntomono1() are:
*     comb()
* 
* 10  Functions that reference berntomono1() are: ParCurv_to_PowCurv1()
* 
******************************************************************************/

double berntomono1(int m, int i, double *c)
{
  int k;
  double C, cmi, minus1, exp;
     
  cmi = comb(m,i);      

  C = 0;
  for (k=0; k<=i; k++) {          
    exp = i - k;
    minus1 = pow(-1.0, exp); 
    C += (minus1*comb(i,k)*c[k]);
  }
  
  C *= cmi;
     
  return C;
}

double comb(int n, int k)
{
  register int i, j;
  register double result;
  int number;
     
  number = n - k;
     
  if (number < k)
    k = number;
     
  j = 1;
     
  result = 1.0;
     
  for (i = n-k+1; i<=n; i++){
    result = result * i;
    result = result / j++;
  }
     
  return result;
}

/********* PowCurv_error() *********
* 1     Purpose
* 
*       Discontinuity errors of parametric spline curve.
* 
* 2     Specification
* 
*       #include "iges.h"
*       void PowCurv_error(PowCurv *pgeom, double errors[2]);
* 
* 3     Description
* 
*       This function calculates the positional and tangent discontinuity
*       errors  of  a  parametric spline curve.
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
*           1.  PowCurv * pgeom
*               On entry:  the address of a structure containing the parametric
* 	      spline curve.
* 
*           2.  double errors[2]
*               On exit:  the address of an array that will contain the
* 	      discontinuity errors.  errors[0] contains the positional
* 	      discontinuity.  errors[1] contains the tangent discontinuity,
* 	      in radians.
* 
* Functions referenced by PowCurv_error() are:
*  add_vect1()
*  distance()
*  dot()
*  mag()
*  scale_vect1()
*
* Functions that reference PowCurv_error() are:
*  ReadIgesPowerCurv()
*/

void PowCurv_error(PowCurv *pgeom, double *errors)
{
  vector v;
  double t;
  int k;
  
  errors[0] = errors[1] = 0.0;

  for (k=0; k<pgeom->nsegmts - 1; k++) { 
    t = pgeom->knots[k+1] - pgeom->knots[k];   /* compute far end of */
    scale_vect1(t, pgeom->contpts[k][3], &v);  /* k-th segment */
    add_vect1(&v, pgeom->contpts[k][2], &v);
    scale_vect1(t, &v, &v);
    add_vect1(&v, pgeom->contpts[k][1], &v);
    scale_vect1(t, &v, &v);                    /* and compare to beginning */
    add_vect1(&v, pgeom->contpts[k][0], &v);   /* of (k+1)-th segment */
    errors[0] = MAX (distance (&v,pgeom->contpts[k+1][0]), errors[0]);

    scale_vect1(1.5*t, pgeom->contpts[k][3], &v);  /* compute 1st deriv of */
    add_vect1(&v, pgeom->contpts[k][2], &v);       /* far end of k-th */
    scale_vect1(2.0*t, &v, &v);                    /* segment and compare to */
    add_vect1(&v, pgeom->contpts[k][1], &v);       /* 1st deriv of beginning */
    t = acos(dot(&v,pgeom->contpts[k+1][1]) /      /* of (k+1)-th segment */
	     (mag(pgeom->contpts[k+1][1])*mag(&v)));
    
    if (t < 1.0/ZERO)
      errors[1] = MAX (t, errors[1]);
  }
}

/********* PowCurv_to_ParCurv() *********
* 1     Purpose
* 
*       Convert parametric spline curve to NURBS curve.
* 
* 2     Specification
* 
*       #include "iges.h"
*       ParCurv  *PowCurv_to_ParCurv(PowCurv  *pgeom,  ParCurv  *egeom,
*                                    double  *maxerr,  int nc);
* 
* 3     Description
* 
*       This  function  converts  a  parametric  spline  curve  into  a  NURBS
*       curve  with  a  specified degree of continuity at the breakpoints of
*       the parametric polynomials.
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
*           1.  PowCurv * pgeom
*               On entry:  the address of a structure containing the parametric
* 	      spline curve.
* 
*           2.  ParCurv * egeom
*               On entry:  the address of a structure that will contain the
* 	      NURBS curve.  If specified as NULL then a new structure will
* 	      be allocated by this function.
* 
*               On  exit:  the  address  of  a  structure  that  contains  the
* 	      NURBS  curve.   This  same address is returned by the function.
* 
*           3.  double * maxerr
*               On exit:  the addrees of a variable that will contain the
* 	      global conversion error calculated as the maximum distance
* 	      between points on the parametric spline and the NURBS curve
* 	      evaluated at the breakpoints and the midpoints of the
* 	      breakpoints.
* 
*           4.  int nc
*               On entry: flag specifying the degree of continuity at the
* 	      breakpoints of the parametric polynomials.
* 
*                 nc   Continuity
*                 --   ----------
*                  0   Curvature
*                  1   Tangent
*                  2   Position
* 
* 6     Return Values, Error Indicators and Warnings
* 
*       The  address  of  the  structure  containing  the  NURBS  curve  is
*       returned.   If  pgeom  is  not NULL on entry, then its address is
*       returned; otherwise, the address of the a newly allocated structure
*       is returned.
* 
* 8     Further Comments
* 
*       The returned structure should be deallocated with free_egeom().
* 
* Functions referenced by PowCurv_to_ParCurv() are:
*  calc_gram_curv()
*  calc_points_curv()
*  copyvector()
*  dbl_array1()
*  dbl_array2()
*  egeomalloc1()
*  free_darray1()
*  free_darray2()
*  knot_normalize()
*  make_knots()
*  nodes_bsp()
*  PowParCurv_compare()
*  solve_gram()
*
* Functions that reference PowCurv_to_ParCurv() are:
*  PowSurf_to_ParSurf_loft()
*  ReadIgesPowerCurv()
*/

ParCurv *PowCurv_to_ParCurv (PowCurv *pgeom, ParCurv *egeom, double *maxerr,
			     int nc)
{
  vector *vect;
  double **basis, *nodes, **xyz;
  int m1, m2;
  int i;

  if (!egeom)			/* allocate nurbs if nil */
    egeom = egeomalloc1(4, (pgeom->nsegmts-1)*nc + 4);
  egeom->type = PCurveOpen;

  nodes = dbl_array1(egeom->ncontpts);
  xyz = dbl_array2(egeom->ncontpts, 4);
  basis = dbl_array2(egeom->ncontpts, egeom->ncontpts);

  make_knots(egeom->knots, pgeom->knots, pgeom->nsegmts, nc);
  nodes_bsp(egeom->knots, egeom->ncontpts, 4, nodes);
  calc_points_curv(nodes, pgeom, xyz, nc);
  knot_normalize(egeom->knots, egeom->order + egeom->ncontpts);
  knot_normalize(nodes, egeom->ncontpts);

  calc_gram_curv(nodes, egeom->knots, egeom->ncontpts, basis, &m1, &m2, 4);
  solve_gram(basis, xyz, egeom->ncontpts, m1, m2);

  for (i=0; i<egeom->ncontpts; i++)
    copyvector((vector *)xyz[i], egeom->contpts[i]);
  
  PowParCurv_compare(egeom, pgeom, maxerr, nc);

  free_darray1(nodes);
  free_darray2(xyz);
  free_darray2(basis);

  return(egeom);
}

/********* PowCurv_to_ParCurv1() *********
* 1     Purpose
* 
*       Convert parametric spline curve to NURBS curve.
* 
* 2     Specification
* 
*       #include "iges.h"
*       ParCurv  *PowCurv_to_ParCurv1(PowCurv  *pgeom,  ParCurv  *egeom,
*                                    double  *maxerr,  int nc);
* 
* 3     Description
* 
*       This  function  converts  a  parametric  spline  curve  into  a  NURBS
*       curve  with  a  specified degree of continuity at the breakpoints of
*       the parametric polynomials.
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
*           1.  PowCurv * pgeom
*               On entry:  the address of a structure containing the parametric
* 	      spline curve.
* 
*           2.  ParCurv * egeom
*               On entry:  the address of a structure that will contain the
* 	      NURBS curve.  If specified as NULL then a new structure will
* 	      be allocated by this function.
* 
*               On  exit:  the  address  of  a  structure  that  contains  the
* 	      NURBS  curve.   This  same address is returned by the function.
* 
*           3.  double * maxerr
*               On exit:  the addrees of a variable that will contain the
* 	      global conversion error calculated as the maximum distance
* 	      between points on the parametric spline and the NURBS curve
* 	      evaluated at the breakpoints and the midpoints of the
* 	      breakpoints.
* 
*           4.  int nc
*               On entry: flag specifying the degree of continuity at the
* 	      breakpoints of the parametric polynomials.
* 
*                 nc   Continuity
*                 --   ----------
*                  0   Curvature
*                  1   Tangent
*                  2   Position
* 
* 6     Return Values, Error Indicators and Warnings
* 
*       The  address  of  the  structure  containing  the  NURBS  curve  is
*       returned.   If  pgeom  is  not NULL on entry, then its address is
*       returned; otherwise, the address of the a newly allocated structure
*       is returned.
* 
* 8     Further Comments
* 
*       The returned structure should be deallocated with free_egeom().
* 
* Functions referenced by PowCurv_to_ParCurv1() are:
*  monotobern1()
*  PowParCurv_compare1()
*
* Functions that reference PowCurv_to_ParCurv1() are:
*  PowSurf_to_ParSurf_loft()
*  ReadIgesPowerCurv()
*/

ParCurv *PowCurv_to_ParCurv1 (PowCurv *pgeom, ParCurv *egeom, double *maxerr,
			     int nc)
{
  ParCurv *egm;
  double *c, coeff;
  int i, j, k, m, n;

  if (!egeom)			/* allocate nurbs if nil */
    egeom = egeomalloc1(pgeom->order, pgeom->order*pgeom->nsegmts);
  egeom->type = PCurveOpen;

  /* convert monomial break points to Bernstein knot vector */

  for (k=0; k<=pgeom->nsegmts; k++)
    for (i=0; i<pgeom->order; i++) 
      egeom->knots[k*pgeom->order + i] = pgeom->knots[k];
 
  c = dbl_array1(pgeom->order);
  m = pgeom->order - 1;   

  /* convert the monomial coefficients to Bernstein control polygon */

  for (k=0; k<pgeom->nsegmts; k++)
    for (j=0; j<3; j++) {
      for (i=0; i<pgeom->order; i++) {
        if (j == 0)
          c[i] = pgeom->contpts[k][i]->x;
        else if (j == 1)
          c[i] = pgeom->contpts[k][i]->y;
        else
          c[i] = pgeom->contpts[k][i]->z;

	/* rescale the parameter space of [0,1] for each segment */

        /*  c[i] /= pow(pgeom->nsegmts, i);    */
	c[i] *= pow(pgeom->knots[k+1]-pgeom->knots[k], i);
      }

      for (i=0; i<pgeom->order; i++) {
        coeff = monotobern1(m, i, c);

        n = k*egeom->order + i;
        if (j == 0)
          egeom->contpts[n]->x = coeff;
        else if (j == 1)
          egeom->contpts[n]->y = coeff;
        else {
          egeom->contpts[n]->z = coeff;
	  egeom->contpts[n]->w = 1.0;
	}
      }
    }

  free_darray1(c);

  /* number of unique internal knot values */

  n = (egeom->ncontpts - egeom->order)/egeom->order;

  /* remove knots so that each unique internal knot value has multiplicity
   * of 1
   */

  for (k=n-1; k>=0; k--) {
    j = (k+2)*egeom->order-1;
    n = RemoveCurveKnot(egeom->ncontpts, egeom->order-1, egeom->knots,
			egeom->contpts, egeom->knots[j], j, egeom->order,
			egeom->order-1);
    egeom->ncontpts -= n;
  }
 
  PowParCurv_compare1(egeom, pgeom, maxerr, nc);

  return(egeom);
}


/*****************************************************************************
*                               monotobern1()
******************************************************************************
* 
* 1   Purpose
*     This function converts a monomial to the the Bernstein form.
*     This function is called once for each i = 0,1,...,m-1, for each
*     spatial coordinate. 
* 
* 2   Specification
*     #include "iges.h"
*     double monotobern1(int m, int i, double *c)
* 
* 3   Description
*     This routine converts monomial basis representation of an algebraic 
*     curve the Bernstein basis.
* 
* 4   References
*     Not applicable
* 
* 5   Parameters
*       1. int m
*         On entry: the maximum degree of the input monomial
*       2. int i
*         On entry: degree for which the Bernstein coefficient is required.
*       3. double * c
*         On entry: array of length m+1 holding the coefficients
* 	   of the monomial
* 
* 6   Return Values, Error Indicators and Warnings
*     Return the Bernstein coefficient 
*     the input monomial polynomial.
* 
* 7   Accuracy
*     Not applicable.
* 
* 8   Further Comments
*     Not applicable.
* 
* 9   Functions referenced by monotobern1() are:
*     Not applicable.
* 
* 10  Functions that reference monotobern1() are:
*     PowCurv_to_ParCurv1()
* 
******************************************************************************/

double monotobern1 (int m, int i, double *c)
{
  double
    C,
    f,
    fi,
    fm;
  int
    k;

  fm = factorial(m);    
  fi = factorial(i);    

  C = 0;
  for (k=0; k<=i; k++) {
    f = factorial(m-k) / factorial(i-k);
    C += (f*c[k]);
  }
     
  C *= (fi/fm);
     
  return C;
}

/*****************************************************************************
*                               RemoveCurveKnot()
******************************************************************************
* 
* 1   Purpose
*     This function removes multiple copies of a single knot value from
*     a B-spline curve knot vector.
* 
* 2   Specification
*     #include "iges.h"
*     double RemoveCurveKnot(int n, int p, double *U, vector **pw,
*                            double u, int r, int s, int num);
* 
* 3   Description
*     This function removes multiple copies of a single knot value from
*     a B-spline curve knot vector.
* 
* 4   References
*     [1] W. Tiller, "Knot-Removal Algorithms for NURBS Curves and Surfaces,"
*         Computer Aided Design, 24(8):445-453, August 1992.
* 
* 5   Parameters
*       1. int n
*         On entry: the largest index of control points.
*       2. int p
*         On entry: the degree (order-1) of the B-spline curve.
*       3. double * U
*         On entry: the initial knot vector.
*         On exit: the knot vector with knots removed.
*       4. vector ** pw
*         On entry: the initial control polygon.
*         On exit: the modified control polygon.
*       5. double u
*         On entry: the knot value to be remove.
*       6. int r
*         On entry: the index of the knot u.
*       7. int s
*         On entry: the multiplicity of knot u.
*       8. int num
*         On entry: the number of instances of u that should be removed.
* 
* 6   Return Values, Error Indicators and Warnings
*     The actual number of instances of u removed is returned. This
*     may be less than the number requested depending upon the 
*     continuity conditions of the curve.
* 
* 7   Accuracy
*     Not applicable.
* 
* 8   Further Comments
*     Not applicable.
* 
* 9   Functions referenced by RemoveCurveKnot() are:
*     addv()
*     copyvector()
*     distance()
*     multv()
*     subv()
*     vect_array1()
*     vectalloc()
* 
* 10  Functions that reference RemoveCurveKnot() are:
*     PowCurv_to_ParCurv1()
* 
******************************************************************************/

int RemoveCurveKnot(int n, int p, double *U, vector **pw, double u, int r,
		    int s, int num)
/**********************************************************************
  This function removes knot u (index r) num times from a B-spline curve
   Input 
       n: the largest index of control points;
       p: the degree of the B spline (order-1);
       U: the knot vector (it changes after this function called);
       pw: homogeneous coordinates of the control points;
       u: the knot which needs to be removed;
       r: the index of u;
       s: the multiplicity of the knot (u);
     num: the times that u is to supposed to be removed;

   Output
       t: returns the actual times u is removed;
       new knots & control points are placed in U & pw
****************************************************************************/

{
  int i,j,k,m,t,ord,fout,last,first,off,ii,jj,remflag, k1;
  double alfi, alfj;
  vector *tmp;
  double Tol=1.0e-9;
  vector **temp;

  m=n+p+1;
  ord=p+1;
  fout=(2*r-s-p)/2;                           /* first control point out  */
  last=r-s;
  first=r-p;

  tmp = vectalloc();
  temp = vec_array1((unsigned)(2*p+1));

  for (t=0; t<num; t++) {
    off = first-1;
    copyvector(pw[off],temp[0]);
    copyvector(pw[last+1],temp[last+1-off]);
    i=first; j=last;
    ii=1;    jj=last-off;
    remflag=0;
    while(j-i>t) {
      alfi=(u-U[i])/(U[i+ord+t]-U[i]);
      alfj=(u-U[j-t])/(U[j+ord]-U[j-t]);
      
      temp[ii] = multv(1.0/alfi,subv(pw[i],multv(1.0-alfi,temp[ii-1])));
      temp[jj] = multv(1.0/(1.0-alfj),subv(pw[j],multv(alfj,temp[jj+1])));

      i=i+1; ii=ii+1;
      j=j-1; jj=jj-1;
    }                                            /* End of while-loop */
    
    if (j-i<t)                                   /* Check if knot removable  */
      {
	if (distance(temp[ii-1],temp[jj+1]) <= Tol) remflag=1;
      }
    
    else {
      alfi=(u-U[i])/(U[i+ord+t]-U[i]);
      tmp = addv(multv(alfi,temp[ii+t+1]),multv(1.0-alfi,temp[ii-1]));
      if (distance(pw[i],tmp)<=Tol)
	remflag=1;
    }
    
    if (remflag==0)                        /* can't remove any more knots*/
      {
	break;
      }
    else                  /* successful removal. save new control points */
      {
	i=first; j=last;
	while(j-i>t) {
	  copyvector(temp[i-off],pw[i]);
	  copyvector(temp[j-off],pw[j]);
	  i=i+1; j=j-1;
	}
      }
    first=first-1; last=last+1;
  }                                                    /*End of for-loop*/

  if(t==0) return t;
  
  for(k=r+1;k<=m;k++) U[k-t]=U[k];                         /*shift knots*/
  
  j=fout; i=j;                         /* pj thru pi will be overwritten*/
  for(k=1;k<t;k++)
    if (k%2==1)                                             /*k modulo 2*/
      i=i+1; else j=j-1;
  for(k=i+1;k<=n;k++)                                            /*shift*/
    {
      copyvector(pw[k],pw[j]);
      j=j+1;
    }
  return t;
}

/* return the sum of two vectors */

vector *addv(vector *a, vector *b) {
  vector *c;
  c = vectalloc();
  c->x = a->x+b->x;
  c->y = a->y+b->y;
  c->z = a->z+b->z;
  c->w = a->w+b->w;
  return (c);
}

/* return the difference between two vectors */
 
vector *subv(vector *a, vector *b) {
  vector *c;
  c = vectalloc();
  c->x = a->x-b->x;
  c->y = a->y-b->y;
  c->z = a->z-b->z;
  c->w = a->w-b->w;
  return (c);
}

/* multiply a vector by a constant value */

vector *multv(double a, vector *b) {
  vector *c;
  c = vectalloc();
  c->x = a*b->x;
  c->y = a*b->y;
  c->z = a*b->z;
  c->w = a*b->w;
  return (c);
}

/********* PowParCurv_compare() *********
* 1     Purpose
* 
*       Position error between parametric spline curve and NURBS curve.
* 
* 2     Specification
* 
*       #include "iges.h"
*       void PowParCurv_compare(ParCurv *egeom, PowCurv *pgeom, double *maxerr,
*                               int nc);
* 
* 3     Description
* 
*       This function calculates the position error between a parametric spline
*       curve and a NURBS curve.
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
*           1.  ParCurv * egeom
*               On entry:  the address of a structure containing the NURBS curve.
* 
* 
*           2.  PowCurv * pgeom
*               On entry:  the address of a structure containing the parametric
* 	      spline curve.
* 
*           3.  double * maxerr
*               On exit:  the addrees of a variable that will contain the global
* 	      conversion error calculated as the maximum distance between
* 	      points on the parametric spline and the NURBS curve evaluated
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
* Functions referenced by PowParCurv_compare() are:
*  distance()
*  PowCurv_eval()
*  rbspeval()
*  vectfree()
*
* Functions that reference PowParCurv_compare() are:
*  ParCurv_to_PowCurv()
*  PowCurv_to_ParCurv()
*/

void PowParCurv_compare(ParCurv *egeom, PowCurv *pgeom, double *maxerr,
			int nc)
{
  double2 u1, u2, u3, u4, dist1, dist2, max;
  struct vector *v1, *v2, *v3, *v4;
  int i, j;
  
  *maxerr = 0.0;		/* init to zero for MAX */
  
  for (i=0; i<pgeom->nsegmts; i++) {
    j = i*nc;
    u1 = egeom->knots[j+3] + 0.5*(egeom->knots[j+4] - egeom->knots[j+3]);
       /* check far end of segment */
    u2 = egeom->knots[j+4] - ZERO;
    u3 = pgeom->knots[i] + 0.5*(pgeom->knots[i+1] - pgeom->knots[i]);
    u4 = pgeom->knots[i+1] - ZERO;
    v1 = rbspeval(egeom, u1, 0);
    v2 = rbspeval(egeom, u2, 0);
    v3 = PowCurv_eval(pgeom, u3, 0);
    v4 = PowCurv_eval(pgeom, u4, 0);
    dist1 = distance(v1, v3);
    dist2 = distance(v2, v4);
    max = MAX(dist1, dist2);
    *maxerr = MAX(max, *maxerr);
    vectfree(v1);
    vectfree(v2);
    vectfree(v3);
    vectfree(v4);
  }
}

/********* PowParCurv_compare1() *********
* 1     Purpose
* 
*       Position error between parametric spline curve and NURBS curve.
*       All internal knots of the NURBS curve must have multiplicity of 1.
* 
* 2     Specification
* 
*       #include "iges.h"
*       void PowParCurv_compare1(ParCurv *egeom, PowCurv *pgeom,
*                                double *maxerr, int nc);
* 
* 3     Description
* 
*       This function calculates the position error between a parametric spline
*       curve and a NURBS curve.
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
*           1.  ParCurv * egeom
*               On entry:  the address of a structure containing the NURBS curve.
* 
* 
*           2.  PowCurv * pgeom
*               On entry:  the address of a structure containing the parametric
* 	      spline curve.
* 
*           3.  double * maxerr
*               On exit:  the address of a variable that will contain the
* 	      global conversion error calculated as the maximum distance
* 	      between points on the parametric spline and the NURBS curve
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
* Functions referenced by PowParCurv_compare1() are:
*  distance()
*  PowCurv_eval()
*  rbspeval()
*  vectfree()
*
* Functions that reference PowParCurv_compare1() are:
*  ParCurv_to_PowCurv()
*  PowCurv_to_ParCurv()
*/

void PowParCurv_compare1(ParCurv *egeom, PowCurv *pgeom, double *maxerr,
			 int nc)
{
  double2 u1, u2, u3, u4, dist1, dist2, max;
  struct vector *v1, *v2, *v3, *v4;
  int i, j;
  
  *maxerr = 0.0;		/* init to zero for MAX */
  
  for (i=0; i<pgeom->nsegmts; i++) {
    j = egeom->order;
    /*    u1 = egeom->knots[j-1] + 0.5*(egeom->knots[j] - egeom->knots[j-1]);  */
       /* check far end of segment */
    /*    u2 = egeom->knots[j] - ZERO;    */

    u1 = egeom->knots[j-1+i] + 0.5*(egeom->knots[j+i] - egeom->knots[j-1+i]);
       /* check far end of segment */
    u2 = egeom->knots[j+i] - ZERO;

    u3 = pgeom->knots[i] + 0.5*(pgeom->knots[i+1] - pgeom->knots[i]);
    u4 = pgeom->knots[i+1] - ZERO;
    
    /*
    v1 = rbspeval(egeom, u1, 0);
    v2 = rbspeval(egeom, u2, 0);
    */

    v1 = rbspeval(egeom, u3, 0);
    v2 = rbspeval(egeom, u4, 0);
    v3 = PowCurv_eval(pgeom, u3, 0);
    v4 = PowCurv_eval(pgeom, u4, 0);
    dist1 = distance(v1, v3);
    dist2 = distance(v2, v4);
    max = MAX(dist1, dist2);
    *maxerr = MAX(max, *maxerr);
    vectfree(v1);
    vectfree(v2);
    vectfree(v3);
    vectfree(v4);
  }
}

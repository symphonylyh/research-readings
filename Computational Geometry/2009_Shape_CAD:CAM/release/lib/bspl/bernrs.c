/* ***************************************************************************
 Copyright (C) 1996 Massachusetts Institute of Technology all rights reserved 
	Programmer: George A. Kriezis
**************************************************************************** */
#include <stdio.h>
#include <math.h>
#include "gen.h"
#include "bspl.h"

/* extern FILE *db; */
static double taorig, tborig;  /* global variables, will be explained later. */

/*****************************************************************************
*                                  bernrs()
******************************************************************************
* 
* 1   Purpose
*     This function is a univariate polynomial root solver for coefficients 
*     expressed in the Bernstein basis.
* 
* 2   Specification
*     #include "bspl.h"
*     void bernrs(double *bcoeff, int order, double ta, double tb, double eps, 
*                 doble *roots, int *nroots)
* 
* 3   Description
*     This routine calculates the roots by first isolating intervals containing 
*     a single root or reducing intervals to a very small value that the root 
*     may be considered to lie in the middle of the interval. The roots are 
*     isolated by recursive binary subdivision of the interval containing the 
*     roots. Roots in isolated intervals are then found by recursive subdivision.
* 
* 4   References
*     The ideas used have been published by Lane and Riesenfeld (1981).
*     [1]J. M. Lane and R. F. Riesenfeld. Bounds on a Polynomial, 
*        BIT: Nordisk Tidskrift for Informations-Behandling, 21(1):112-117, 1981.
* 
* 5   Parameters
*       1.double * bcoeff
*         On entry: one dimensional array with bcoeff[0] containing the 
* 	          coefficient of the first Bernstein polynomial
*       2.int order
*         On entry: the degree of the polynomial curve in the Bernstein basis 
* 	          plus one
*       3.double ta
*         On entry: parameter value of the left end of the interval containing 
* 	          the roots
*       4.double tb
*         On entry: parameter value of the right end of the interval containing 
* 	          the roots
*       5.double eps
*         On entry: accuracy to which the roots are desired
*       6.double * roots
*         On exit: one dimensional array containing the roots found in the 
* 	         interval
*       7.int * nroots
*         On exit: the number of roots found in the interval
* 
* 6   Return Values, Error Indicators and Warnings
*     The presence of multiple roots leads to extreme subdivision.  The point 
*     at which this is stopped depends on the degree of the polynomial. The 
*     present stopping criteria used may fail at some high degree cases
* 
* 7   Accuracy
*     After the roots have been isolated, other root finding techniques such 
*     as Newton Raphson may be used for faster convergence.
* 
* 8   Further Comments
*     See also bernrs1 (bspl _ bernrs.c).
* 
* 9   Functions referenced by bernrs() are:
*     bernrs1()
*     lf2comp()
* 
* 10  Functions that reference bernrs() are: None
* 
******************************************************************************/

void bernrs(double bcoeff[], int order, double ta, double tb, double eps, 
       double roots[], int *nroots)
{
  int i,iedge;

  /* check if the polynomial is 0. */
  /* if so, there're infinite number of roots? (glshen) */
  iedge=0;
  i=0;
  while (iedge==0 && i<order){
        if (lf2comp(bcoeff[i],0.0,eps) == 0)
           i++;
        else
           iedge=1;
        };

  /* if not, call bernrs1() to find the roots. */
  if (iedge==1){
      bernrs1(bcoeff,order,eps,roots,nroots);
      for (i=0; i<*nroots; i++)
	  /* map the roots from [0,1] back to [ta,tb]. */
          roots[i] = ta + (tb-ta)*roots[i];
      }
  /* if so, give tho two roots ta and tb. */
  else{
      *nroots=2;
      roots[0]=ta;
      roots[1]=tb;
      }
}

/*****************************************************************************
*                                 bernrs1()
******************************************************************************
* 1   Purpose
*     This function is a univariate polynomial root solver for coefficients 
*     expressed in the Bernstein basis.
* 
* 2   Specification
*     #include "bspl.h"
*     void bernrs1(double bcoeff[], int order, double eps, double roots[],
* 	         int *nroots)
* 
* 3   Description
*     This routine calculates the roots by first isolating intervals 
*     containing a single root or reducing intervals to a very small value 
*     that the root may be considered to lie in the middle of the interval.
*     The roots are isolated by recursive binary subdivision of the interval 
*     containing the roots. Roots in isolated intervals are then found by
*     recursive subdivision.  In this function, different from bernrs(), the
*     lower and upper bounds of the interval considered is [0.0 1.0].
* 
* 4   References
*     The ideas used have been published by Lane and Riesenfeld (1981).
*     [1]J. M. Lane and R. F. Riesenfeld. Bounds on a Polynomial, 
*        BIT: Nordisk Tidskrift for Informations-Behandling, 21(1):112-117,
*        1981.
* 
* 5   Parameters
*       1.double * bcoeff
*         On entry: one dimensional array with bcoeff[0] containing the 
*                   coefficient of the first Bernstein polynomial
*       2.int order
*         On entry: the degree of the polynomial curve in the Bernstein basis 
*                   plus one
*       3.double eps
*         On entry: accuracy to which the roots are desired
*       4.double roots[]
*         On exit: one dimensional array containing the roots found in the 
*                  interval
*       7.int * nroots
*         On exit: the number of roots found in the interval
* 
* 6   Return Values, Error Indicators and Warnings
*     The presence of multiple roots leads to extreme subdivision.  The point 
*     at which this is stopped depends on the degree of the polynomial. The 
*     present stopping criteria used may fail at some high degree cases
* 
* 7   Accuracy
*     After the roots have been isolated, other root finding techniques such 
*     as Newton Raphson may be used for faster convergence.
* 
* 8   Further Comments
*     See also bernrs (bspl _ bernrs.c).
* 
* 9   Functions referenced by bernrs1() are:
*     dbl_array1()
*     dbl_array2()
*     egeomalloc1()
*     findroot()
*     free_darray1()
*     free_darray2()
*     free_egeom()
*     rootisolate()
*     splitpoly()
* 
* 10  Functions that reference bernrs1() are:
*     bernrs()
* 
****************************************************************************/

void bernrs1(double bcoeff[], int order, double eps, double roots[],
	     int *nroots)
/* bernstein root solver assuming bcoeff defined in 0. to 1. knot range  */
{
  int i,ninterval;    /* ninterval -- the number of intervals each of which
		       *              contains only one root. */
  double u,max;       /* u   -- working variable, to be used as the parametric
		       *        value of the splitting point. 
		       * max -- no longer used. (glshen)*/
  double **interval, *rinterval;   /* interval -- 2d array containing lower and
				    *            upper bounds of the intervals.
				    * rinterval -- no longer used. (glshen) */
  double *c1,*c2,ta,tb;         /* c1 -- the array containing the coefficients
				 *       of the lower part after splitting.
				 * c2 -- the array containing the coefficients
				 *       of the upper part after splitting.
				 * ta -- the lower bound of the interval
				 * tb -- the upper bound of the interval   */
  struct vector *v;      /* no longer used (glshen) */
  struct egeom *egeom;   /* no longer used (glshen) */

  /* allocate memory */
  interval = dbl_array2((unsigned)order,3);
  rinterval = dbl_array1((unsigned)order);
  c1 = dbl_array1((unsigned)order);
  c2 = dbl_array1((unsigned)order);

/*
 * max = absmaxarray(order, 1, bcoeff, 1);
 *
 * fprintf(db, "\nThe Bernstein coefficients");
 *
 * for(i=0; i<order; i++)
 *   bcoeff[i] = bcoeff[i]/max;
 */

  /* initialize the number of intervals where exists one and only one root
   * each. */
  ninterval = 0;
  /* obtain the intervals containing roots by calling rootisolate(). */
  rootisolate(bcoeff, bcoeff, order,0.0,1.0, interval, &ninterval);

  /* set the number of roots as 0 initially. */
  *nroots = 0;

  /* allocate memory for Bezier form of the polynomial. (glshen) */
  egeom = egeomalloc1(order,order);

/* this is no longer used, so, what's it for? (glshen) */
  /* set the values of the polynomial as Bezier curve. */
  for (i=0; i<order; i++)
      {
       /* knot vector */
       egeom->knots[i] = 0.0;
       egeom->knots[i+order] = 1.0;
       /* control points */
       egeom->contpts[i]->x = bcoeff[i];
       egeom->contpts[i]->y = 0.0;
       egeom->contpts[i]->z = 0.0;
       egeom->contpts[i]->w = 1.0;
      }

  /* get the root in each interval. */
  for (i=0; i< ninterval; i++)
      {
       /* get the lower and upper bounds of the i-th interval */
       ta = interval[i][0];
       tb = interval[i][1];

/*    
 *   root[i] = 0.5*(ta+tb);
 *
 *   v = evalbsp(&egeom, root[i]);
 */

       /* split the curve to get the part in [ta,tb] */    
       u = (tb-ta)/(1.0-ta);    /* parametric value of tb after mapping [ta,1]
				   to [0,1]. */
       /* split at ta first. */
       splitpoly(bcoeff,order,ta,c1,c2);
       /* split the upper part [ta,1] at tb. */
       splitpoly(c2,order,u,c1,c2);  /* now c1 is for the part [ta,tb]. */

       /* these are global variable declared at the beginning of the file,
        * which remember the original bounds of an interval containing only
	* one root, and will be used at findroot(). */
       taorig = ta;
       tborig = tb;

/* what is the meaning of interval[i][2]. it seems here it makes no
 * differences. (glshen) */
       /* find the root in the current interval. */
       if (interval[i][2] == 1.0)
          roots[i] = findroot(c1, c1, order, ta, tb);
       else if (interval[i][2] == 2.0)
               {
        	roots[i] = findroot(c1, c1, order, ta, tb);
               }
       *nroots = *nroots + 1;
      }

  /* free memory */
  free_egeom(egeom);
  free_darray2(interval);
  free_darray1(rinterval);
  free_darray1(c1);
  free_darray1(c2);
}

/*****************************************************************************
*                                   rootisolate()
******************************************************************************
* 
* 1   Purpose
*     This function is a subroutine of the routine bernrs(), which finds the
*     roots of a polynomial in Bernstein basis.
* 
* 2   Specification
*     #include "bspl.h"
*     int rootisolate(double c[], double cpure[], int order, double ua,
* 	            double ub, double **interval, int *ninterval)
* 
* 3   Description
*     Given a polynomial in Bernstein basis, this function computes intervals
*     such that there exists one and only one root in each interval.
* 
* 4   References
*     The ideas used have been published by Lane and Riesenfeld (1981).
*     [1]J. M. Lane and R. F. Riesenfeld. Bounds on a Polynomial, 
*        BIT: Nordisk Tidskrift for Informations-Behandling, 21(1):112-117, 1981.
* 
* 5   Parameters
*       1.double c[]
*         On entry: the coefficients of a part of the given polynomial in
*                 Bernstein basis. the initial call of this routine uses the
*                 whole polynomial.
*       2.double cpure[]
*         On entry: the coefficients of the given polynomial in Bernstein basis
*       3.int order
*         On entry: the order of the polynomial
*       4.double ua
*         On entry: the lower bound of the interval to be considered.
*       5.double ub
*         On entry: the upper bound of the interval to be considered.
*       6.double **interval
*         On entry: a 2D array
*         On exit:  the matrix containing information such as the lower and upper
* 	          bounds of the intervals containing one root.
*       7.int *ninterval
*         On entry: the address of an integer
*         On exit:  the number of intervals which contain roots.
* 
* 6   Return Values, Error Indicators and Warnings
*     Not applicable  
* 
* 7   Accuracy
*     Not applicable
* 
* 8   Further Comments
*     Not applicable
* 
* 9   Functions referenced by rootisolate() are:
*     abssumarray()
*     dbl_array1()
*     free_darray1()
*     signchange()
*     splitpoly2()
*     sum()
* 
* 10  Functions that reference rootisolate() are: 
*     bernrs1()
*     rootisolate()
* 
*****************************************************************************/

int rootisolate(double c[], double cpure[], int order, double ua,
	    double ub, double **interval, int *ninterval)
{
  double um;         /* to be used as the midpoint of an interval */
  double *c1,*c2;    /* containing the coefficients of the lower and upper
		      * parts everytime sudivision is done. */
  double rmin, Cm;   /* define a small number as 0 for comparison, which is 
		      * depends on the degree of the polynomial and will 
		      * determine the accuracy of the roots. */
  int nsign;         /* number of sign changes */

  /* allocate memory */
  c1=dbl_array1((unsigned)order);
  c2=dbl_array1((unsigned)order);

  /* define the samll number. */
  Cm = 1.0;
  rmin = 2*(order-1)*Cm*PRECISION;

  /* find how many times the coefficients change sign. */
  nsign = signchange(c,order);

  /* if the coefficients do not change sign, no root exists. */
  if (nsign == 0 )
     {
      free_darray1(c1);
      free_darray1(c2);
      return(0);
     }
  /* otherwise, there exist roots. */
  else if (nsign == 1)   /* if only changes sign once, only one root. */
          {
           /* set the interval for the root. */
           interval[*ninterval][0] = ua;   /* the lower bound */
           interval[*ninterval][1] = ub;   /* the upper bound */
           interval[*ninterval][2] = 1.0;  /* what is this ? (glshen) */
           *ninterval = *ninterval + 1;    /* the number of roots */
	   /* free memory */
           free_darray1(c1);
           free_darray1(c2);
           return(1);
          }
       else if (nsign >= 2)   /* if there are more than 1 root. */
               {
		/* if the curve is extremely close to 0-axis, (glshen) */
                if (abssumarray(c,order)/order < rmin)
                   {
                    interval[*ninterval][0] = ua;
                    interval[*ninterval][1] = ub;
	            interval[*ninterval][2] = 2.0;
	            *ninterval = *ninterval + 1;
	            free_darray1(c1);
	            free_darray1(c2);
	            return(2);
                   }

		/* subdivide at the midpoint. */
                um = 0.5*sum(ua,ub);
                splitpoly2(cpure, order, ua, um, ub, c1, c2);
		/* recursively call rootisolate() with the two parts. */
                rootisolate(c1, cpure, order, ua, um, interval, ninterval);
                rootisolate(c2, cpure, order, um, ub, interval, ninterval);
                /* free memory */
                free_darray1(c1);
                free_darray1(c2);
               }
}

/*************************************************************************
*                               abssumarray()
**************************************************************************
* 
* 1   Purpose
*     Calculate the absolute value sum of an array.
* 
* 2   Specification
*     #include "bspl.h"
*     double abssumarray(double *array, int n);
* 
* 3   Description
*     This routine calculates the absolute value sum of an array, 
*     sum(i=0..n-1)array[i].
* 
* 4   References
*     Not applicable.
* 
* 5   Parameters
*        1.double * array
*          On entry: the address of the array.
*        2.int n
*          On entry: the number of elements in the array.
* 
* 6   Return Values, Error Indicators and Warnings
*     The absolute value sum of the array, sum(i=0..n-1)array[i], is returned.
* 
* 7   Accuracy
*     Not applicable.
* 
* 8   Further Comments
*     Not applicable.
* 
* 9   Functions referenced by abssumarray() are: None
* 
* 10  Functions that reference abssumarray() are:
*     rootisolate()
* 
**************************************************************************/

double abssumarray(double array[],int size)
{
  int i;
  double sm;

  /* initialization */
  sm = 0.0;

  /* summation */
  for (i=0; i<size; i++)
      sm += FABS(array[i]);

  return(sm);
}

/***************************************************************************
*                             splitpoly2()
****************************************************************************
* 
* 1   Purpose
*     This function is a subroutine of the routine which finds the roots of
*     a polynomial in Berstein basis. It splits the polynomial at given
*     positions. 
* 
* 2   Specification
*     #include "bspl.h"
*     void splitpoly2(double c[], int order, double ua, 
* 	            double um, double ub, double c1[], double c2[])
* 
* 3   Description
*     This function calculates the coefficients of two polynomials in Bernstein
*     basis which are obtained by splitting a polynomial also in Bernstein
*     basis at certain parametric values. It first maps the parametric to
*     interval [0,1] and then call splitpoly() which uses De Castgeljau
*     algorithm.
*  
* 4   References
*     Not applicable
* 
* 5   Parameters
*     1.double c[]
*       On entry: the coefficients of the polynomial to be split.
*     2.int order
*       On entry: the order of the polynomial
*     3.double
*       On entry: the lower bound of the interval
*     4.double um
*       On entry: the parametric value where the polynomial is to be split.
*     5.double ub
*       On entry: the upper bound of the interval
*     6.double c1[]
*       On entry: an array of length = order of the polynomial
*       On exit:  the coefficients of the new polynomial which is the lower
*                 part of the old one.
*     7.double c2[]
*       On entry: an array of length = order of the polynomial
*       On exit:  the coefficients of the new polynomial which is the upper
*                 part of the old one.
* 
* 6   Return Values, Error Indicators and Warnings
*     Not applicable
* 
* 7   Accuracy
*     Not applicable
* 
* 8   Further Comments
*     Not applicable
* 
* 9   Functions referenced by splitpoly2() are:
*     dbl_array1()
*     free_darray1()
*     splitpoly()
* 
* 10  Functions that reference splitpoly2() are: 
*     rootisolate()
* 
****************************************************************************/

void splitpoly2(double c[], int order, double ua, 
	   double um, double ub, double c1[], double c2[])
{
  double *ctemp;  /* temporarily used for containing the coefficients of a
		     split part. */
  double u;       /* to be used as parametric value. */

  /* allocate memory */
  ctemp = dbl_array1((unsigned)order);

  /* if the interval is [0,1] and the splitting point is 0.5. */ 
  if (ua == 0.0 && um == 0.5 && ub == 1.0)
     {
      splitpoly(c, order, 0.5, c1, c2);
     }
  /* otherwise, first subdivide at ua and then subdivide the upper part at
   * um to get the coefficients for the part from ua to um.
   * then, subdivide at ub and then subdivided at um to get the coefficients
   * for the part from um to ub. */
  else
     {
      u = (um-ua)/(1.0-ua);          /* the splitting point after map [ua, 1]
				      * to [0,1]. */ 
      splitpoly(c,order,ua,c1,c2);   /* split at ua to get the part from ua to
				      * 1. */
      splitpoly(c2,order,u,c1,c2);   /* split at um to get the part from ua
				      * to um */

      u = um/ub;   /* splitting point after map [0, ub] to [0,1]. */    
      splitpoly(c,order,ub,ctemp,c2);     /* split at ub */
      splitpoly(ctemp,order,u,ctemp,c2);  /* split at um. */
    
  }
  
  /* free memory */
  free_darray1(ctemp);
}

/***************************************************************************
*                              splitpoly()
****************************************************************************
* 
* 1   Purpose
*     This function is a subroutine of the routine which finds the roots of
*     a polynomial in Berstein basis. It splits a Bezier curve at a given
*     position.
* 
* 2   Specification
*     #include "bspl.h"
*     void splitpoly(double c[], int order, double u, double c1[], double c2[])
* 
* 3   Description
*     This function calculates the control polygons of two Bezier curves which
*     are obtained by splitting a Bezier curve at certain parametric value. It 
*     uses De Castgeljau algorithm.
* 
* 4   References
*     Not applicable
* 
* 5   Parameters
*     1.double c[]
*       On entry: the control points of the Bezier curve to be split.
*     2.int order
*       On entry: the order of the Bezier curve
*     3.double u
*       On entry: the parametric value at which the curve is to be subdivided.
*     4.double c1[]
*       On entry: the array of certain length
*       On exit:  the array containing the control points of the new Bezier
*                 curve which is the lower part of the old one.
*     5.double c2[]
*       On entry: the array of certain length
*       On exit:  the array containing the control points of the new Bezier
*                 curve which is the upper part of the old one.
* 
* 6   Return Values, Error Indicators and Warnings
*     Not applicable
* 
* 7   Accuracy
*     Not applicable
* 
* 8   Further Comments
*     Not appicable
* 
* 9   Functions referenced by splitpoly() are:
*     dbl_array2()
*     free_darray2()
*     sum()
* 
* 10  Functions that reference splitpoly() are:
*     bernrs1()
*     findroot()
*     splitpoly2()
* 
****************************************************************************/

void splitpoly(double c[], int order, double u, double c1[], double c2[])
{
  double u1,u2;    /* two linear interpolation factors */
  double **cm;     /* the matrix of control points. each column corresponds
		      to the newly added control points of the j-th evaluation 
		      in De Casteljau algorithm. */
  int i,j,k,km1;   /* indices */

  /* allocate memory */
  cm = dbl_array2((unsigned)order,(unsigned)order);
  
  /* evaluate u1 and u2 */
  u1 = 1-u;
  u2 = u;

  /* initialize the 0-th column of the matrix of control points, which actually
   * is the control points before subdivision, i.e., c[]. */
  for (j=0; j<order; j++)
      cm[j][0] = c[j];

  /* if the parametric value at which the subdivision is to be proceded, is
   * 0.5,  the new control points are the midpoints of edges of the old
   * control polygon. */
  if (u == 0.5)
     {
      for (k=1; k<order; k++)
          {
	   km1 = k-1;
	   for (j=k; j<order; j++)
	       cm[j][k] = 0.5*sum(cm[j-1][km1],cm[j][km1]);
          }
     }
  else 
  /* otherwise, each control point is linear interpolations of the two vertices
   * on each edge of the old control polygon. */ 
     {
      for (k=1; k<order; k++)
          {
	   km1 = k-1;
	   for (j=k; j<order; j++)
	       cm[j][k] = sum(u1*cm[j-1][km1],u2*cm[j][km1]);
          }
     }

  /* assign the control points to two arrays. */
  for (i=0; i<order; i++)
      {
       c1[i] = cm[i][i];
       c2[i] = cm[order-1][order-1-i];
      }

  /* free memory */
  free_darray2(cm);
}

/************************************************************************
*                                  sum()
*************************************************************************
* 
* 1   Purpose
*     Add to double precision values.
* 
* 2   Specification
*     #include "bspl.h"
*     double sum(double a, double b);
* 
* 3   Description
*     This function adds two double precision values.
* 
* 4   References
*     Not applicable.
* 
* 5   Parameters
*       1.double a
*         On entry: the first value.
*       2.double b
*         On entry: the second value.
* 
* 6   Return Values, Error Indicators and Warnings
*     The sum of the two values a + b is returned.
* 
* 7   Accuracy
*     Not applicable.
* 
* 8   Further Comments
*     Not applicable.
* 
* 9   Functions referenced by sum() are:
* 
* 10  Functions that reference sum() are:
*     rootisolate()
*     splitpoly()
* 
************************************************************************/

double sum(double a, double b)
{
  return(a+b);
}

/*************************************************************************
*                                 signchange()
**************************************************************************
* 
* 1   Purpose
*     This function is a subroutine of the routine which finds the roots
*     of a polynomial in Berstein basis. It checks how many times the 
*     elements of an array change sign.
* 
* 2   Specification
*     #include "bspl.h"
*     int signchange(double a[], int n)
* 
* 3   Description
*     This function counts the number of times that the elements of an array 
*     of length n changes sign from the beginning to the end. 
* 
* 4   References
*     Not applicable
* 
* 5   Parameters
*     1.double a[]
*       On entry: an array of length n.
*     2.int n
*       On entry: the length of the array
* 
* 6   Return Values, Error Indicators and Warnings
*     This function returns an integer which is the number of times that the
*     elements change sign. 
* 
* 7   Accuracy
*     Not applicable
* 
* 8   Further Comments
*     Not applicable
* 
* 9   Functions referenced by signchange() are: None
* 
* 10  Functions that reference signchange() are: 
*     findroot()
*     rootisolate()
* 
**************************************************************************/

int signchange(double a[], int n)
{
  int nsign,osign,sign,i;   /* nsign -- the number of sign change 
			     * osign -- the sign of the preceeding element
			     * sign  -- 1  positive
			     *          0  zero
			     *         -1 negative                    */
  double rmin,cm;           /* rmin is to be used as 0 for comparison. */

  cm = 1.0;
  nsign = 0;
  /* rmin defines an extremely small number which is to be used as 0. it 
   * depends on the order of the polynomial. numbers within [-rmin,rmin]
   * are considered as 0. */
  rmin = 2.0*((double)(n-1))*cm*PRECISION;

  /* starts from the first element in the array */
  i = 0;
  if (a[i] < -rmin)  /* if it is less than 0 */
     sign = -1;
  else if (a[i] > rmin)
          sign = 1;  /* if it is greater than 0 */
       else
          sign = 0;  /* if it is 0 */

  /* set the sign of the preceeding element as the sign of the first element */
  osign = sign;

  /* check through the whole array for sign changes. */
  for (i=1; i<n; i++)
      {
       if (a[i] < -rmin)  
          sign = -1;
       else if (a[i] > rmin)
               sign = 1;
            else
               sign = 0;
       /* if the sign of current element is different from that of its
	  preceeding element, increase nsign by 1. */
       if (sign != osign)
          nsign++;
       /* set the old sign as that of the current element and go to next
	* element. */
       osign = sign;
      }

  return(nsign);
}

/*****************************************************************************
*                                   findroot()
******************************************************************************
* 
* 1   Purpose
*     This function is a subtoutine of the routine which finds the roots of a 
*     polynomial expressed in Bernstein basis. It finally finds the root in an
*     given interval where exists only one root.
* 
* 2   Specification
*     #include "bspl.h"
*     double findroot(double c[], double cpure[], int order, double ua,
*                     double ub)
* 
* 3   Description
*     This function calculates a position to split the interval further down
*     until the curve in the smaller interval is very close to 0 or one of
*     its ends is 0, and then return the midpoint or one of the ends as root.
*     The algorithm used here is from the following reference.
* 
* 4   References
*     The ideas used have been published by Lane and Riesenfeld (1981).
*     [1]J. M. Lane and R. F. Riesenfeld. Bounds on a Polynomial, 
*        BIT: Nordisk Tidskrift for Informations-Behandling, 21(1):112-117, 1981.
* 
* 5   Parameters
*     1.double c[]
*       On entry: the array containing the coefficients of a subdivided part
*                 of the original polynomial.
*     2.double cpure[]
*       On entry: the array containing the coefficients of the original
*                 polynomial.
*     3.int order
*       On entry: the order of the polynomial
*     4.double ua
*       On entry: the lower bound of the considered interval
*     5.double ub
*       On entry: the upper bound of the considered interval
* 
* 6   Return Values, Error Indicators and Warnings
*     This function returns the root in the considered interval.
* 
* 7   Accuracy
*     Not applicable
* 
* 8   Further Comments
*     Not applicable
* 
* 9   Functions referenced by findroot() are:
*     dbl_array1()
*     findichange()
*     free_darray1()
*     signchange()
*     splitpoly()
* 
* 10  Functions that reference findroot() are:
*     bernrs1()
*     findroot()
* 
******************************************************************************/

double findroot(double c[], double cpure[], int order, double ua, double ub)
{
  double *c1,*c2;    /* array containing the coefficents of the lower and upper
		      * part everytime when subdivision is done. */
  double utemp,us,dc,um;  /* to be used as parametric values */
  double root;            /* root of the polynomial in the interval
			   * considered. */
  double cm,zero;         /* define a small number as 0. */
  int ic,icp1,flag,i;     /* indices and flag */

  /* allocate memory */
  c1 = dbl_array1((unsigned)order);
  c2 = dbl_array1((unsigned)order);

  /* define zero */
  cm = 1.0;
  zero = 1000.0*2*(order-1)*cm*PRECISION;

  /* initially, set flag as TRUE. */
  flag = TRUE;

  /* check if all the coefficients are zero. */
  for (i=0; i<order; i++)
      flag = flag && (FABS(c[i]) < zero);
  /* if so, return the midpoint as root. */
  if (flag == TRUE)
     {
/*    fprintf(db, "\n All points within round off error\n"); */
      free_darray1(c1);
      free_darray1(c2);
      return(0.5*(ua+ub));
     }

  /* if the first coefficient is zero, then the root is the lower bound of the 
   * interval. */
  if (FABS(c[0]) < zero)
     {
/*    fprintf(db, "\nReturnl %+.16le %+.16le %+.16le",ua,c[0],c[order-1]); */
      free_darray1(c1);
      free_darray1(c2);
      return(ua);
     }
  /* if the last coefficient is zero, then the root is the upper bound of the
   * interval. */
  else if(FABS(c[order-1]) < zero)
     {
/*    fprintf(db, "\nReturnr %+.16le %+.16le %+.16le",ub,c[0],c[order-1]); */
      free_darray1(c1);
      free_darray1(c2);
      return(ub);
     }
/*
 * else if(abssumarray(c,order)/order < zero)
 * {
 *   um = 0.5*(ua+ub);
 *
 *   fprintf(db, "\nReturna %+.16le %+.16le %+.16le",um,c[0],c[order-1]);
 *   free_darray1(c1);
 *   free_darray1(c2);
 *   return(um);
 * }
 */

  /* find the index of the element whose sign is different from that of the
   * next. */
  ic = findichange(c,order);
  icp1 = ic+1;
  dc = c[icp1] - c[ic];    /* the width of interval [c[ic],c[ic+1]]. */

  /* calculate the point to subdivide the interval, see references. */
  us = (ic*dc - c[ic])/((order-1)*dc);
  um = (1.0 - us)*ua + us*ub;
  utemp = (um-taorig)/(tborig - taorig);   /* here taorig and tborig are the
					    * lower and upper bounds of the
					    * interval under consideration.
					    * see bernrs1(). */
  /* subdivide the interval at utemp. */
  splitpoly(cpure, order, utemp, c1,c2);

  /* obtain the root in this interval by recursively calling findroot(). */
  if (signchange(c1,order) > 0)
     /* if the root is in the lower part, only deal with that part. */
     {
      root = findroot(c1,cpure,order,ua,um);
      free_darray1(c1);
      free_darray1(c2);
      return(root);
     }
  else if (signchange(c2,order) > 0)
          /* if the root is in the upper part, only deal with that part. */
          {
           root = findroot(c2,cpure,order,um,ub);
           free_darray1(c1);
           free_darray1(c2);
           return(root);
          }
       else
	  /* if no sign change in both parts, no root. */
          printf("There is something wrong()\n");

  /* free memory */
  free_darray1(c1);
  free_darray1(c2);
}
	     
/****************************************************************************
*                                    findichange()
*****************************************************************************
* 
* 1   Purpose
*     This function is a subroutine of the routine which finds the roots of a 
*     polynomial expressed in Bernstein basis. It finds where the elements of an
*     array change sign.
* 
* 2   Specification
*     #include "bspl.h"
*     int findichange(double c[], int order)
* 
* 3   Description
*     This function finds the index in an array where the sign change of the 
*     elements will occur when going from this element to the next.
* 
* 4   References
*     Not applicable
* 
* 5   Parameters
*     1.double c[]
*       On entry: the array to be checked for sign change.
*     2.int order
*       On entry: the length of the array
* 
* 6   Return Values, Error Indicators and Warnings
*     This function returns an integer which is the index of the element where
*     sign change is to occur. 
* 
* 7   Accuracy
*     When checking sign change, an extremely small number is used instead of 0.
* 
* 8   Further Comments
*     Not applicable
* 
* 9   Functions referenced by findichange() are: None
* 
* 10  Functions that reference findichange() are:
*     findroot()
*   
****************************************************************************/

int findichange(double c[], int order)
{
  int sign,osign,i;    /* sign -- sign of the current element
			* osign -- sign of the preceeding element */
  double zero;         /* defines a small number as zero */
  zero = 2*(order-1)*1.0*PRECISION;

  osign = 0;   /* initialization */

  for (i=0; i<order; i++)
      {
       if (c[i] < -zero)     /* if the current element is negative. */
          sign = -1;
       else if(c[i] > zero)  /* if the current element is positive. */
          sign = 1;
       /* if the change of sign occurs, return the index of the preceeding
	* element. */
       if (sign*osign == -1)
          return(i-1);
       /* reset the sign of the preceeding element and goto next element. */
       osign = sign;
      }
}

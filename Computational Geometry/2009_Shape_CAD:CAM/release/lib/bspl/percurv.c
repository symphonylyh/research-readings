/************************************************************************
 *									*
			Copyright (C) 1996 by
	Massachusetts Institute of Technology, Cambridge, MA
			 All rights reserved

 *									*
 ************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <malloc.h>
#include "bspl.h"

/***************************************************************************
*                                    nbasisd_per()
* ****************************************************************************
* 
* 1   Purpose
*     This function evaluates the periodic B-spline basis at a given
*     parameter value.
* 
* 2   Specification
*     #include "bspl.h"
* 
* 3   Description
*     This function calculates all the non-zero periodic B-spline basis 
*     functions for a given parameter value using the Cox-Deboor algorithm.
*     These values are all calculated recursively and are placed in the matrix
*     N.  The indexing of the matrix N has been adjusted to keep the size of
*     the matrix independent of the index corresponding to B-spline control
*     points and only dependent on the order of the curve.  
*     Here we have set Ni,k(u) = N[i + k - tl][k](u) 
*     where Ni,k(u) is the the ith periodic B-spline basis function of order k.
*     The value of tl provides the neccessary shift in the index. This 
*     complete evaluation of all non-zero periodic B-spline basis functions
*     allows the calculation of the  derivatives with no additional
*     computational expense. It first adds (k-1) spans to the knot vector of
*     the periodic B-spline curve to make the knot vector repeat itself as a
*     ring structure, and then evaluate the B-spline basis using the general
*     recursive formula of B-spline basis. Note that due to the shift of the
*     index of knot vector, it is also necessary to remember the indices
*     of control points which correspond to those nonzero basis functions.
* 
* 4   References
*     [1] C. De Boor. A Practical Guide to Splines, Springer, New York, 1978.
*  
* 5   Parameters
*     1.int k
*       On entry: order of the periodic B-spline basis function.
*     2.int n
*       On entry: number of control points of the periodic B-spline curve.
*     3.double u
*       On entry: parameter value at which the basis function is to be
*                 evaluated.
*     4.double t[]
*       On entry: knot vector of the periodic B-spline basis, with length n+1.
*     5.double **N
*       On exit:  elements are the non-zero evaluated basis functions.
*     6.int *P
*       On exit:  array containing the indices of those control points whose
*                 associated basis functions have nonzero values.
* 
* 6   Return Values, Error Indicators and Warnings
*     Not applicable.
* 
* 7   Accuracy
*     Not applicable.
* 
* 8   Further Comments
*     Not applicable
* 
* 9   Functions referenced by nbasisd_per() are:
*     dbl_array1()
*     find()
*     free_darray1()
*     nbasisd()
* 
* 10  Functions that reference nbasisd_per() are:
*     calc_gram_percurv()
*     evalbsp_per()
*     evalderivbsp_per()
*     evalderivsurf_per()
*     evalrsurf_per()
*     evalsurf_per()
*     n_matrix_per()
* 
****************************************************************************/
void nbasisd_per(int k, int n, double u, double t[], double **N, int *P)
{
  double *t0;          /* knot vector in ring structure */
  int i,j,l,m,tl;      /* i,j -- loop index
			* l,m,tl -- index   */
			  
  /* allocate memory */
  t0 = dbl_array1((unsigned) k+n);

  j = n+1 ;            /* here, j is the length of the knot vector */
  tl = find(j,t,u) ;   /* find index tl such that t[tl] <= u <= t[tl+1]. */
  m = k-tl-1 ; 
  /* if 0 <= tl < k-1, i.e., tl is in the first (k-1) spans, repeat the last
   * (m) spans of t[] at the beginning of t0[], and repeat the first tl spans
   * of t[] at the end of t0[]. */
  if (m>0) {
     j = m ; 
     for (i=0;i<=n;i++) 
         t0[j++] = t[i] ;
     l = 0 ; 
     for (j=m;j<k;j++) 
         P[j] = l++ ;
     j = m-1 ; 
     l = n ; 
     for (i=0;i<m;i++) {
         t0[j] = t0[j+1]-t[l]+t[l-1] ; 
	 P[j] = l-1 ; 
	 l-- ; j-- ;
         }
     if (tl) {
        j = n+m+1 ; l = 0 ; 
        for (i=0;i<tl;i++) {
	    t0[j] = t0[j-1]+t[l+1]-t[l] ; 
	    j++ ; l++ ;
            }
       }
     }
  /* if tl >= k-1, simply repeat the first (k-1) spans of t[] at the end of 
   * the new knot vector t0[]. */
  else {
     j = 0 ; 
     for (i=0;i<=n;i++) 
         t0[j++] = t[i] ;
     for (i=0;i<k;i++) 
         P[i] = (tl+i-k+1)%n ;
     for (i=1;i<k;i++) {
         t0[j] = t0[j-1]+t[i]-t[i-1] ; 
	 j++ ;
        }
     }

  j = n+k ;            /* j now is the length of the new knot vector t0. */
  tl = find(j,t0,u);   /* find tl such that t0[tl] <= u <= t0[tl+1]. */
  /* evaluate the periodic B-spline basis with knot vector t[] by evaluating
   * the non-periodic B-spline basis with knot vector t0[]. */
  nbasisd(tl,k,u,t0,N) ;

  /* free memory */
  free_darray1(t0);
}

/****************************************************************************
*                                   n_matrix_per()
*****************************************************************************
* 
* 1   Purpose
*     This function calculates the Gram matrix of a periodic B-spline curve.
* 
* 2   Specification
*     #include "bspl.h"
*     void n_matrix_per(double **nij, double *knot, double *points,
* 		        short order, short npoints)
* 
* 3   Description
*     The element nij of the Gram matrix of a periodic B-spline, is the j-th
*     B-spline basis evaluated at i-th node.
* 
* 4   References
*     Not applicable
* 
* 5   Parameters
*     1.double **nij
*       On exit:  the Gram matrix
*     2.double *knot
*       On entry: the knot vector of the periodic B-spline
*     3.double *point
*       On entry: the nodes of the periodic B-spline
*     4.short order
*       On entry: order of the periodic B-spline
*     5.short npoints
*       On entry: length of the knot vector
* 
* 6   Return Values, Error Indicators and Warnings
*     Not applicable.
* 
* 7   Accuracy
*     Not applicable.
* 
* 8   Further Comments
*     Not applicable
* 
* 9   Functions referenced by n_matrix_per() are:
*     dbl_array2()
*     free_darray2()
*     free_iarray1()
*     int_array1()
*     nbasisd_per()
* 
* 10  Functions that reference n_matrix_per() are:
*     interp_points_per()
* 
*****************************************************************************/

void n_matrix_per(double **nij, double *knot, double *points, short order,
		  short npoints)
{
  double **N;      /* the periodic B-spline basis */
  int i,j,k,*P ;   /* i,j -- loop indices 
		    * P -- the indices of those control points whose associated
		    *      B-spline basis functions are nonzero    */

  /* allocate memory */
  P = int_array1(order+1);
  N = dbl_array2(order+2, order+1);

  /* initializtion */
  for (i=0; i<npoints; i++) 
      for (j=0; j<npoints;j++)
          nij[i][j] = 0.0 ;

  /* calculate the Gram matrix */
  for (i=0; i<npoints;i++) {
      nbasisd_per(order, npoints, points[i], knot, N, P) ;
      for (j=0; j<order; j++)
          nij[P[j]][i] = N[j+1][order] ;
  }

  /* free memory */
  free_darray2(N);
  free_iarray1(P);
}

/*****************************************************************************
*                               nodes_bsp_per()
******************************************************************************
* 
* 1   Purpose
*     This function computes the nodes of a periodic B-spline basis.
* 
* 2   Specification
*     #include "bspl.h"
*     void nodes_bsp_per(double knot[], int number, int order,double node[])
* 
* 3   Description
*     The nodes of a periodic B-spline basis of order M over the knot vector 
*     T = { t[0], t[1], ... , t[n] } are given by
*                      1
*          node[i] = ----- ( t[i+1] + t[i+2] + ... + t[i+m-1] )
*                     M-1
*     where 0 <= i <= n-1.
* 
* 4   References
*     Not applicable
* 
* 5   Parameters
*     1.double knot[]
*       On entry: knot vector of length (n+1).
*     2.int number
*       On entry: number of control points of the periodic B-spline curve.
*                 the length of knot vector is number+1.
*     3.int order
*       On entry: order of the periodic B-spline basis.
*     4.double node[]
*       On exit:  array containing the nodes.
* 
* 6   Return Values, Error Indicators and Warnings
*     Not applicable.
* 
* 7   Accuracy
*     Not applicable.
* 
* 8   Further Comments
*     The value of tl is determined by a call the the routine find 
*     (bspl _ mscel.c).
* 
* 9   Functions referenced by nodes_bsp_per() are: None
* 
* 10  Functions that reference nodes_bsp_per() are:
*     opt_param_per()
* 
**************************************************************************/

void nodes_bsp_per(double knot[], int number, int order,double node[])
{
  int i,j,k,l ;      /* indices */
  double s,t,u,v ;   /* s -- sum of knots
		      *	t -- values of knots
		      *	u -- order - 1
		      *	v -- node      */

  u = (double)(order-1) ; 
  /* calculate node[i], where i is from 0 to n-1. */
  for (i=0;i<number;i++) {
      s = 0 ; /* initialize the sum */
      /* do the summation */
      for (j=1;j<order;j++) {
          k = i+j ; 
	  /* if the index of the knot is over the maximum, start from the 
	   * beginning again. */
	  if (k>number) {
             l = k%number ; 
	     t = knot[number]+knot[l]-knot[0] ;
             }
          else t = knot[k] ;
          s += t ;
          }
      v = s/u ; 
      node[i] = (v<1.0)? v : v-1.0 ;
      }
}

/*****************************************************************************
*                                evalbsp_per()
******************************************************************************
* 
* 1   Purpose
*     This function evaluates a periodic B-spline curve.
* 
* 2   Specification
*     #include "bspl.h"
*     vector *evalbsp_per(ParCurv *egeom, double u)
* 
* 3   Description
*     This routine calculates and returns a vector data structure that 
*     contains the point value of the given periodic B-spline curve at the 
*     specific parameter value.
* 
* 4   References
*     [1] C. De Boor. A Practical Guide to Splines, Springer, New York, 1978.
* 
* 5   Parameters
*        1.ParCurv * egeom
*          On entry: NURBS data structure containing geometry of curve to be 
* 	           evaluated.
*          On exit:
*        2.double u
*          On entry: the parameter value at which the periodic B-spline curve 
* 	           is to be evaluated.
* 
* 6   Return Values, Error Indicators and Warnings
*     The value of u should lie between egeom -> knots[0] and egeom ->
*     knots[k]  where k = (egeom -> ncontpts) + (egeom -> order) - 1.  If not,
*     it will print an error message and stop.
* 
* 7   Accuracy
*     Not applicable.
* 
* 8   Further Comments
*     This basis function evaluation is performed using routine nbasisd_per
*     (bspl).
*     This routine should not be used to evaluate derivatives of periodic 
*     NURBS curves. For the evaluation of periodic NURBS curves, see
*     rbspeval_per (bspl) for details.
* 
* 9   Functions referenced by evalbsp_per() are:
*     dbl_array2()
*     errormsg()
*     free_darray2()
*     free_iarray1()
*     int_array1()
*     nbasisd_per()
*     vectalloc()
* 
* 10  Functions that reference evalbsp_per() are:
*     DrawCosEvaluate()
*     DrawCosKnots()
*     DrawCurvKnots()
*     EvaluateCos()
*     EvaluateCurv()
*     GetCosValues()
*     GetCurvValues()
*     normalcurv()
*     PostScriptCosEvaluate()
*     PostScriptCosKnots()
*     PostScriptCurvKnots()
*     rbspeval_per()
*     sample_per()
*     SaveIgesCos()
*     SaveIgesCurv()
* 
*****************************************************************************/

vector *evalbsp_per(ParCurv *egeom, double u)
/* evaluate a periodic integral B-spline at u */
{
  int ncontpts,order,i,*P;    /* ncontpts -- number of control points 
			       * order -- order of the B-spline basis
			       * i -- loop index
			       * P -- array containing the indices of control
			       *      points associated with nonzero basis at 
			       *      a given parameter value                */
  vector *p,*eval;            /* p -- vector containing a control point
			       * q -- vector containing the value of the
			       *      evaluated point              */
  double q,**N;               /* q -- a B-spline basis
			       * N -- the matrix containing B-spline basis   */

  /* get the number of control points and the order */
  ncontpts = egeom->ncontpts; order = egeom->order;

  /* allocate memory */
  P = int_array1((unsigned) order+1);
  N = dbl_array2((unsigned) order+1,(unsigned) order+1);

  /* if u is not inside the range of the curve. */
  if ( u > egeom->knots[ncontpts] || u < egeom->knots[0] )
     {
      errormsg(3,"outside  the knot range in function evalbsp_per\n");
      getchar();
     }
  /* if u is inside the range. */
  else  
     {
      /* calculate the matrix containing periodic B-spline basis. */
      nbasisd_per(order,ncontpts,u,egeom->knots,N,P); 
      /* initialize the evaluated value. */
      eval = vectalloc();
      eval->x = eval->y = eval->z = eval->w = 0.0;
      /* evaluate the point. summation is only applied on those control points 
       * which have nonzero coefficients/basis. */
      for(i=0; i<order; i++)
	{
	  p = egeom->contpts[P[i]] ; q = N[i+1][order];
	  eval->x += q*p->x ; eval->y += q*p->y;	
	  eval->z += q*p->z ; eval->w += q*p->w;	
	}
    }
  
  /* free memory */
  free_darray2(N);
  free_iarray1(P);

  return(eval);
}

/**************************************************************************
*                             evalderivbsp_per()
***************************************************************************
* 
* 1   Purpose
*     This function evaluates a periodic non-uniform integral B-spline curve 
*     or derivative.
* 
* 2   Specification
*     #include "bspl.h"
*     vector *evalderivbsp_per(ParCurv *egeom, double u, int deriv)
* 
* 3   Description
*     This routine calculates and returns a vector data structure that 
*     contains the point or derivative value of the given periodic non-uniform 
*     integral B-spline curve at the specific parameter value.
* 
* 4   References
*     [1] C. De Boor. A Practical Guide to Splines, Springer, New York, 1978.
* 
* 5   Parameters
*        1.ParCurv * egeom
*          On entry: NURBS data structure containing geometry of curve to be 
* 	           evaluated.
*        2.double u
*          On entry: the parameter value at which the non-uniform integral 
* 	           B-spline curve is to be evaluated.
*        3.int deriv
*          On entry: an index referring to the required derivative to be 
* 	           calculated such that if the curve is defined as R(u) then 
* 		   the routine will evaluate
*                                        derive
*                                       @      R(u)
*                                      -------------
*                                            deriv
*                                          @u
* 		   ( @ represents partial derivative. )
* 
* 6   Return Values, Error Indicators and Warnings
*     The index deriv <= the order of the perodic non-uniform integral
*     B-spline curve.
*     The value of u should lie between egeom -> knots[0] and egeom ->
*     knots[k]  where k = (egeom -> ncontpts) + (egeom -> order) - 1.  If not,
*     it will print an error message and stop.
* 
* 7   Accuracy
*     Not applicable.
* 
* 8   Further Comments
*     This basis function evaluation is performed using routine nbasisd_per
*     (bspl).
*     This routine should not be used to evaluate derivatives of periodic 
*     NURBS curves. 
*     For the evaluation of periodic NURBS curves, see rbspeval_per (bspl)
*     for details.
* 
* 9   Functions referenced by evalderivbsp_per() are:
*     Aij_per()
*     dbl_array2()
*     dbl_array3()
*     free_darray2()
*     free_darray3()
*     free_iarray1()
*     int_array1()
*     nbasisd_per()
*     vectalloc()
* 
* 10  Functions that reference evalderivbsp_per() are:
*     rbspeval_per()
* 
****************************************************************************/

vector *evalderivbsp_per(ParCurv *egeom, double u, int deriv)

/* compute the parametric derivatives of a periodic integral B-spline at
 * u. deriv is the order of the derivative. */
{
  int prod,ncontpts,order,i,*P; 
                /* prod -- the scalor factor in the derivative of B-spline
                 *         curve
		 * ncontpts -- number of control points
		 * order -- order of the B-spline basis
		 * i -- loop index
		 * P -- array containing the indices of the control points 
		 *      associated with nonzero basis functions      */
  vector *p,*eval;
                /* p -- never used (glshen) 
		 * eval -- the evaluated value */
  double q,**N,***derivpts;
                /* q -- value of periodic B-spline basis 
		 * N -- the periodic B-spline basis matrix 
		 * deripts -- control points of the derivatives  */

  prod = 1; 

  /* get the number of control points and the order of the curve. */
  ncontpts = egeom->ncontpts; order = egeom->order;

  /* allocate memory */
  P = int_array1((unsigned) order+1);
  N = dbl_array2((unsigned) order+1,(unsigned) order+1);
  derivpts = dbl_array3((unsigned) order,(unsigned) order,4);

  /* if u is not inside the range of the curve. */
  if( u > egeom->knots[ncontpts] || u < egeom->knots[0] )
    {
      printf("u = %f outside the knot range in fn. evalderivbsp\n", u);
      exit(0);
    } 
	
  /* calculate the factor to be mulitplied. */
  for (i=1; i<=deriv; i++) 
      prod *= order-i;

  /* evaluate the periodic B-spline basis */
  nbasisd_per(order,ncontpts,u,egeom->knots,N,P);

  /* calculate the control points of the derivative with no factor */
  Aij_per(ncontpts,deriv,order,P,egeom->knots,egeom->contpts,derivpts);

  /* initialize the value of evaluated point */
  eval = vectalloc(); eval->x = eval->y = eval->z = eval->w = 0.0;

  /* evaluate the point at the given parameter value */
  for (i=0;i<order-deriv;i++) {
      q = N[i+1][order-deriv];
      eval->x  += q*derivpts[i+deriv][deriv][0];
      eval->y  += q*derivpts[i+deriv][deriv][1];
      eval->z  += q*derivpts[i+deriv][deriv][2];
      eval->w  += q*derivpts[i+deriv][deriv][3];
     }

  /* multiply the factor */
  eval->x = prod*eval->x ; eval->y = prod*eval->y;
  eval->z = prod*eval->z ; eval->w = prod*eval->w;

  /* free memory */
  free_darray2(N);
  free_iarray1(P);
  free_darray3(derivpts);

  return(eval);
}

/****************************************************************************
*                               rbspeval_per()
*****************************************************************************
* 
* 1   Purpose
*     This function evaluates a periodic NURBS curve and computes derivatives.
* 
* 2   Specification
*     #include "bspl.h"
*     vector *rbspeval_per(ParCurv *egeom, double u, int deriv)
* 
* 3   Description
*     This routine calculates and returns a vector data structure that 
*     contains the point or derivative of a point value of the given periodic
*     NURBS curve at the specific parameter value. The returned vector data
*     structure contains values for all four homogeneous coordinates.
* 
* 4   References
*     [1] C. De Boor. A Practical Guide to Splines, Springer, New York, 1978.
* 
* 5   Parameters
*        1.ParCurv * egeom
*          On entry: NURBS data structure containing geometry of curve to be 
* 	           evaluated.
*        2.double u
*          On entry: the parameter value at which the periodic NURBS curve is 
* 	           to be evaluated.
*        3.int deriv
*          On entry: index denoting required derivative where if the curve is 
* 	           defined as R(u) then the routine will return
*                                        deriv
*                                       d     R(u)
*                                      ------------
*                                           deriv
*                                         du
* 
* 6   Return Values, Error Indicators and Warnings
*     The value of u should lie between egeom -> knots[0] and egeom ->knots[k] 
*     where k = (egeom -> ncontpts)+(egeom -> order)-1. If not, it will print 
*     an error message and stop.
*     The value for deriv <= 4.
* 
* 7   Accuracy
*     Not applicable.
* 
* 8   Further Comments
*     This basis function evaluation is performed using routine nbasisd_per() 
* 
* 9   Functions referenced by rbspeval_per() are:
*     evalbsp_per()
*     evalderivbsp_per()
* 
* 10  Functions that reference rbspeval_per() are:
*     calc_disc()
*     cos_conv_nurbs()
*     deriv_dev_per()
*     DrawCurvTorsion()
*     EvaluateCurv()
*     find_error_percurv()
*     find_points()
*     GetCurvTang()
*     GetCurvValues()
*     max_curve_curvature()
*     normalcurv()
*     offset_curv2d()
*     offset_geod()
*     offset_geod_par()
*     offset_normal()
*     PostScriptCurvTorsion()
* 
******************************************************************************/

vector *rbspeval_per(ParCurv *egeom, double u, int deriv)
/* evaluate a periodic, rational B-spline */
{
  vector *eval0,*eval1,*eval2; /* the 1st, 2nd, and 3rd derivative. */
  double w,wp,wpp;             /* w-coordinates of homogeneous system */

  /* evaluate the periodic B-spline curve at parametric value u */
  eval0 = evalbsp_per(egeom,u);
  w = eval0->w;

  /* if only the value at the given point is wanted. */
  if (deriv == 0)
     return(eval0);

  /* if the order of derivative is no less than 1. */
  if (deriv >= 1)
     {
      /* evaluate the 1st derivative */
      eval1 = evalderivbsp_per(egeom,u,1);
      wp = eval1->w;
      /* calculate the homogeneous coordinates */
      eval1->x = (eval1->x - (eval0->x/w)*wp)/w;
      eval1->y = (eval1->y - (eval0->y/w)*wp)/w;
      eval1->z = (eval1->z - (eval0->z/w)*wp)/w;
      eval1->w = 1.0;
      /* if the 1st derivative is wanted. */
      if (deriv == 1)
	 {
	  free((char *) eval0);
	  return(eval1);
	 }
      }

  /* if the order of derivative is no less than 2. */
  if (deriv >= 2)
     {
      /* evaluate the 2nd derivative */
      eval2 = evalderivbsp_per(egeom,u,2);
      wpp = eval2->w;
      /* calculate the homogeneous coordinates. */
      eval2->x = (eval2->x - 2.0*eval1->x*wp - (wpp*eval0->x/w))/w;
      eval2->y = (eval2->y - 2.0*eval1->y*wp - (wpp*eval0->y/w))/w;
      eval2->z = (eval2->z - 2.0*eval2->z*wp - (wpp*eval0->z/w))/w;
      eval2->w = 1.0;
      /* if the 2nd derivative is wanted. */
      if (deriv == 2)
	 {
	  free((char *) eval0);
          free((char *) eval1);
          return(eval2);
	  }
      }
}

/***************************************************************************
*                                  Aij_per()
****************************************************************************
* 
* 1   Purpose
*     This is a subroutine of function evalderivbsp_per(), which evaluates the 
*     derivatives of a periodic B-spline curve, see evalderivbsp_per(). 
* 
* 2   Specification
*     #include "bspl.h"
*     void Aij_per(int n, int q, int k, int P[], double t[], vector *contpts[],
* 	         double ***derivpts)
* 
* 3   Description
*     This function calculates the control points of the q-th derivative curve
*     of a periodic B-spline curve. The derivative curve itself is also a
*     periodic B-spline curve with lower order (order of the curve - q ). The
*     control points  of the derivative curve are recursively determined from
*     the control points of  the original curve.
* 
* 4   References
*     [1]C. De Boor. A Practical Guide to Splines, Springer, New York, 1978.
* 
* 5   Parameters
*       1.int n
*         On entry: the number of control points, and also (n+1) is the
* 	          length of the  knot vector.
*       2.int q
*         On entry: the order of the derivative.
*       3.int k
*         On entry: the order of the curve to be evaluated.
*       4.int P[]
*         On entry: the array containing the indices of the control points with
* 	          nonzero associated B-spline basis function.
*       5.double t[]
*         On entry: the knot vector of the curve to be evaluated.
*       6.vector *contpts[]
*         On entry: the control points of the curve to be evaluated.
*       7.double ***derivpts
*         On entry: a 3D array
*         On exit:  the 3D array containing the control points of the
*                   derivative curves with derivative order varying from 0 to
*                   deriv.
* 
* 6   Return Values, Error Indicators and Warnings
*     Only derivpts[i][deriv][], deriv <= i <= n, are used by function
*     evalderivbsp_per(). Here n is # of control points of the original curve.
* 
* 7   Accuracy
*     Not applicable.
* 
* 8   Further comments
*     Not appicable
* 
* 9   Functions referenced by Aij_per() are: None
* 
* 10  Functions that reference Aij_per() are:
*     evalderivbsp_per()
* 
*****************************************************************************/

void Aij_per(int n, int q, int k, int P[], double t[], vector *contpts[],
	     double ***derivpts)
/************************************************************************ 
* The Aij_per are such that A sub(i) sup(j) = A[i-il+(k-1)][j]
* Although the variable il is not available here, this value is useful
* in indexing A[][] compactly and to calculate only Aij required for a 
* a given il, see program evalderivbsp() where Aij is used as such. 
**********************************************************************/
{
  vector *v1, *v2, *v3; /* v2,v3 -- never used (glshen) */ 
  double dt;            /* difference between two knots */
  int i,j,ireal;        /* indices */

  /* set the control points of 0-th derivative curve as that of the original
   * curve. */
  for (i=0; i<k; i++) {
      v2 = contpts[P[i]];
      derivpts[i][0][0] = v2->x; derivpts[i][0][1] = v2->y;
      derivpts[i][0][2] = v2->z; derivpts[i][0][3] = v2->w;
      }

  /* calculate the control points of j-th derivative curves, where j varies
   * from 1 to q. */
  for (j=1; j<=q; j++) 
      for (i=j; i<k; i++) {
          ireal = P[i]+k-j ;      /* index translation */
	  /* if the index ireal is beyond the maximum, repeat from the
           * beginning.
	   * in this case, dt = t[ireal]-t[P[i]] < 0, and should be adjusted
           * as 1.0-dt. */
          if (ireal>n) { ireal = ireal%n ; dt = 1.0 ; } 
          else dt = 0.0 ;
          dt += t[ireal] - t[P[i]]; 
          if (dt < ZERO)
	     printf("dt is zero in routine Aij()\n");	
          derivpts[i][j][0] = (derivpts[i][j-1][0] - derivpts[i-1][j-1][0])/dt;
          derivpts[i][j][1] = (derivpts[i][j-1][1] - derivpts[i-1][j-1][1])/dt;
          derivpts[i][j][2] = (derivpts[i][j-1][2] - derivpts[i-1][j-1][2])/dt;
          derivpts[i][j][3] = (derivpts[i][j-1][3] - derivpts[i-1][j-1][3])/dt;
	 }
}

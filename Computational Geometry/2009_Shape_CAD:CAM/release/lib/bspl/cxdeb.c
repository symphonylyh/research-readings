/* ***************************************************************************
 Copyright (C) 1996 Massachusetts Institute of Technology all rights reserved 
	Programmer: George A. Kriezis
**************************************************************************** */
# include "gen.h"
# include "bspl.h"

/*****************************************************************************
*                                nbasisd()
******************************************************************************
* 
* 1   Purpose
*     This function computes the B-spline basis function for the given
*     parameters.
* 
* 2   Specification
*     #include "bspl.h"
*     void nbasisd(int tl, int k, double u, double *t, double **N)
* 
* 3   Description
*     This routine calculates all the non-zero B-spline basis functions for a 
*     given parameter value using the Cox-Deboor algorithm.  These values are 
*     all calculated recursively and are placed in the matrix N.  The indexing 
*     of the matrix N has been adjusted to keep the size of the matrix 
*     independent of the index corresponding to B-spline control points and 
*     only dependent on the order of the curve.  Here we have set Ni,k(u) =
*     N[i + k - tl][k](u) 
*     where Ni,k(u) is the the ith B-spline basis function of order k. The 
*     value of tl provides the neccessary shift in the index. This complete 
*     evaluation of all non-zero B-spline basis functions allows the 
*     calculation of the derivatives with no additional computational expense.
* 
* 4   References
*     [1] C. De Boor. A Practical Guide to Splines, Springer, New York, 1978.
* 
* 5   Parameters
*        1.double u
*          On entry: parameter value at which the basis function is evaluated.
*        2.int k
*          On entry: order of the B-spline basis function.
*        3.int tl
*          On entry: is the index such that t[tl] <= u <= t[tl + 1].
*        4.double t
*          On entry: array of length >= tl that contains the knot vector to be 
* 	           used by the B-spline basis function.
*        5.double ** N
*          On entry: matrix of size order x order.
*          On exit:  elements are the non-zero evaluated basis functions.
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
* 9   Functions referenced by nbasisd() are:
*     dbl_array1()
*     errormsg()
*     free_darray1()
* 
* 10  Functions that reference nbasisd() are:
*     build_offset()
*     build_offsets_int()
*     calc_gram_curv()
*     calc_gram_surf()
*     evalbsp()
*     evalderivbsp()
*     evalderivsurf()
*     evalderivsurf_per()
*     evalrsurf()
*     evalrsurf_per()
*     evalsurf()
*     evalsurf_per()
*     integral_build_offset()
*     loft_integral()
*     loft_rational()
*     nbasisd_per()
*     ParCurv_iso()
*     tol_loft_rational()
*     zevalderivsurf()
* 
*****************************************************************************/

void nbasisd(int tl, int k, double u, double t[], double **N)
{
  int r,s;      /* r -- the index of r-th basis. 
		 * s -- the index of order. */
  double  M, *dp, *dm, test;   /* working variables and matices. */

  /************************************************************************
   * Here we set N sub(i,k)(u) = N[i+k-tl][k](u)
   * So use it in the above form - this reassignment is to keep the 
   * order of the matrix sizes to known limits - here tl is such that  
   * t(tl) <= u < t(tl+1) 
   **********************************************************************/

  /* allocate memory */
  dp = dbl_array1((unsigned)k);
  dm = dbl_array1((unsigned)k);

  /* initialize the B-spline basis matrix to be evaluated */
  for (r=0; r<k+1; r++)
      for (s=0; s<k+1; s++)
	  N[r][s] = 0.0;

  /* N sub(tl,1) = 1, because t[tl] <= u <= t[tl+1] */
  N[1][1] = 1.0;     

  /* evaluate the basis matrix by recursive definition of B-spline basis */
  for (s=1; s<=k-1; s++)   /* order varies from 1 to k-1 */
      {
       dp[s] = t[tl+s] - u; 
       dm[s] = u - t[tl-s+1];
       N[1][s+1] = 0.0;
       for (r=1; r<=s; r++)  /* r-th basis of order s */
	   {
	    test = dp[r] + dm[s -r + 1];
	    if (test == 0.0){
	       errormsg(1,"error in nbasisd, divide by zero");
	       M = 0.0;
	       }
	    else
	       M = N[r][s]/test;
	    N[r][s+1] = N[r][s+1] + dp[r]*M;  
	    N[r+1][s+1] = dm[s - r + 1]*M;
	   }	
      }

  /* free memory */
  free_darray1(dp);
  free_darray1(dm);
}

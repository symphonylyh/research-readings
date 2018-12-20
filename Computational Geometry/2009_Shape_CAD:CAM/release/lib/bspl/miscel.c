/* ***************************************************************************
 Copyright (C) 1996 Massachusetts Institute of Technology all rights reserved 
	Programmer: George A. Kriezis
**************************************************************************** */
#include <math.h>
#include "gen.h"
#include "bspl.h"

#define MACHPREC 1.1102230246251567E-13	  /* changed on 25 JUN 90 */

/******************************************************************************
*                                    lf2comp()
*******************************************************************************
* 
* 1   Purpose
*     Compare double precision values with tolerance.
* 
* 2   Specification
*     #include "bspl.h"
*     int lf2comp(double a, double b, double eps);
* 
* 3   Description
*     This function compares two double precision values with a specified 
*     tolerance.
* 
* 4   References
*     Not applicable.
* 
* 5   Parameters
*       1.double a
*         On entry: the first value to be compared.
*       2.double b
*         On entry: the second value to be compared.
*       3.double eps
*         On entry: the tolerance value.
* 
* 6   Return Values, Error Indicators and Warnings
*     If a >= b + eps then return 1, if a <= b - eps then return -1, otherwise 
*     return 0.
* 
* 7   Accuracy
*     Not applicable.
* 
* 8   Further Comments
*     Not applicable.
* 
* 9   Functions referenced by lf2comp() are:
*     errormsg()
* 
* 10  Functions that reference lf2comp() are:
*     bernrs()
*     calc_intsct()
*     clip_bezier_curve()
*     common_interval()
*     cvx_hull_axis_intrsct()
*     find()
*     lowindex()
*     min_bez_curve_to_axis_intrsct()
*     subbezier()
* 
******************************************************************************/

int lf2comp(double a, double b, double eps)
{
double c;

c = a-b;

if	( c >= eps)
	return(1);
else if	( -eps < c && c < eps )
	return(0);
else if	( c <= -eps )
	return(-1);
else 
	errormsg(1,"error in lf2comp()");
return (0);
} 

/*****************************************************************************
*                                bernmat()
******************************************************************************
* 
* 1   Purpose
*     Calculate the transformation matrix from Bernstein polynomials to the 
*     monomial basis with degree n.
* 
* 2   Specification
*     #include "bspl.h"
*     void bernmat(int order, double mat[], int dim2)
* 
* 3   Description
*     This function calculates the matrix M = {m(i,j)} which describes the
*     basis transformation between Bernstein polynomials and the monomial
*     basis.                     j-i
*                            (-1)   n!
*                m(i,j) = ----------------
*                          i!(n-j)!(j-i)!
*     where degree n = order-1
* 
* 4   References
*     Not applicable
* 
* 5   Parameters
*       1.int order
*         On entry: the order of the Bezier curve.
*       2.double mat[]
*         On exit:  the transformation matrix.
*       3.int dim2
*         On entry: the dimension of each row in the matrix.
* 
* 6   Return Values, Error Indicators and Warnings
*     Because matrix M is symmetric, only the upper triangular metrix
*     is calculated.
* 
* 7   Accuracy
*     Not applicable
* 
* 8   Further comments
*     Not applicable
* 
* 9   Functions referenced by bernmat() are:
*     factorial()
* 
* 10  Functions that reference bernmat() are:
*     berntomono()
*     monotobern()
* 
*****************************************************************************/

void bernmat(int order, double mat[], int dim2)
{
int i,j,m;

m = order;

for(i=0; i<m; i++)
  for(j=0; j<m; j++)
    {
      if (j<i) 
	mat[i*dim2+j] = 0.0;  /* set the element 0 if i<j */
      else
	mat[i*dim2+j] = UP_I((double)-1.0,(double)j-i)*factorial(m-1)/
	  (factorial(i)*factorial(m-1-j)*factorial(j-i));
    }
}

/*****************************************************************************
*                                combine()
******************************************************************************
* 
* 1   Purpose
*     Calculate ratio of factorials.
* 
* 2   Specification
*     #include "bspl.h"
*     double combine(int n, int r);
* 
* 3   Description
*     This function calculates the ratio of factorials, n!/[r!(n-r)!].
* 
* 4   References
*     Not applicable.
* 
* 5   Parameters
*       1.int n
*         On entry: the first value.
*       2.int r
*         On entry: the second value.
* 
* 6   Return Values, Error Indicators and Warnings
*     The ratio n!/[r!(n-r)!] is returned.
* 
* 7   Accuracy
*     Not applicable.
* 
* 8   Further Comments
*     The function factorial() is used to calculate the factorials.
* 
* 9   Functions referenced by combine() are:
*     factorial()
* 
* 10  Functions that reference combine() are:
*     linear_trans()
* 
*****************************************************************************/

double combine(int n, int r)
{
  return(factorial(n)/(factorial(r)*factorial(n-r)) );
}

/******************************************************************************
*                                factorial()
*******************************************************************************
*
* 1   Purpose
*     Calculate factorial.
*
* 2   Specification
*     #include "bspl.h"
*     double factorial(int r);
* 
* 3   Description
*     This function calculates the factorial, r! = r(r - 1)(r - 2)..
* 
* 4   References
*     Not applicable.
* 
* 5   Parameters
*        1.int n
*          On entry: the value.
* 
* 6   Return Values, Error Indicators and Warnings
*     The value of the factorial r! = r(r - 1)(r - 2). is returned.
* 
* 7   Accuracy
*     Not applicable.
* 
* 8   Further Comments
*     Not applicable.
* 
* 9   Functions referenced by factorial() are: None
* 
* 10  Functions that reference factorial() are:
*     bernmat()
*     combine()
*     knot_holzle()
*     knot_holzle_fn()
* 
******************************************************************************/

double factorial(int r)
{
int i;
double prod;

prod = 1.0;

if(r==0)
	return(prod);
else
	{
	for(i=0; i<r; i++)
		prod *= r-i;
	return(prod);
	}
}

/***************************************************************************
*                             absmaxarray()
****************************************************************************
*
* 1   Purpose
*     Find the absolute maximum value in a 2D array.
* 
* 2   Specification
*     #include "bspl.h"
*     double absmaxarray(int m, int n, double w[], int dim2);
* 
* 3   Description
*     This routine finds the maximum absolute value (maximum magnitude) in a 
*     2D array.
* 
* 4   References
*     Not applicable.
* 
* 5   Parameters
*       1.int m
*         On entry: the number of rows in the array
*       2.int n
*         On entry: the number of columns in the array
*       3.double 
* 
* 6   Return Values, Error Indicators and Warnings
*     The value of the maximum absolute value (maximum magnitude) is returned.
* 
* 7   Accuracy
*     Not applicable.
* 
* 8   Further Comments
*     Not applicable.
* 
* 9   Functions referenced by absmaxarray() are: None
* 
* 10  Functions that reference absmaxarray() are: monotobern()
* 
***************************************************************************/

    

double absmaxarray(int m, int n, double w[], int dim2)
{
int i,j;
double max,wij;

max = 0.0;

for(i=0; i<m; i++)
	for(j=0; j<n; j++)
		{
		wij = FABS(w[i*dim2 + j]);
		if(wij > max) 
			max = wij;
		}
return(max);
}	

/***************************************************************************
*                              find()
****************************************************************************
* 
* 1   Purpose
*     This function finds the index in an array such that a given value is 
*     bounded by adjacent elements in the array.
* 
* 2   Specification
*     #include "gen.h"
*     int find(int kn, double *tau, double t1)
* 
* 3   Description
*     This routine returns a value i such that it is the last index of the 
*     array tau where the relation tau[i] <=  t1 is valid.
* 
* 5   Parameters
*       1.int kn
*         On entry: number of elements in the array tau or the number of 
* 	          elements of tau that need to be offset from tau[0].
*       2.double * tau
*         On entry: array of length>=kn whose elements form the bound in which 
* 	          t1 is to be located and are such that tau[i]<=tau[i+1].
*       3.double t1
*         On entry: the value whose position in the array is to be found such 
* 	          that tau[i] = t1.
* 
* 6   Return Values, Error Indicators and Warnings
*     The value returned is equal to the number of elements in the array if t1 
*     is greater than the range in tau[] and is equal to 0 if less than the 
*     range in tau[].
* 
* 
* 7   Accuracy
*     The accuracy is determined by the machine precision used in comparing 
*     two double precision numbers.
* 
* 8   Further Comments
*     Not applicable
* 
* 9   Functions referenced by find() are:
*     lf2comp()
* 
* 10  Functions that reference find() are:
*     build_offset()
*     build_offsets_int()
*     calc_gram_curv()
*     calc_gram_surf()
*     convexhull_test()
*     evalbsp()
*     evalderivbsp()
*     evalderivsurf()
*     evalderivsurf_per()
*     evalrsurf()
*     evalrsurf_per()
*     evalsurf()
*     evalsurf_per()
*     fn_camber()
*     fn_camber3()
*     fn_thckns()
*     fn_thckns3()
*     integral_build_offset()
*     integral_sample_offset()
*     int_sample_gencyl()
*     loft_integral()
*     loft_rational()
*     nbasisd_per()
*     ParCurv_iso()
*     PowCurv_eval()
*     PowSurf_eval()
*     PowSurf_iso()
*     sample_offset()
*     soslo()
*     subbezier()
*     subdivbs()
*     tol_loft_rational()
*     tol_sample_offsets()
*     transform()
*     trans_per()
*     zevalderivsurf()
* 
*****************************************************************************/

int find ( int kn, double tau[], double t1)
{
  int i,mu;
  mu = 0;

  for(i=0; i < kn; i++)
     {
      if(lf2comp(t1 , tau[i], MACHPREC) >= 0 )
         mu = i;
     }

  return(mu);
} 

/****************************************************************************
*                               lowindex()
*****************************************************************************
* 
* 1   Purpose
* 
* This function always returns the last index of tau[] if one of the elements
* is larger than tl. (glshen )
* 
****************************************************************************/

int lowindex( int kn, double tau[], double t1)
{
  int i, mu;

  for (i=0; i<kn; i++) {
    if (lf2comp(t1, tau[i], MACHPREC) >= 0)
      mu = i;
  }

  return(mu);
} 

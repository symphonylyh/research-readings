/* ***************************************************************************
 Copyright (C) 1996 Massachusetts Institute of Technology all rights reserved 
	Programmer: George A. Kriezis
**************************************************************************** */
#include "gen.h"
#include "bspl.h"

/*****************************************************************************
*                               berntomono()
******************************************************************************
* 
* 1   Purpose
*     This function converts a bivariate polynomial in Bernstein form to one 
*     in the monomial form.
* 
* 2   Specification
*     #include "bspl.h"
*     void berntomono(ParSurf fgeom, double **fuv, int *m, int *n)
* 
* 3   Description
*     This routine converts Bernstein basis representation of an algebraic 
*     curve (stored in a NURBS surface data format) to the monomial basis.
* 
* 4   References
*     Not applicable
* 
* 5   Parameters
*       1.ParSurf * fgeom
*         On entry: NURBS surface data structure containing the Bernstein 
* 	          coefficients as z values of the control points and 
* 		  appropriate open knot vectors in the u and v directions
* 		  along with uniformly spaced u and v values for the x and 
* 		  y values of the control points.
*       2.double ** fuv
*         On exit: Two dimensional array containing the coefficients of the 
* 	         algebraic curve in the monomial basis, the indices 
* 		 corresponding to the degree of the two variables.
*       3.int * m
*         On entry: Value of the order (degree + 1) in the u direction of 
* 	          the algebraic curve represented by fgeom.
*       4.int * n
*         On entry: Value of the order (degree + 1) in the v direction of the 
* 	          algebraic curve represented by fgeom.
* 
* 6   Return Values, Error Indicators and Warnings
*     Not applicable.
* 
* 7   Accuracy
*     Not applicable.
* 
* 8   Further Comments
*     Not applicable.
* 
* 9   Functions referenced by berntomono() are:
*     bernmat()
*     dbl_array2()
*     free_darray2()
*     matrixmu()
*     matrix_transpose()
* 
* 10  Functions that reference berntomono() are: None
* 
******************************************************************************/

void berntomono(ParSurf *wfgeom,  double **fuv, int *m, int *n)
{
  double **l1,**l1t;    /* transformation matrix and its transpose. */
  double **temp, **wf;  /* working matrix and control points matrix */
  int i,j,size;         /* i,j -- indices
			 * size -- size of matrices                 */

  /* determine the working space */
  *m = wfgeom->uorder;
  *n = wfgeom->vorder;
  if (*m >= *n)
     size= *m;
  else
     size= *n;

  /* allocate memory */
  l1 = dbl_array2((unsigned)size,(unsigned)size);
  l1t = dbl_array2((unsigned)size,(unsigned)size);
  temp = dbl_array2((unsigned)size,(unsigned)size);
  wf = dbl_array2((unsigned)size,(unsigned)size);

  /* matrix containing z-coordinates of the control points */
  for (i=0; i<*m; i++)
      for (j=0; j<*n; j++)
          wf[i][j] = wfgeom->contpts[i][j]->z;

  /* transformation matrix from Bernstein basis with order m to the monomial */
  /* basis do we need to calculate another transformation matrix for v
   * direction? (glshen) */
  bernmat(*m, &l1[0][0], size);
  /* its transpose */
  matrix_transpose(*m,*m,&l1[0][0],size,&l1t[0][0],size);
    
  /* matrices multiplication: fuv = l1t(mxm)*wf(mxn)*l1(nxn) */
  matrixmu(*m,*m,*n,size,size,size,&l1t[0][0], &wf[0][0], &temp[0][0]);
  matrixmu(*m,*n,*n,size,size,size,&temp[0][0],&l1[0][0], &fuv[0][0]);

  /* free memory */
  free_darray2(l1);
  free_darray2(l1t);
  free_darray2(temp);
  free_darray2(wf);
}

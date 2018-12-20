/* ***************************************************************************
 Copyright (C) 1996 Massachusetts Institute of Technology all rights reserved 
	Programmer: George A. Kriezis
**************************************************************************** */
#include <stdio.h>
#include <math.h>
#include "gen.h"
#include "bspl.h"

void simul_(int*,double*,double*,double*,int*,int*,int*,double*);
/* FORTRAN function, solving a set of simultaneous equations. */

/******************************************************************************
*                               monotobern()
*******************************************************************************
* 
* 1   Purpose
*     This routine converts a polynomial surface expressed in power basis to
*     bernstein basis.
* 
* 2   Specification
*     #include "bspl.h"
*     void monotobern(double **a, int m, int n, double xa, double xb,
*                     double ya, double yb, ParSurf *fgeom, int nondim)
* 
* 3   Description                                    
*     This function convert a polynomial surface z = f(x,y) in power basis to
*     a vector-valued patch {x,y,z} = sum{i=0..m}sum{j=0..n}ctrl_pts(i,j)
*     *Bi,m(u)Bj,n(v), where ctrl_pts are Bernstein coefficients, Bi,m(u) and
*     B(j,n) are Bernstein basis.
* 
* 4   References
*     [1] R. T. Farouke, V. T. Rajan,  Polynomials in Bernstein form
* 
* 5   Parameters
*       1.double **a                                                      i j
*         On entry: the coefficient matrix. a(i,j) is the coefficient of x y .
*       2.int m
*         On entry: the degree of x
*       3.int n
*         On entry: the degree of y
*       4.double xa
*         On entry: the lower bound of x
*       5.double xb
*         On entry: the upper bound of x
*       6.double ya
*         On entry: the lower bound of y
*       7.double yb
*         On entry: the upper bound of y
*       8.ParSurf *fgeom
*         On entry: a NURBS data structure
*         On exit:  the Bezier form of the input polynomial surface 
*       9.int nondim
*         On entry: an integer indicating whether to dimensionize, 1 for not,
*                   0 for yes.
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
* 9   Functions referenced by monotobern() are:
*     absmaxarray()
*     bernmat()
*     dbl_array1()
*     dbl_array2()
*     free_darray1()
*     free_darray2()
*     free_iarray2()
*     int_array2()
*     linear_trans()
*     matrixmu()
*     matrix_transpose()
*     simul_()
*     vectalloc()
* 
* 10 Functions that reference monotobern() are: None
* 
*****************************************************************************/

void monotobern(double **a, int m, int n, double xa, double xb, double ya,
		double yb, ParSurf *fgeom, int nondim) 
{
  vector *v;               /* working vector */
  int size,i,j,**iwks;     /* size -- the size of working space.
			    * i,j -- indices
			    * iwks -- an integer working space for simul_(). */
  double xk,yk,max;
  double eps,one;          /* tolerance and 1.0 (double) */
  double *x,*y;            /* working array for simul_() */
  double **acopy;          /* copy of matrix a */
  double **bx,**bxt,**by;  /* transformation matrices */
  double **l1,**l2,**l1t;  /* transformation matrices */ 
  double **binv,**temp1,**temp2;  /* working space */
  double **w;              /* matrix containing control points of Bezier
			    * patch. */

  if (m >= n)
     size= m;
  else
     size= n;

  /* allocate vectors */
  bx = dbl_array2((unsigned)size,(unsigned)size);
  bxt = dbl_array2((unsigned)size,(unsigned)size);
  by = dbl_array2((unsigned)size,(unsigned)size);
  l1 = dbl_array2((unsigned)size,(unsigned)size);
  l2 = dbl_array2((unsigned)size,(unsigned)size);
  l1t = dbl_array2((unsigned)size,(unsigned)size);
  binv = dbl_array2((unsigned)size,(unsigned)size);
  temp1 = dbl_array2((unsigned)size,(unsigned)size);
  w = dbl_array2((unsigned)size,(unsigned)size);
  temp2 = dbl_array2((unsigned)size,(unsigned)size);
  acopy = dbl_array2((unsigned)size,(unsigned)size);
  x = dbl_array1((unsigned)size);

  eps = 1.0e-15;
  one = 1.0;

  /* copy a to acopy such that a can remain unchanged. */
  for (i=0; i<m; i++)
      for(j=0; j<n; j++)
         acopy[i][j] = a[i][j];

  /* transformation matrices from Bernstein basis to monomial, in x,y
   * directions, respectively. */
  bernmat(m,&bx[0][0],size);
  bernmat(n,&by[0][0],size);

  /* x transformation matrix needs to be transposed for future use. */
  matrix_transpose(m,m,&bx[0][0],size,&bxt[0][0],size);

  /* obtain the transformation matrices from monomial to Bernstein basis
   * by inversing the transformation matrices from Bernstein to monomial. */
  i = -1;    /* indicate the simul function returns an inverse matrix. */
  j = size;  /* maximum size of the matrices. */
  /* allocate memory for working space */
  y = dbl_array1((unsigned)size);  
  iwks = int_array2((unsigned)size,3);
  /* inverse the matrices by solving a set of simultaneous equations. */
  simul_(&m,&bxt[0][0],x,&eps,&i,&j,&iwks[0][0],y);
  simul_(&n,&by[0][0],x,&eps,&i,&j,&iwks[0][0],y);
  /* free the working space for simul function. */
  free_iarray2(iwks);
  free_darray1(y);

  /* calculate the transformation matrices for linearization in x,y
   * directions. */
  linear_trans(m,l1,xb-xa,xa);
  linear_trans(n,l2,yb-ya,ya);
  /* transpose the x transformation matrix for pre-multiplication. */
  matrix_transpose(m,m,&l1[0][0],size,&l1t[0][0],size);

  /* calculate the final control points in z direction,
   * w = bxt*l1t*acopy*l2*by */
  matrixmu(m,m,m,size,size,size,&bxt[0][0],&l1t[0][0],&temp1[0][0]);
  matrixmu(m,m,n,size,size,size,&temp1[0][0],&acopy[0][0],&temp2[0][0]);
  matrixmu(m,n,n,size,size,size,&temp2[0][0],&l2[0][0],&temp1[0][0]);
  matrixmu(m,n,n,size,size,size,&temp1[0][0],&by[0][0],&w[0][0]);

  /* set the orders in u,v directions */
  fgeom->uorder = m;
  fgeom->vorder = n;
  fgeom->ucontpts = m;
  fgeom->vcontpts = n;
	
  /* set the knot vector in u direction */
  for (i=0; i<m; i++)
      {
       fgeom->uknots[i]   = xa; 
       fgeom->uknots[m+i] = xb;
      }
  /* set the knot vector in v direction */
  for(i=0; i<n; i++)
      {
       fgeom->vknots[i]   = ya; 
       fgeom->vknots[n+i] = yb;
      }
	
  /* nondimension is wanted or not. */
  if (nondim == TRUE)
     max = absmaxarray(m,n,&w[0][0],size);
  else
     max = 1.0;

  for (i=0; i<m; i++)
      {
       /* calculate x-coordinates of the control points by linear precision
	* property */
       xk = xa + ((double) i/(m-1))*(xb-xa);
       for (j=0; j<n; j++)
	   {
	    /* calculate y-coordinates of the control points by linear
	     * precision property */
	    yk = ya + ((double) j/(n-1))*(yb-ya);
	    if ( (v = fgeom->contpts[i][j]) == NULL)	
	      	v = fgeom->contpts[i][j] = vectalloc();
	    /* assign the coordinates to the control point */
	    v->x = xk;
	    v->y = yk;
	    v->z = w[i][j]/max;   
	    v->w = one;
	   }
      }

  /* free memory */
  free_darray2(bx);
  free_darray2(bxt);
  free_darray2(by);
  free_darray2(l1);
  free_darray2(l2);
  free_darray2(l1t);
  free_darray2(binv);
  free_darray2(temp1);
  free_darray2(w);
  free_darray2(temp2);
  free_darray2(acopy);
  free_darray1(x);
}

/***************************************************************************
*                                  linear_trans()
****************************************************************************
* 
* 1   Purpose
*     This function is a subroutine of the function which converts a polynomial
*     with power basis to Bernstein basis. It calculates a transformation
*     matrix for linearization.
* 
* 2   Specification
*     #include "bspl.h"
*     void linear_trans(int n, double **t, double a, double b)
*     
* 3   Description
*     This function calculates the matrix which transform a monomial with its 
*     variable x in interval [a,b] to another monomial with its variable t in 
*     interval [0,1], where x = a + (b-a) * t.
* 
*                         i!        j        (i-j)
*            t(i,j) = ---------- * a  * (b-a)
*                      j!(i-j)!
* 
* 4   References
*     Not applicable
* 
* 5   Parameters
*       1.int n
*         On entry: the order of the monomial.
*       2.double **t
*         On entry: a 2D array
* 	On exit:  the transformation matrix
*       3.double a
*         On entry: the lower bound of x
*       4.double b
*         On entry: the upper bound of x
* 
* 6   Return Values, Error Indicators and Warnings
*     Not applicable
* 
* 7   Accuracy
*     Not applicable.
* 
* 8   Further Comments
*     Not applicable
* 
* 9   Functions referenced by linear_trans() are:
*     combine()
* 
* 10  Functions that reference linear_trans() are:
*     monotobern()
* 
****************************************************************************/

linear_trans(int n, double **t, double a, double b)
{ 
  int i,j;       /* loop indices */
  int x;         /* working variable */
  double y,z;    /* working variables */

  /* the transformation matrix has dimension nxn */
  for (i=0; i<n; i++)
      for (j=0; j<n; j++)
	  {
	   if (i<j)   /* it is an upper triangluar matrix. */
	      t[i][j] = 0.0;
	   else
	      {
	       x = combine(i,j);  
	       /* calculate the element of the matrix. */
	       if ( a==0.0 && j==0)  /* handle the boundary */
		  y = 1.0;
	       else
		  y = LPOW(a, (double) j);
	       if ( b==0.0 && (i-j) == 0)   /* handle the boundary */
		  z = 1.0;
	       else
		  z = LPOW(b, (double) (i-j));
	       t[i][j] = x*y*z;
	      }
	  }
}




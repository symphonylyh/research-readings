/***************************************************************************
 *                                                                         *
                           Copyright (C) 1992 by
            Massachusetts Institute of Technology, Cambridge, MA
                             All Rights Reserved
 *                                                                         *
 **************************************************************************/
/*  Filename:    partan.c

    Written by:  Peter C. Filkins
    Date:        15 February 1991
===========================================================================
    File Modification History:
===========================================================================
    Description:  Solves the system Ax = b using minimal least squares
                  technique.
===========================================================================
    Subroutines called:    NAG routine FO4JAF
    Library dependencies:  libgen.a, lnag
===========================================================================
    Arguments:  vector *r       corresponding to a known direction
                vector *ru,*rv  first parametric derivatives of surface
===========================================================================
    NAG routine variables:
           int m       number of rows of matrix A
	   int n       number of columns of matrix A
	   double a    an array of (nra,t) dimension.  The matrix A in the
                       system Ax = b.
	   int nra     first dimension of A 
	   double b    array of dimension at least m; contains the vector
	               b of the system; on exit contains the solution
	   double tol  relative tolerance to determine the rank of A.  See
	               NAG manual.
	   double sigma  standard error
	   int rank    rank of matrix A
	   double work dimension (lw) workspace
	   int lw      length of workspace; lw >= 4n
	   int fail    error flag; 0 on entry; 0 on normal exit
=========================================================================*/
#include <stdio.h>
#include <math.h>
#include "gen.h"
#include "editor.h"

#define TOL  1e-05                          /* relative tolerance level */

void f04jaf_(int *, int *, double *, int *, double *, double *, double *,
	     int *, double *, int *, int *);

double *partan_dir(vector *r, vector *ru, vector *rv)
{
  double *dir;
  double **a, *b, tol, *w, err;
  int m, n, nra, rank, lw, fail;

             /* set initial values.  These are hardcoded to reflect this
		specific problem.  A more  general routine would pass these
		dimensions as arguments from the calling function */
  m = nra = 3;
  n = 2;
  lw = 4 * n;
  tol = TOL;
  fail = 0;

  a = dbl_array2(n, m);
  b = dbl_array1(m);
  w = dbl_array1(lw);

     /* initialize A, b matrices */
  a[0][0] = ru->x;
  a[1][0] = rv->x;
  a[0][1] = ru->y;
  a[1][1] = rv->y;
  a[0][2] = ru->z;
  a[1][2] = rv->z;

  b[0] = r->x;
  b[1] = r->y;
  b[2] = r->z;
     
  f04jaf_(&m, &n, &a[0][0], &nra, b, &tol, &err, &rank, w, &lw, &fail);
  if (fail)
    errormsg(1, "failure in f04jaf - called from partan_dir");

  dir = dbl_array1(2);               /* dir[0] = du/dt; dir[1] = dvdt */

  dir[0] = b[0];
  dir[1] = b[1];

  free_darray2(a);
  free_darray1(b);
  free_darray1(w);

  return (dir);
}
/*=========================================================================*/

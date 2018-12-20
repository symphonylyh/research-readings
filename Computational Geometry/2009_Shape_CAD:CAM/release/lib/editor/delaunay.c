/* Copyright (C) 1993 Massachusetts Institute of Technology, Cambridge, MA
   All rights reserved
*/

/* delaunay.c */

/* Delaunay()
   DelaunayCompare()
*/

#include <stdlib.h>
#include "gen.h"
#include "editor.h"

void e01saf_(int *, double *, double *, double *, int *, double *, int *);

short Delaunay(int nPts, double **pts, int *iEnd, int *nAdj,
	       int **iAdj, int sort)
{
  double **grads, *x, *y, *z;
  int ifail = 0, n, *triang;
  int i;
  
  n = nPts;

  if (sort)
    qsort(&pts[0][0], n, 3*sizeof(double), DelaunayCompare);

  x = dbl_array1(n);
  y = dbl_array1(n);
  z = dbl_array1(n);
  for (i=0; i<n; i++) {
    x[i] = pts[i][0];
    y[i] = pts[i][1];
    z[i] = pts[i][2];
  }

  triang = int_array1(n*7);
  grads = dbl_array2(n, 2);

  e01saf_(&n, x, y, z, triang, &grads[0][0], &ifail);
  if (!ifail) {
    for (i=6*n; i<7*n; i++)
      iEnd[i - 6*n] = triang[i] - 1;
    *nAdj = iEnd[nPts-1] + 1;
    (*iAdj) = int_array1(*nAdj);
    for (i=0; i < *nAdj; i++)
      (*iAdj)[i] = triang[i] - 1;
  }

  free_darray1(x);
  free_darray1(y);
  free_darray1(z);
  free_iarray1(triang);
  free_darray2(grads);

  return (ifail);
}

int DelaunayCompare(const void *xx1, const void *xx2)
{
  double *x1, *x2;

  x1 = (double *)xx1;
  x2 = (double *)xx2;
  if (x1[0] < x2[0])         /* primary sort on x */
    return -1;
  else if (x1[0] > x2[0])
    return 1;
  else {                     /* secondary sort on y */
    if (x1[1] < x2[1])
      return -1;
    else if (x1[1] > x2[1])
      return 1;
  }
  return 0;
}





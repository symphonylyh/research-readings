/* Copyright (C) 1992 Massachusetts Institute of Technology, Cambridge, MA
   All rights reserved
*/

/* farray.c */

/* flt_array3()
   free_farray3()
*/

#include <malloc.h>
#include "gen.h"

float ***flt_array3(unsigned nrow, unsigned ncol, unsigned depth)
/* nrow by ncol by depth array of doubles */
{
  float ***t;
  short i, j;
  char *c1, *c2;

  t = (float ***)malloc(nrow*sizeof(float **));
  if (!t)
    errormsg(14, "allocation failure 1 in Flt_Array3()");
  c1 = (char *)malloc(nrow*ncol*sizeof(float *));
  c2 = (char *)calloc(nrow*ncol*depth, sizeof(float));
  if (!c1)
    errormsg(15, "allocation failure 2 in Flt_Array3()");
  if (!c2)
    errormsg(16, "allocation failure 3 in Flt_Array3()");
  for (i = 0; i < nrow; i++) {  /* initialize the row pointers */
    t[i] = (float **)(c1 + i*ncol*sizeof(float *));
    for (j = 0; j < ncol; j++)  /* initialize the row and col pointers */
      t[i][j] = (float *)(c2 + (i*ncol + j)*depth*sizeof(float));
  }
  return (t);
}

void free_farray3(float ***t)
{
  free(&t[0][0][0]);
  free(&t[0][0]);
  free(t);
}

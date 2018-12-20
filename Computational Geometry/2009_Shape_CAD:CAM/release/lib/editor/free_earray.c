/* Copyright (C) 1992 Massachusetts Institute of Technology, Cambridge, MA
   All rights reserved
*/

/* free_earray.c */

/* free_earray1()
   free_earray2()
*/

#include <malloc.h>
#include "gen.h"

void free_earray1(ParCurv **egm, unsigned nel)
{
  short i;

  for (i=0; i<nel; i++)
    free_egeom(egm[i]);
  free(egm);
}

void free_earray2(ParCurv ***egm, short nrow, short ncol)
{
  short i, j;

  for (i=nrow-1; i>=0; i--)
    for (j=ncol-1; j>=0; j--)
    if (egm[i][j])
      free_egeom(egm[i][j]);
  free(egm[0]);
  free(egm);
}

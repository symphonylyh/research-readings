/* Copyright (C) 1993 Massachusetts Institute of Technology, Cambridge, MA
   All rights reserved
*/

/* efg_array.c */

/* fgeom_array2()
   free_egeom_array1()
   free_fgeom_array2()
*/

#include "gen.h"

void free_egeom_array1(ParCurv **egeoms, int np)
{
  short i;

  for (i=0; i<np; i++)
    free_egeom(egeoms[i]);

  free_garray1((char *)egeoms);
}

ParSurf ***fgeom_array2(int nps, int npt)
{
  ParSurf ***fgeoms;

  fgeoms = (ParSurf ***)gen_array2(npt, nps, sizeof(ParSurf *));

  return fgeoms;
}

void free_fgeom_array2(ParSurf ***fgeoms, int nps, int npt)
{
  short i, j;

  for (i=0; i<nps; i++)
    for (j=0; j<npt; j++)
      free_fgeom(fgeoms[j][i]);

  free_garray2((char **)fgeoms);
}

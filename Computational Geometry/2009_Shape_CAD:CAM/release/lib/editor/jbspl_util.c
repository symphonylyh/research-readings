/* Copyright (C) 1993 Massachusetts Institute of Technology, Cambridge, MA
   All rights reserved
*/

/* jbspl_util.c */

/* arc_length()
   avg_arc_length()
   free_egeom_array1()
   free_fgeom_array2()
   swapsurf()
   vect_colinear()
*/

#include <math.h>
#include "gen.h"
#include "bspl.h"
#include "appr.h"
#include "editor.h"

double arc_length(ParCurv *egeom)
{
  vector *v1, *v2, *v3;
  double s = 0.0, dt = 0.01;
  short i;

  v3 = vectalloc();
  for (i=0; i<100; i++) {
    v1 = evalbsp(egeom, dt*i);
    v2 = evalbsp(egeom, dt*(i+1));

    sub_vect1(v2, v1, v3);
    s += mag(v3);

    vectfree(v1);
    vectfree(v2);
  }
  vectfree(v3);

  return s;
}

double avg_arc_length(ParSurf *fgeom)
{
  ParCurv *egeom;
  double *nodes, s = 0.0;
  short i;
  
  nodes = dbl_array1(fgeom->vcontpts);
  nodes_bsp(fgeom->vknots, fgeom->vcontpts, fgeom->vorder, nodes);

  for (i=0; i<fgeom->vcontpts; i++) {
    egeom = ParCurv_iso(fgeom, nodes[i], 0, NULL);
    s += arc_length(egeom);
    free_egeom(egeom);
  }

  return s/fgeom->vcontpts;
}

ParSurf *swapsurf(ParSurf *fgeom)
{
  ParSurf *newsurf;
  short i, j;
  
  newsurf = fgeomalloc1(fgeom->vorder, fgeom->uorder, 
			fgeom->vcontpts, fgeom->ucontpts);

  for (i=0; i<fgeom->vorder+fgeom->vcontpts; i++)
    newsurf->uknots[i] = fgeom->vknots[i];
  for (i=0; i<fgeom->uorder+fgeom->ucontpts; i++)
    newsurf->vknots[i] = fgeom->uknots[i];

  for (i=0; i<fgeom->vcontpts; i++)
    for (j=0; j<fgeom->ucontpts; j++)
      copyvector(fgeom->contpts[j][i], newsurf->contpts[i][j]);

  return newsurf;
}

int vect_colinear(vector *v1, vector *v2)
{
  double a, b, c;

  a = dot(v1, v2);
  b = mag(v1);
  c = mag(v2);

  if (fabs(a/(b*c) + 1) < EDITOR_JBSPL_TOL)
    return 1;

  return 0;
}

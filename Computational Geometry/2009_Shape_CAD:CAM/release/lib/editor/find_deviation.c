/************************************************************************
 *									*
			Copyright (C) 1992 by
	Massachusetts Institute of Technology, Cambridge, MA
			 All rights reserved
 *									*
 ************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include "gen.h"
#include "bspl.h"
#include "editor.h"

void find_deviation(ParCurv *egeom1, ParCurv *egeom2, double *maxd)
{
  vector *v1, *v2;
  double dist;
  short i;
  
  *maxd = 0.0;
  
  for (i=0; i<egeom1->ncontpts; i++) {
    v1 = egeom1->contpts[i];
    v2 = egeom2->contpts[i];
    dist = distance(v1, v2);
    *maxd = MAX(*maxd, dist);
  }
}

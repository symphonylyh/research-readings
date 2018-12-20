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
#include <malloc.h>
#include "gen.h"
#include "bspl.h"
#include "editor.h"

void fair_per_knot(ParCurv *egeom, short i)
{
  vector *v1, *v2, *v3, *v4;
  double2 a, b, c, d, e, f, g, h, q;
  double2 u0, u1, u2, u3, d1, d2, d3;
  short j, n;
  
  v1 = vectalloc();
  v2 = vectalloc();
  v3 = vectalloc();
  v4 = vectalloc();
  
  n = egeom->ncontpts;
  j = i + n;
  
  u0 = egeom->knots[i];
  u1 = egeom->knots[i+1];
  u2 = egeom->knots[i+2];
  u3 = egeom->knots[i+3];

  if (i == n) {
    u1 = egeom->knots[1] + egeom->knots[n];
    u2 = egeom->knots[2] + egeom->knots[n];
    u3 = egeom->knots[3] + egeom->knots[n];
  }

  if (i == n-1) {
    u2 = egeom->knots[1] + egeom->knots[n];
    u3 = egeom->knots[2] + egeom->knots[n];
  }

  if (i == n-2)
    u3 = egeom->knots[1] + egeom->knots[n];

  d1 = egeom->knots[i-1];
  d2 = egeom->knots[i-2];
  d3 = egeom->knots[i-3];

  if (i == 0) {
    d1 = egeom->knots[n-1] - egeom->knots[n];
    d2 = egeom->knots[n-2] - egeom->knots[n];
    d3 = egeom->knots[n-3] - egeom->knots[n];
  }
  
  if (i == 1) {
    d2 = egeom->knots[n-1] - egeom->knots[n];
    d3 = egeom->knots[n-2] - egeom->knots[n];
  }
  
  if (i == 2)
    d3 = egeom->knots[n-1] - egeom->knots[n];

  a = u1 - d3;
  b = u0 - d3;
  c = u1 - u0;
  d = u3 - d1;
  e = u3 - u0;
  f = u0 - d1;
  g = u2 - u0;
  h = u2 - d2;
  q = u0 - d2;
  
  copyvector(egeom->contpts[(j-4)%n], v1);
  copyvector(egeom->contpts[(j-3)%n], v2);
  copyvector(egeom->contpts[(j-1)%n], v3);
  copyvector(egeom->contpts[j%n],     v4);
  
  scale_vect1(c/b, v1, v1);
  scale_vect1(a/b, v2, v2);
  scale_vect1(d/e, v3, v3);
  scale_vect1(f/e, v4, v4);
  
  sub_vect1(v2, v1, v1);
  sub_vect1(v3, v4, v2);
  
  scale_vect1(g/h, v1, v1);
  scale_vect1(q/h, v2, v2);
  
  add_vect1(v1, v2, egeom->contpts[(j-2)%n]);
  
  vectfree(v1);
  vectfree(v2);
  vectfree(v3);
  vectfree(v4);
}

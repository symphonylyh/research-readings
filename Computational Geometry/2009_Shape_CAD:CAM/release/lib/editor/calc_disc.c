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

#define EPS 1.0e-8

/* function for calculating local and global k' discontinuity of a cubic 
   B-spline curve.
   input: egeom, pointer to a structure that stores the B-spline geometry
   output: zglob, the sum of the local discontinuities    */

void calc_disc(ParCurv *egeom, double *zglob,int *ik)
{
  vector *v1, *v2, *v3;
  double g, k0, max, p, q, r, s, t;
  double k1[2], u[2], *zloc;
  short i, j, low_end;
  
  *zglob = 0.0;
  zloc = dbl_array1(egeom->ncontpts);
  
  if (egeom->type == PCurvePer)
    low_end = 1;
  else
    low_end = egeom->order;

  for (i=low_end; i<egeom->ncontpts; i++) {
    u[0] = egeom->knots[i] - EPS;
    u[1] = egeom->knots[i] + EPS;
      
    for (j=0; j<2; j++) {

	  /*  evalrbsp has been set to correct function,  */
	  /*  periodic or non-periodic.                   */

      if (egeom->type == PCurvePer) {
	v1 = rbspeval_per(egeom, u[j], 1);
	v2 = rbspeval_per(egeom, u[j], 2);
	v3 = rbspeval_per(egeom, u[j], 3);
      }
      else {
	v1 = rbspeval(egeom, u[j], 1);
	v2 = rbspeval(egeom, u[j], 2);
	v3 = rbspeval(egeom, u[j], 3);
      }
	  
      s = pow(v1->x, 2.0) + pow(v1->y, 2.0) + pow(v1->z, 2.0);
      p = sqrt(pow(((v1->x)*(v2->y)-(v1->y)*(v2->x)), 2.0) +
	       pow(((v1->y)*(v2->z)-(v1->z)*(v2->y)), 2.0) +
	       pow(((v1->z)*(v2->x)-(v1->x)*(v2->z)), 2.0));
	  
      r = ((v1->x)*(v2->y) - (v1->y)*(v2->x)) *
	  ((v1->x)*(v3->y) - (v1->y)*(v3->x)) +
	  ((v1->y)*(v2->z) - (v1->z)*(v2->y)) *
	  ((v1->y)*(v3->z) - (v1->z)*(v3->y)) +
	  ((v1->z)*(v2->x) - (v1->x)*(v2->z)) *
	  ((v1->z)*(v3->x) - (v1->x)*(v3->z));
	  
      t = (v1->x)*(v2->x) + (v1->y)*(v2->y) + (v1->z)*(v2->z);
	  
      g = sqrt(s);
      q = pow(g, 3.0);
      k0 = p/q;
      k1[j] = (r/(q*p) - (3.0*k0*t)/s)/g;

      vectfree(v1);
      vectfree(v2);
      vectfree(v3);
    }
    zloc[i] = fabs(k1[0] - k1[1]);
    *zglob = *zglob + zloc[i];
  } 
  
  *ik = low_end;
  max = zloc[low_end];
  for (i=low_end; i<egeom->ncontpts; i++ ) {
    if (max <= zloc[i]) {
      *ik = i;
      max = zloc[i];
    }
  }

  free_darray1(zloc);
}

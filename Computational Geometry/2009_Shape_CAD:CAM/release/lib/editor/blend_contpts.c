/***************************************************************************
 *                                                                         *
                           Copyright (C) 1992 by
            Massachusetts Institute of Technology, Cambridge, MA
                             All Rights Reserved
 *                                                                         *
 **************************************************************************/
/*  Filename:    blend_contpts.c

    Written by:  Peter C. Filkins
    Date:        23 January 1991
===========================================================================
    File Modification History:
         25 January 1991 - Made the functions more general (to allow for
	                   more than 6 control points)
===========================================================================
    Description:  computes the control points P(2) and P(m-3) (m vertices
                  in the control polygon) based on second derivative
		  and normal curvature information for use in developing 
		  the blend crosslink curve
===========================================================================
    Subroutines called: sub_vect(), scale_vect1(), add_vect(), add_vect1()
    Library dependencies:  libgen.a
===========================================================================
    Arguments:    ParCurv *egm         cross link curve
                  vector *d2[0]        second derivative at cross link 
		                       curve parameter t = 0
		  vector *d2[1]        second derivative at t = 1
===========================================================================*/

#include <stdio.h>
#include <malloc.h>
#include "gen.h"
#include "editor.h"

void curv_contpts0(ParCurv *egm, vector *d2)
{
  double *s0;

  s0 = knot_factor0(egm->order, egm->knots);
  egm->contpts[2] = eval_vert(egm->contpts[0], egm->contpts[1], s0, d2);

  return;
}

void curv_contpts1(ParCurv *egm, vector *d2, short m)
{
  double *s1;

  s1 = knot_factor1(egm->order, egm->ncontpts, egm->knots);
  egm->contpts[m-3] = eval_vert(egm->contpts[m-1], egm->contpts[m-2], s1, d2);

  return;
}

double *knot_factor0(short order, double *s)
/* short order  order of cross link curve */
/* double *s    knot vector */
{
  double *factor;

  factor = dbl_array1(3);

  factor[0] = (order - 1) / (s[order] - s[1]);
  factor[1] = (s[order]+s[order+1] - s[1]-s[2]) / (order-1);
  factor[2] = (s[order]-s[2]) * (s[order+1]-s[2]) / ((order-1)*(order-2));
  
  return (factor);
}

double *knot_factor1(short order, short m, double *s)

/* short order, m  order and numcontpts */
/* double *s       knot vector */
{
  double *f;

  f = dbl_array1(3);

  f[0] = (order-1) / (s[m+order-2]-s[m-1]);
  f[1] = (s[m+order-2]+s[m+order-3]-s[m-1]-s[m-2]) / (order-1);
  f[2] = (s[m+order-3]-s[m-1])*(s[m+order-3]-s[m-2])/((order-1)*(order-2));

  return (f);
}

vector *eval_vert(vector *p0, vector *p1, double *s, vector *d2)
/* vector *p0  control point P(0) or P(m-1) */
/* vector *p1  control point P(1) or P(m-2) */
/* double *s   array of knot products */
/* vector *d2  second derivative */
{
  vector *d1, *sum;

  d1 = sub_vect(p1, p0);
  scale_vect1(s[0], d1, d1);
  scale_vect1(s[1], d1, d1);
  scale_vect1(s[2], d2, d2);

  sum = add_vect(d1, d2);
  add_vect1(p0, sum, sum);

  vectfree(d1);
     
  return (sum);
}

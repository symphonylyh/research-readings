/************************************************************************
 *									*
			Copyright (C) 1992 by
	Massachusetts Institute of Technology, Cambridge, MA
			 All rights reserved
 *									*
 ************************************************************************/

#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include "gen.h"
#include "bspl.h"
#include "appr.h"
#include "editor.h"

static short m;

ParCurv *cross_link(double u, ParSurf **fgm, ParCurv **egm, short bc1,
		    short bc2, double surf1, double surf2)
{
  ParCurv *clink;
  vector *v0[6], *v1[6], *n0, *n1, *p;
  vector *tanv0[2], *tanv1[2];
  vector *qtan0, *qtan1, *sec_deriv0, *sec_deriv1;
  vector *d_bi;
  double b, dist;
  short j;

  clink = link_alloc(bc1, bc2, &m);

  n0 = surf_data(egm[0], fgm[0], u, v0);   /*  Surface data, normals  */
  n1 = surf_data(egm[1], fgm[1], u, v1);

  p = sub_vect(v1[0], v0[0]);
  unitvector1(p, p);
  d_bi = rbspeval(egm[2], u, 0);
  sethomogeq(d_bi, 1.0/d_bi->w);

  tanv0[0] = sub_vect(d_bi, v0[0]);
  tanv0[1] = cross(p, tanv0[0]);
  unitvector1(tanv0[1], tanv0[1]);
  cross1(tanv0[1], n0, tanv0[0]);
  unitvector1(tanv0[0], tanv0[0]);

  tanv1[0] = sub_vect(d_bi, v1[0]);
  tanv1[1] = cross(p, tanv1[0]);
  unitvector1(tanv1[1], tanv1[1]);
  cross1(tanv1[1], n1, tanv1[0]);
  unitvector1(tanv1[0], tanv1[0]);

  vectfree(p);
  dist = distance(v1[0], v0[0]); 
                            
  if (bc1 > 0) {      /*  scaled tan vec for contpt,curvature  */
    b = bias_factor1(0, dist, surf1, surf2);
    qtan0 = bias(tanv0[0], v0[0], d_bi, b);
    if (bc1 == 2)
      sec_deriv0 = curv_cont(v0, n0, tanv0, qtan0);
  }
  if (bc2 > 0) {
    b = bias_factor1(1, dist, surf1, surf2);
    qtan1 = bias(tanv1[0], v1[0], d_bi, b);
    if (bc2 == 2) 
      sec_deriv1 = curv_cont(v1, n1, tanv1, qtan1);
  }
  vectfree(tanv0[0]);
  vectfree(tanv0[1]);
  vectfree(tanv1[0]);
  vectfree(tanv1[1]);
  vectfree(n0);
  vectfree(n1);
  vectfree(d_bi);

/*  control points  */

  copyvector(v0[0], clink->contpts[0]);    /*  Position  */
  copyvector(v1[0], clink->contpts[m-1]);

  if (bc1 > 0) {                           /*  Tangent  */
    clink->contpts[1] = add_vect(v0[0], qtan0);
    vectfree(qtan0);
  }
  if (bc2 > 0) {
    clink->contpts[m-2] = add_vect(v1[0], qtan1);
    vectfree(qtan1);
  }

  if (bc1 == 2) {                          /*  Curvature  */
    curv_contpts0(clink, sec_deriv0);
    vectfree(sec_deriv0);
  }
  if (bc2 == 2) {
    curv_contpts1(clink, sec_deriv1, m);
    vectfree(sec_deriv1);
  }
  
  free_varray1(v0,6);
  free_varray1(v1,6);

  return(clink);
}

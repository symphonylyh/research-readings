/************************************************************************
 *									*
			Copyright (C) 1992 by
	Massachusetts Institute of Technology, Cambridge, MA
			 All rights reserved
 *									*
 ************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <malloc.h>
#include "gen.h"
#include "bspl.h"
#include "editor.h"

#define MACHPREC 1.1102230246251567E-13	  /* changed on 25 JUN 90 */

static int trans_key = 0;

ParCurv *cyl_generatrix(ParCurv *directrix, double t, double tz, double theta)

/* This routine returns the generatrix at the parameter value t of the
   spine. The spine is the line directrix.
   The generatrix is trimmed in position using the Frenet's
   trihedron at t, so that the local y axis points along the normal and the
   local x axis along the binormal. The routine generatrix_lk is called,
   which returns the generatrix at point t, defined as a rational B-spline
   in the local coordinate system. */

{
  ParCurv *eg;
  vector *v1,**v0;
  double **tt,u;
  short i;

/* Meaning of the most important local variables (in alphabetical order)

eg       =  the generatrix expressed as a rational B-spline in its local
            system. eg is transformed to the global system and returned by this
            routine.
tt       =  transformation matrix (4 x 4). The control polygon
            coordinates of the generatrix at its final position are obtained by
            multiplying [t] with [x y 0 w]T, where (x,y,w) are the homogeneous
            coordinates in the local system.
v0       =  the Frenet trihedron.

*/

  v1 = rbspeval(directrix, t, 0);
  v0 = cyl_frenet_tr(directrix, t, theta);
  eg = cyl_generatrix_lk(t, tz);
  tt = dbl_array2(4, 4);
  for (i=0; i<3; i++) {
    tt[0][i] = v0[i]->x;
    tt[1][i] = v0[i]->y;
    tt[2][i] = v0[i]->z;
  }
  u = v1->w;
  tt[0][3] = v1->x/u;
  tt[1][3] = v1->y/u;
  tt[2][3] = v1->z/u;
  tt[3][0] = tt[3][1] = tt[3][2] = 0.0;
  tt[3][3] = 1.0 ;
  free_varray1(v0, 3);
  vectfree(v1);
  for (i=0; i<eg->ncontpts; i++)
    mult4x4(tt, eg->contpts[i], eg->contpts[i]);
  free_darray2(tt);
  
  return(eg);
}

ParCurv *cyl_generatrix_lk(double t, double sg)

/* this routine returns the generatrix at the parameter value t of the
   spine, expressed as a rational B-spline in its local system. In its
   current form the routine returns a circle.
   sg is the radius of the circle. */

{
  ParCurv *eg;
  double d,f;
  short i;
  
  eg = egeomalloc1(3, 9);
  for (i=0; i<3; i++) {
    eg->knots[i] = 0.0;
    eg->knots[i+9] = 1.0;
  }
  eg->knots[3] = eg->knots[4] = 0.25;
  eg->knots[5] = eg->knots[6] = 0.5;
  eg->knots[7] = eg->knots[8] = 0.75;
  d = 1.0/sqrt((double)2.0);
  f = d*sg;
  for (i=0; i<9; i++)
    eg->contpts[i]->z = 0.0;
  eg->contpts[0]->x = sg;
  eg->contpts[0]->y = 0.0;
  eg->contpts[0]->w = 1.0;
  eg->contpts[1]->x = eg->contpts[1]->y = f ; eg->contpts[1]->w = d ;
  eg->contpts[2]->x = 0.0;
  eg->contpts[2]->y = sg;
  eg->contpts[2]->w = 1.0;
  eg->contpts[3]->x = -f;
  eg->contpts[3]->y = f;
  eg->contpts[3]->w = d;
  eg->contpts[4]->x = -sg;
  eg->contpts[4]->y = 0.0;
  eg->contpts[4]->w = 1.0;
  eg->contpts[5]->x = -f;
  eg->contpts[5]->y = -f;
  eg->contpts[5]->w = d;
  eg->contpts[6]->x = 0.0;
  eg->contpts[6]->y = -sg;
  eg->contpts[6]->w = 1.0;
  eg->contpts[7]->x = f;
  eg->contpts[7]->y = -f;
  eg->contpts[7]->w = d;
  eg->contpts[8]->x = sg;
  eg->contpts[8]->y = 0.0;
  eg->contpts[8]->w = 1.0 ;
  
  return(eg);
}

void rotate_section(vector **v0, double theta)
{
  double **tt,**tr,th;
  short i,j;

  tt = dbl_array2(4, 4);      /* transform global to local */
  for (i=0; i<3; i++) {
    tt[i][0] = v0[i]->x;
    tt[i][1] = v0[i]->y;
    tt[i][2] = v0[i]->z;
  }
  tt[3][3] = 1.0;
  mult4x4(tt,v0[0],v0[0]);
  mult4x4(tt,v0[1],v0[1]);

  tr = dbl_array2(4,4);      /* rotate by theta degrees */
  th = theta*DEG_TO_RAD;
  tr[0][0] = cos(th);
  tr[0][1] = sin(th);
  tr[1][0] = -sin(th);
  tr[1][1] = cos(th);
  tr[2][2] = 1.0;
  tr[3][3] = 1.0;
  mult4x4(tr,v0[0],v0[0]);
  mult4x4(tr,v0[1],v0[1]);
  free_darray2(tr);

  for (j=0; j<3; j++)        /* transform local to global */
    for (i=j+1; i<4; i++) {  /* transpose of tt */
      th       = tt[j][i];
      tt[j][i] = tt[i][j];
      tt[i][j] = th;
    }
  mult4x4(tt,v0[0],v0[0]);
  mult4x4(tt,v0[1],v0[1]);
  free_darray2(tt);
}

vector **cyl_frenet_tr(ParCurv *eg, double t, double theta)

/* Compute the components of the Frenet trihedron at the parameter t
   Output is th Frenet trihedron in the order (X(s),Y(s),Z(s)). */

{
  static vector t0,n0,b0;
  static double d0;
  vector **v0,*v1,*v2,tu;
  double **rt,d,du;
  short i;

  v0 = vec_array1(3);

  v1 = rbspeval(eg, t, 1);
  d = mag(v1);
  if (fabs(d) < MACHPREC) {
    printf("cyl_frenet_tr: First derivative = 0.0 at t = %f\n", t);
    exit(1);
  }
  v0[2] = unitvector(v1);      /* tangent vector */

  v2 = rbspeval(eg, t, 2);
  d = dot(v1,v2)/(d*d);
  scale_vect1(d, v1, v1);
  sub_vect1(v2, v1, v2);
  v0[0] = unitvector(v2);       /* normal vector */
  v0[1] = cross(v0[2], v0[0]);  /* binormal vector */

  vectfree(v2);
  vectfree(v1);

  if (trans_key) {
    if (t == 0.0) {
      copyvector(v0[0], &n0);     /* save normal */
      copyvector(v0[1], &b0);     /* save binormal */
      copyvector(v0[2], &t0);     /* save tangent */
      d0 = sqrt(t0.x*t0.x + t0.z*t0.z);
      if (d0 == 0.0)
	d0 = sqrt(1.0 - t0.y*t0.y);
    }
    else {
      copyvector(v0[2], &tu);
      du = sqrt(tu.x*tu.x + tu.z*tu.z);
      if (du == 0.0)
	du = sqrt(1.0 - tu.y*tu.y);

      rt = dbl_array2(4,4);

      rt[0][0] =  tu.z*t0.z/(du*d0) + tu.y*tu.x*t0.y*t0.x/(du*d0) + tu.x*t0.x;
      rt[0][1] = -tu.y*tu.x*d0/du + tu.x*t0.y;
      rt[0][2] = -tu.z*t0.x/(du*d0) + tu.y*tu.x*t0.y*t0.z/(du*d0) + tu.x*t0.z;
      rt[1][0] = -du*t0.y*t0.x/d0 + tu.y*t0.x;
      rt[1][1] =  du*d0 + tu.y*t0.y;
      rt[1][2] = -du*t0.y*t0.z/d0 + tu.y*t0.z;
      rt[2][0] = -tu.x*t0.z/(du*d0) + tu.y*tu.z*t0.y*t0.x/(du*d0) + tu.z*t0.x;
      rt[2][1] = -tu.y*tu.z*d0/du + tu.z*t0.y;
      rt[2][2] =  tu.x*t0.x/(du*d0) + tu.y*tu.z*t0.y*t0.z/(du*d0) + tu.z*t0.z;

      for (i=0; i<3; i++) {
	rt[i][3] = 0.0;
	rt[3][i] = 0.0;
      }
      rt[3][3] = 1.0;

      mult4x4(rt, &n0, v0[0]);
      v0[0] = unitvector(v0[0]);

      v0[1] = cross(v0[2], v0[0]);
      v0[1] = unitvector(v0[1]);

      if (dot(&t0, &tu) < 0.766) {
	copyvector(v0[0], &n0);
	copyvector(v0[1], &b0);
	copyvector(v0[2], &t0);
	d0 = sqrt(t0.x*t0.x + t0.z*t0.z);
	if (d0 == 0.0)
	  d0 = sqrt(1.0 - t0.y*t0.y);
      }
    }
  }

  if (theta != 0.0)
    rotate_section(v0, theta);

  return(v0);
}

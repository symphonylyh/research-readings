/************************************************************************
 *									*
			Copyright (C) 1992 by
	Massachusetts Institute of Technology, Cambridge, MA
			 All rights reserved
 *									*
 ************************************************************************/

#include <math.h>
#include "gen.h"
#include "editor.h"

double **compute_transfm(double *t)

/* Compute the homogeneous coordinate transformation matrix (4 x 4) from
   the components of the translation vector and the Euler angles of
   rotation
   t[0] = x-coordinate of the translation vector.
   t[1] = y-coordinate of the translation vector.
   t[2] = z-coordinate of the translation vector.
   t[3] = phi (angle of rotation with respect to Oz)
   t[4] = theta (angle of rotation with respect to Oy)
   t[5] = psi (angle of rotation with respect to Ox)
   The routine returns the coordinate transformation matrix. */
{
  double **c,cf,ct,cp,sf,st,sp ;
  short i,j,k,m,n;

/* Meaning of the most important local variables (in alphabetical order)
   c  = homogeneous coordinate transformation matrix.
   cf = cos(phi)
   cp = cos(psi)
   ct = cos(theta)
   sf = sin(phi)
   sp = sin(psi)
   st = sin(theta) */

  c = dbl_array2(4, 4);

  cf = cos(t[3]);
  sf = sin(t[3]);
  ct = cos(t[4]);
  st = sin(t[4]);
  cp = cos(t[5]);
  sp = sin(t[5]);

  c[0][0] = ct*cf;
  c[0][1] = -ct*sf;
  c[0][2] = st;
  c[0][3] = t[0];

  c[1][0] = st*sp*cf+cp*sf;
  c[1][1] = cp*cf-st*sf*sp;
  c[1][2] = -sp*ct;
  c[1][3] = t[1];

  c[2][0] = sp*sf-st*cf*cp;
  c[2][1] = st*sf*cp+sp*cf;
  c[2][2] = ct*cp;
  c[2][3] = t[2];

  c[3][0] = c[3][1] = c[3][2] = 0.0;
  c[3][3] = 1.0;

  return(c);
}

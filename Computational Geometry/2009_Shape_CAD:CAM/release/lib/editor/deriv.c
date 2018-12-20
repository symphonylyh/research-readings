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

double **partial(short nd, double *t)

/* Compute the matrix (4 x 4) of partial derivatives from
   the components of the translation vector and the Euler angles of
   rotation
   t[0] = x-coordinate of the translation vector.
   t[1] = y-coordinate of the translation vector.
   t[2] = z-coordinate of the translation vector.
   t[3] = phi (angle of rotation with respect to Oz)
   t[4] = theta (angle of rotation with respect to Oy)
   t[5] = psi (angle of rotation with respect to Ox)
   The routine returns the partial derivative matrix.
*/
{
  double **c,cf,ct,cp,sf,st,sp;

/* Meaning of the most important local variables (in alphabetical order)
   c        =  homogeneous coordinate transformation matrix.
   cf       =  cos(phi)
   cp       =  cos(psi)
   ct       =  cos(theta)
   sf       =  sin(phi)
   sp       =  sin(psi)
   st       =  sin(theta)
*/

  c = dbl_array2(4,4) ; cf = cos(t[3]) ; sf = sin(t[3]) ; 
  ct = cos(t[4]) ; st = sin(t[4]) ; cp = cos(t[5]) ; sp = sin(t[5]) ;

  if (nd == 3) {
    c[0][0] = -ct*sf; c[0][1] = -ct*cf; c[0][2] = c[0][3] = 0.0;
    c[1][0] = cp*cf-st*sp*sf; c[1][1] = -cp*sf-st*cf*sp;
    c[1][2] = c[1][3] = 0.0;
    c[2][0] = sp*cf+st*sf*cp; c[2][1] = -sp*sf+st*cf*cp;
    c[2][2] = c[2][3] = 0.0;
    c[3][3] = 1.0; c[3][0] = c[3][1] = c[3][2] = 0.0;
  }
  else if (nd == 4) {
    c[0][0] = -st*cf; c[0][1] = st*sf; c[0][2] = ct; c[0][3] = 0.0;
    c[1][0] = ct*sp*cf; c[1][1] = -ct*sf*sp;
    c[1][2] = sp*st; c[1][3] = 0.0;
    c[2][0] = -ct*cf*cp; c[2][1] = ct*sf*cp;
    c[2][2] = -st*cp; c[2][3] = 0.0;
    c[3][3] = 1.0; c[3][0] = c[3][1] = c[3][2] = 0.0;
  }
  else if (nd == 5) {
    c[0][0] = c[0][1] = c[0][2] = c[0][3] = 0.0;
    c[1][0] = -sp*sf+st*cp*cf; c[1][1] = -sp*cf-st*sf*cp;
    c[1][2] = -cp*ct; c[1][3] = 0.0;
    c[2][0] = cp*sf+st*cf*sp; c[2][1] = cp*cf-st*sf*sp; 
    c[2][2] = -ct*sp; c[2][3] = 0.0;
    c[3][3] = 1.0; c[3][0] = c[3][1] = c[3][2] = 0.0;
  }

  return (c);
}

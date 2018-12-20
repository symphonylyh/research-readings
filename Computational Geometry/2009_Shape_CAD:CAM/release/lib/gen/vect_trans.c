/************************************************************************
 *									*
			Copyright (C) 1991 by
	Massachusetts Institute of Technology, Cambridge, MA
			 All rights reserved
 *									*
 ************************************************************************/

/* Some simple rotations and translations of vector entities. Written
   by Seamus Tuohy. */

#include <math.h>
#include "gen.h"

#define M_PI		3.14159265358979323846

/* Function: rotx1()
 * Purpose: Rotate a vector
 * Method: The input vector is rotated about the x-axis
 *         The input vector remains unchanged.
 * Arguments:
 *  angle - the angle of rotation in degrees
 *  in - the address of the input vector
 *  out - the addresss of the rotated vector
 */

void rotx1(double2 angle, vector *in, vector *out)

{
  double2 y,z,rad,sine,cosine;

  rad = (angle * M_PI) / 180;           /* convert degrees to radians */
  sine = sin(rad), cosine = cos(rad);

  y = in->y*cosine - in->z*sine, z = in->y*sine + in->z*cosine;

  out->x = in->x;
  out->y = y;
  out->z = z;
  out->w = in->w;
}

/* Function: roty1()
 * Purpose: Rotate a vector
 * Method: the address of the vector is rotated about the y-axis
 *         the address of the vector remains unchanged.
 * Arguments:
 *  angle - the angle of rotation in degrees
 *  in - the address of the input vector
 *  out - the addresss of the rotated vector
 */

void roty1(double2 angle, vector *in, vector *out)

{
  double2 x,z,rad,sine,cosine;

  rad = (angle * M_PI) / 180;           /* convert degrees to radians */
  sine = sin(rad), cosine = cos(rad);

  x = in->x*cosine + in->z*sine, z = -in->x*sine + in->z*cosine;

  out->x = x;
  out->y = in->y;
  out->z = z;
  out->w = in->w;
}

/* Function: rotz1()
 * Purpose: Rotate a vector
 * Method: the address of the vector is rotated about the z-axis.
 *         the address of the vector remains unchanged.
 * Arguments:
 *  angle - the angle of rotation in degrees
 *  in - the address of the input vector
 *  out - the addresss of the rotated vector
 */

void rotz1(double2 angle, vector *in, vector *out)

{
  double2 x,y,rad,sine,cosine;

  rad = (angle * M_PI) / 180;           /* convert degrees to radians */
  sine = sin(rad), cosine = cos(rad);

  x = in->x*cosine - in->y*sine, y = in->x*sine + in->y*cosine;

  out->x = x;
  out->y = y;
  out->z = in->z;
  out->w = in->w;
}

/* Function: translate1()
 * Purpose: Translate a vector
 * Method: the address of the vector remains unchanged.
 * Arguments:
 *  tx - the translation in the x direction
 *  ty - the translation in the y direction
 *  tz - the translation in the y direction
 *  in - the address of the input vector
 *  out - the addresss of the rotated vector
 */

void translate1(double2 tx, double2 ty, double2 tz, vector *in, vector *out)

{
  out->x = in->x + in->w*tx;
  out->y = in->y + in->w*ty;
  out->z = in->z + in->w*tz;
  out->w = in->w;
}

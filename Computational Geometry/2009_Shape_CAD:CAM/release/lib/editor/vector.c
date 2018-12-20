/* Copyright (C) 1992 Massachusetts Institute of Technology, Cambridge, MA
   All rights reserved
*/

/* vector.c */

/* Add()
   Copy()
   Cross()
   Dot()
   IsCollinear2()
   Normalize()
   Rotate()
   Scale()
   Sub()
*/

#include <math.h>
#include "editor.h"

void Add(short n, float *a, float *b, float *c)
{
  short i;

  for (i=0; i<n; i++)
    c[i] = a[i] + b[i];
}

void Copy(short n, float *a, float *b)
{
  short i;

  for (i=0; i<n; i++)
    b[i] = a[i];
}

void Cross(float *a, float *b, float *c)
{
  c[0] = a[1]*b[2] - a[2]*b[1];
  c[1] = a[2]*b[0] - a[0]*b[2];
  c[2] = a[0]*b[1] - a[1]*b[0];
}

float Dot(float *a, float *b)
{
  return (a[0]*b[0] + a[1]*b[1] + a[2]*b[2]);
}

int IsCollinear2(float x0, float y0, float x1, float y1, float x2, float y2)
{
  float cth, dx1, dy1, dx2, dy2, len;

  dx1 = x1 - x0;
  dy1 = y1 - y0;
  len = sqrt(dx1*dx1 + dy1*dy1);
  if (len != 0.0) {
    dx1 /= len;
    dy1 /= len;
  }
  dx2 = x2 - x0;
  dy2 = y2 - y0;
  len = sqrt(dx2*dx2 + dy2*dy2);
  if (len != 0.0) {
    dx2 /= len;
    dy2 /= len;
  }

  cth = dx1*dx2 + dy1*dy2;
  if (fabs(cth - 1.0) < 1.0e-10)
    return 1;

  return 0;
}

void Normalize(short n, float *v)
{
  float len = 0.0;
  short i;

  for (i=0; i<n; i++)
    len += v[i]*v[i];

  len = sqrt(len);
  if (len != 0.0)
    for (i=0; i<n; i++)
      v[i] /= len;
}

void Rotate(float *a, float *axis, float angle, float *b)
{
/* rotate point a about axis by angle degrees */
/* from Rogers, "Mathematical Elements for Computer Graphics" */

  float ctheta, stheta, theta, n[3], r[3];

  theta = angle*DEG_TO_RAD;
  ctheta = cos(theta);
  stheta = sin(theta);

  Copy(3, axis, n);
  Normalize(3, n);

  r[0] = a[0] * (n[0]*n[0] + (1.0 - n[0]*n[0])*ctheta) +
         a[1] * (n[0]*n[1]*(1.0 - ctheta) - n[2]*stheta) +
	 a[2] * (n[0]*n[2]*(1.0 - ctheta) + n[1]*stheta);
  r[1] = a[0] * (n[0]*n[1]*(1.0 - ctheta) + n[2]*stheta) +
         a[1] * (n[1]*n[1] + (1.0 - n[1]*n[1])*ctheta) +
	 a[2] * (n[1]*n[2]*(1.0 - ctheta) - n[0]*stheta);
  r[2] = a[0] * (n[0]*n[2]*(1.0 - ctheta) - n[1]*stheta) +
         a[1] * (n[1]*n[2]*(1.0 - ctheta) + n[0]*stheta) +
	 a[2] * (n[2]*n[2] + (1.0 - n[2]*n[2])*ctheta);
  Copy(3, r, b);
}

void Scale(short n, float *a, float s, float *b)
{
  short i;

  for (i=0; i<n; i++)
    b[i] = a[i]*s;
}

void Sub(short n, float *a, float *b, float *c)
{
  short i;

  for (i=0; i<n; i++)
    c[i] = a[i] - b[i];
}

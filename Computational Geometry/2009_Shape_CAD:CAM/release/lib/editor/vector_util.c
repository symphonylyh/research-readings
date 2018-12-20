/************************************************************************
 *                                                                      *
                        Copyright (C) 1992 by
        Massachusetts Institute of Technology, Cambridge, MA
                         All rights reserved
 *                                                                      *
 ************************************************************************/

#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include "gen.h"
#include "editor.h"

/************************************************************************
 *                         midpoint()                                   *
 *                                                                      *
 *     The function returns a vector which is the midpoint              *
 *     of two homogeneous vectors                                       *
 ************************************************************************
 *                        INPUT                                         *
 *            type   name    description                                *
 *           vector   v1     point 1                                    *
 *           vector   v2     point 2                                    *
 ************************************************************************
 *                        OUTPUT                                        *
 *            type   name    description                                *
 *           vector   vm     midpoint, allocated and returned if vm     *
 *                           is not allocated at the time of calling    *
 ************************************************************************/

vector *midpoint(vector *v1, vector *v2, vector *vm)
     /* return a pointer to vector if not allocated*/
{
  if (!vm)                    /* allocate target if nil */
    vm = vectalloc();

  vm->x = (v1->x/v1->w + v2->x/v2->w)/2.0;
  vm->y = (v1->y/v1->w + v2->y/v2->w)/2.0;
  vm->z = (v1->z/v1->w + v2->z/v2->w)/2.0;
  vm->w = 1.0;

  return(vm);
}

/************************************************************************
 *                                                                      *
 *     Function vector_dircos() returns a vector containing the         *
 *     value of the direction cosines of vector v.                      *
 ************************************************************************
 *                        INPUT                                         *
 *            type   name    description                                *
 *           vector   v1     point 1                                    *
 ************************************************************************
 *                        OUTPUT                                        *
 *            type   name    description                                *
 *           vector   vm    direction cosines, allocated and returned   *
 *                          if a is not allocated at the time of calling*
 ************************************************************************/

vector *vector_dircos(vector *v, vector *a)
{
  double2 m;

  if(!a)
    a = vectalloc();
  m = mag(v); 
  a->x = (v->x/m);
  a->y = (v->y/m);
  a->z = (v->z/m);
  a->w = 1.0;

  return(a);
}

/************************************************************************
 *                                                                      *
 *      Function vector_rotate() returns a vector containing            *
 *      the coordinates of the rotated vector v, rotated by             *
 *      the angles g,b,t from Mortensen pg. 351                         *
 ************************************************************************
 *                        INPUT                                         *
 *            type   name    description                                *
 *           vector   v      vector to be rotated                       *
 *           double2  g      angle about the x-axis to rotate           *
 *           double2  b      angle about the y-axis to rotate           *
 *           double2  t      angle about the z-axis to rotate           *
 ************************************************************************
 *                        OUTPUT                                        *
 *            type   name    description                                *
 *           vector   f      rotated vector, allocated and returned if  *
 *                           f is not allocated at the time of calling  *
 ************************************************************************/

vector *vector_rotate(vector *v, double g, double b, double t, vector *f)
{
  double2 r[3][3];

  if(!f)
    f = vectalloc(); 

  r[0][0] = cos(t)*cos(b);
  r[0][1] = sin(t)*cos(g) + cos(t)*sin(b)*sin(g);
  r[0][2] = sin(t)*sin(g) - cos(t)*sin(b)*cos(g);
  r[1][0] = -sin(t)*cos(b);
  r[1][1] = cos(t)*cos(g) - sin(t)*sin(b)*sin(g);
  r[1][2] = cos(t)*sin(g) + sin(t)*sin(b)*cos(g);
  r[2][0] = sin(b);
  r[2][1] = -cos(b)*sin(g);
  r[2][2] = cos(b)*cos(g);

  f->x = (v->x * r[0][0] + v->y * r[1][0] + v->z * r[2][0])/v->w;
  f->y = (v->x * r[0][1] + v->y * r[1][1] + v->z * r[2][1])/v->w;
  f->z = (v->x * r[0][2] + v->y * r[1][2] + v->z * r[2][2])/v->w;
  f->w = 1.0;

  return(f);
}
/************************************************************************
 *              This routine returns the angle between                  *
 *                 two vectors v1, and v2 in radians                    *
 ************************************************************************
 *                        INPUT                                         *
 *            type   name    description                                *
 *           vector   v1     first vector                               *
 *           vector   v2     second vector                              *
 ************************************************************************
 *                        OUTPUT                                        *
 *           type    name    description                                *
 *         double2  return   angle between v1 and v2 in radians         *
 ************************************************************************/

double2 vectang(vector *v1, vector *v2)
{
  vector *v;
  double2 angle;

  angle = atan2(mag(v = cross(v1,v2)),dot(v1,v2));
  vectfree(v);

  return angle;
}

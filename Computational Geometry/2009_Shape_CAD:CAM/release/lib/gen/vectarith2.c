/************************************************************************
 *                                                                      *
                        Copyright (C) 1989 by
        Massachusetts Institute of Technology, Cambridge, MA
                         All rights reserved
 *                                                                      *
 ************************************************************************/

#include <stdio.h>
#include <math.h>

#include "gen.h"

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
/* Functions referenced by midpoint():
 *  vectalloc()
 */

vector *midpoint(vector *v1, vector *v2, vector *vm)

/* return a pointer to vector if not allocated */

{
  if (!vm)                    /* allocate target if nil */
    vm = vectalloc();

  vm->x = (v1->x/v1->w + v2->x/v2->w)/2.0;
  vm->y = (v1->y/v1->w + v2->y/v2->w)/2.0;
  vm->z = (v1->z/v1->w + v2->z/v2->w)/2.0;
  vm->w = 1.0;

  return(vm);
}

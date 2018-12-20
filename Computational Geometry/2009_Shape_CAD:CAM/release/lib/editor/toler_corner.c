/************************************************************************
 *									*
			Copyright (C) 1992 by
	Massachusetts Institute of Technology, Cambridge, MA
			 All rights reserved
 *									*
 ************************************************************************/

#include <math.h>
#include <malloc.h>
#include "gen.h"
#include "bspl.h"
#include "editor.h"

ParSurf *corner_offset(ParSurf **canals, int cm, double tz)

/* Express a spherical segment corresponding to the offset of one corner
   point of plate as a rational B-spline.
   canals is an array of ParSurf containing the canal surfaces (offsets
   to the edges of the deformed plate geometry).
   canals[0] is the canal surface with spine the line v=0 of the
             deformed plate geometry.
   canals[1] is the canal surface with spine the line v=1 of the
             deformed plate geometry.
   canals[2] is the canal surface with spine the line u=0 of the
             deformed plate geometry.
   canals[3] is the canal surface with spine the line u=1 of the
             deformed plate geometry.
   cm characterizes one of the four corners.
   cm = 0 corresponds to the corner u = v = 0 .
   cm = 1 corresponds to the corner u = 1 , v = 0 .
   cm = 2 corresponds to the corner u = 0 , v = 1 .
   cm = 3 corresponds to the corner u = v = 1 .
   tz is the radius of the spherical segment (tolerance).
   The routine returns the spherical segment. */

{
  ParSurf *corner;
  ParCurv *eg,*ef;
  vector *v,*v0,*v1,*v2;
  double a,b;
  int i,j,nu,nv;

/* Meaning of the most important local variables (in alphabetical order)

a        =  space angle between the bounding great semi-circles of the
            spherical segment
corner   =  the spherical segment geometry, returned by this routine.
eg       =  geometry of the first edge of corner.
ef       =  geometry of the second edge of corner.
nu       =  number of vertices along the v=ct. edges.
nv       =  number of vertices along the u=ct. edges.
v        =  position vector to the mid-vertices of corner.

*/

/* Allocate memory */

  corner = fgeomalloc1(3, 3, 3, 5);
  eg = egeomalloc1(canals[0]->vorder, canals[0]->vcontpts);
  ef = egeomalloc1(canals[0]->vorder, canals[0]->vcontpts);
  v = vectalloc();
  v0 = vectalloc();
  v1 = vectalloc();
  v2 = vectalloc();
  nu = canals[0]->ucontpts-1;
  nv = canals[2]->ucontpts-1;
  
/* Extract the appropriate edges from canals */

  switch (cm) {
  case 0 :
    extract_edge(eg, canals[0], 0, 0);
    extract_edge(ef, canals[2], 0, 0);
    break;
  case 1 :
    extract_edge(eg, canals[0], 0, nu);
    extract_edge(ef, canals[3], 0, 0);
    break;
  case 2 :
    extract_edge(eg, canals[1], 0, 0);
    extract_edge(ef, canals[2], 0, nv);
    break;
  case 3 :
    extract_edge(eg, canals[1], 0, nu);
    extract_edge(ef, canals[3], 0, nv);
    break;
  }
  for (i=0; i<3; i++) {
    corner->uknots[i] = 0.0;
    corner->vknots[i] = 0.0;
    corner->uknots[i+3] = 1.0;
    corner->vknots[i+5] = 1.0;
  }
  corner->vknots[3] = corner->vknots[4] = 0.5;
  
/* Compute the control polygon vertices of corner */

  for (i=0; i<3; i++) {
    corner->contpts[i][0] = copyvector(eg->contpts[0], 0);
    corner->contpts[i][4] = copyvector(eg->contpts[4], 0);
    j = i+1;
    corner->contpts[0][j] = copyvector(eg->contpts[j], 0);
    corner->contpts[2][j] = copyvector(ef->contpts[j], 0);
  }
  sub_vect1(eg->contpts[1], eg->contpts[0], v1);
  sub_vect1(ef->contpts[1], ef->contpts[0], v2);
  add_vect1(v1, v2, v);
  unitvector1(v, v);
  scale_vect1(tz, v, v);
  add_vect1(v,eg->contpts[0], v);
  corner->contpts[1][1] = copyvector(v, 0);
  sub_vect1(eg->contpts[3], eg->contpts[4], v1);
  sub_vect1(ef->contpts[3], ef->contpts[4], v2);
  add_vect1(v1, v2, v);
  unitvector1(v, v);
  scale_vect1(tz, v, v);
  add_vect1(v, eg->contpts[3],v);
  corner->contpts[1][3] = copyvector(v, 0);
  add_vect1(eg->contpts[0], eg->contpts[4], v0);
  scale_vect1((double)0.5, v0, v0);
  sub_vect1(eg->contpts[2], v0, v1);
  sub_vect1(ef->contpts[2], v0, v2);
  
/**************    Revised piece of code !! *****************/

  add_vect1(v1, v2, v);
  unitvector1(v, v);
  a = acos(dot(v1, v2)/(tz*tz));
  b = tz/cos(a/2.0);
  scale_vect1(b, v, v);
  add_vect1(v, v0, v);
  corner->contpts[1][2] = copyvector(v, 0);
  
/* Free auxiliary memory */

  vectfree(v2);
  vectfree(v1);
  vectfree(v0); vectfree(v);
  free_egeom(ef);
  free_egeom(eg); 

  return (corner);
}


/************************************************************************
 *                                                                      *
                        Copyright (C) 1996 by
        Massachusetts Institute of Technology, Cambridge, MA
                         All rights reserved

 *                                                                      *
 ************************************************************************/
#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include "gen.h"
#include "bspl.h"


#define M_2_PI		6.28318530717958647692
#define M_PI		3.14159265358979323846
#define M_PI_2		1.57079632679489661923
#define M_SQRT2		1.41421356237309504880
#define M_SQRT1_2	0.70710678118654752440

/* knot vector for the B-spline representation of a circle */  
static double U[12] = {0.0,0.0,0.0,0.25,0.25,0.5,0.5,0.75,0.75,1.0,1.0,1.0};

/****************************************************************************
*                               line_gen()
*****************************************************************************
* 
* 1   Purpose
*     Given two vectors, this function generates the geometry for a straight 
*     line segment in the form of a NURBS curve.
* 
* 2   Specification
*     #include "bspl.h"
*     ParCurv *line_gen(vector **v)
* 
* 3   Description
*     This function returns a NURBS data structure that contains the geometry 
*     for a linear NURBS curve that represents a line that is defined by the 
*     given two points V0, and V1.
* 
* 4   References
*     [1] S. T. Tuohy. Sculptured shape creation, approximation, and 
*         interrogation. Engineer's Thesis, Massachusetts Institute of, 
* 	  Technology Department of Ocean Engineering, Cambridge, Massachusetts,
*         1992.
* 
* 5   Parameters
*        1.vector ** v
*          On entry: an array of length two of vector data structures that 
* 	           contain the points:
* 		   v[0] = V0, and v[1] = V1.
* 
* 6   Return Values, Error Indicators and Warnings
*     The address of the newly allocated structure containing the NURBS curve 
*     is returned.
* 
* 7   Accuracy
*     Not applicable.
* 
* 8   Further Comments
*     Deallocate the NURBS curve structure with free_egeom().
*     For the generation of other primitive curves, see functions circle_gen 
*     (bspl _ prim.c).
* 
* 9   Functions referenced by line_gen() are:
*     copyvector()
*     egeomalloc1()
* 
* 10  Functions that reference line_gen() are:
*     PrimitiveCurvCB()
*     ReadIgesLine()
* 
******************************************************************************/

ParCurv *line_gen(vector **v)
{
  ParCurv *egm;   /* NURBS data structure which will the NURBS curve form
		   * of the straight line. */
  int i;

  /* allocate memory for egm with order 2 and 2 control points */
  egm = egeomalloc1(2, 2);

  /* set the knot vector as {0.0, 0.0, 1.0, 1.0} */
  for(i=0; i<2; i++) {
    egm->knots[i] = 0.0;
    egm->knots[i+2] = 1.0;
  }

  /* set two control points to two given points */
  copyvector(v[0], egm->contpts[0]);
  copyvector(v[1], egm->contpts[1]);

  return egm;
}

/****************************************************************************
*                                circle_gen()
*****************************************************************************
* 
* 1   Purpose
*     Given two vectors and a radius, this function generates the geometry for 
*     a circle in the form of a NURBS curve.
* 
* 2   Specification
*     #include "bspl.h"
*     ParCurv *circle_gen(vector **v, double r)
* 
* 3   Description
*     This function returns a NURBS data structure that contains the geometry 
*     for a quadratic NURBS curve that represents a line that is defined by 
*     the given two points V0, and V1 and the radius r.
*     The circle is generated in the plane through the first vector and normal 
*     to the line defined by the two vectors.
* 
* 4   References
*     [1] S. T. Tuohy. Sculptured shape creation, approximation, and
*     interrogation.
*     Engineer's Thesis, Massachusetts Institute of Technology, Department of
*     Ocean Engineering, Cambridge, Massachusetts, 1992.
* 
* 5   Parameters
*        1.vector ** v
*          On entry: an array of length two of vector data structures that
* 	 contain the points: v[0] = V0, and v[1] = V1.
*        2.double r
*          On entry: the radius of the circle.
* 
* 6   Return Values, Error Indicators and Warnings
*     The address of the newly allocated structure containing the NURBS curve 
*     is returned.
* 
* 7   Accuracy
*     Not applicable.
* 
* 8   Further Comments
*     Deallocate the NURBS curve structure with free_egeom().
*     For the generation of other primitive curves, see functions 
*     line_gen ( in prim.c ).
* 
* 9   Functions referenced by circle_gen() are:
*     copyvector()
*     cross1()
*     egeomalloc1()
*     mag()
*     midpoint1()
*     scale_vect1()
*     sub_vect()
*     unitvector1()
* 
* 10  Functions that reference circle_gen() are:
*     PrimitiveCurvCB()
* 
*****************************************************************************/

ParCurv *circle_gen(vector **v, double r)
{
  ParCurv *egm;   /* NURBS data structure which will contain the NURBS curve
		   * form of the given circle. */
  vector *n;      /* normal vector of the plane on which the circle will be */
  short i;        /* loop index */

  /* allocate memory for the BURBS curve with order 3 and 9 control points */
  egm = egeomalloc1(3, 9);

  /* set the knot vector as the vector given at the beginning of this file */
  for (i=0; i<12; i++)
      egm->knots[i] = U[i];

  n = sub_vect(v[1], v[0]);   /* get the normal by vector substraction
			       * v[1]-v[0] */

  /* initialize the first control point as (1.0, 0.0, 0.0, 1.0) */
  egm->contpts[0]->x = 1.0; egm->contpts[0]->y = 0.0; 
  egm->contpts[0]->z = 0.0; egm->contpts[0]->w = 1.0;

  /* initialize the third control point as the cross product of the vector 
   * (1.0,0.0,0.0) and n. it is orthogonal to both the normal and the 
   * vector (1.0,0.0,0.0). */
  cross1(egm->contpts[0],n,egm->contpts[2]);

  /* if the normal is parallel to the direction (1.0, 0.0, 0.0), then the
   * first control point is set as (0.0, 1.0, 0.0) */
  if (mag(egm->contpts[2]) < ZERO) {
     egm->contpts[0]->x = 0.0; egm->contpts[0]->y = 1.0;
     egm->contpts[0]->z = 0.0; egm->contpts[0]->w = 1.0;
     /* set the third control point as the cross product of the normal and
      *the vector (0.0,1.0,0.0). */
     cross1(egm->contpts[0],n,egm->contpts[2]);
     }

  /* if the normal now is parallel to the direction (0.0,1.0,0.0), then the
   * 1st control point is set as (0.0,0.0,1.0). */
  if (mag(egm->contpts[2]) < ZERO) {
     egm->contpts[0]->x = 0.0; egm->contpts[0]->y = 0.0; 
     egm->contpts[0]->z = 1.0; egm->contpts[0]->w = 1.0;
     cross1(egm->contpts[0],n,egm->contpts[2]);
     }

  /* set the 3rd control point so that it is a control point of a unit circle
   * on the specific plane. */
  unitvector1(egm->contpts[2],egm->contpts[2]);
  /* the 7th control point is symmetrical to the 3rd w.r.t. the origin. */
  scale_vect1(-1.0,egm->contpts[2],egm->contpts[6]);

  /* the 1st control point is 90 degrees clockwise from the 2nd one. */
  cross1(egm->contpts[2],n,egm->contpts[0]);
  unitvector1(egm->contpts[0],egm->contpts[0]);
  /* the 9th control point is at the same position as the 1st one */
  copyvector(egm->contpts[0],egm->contpts[8]);
  /* the 5th control point is symmetrical to the 1st one w.r.t. the origin. */
  scale_vect1(-1.0,egm->contpts[0],egm->contpts[4]);
  
  /* the 2nd control point is 45 degrees clockwise from the 3rd one, but
   * with its magnitude as sqrt(2). */
  midpoint1(egm->contpts[0],egm->contpts[2],egm->contpts[1]);
  unitvector1(egm->contpts[1],egm->contpts[1]);
  scale_vect1(M_SQRT2,egm->contpts[1],egm->contpts[1]);
  /* the 6th control point is symmetrical to the 2nd one w.r.t. the origin. */
  scale_vect1(-1.0,egm->contpts[1],egm->contpts[5]);

  /* the 4th control point is 45 degrees count-clockwise from the 3rd one with
   * its magnitude as sqrt(2). */
  midpoint1(egm->contpts[4],egm->contpts[2],egm->contpts[3]);
  unitvector1(egm->contpts[3],egm->contpts[3]);
  scale_vect1(M_SQRT2,egm->contpts[3],egm->contpts[3]);
  /* the 8th control point is symmetrical to the 4th one w.r.t. the origin. */
  scale_vect1(-1.0,egm->contpts[3],egm->contpts[7]);

  /* scale the control points such that the resulting circle has its
   * radius = r. */
  for (i=0;i<9;i++) {
      scale_vect1(r,egm->contpts[i],egm->contpts[i]);
      if (i%2) {
         scale_vect1(M_SQRT1_2,egm->contpts[i],egm->contpts[i]);
         egm->contpts[i]->w = M_SQRT1_2;
         }
      }

  /* free memory */
  free((char *)n);

  return egm;
}

/****************************************************************************
*                                 plane_gen()
*****************************************************************************
* 
* 1   Purpose
*     Given three vectors, this function generates the geometry for a plane in 
*     the form of a NURBS surface.
* 
* 2   Specification
*     #include "bspl.h"
*     ParSurf *plane_gen(vector **v)
* 
* 3   Description
*     This function returns a NURBS data structure that contains the geometry 
*     for a bilinear NURBS surface that represents a plane that is defined by 
*     the given three points V0, V1, and V2. The plane is such that V1 is a 
*     corner point and V1 and V2 form the adjacent extent of the edges of the 
*     plane.
* 
* 5   Parameters
*        1.vector ** v
*          On entry: an array of length three of vector data structures that 
* 	           contain the points:
* 		   v[0] = V0, v[1] = V1, and v[2] = V2.
* 
* 6   Return Values, Error Indicators and Warnings
*     The address of the NURBS surface structure containing the planar surface 
*     is returned.
* 
* 7   Accuracy
*     Not applicable.
* 
* 8   Further Comments
*     Deallocate the surface with free_fgeom().
*     For the generation of other primitive surfaces, see functions
*     cylinder_gen (bspl-prim.c), cone_gen (bspl-prim.c), sphere_gen (bspl-
*     prim.c), and torus_gen (bspl-prim.c).
* 
* 9   Functions referenced by plane_gen() are:
*     add_vect1()
*     copyvector()
*     fgeomalloc1()
*     sub_vect1()
*  
* 10  Functions that reference plane_gen() are:
*     PrimitiveSurfCB()
* 
******************************************************************************/

ParSurf *plane_gen(vector **v)
{
  ParSurf *fgm;    /* NURBS data structure which will contain the plane
		    * defined by the points. */
  int i;

  /* allocate memory */
  fgm = fgeomalloc1(2,2,2,2);

  /* set the knot vectors in u,v directions. */
  for(i=0;i<2;i++) {
    fgm->uknots[i] = 0.0;
    fgm->vknots[i] = 0.0;
    fgm->uknots[i+2] = 1.0;
    fgm->vknots[i+2] = 1.0;
  }

  /* the given points are the corners of the control polygon. */
  copyvector(v[1],fgm->contpts[0][0]);
  copyvector(v[0],fgm->contpts[0][1]);
  copyvector(v[2],fgm->contpts[1][0]);

  /* calculate another control point */
  sub_vect1(v[2],v[1],fgm->contpts[1][1]);
  add_vect1(v[0],fgm->contpts[1][1],fgm->contpts[1][1]);

  return(fgm);
}

/****************************************************************************
*                                cylinder_gen()
*****************************************************************************
*
* 1   Purpose
*     Given two vectors and a radius, this function generates the geometry for 
*     a cylinder in the form of a NURBS surface.
* 
* 2   Specification
*     #include "bspl.h"
*     ParSurf *cylinder_gen(vector **v, double r)
* 
* 3   Description
*     This function returns a NURBS data structure that contains the geometry 
*     for a linear by quadratic NURBS surface that represents a cylinder that 
*     is defined by the given two points V0, and V1, and the radius r.  The 
*     cylinder is oriented such that the axis connects V0 and V1, and that 
*     this axis is expanded by a radius of value r. The u parametric direction
*     corresponds to the axial direction. The v parametric direction 
*     corresponds to the radial direction,
* 
* 4   References
*     [1] S. T. Tuohy. Sculptured shape creation, approximation, and 
*         interrogation. Engineer's Thesis, Massachusetts Institute of, 
* 	  Technology Department of Ocean Engineering, Cambridge, Massachusetts,
*         1992.
* 
* 5   Parameters
*        1.vector ** v
*          On entry: an array of length two of vector data structures that 
* 	           contain the points: v[0] = V0, and v[1] = V1.
*        2.double r
*          On entry: the radius of the cylinder
* 
* 6   Return Values, Error Indicators and Warnings
*     Not applicable.
* 
* 7   Accuracy
*     Not applicable.
* 
* 8   Further Comments
*     Deallocate the surface with free_fgeom().
*     For the generation of other primitive surfaces, see functions plane_gen 
*     (bspl _ prim.c), cone_gen (bspl _ prim.c), sphere_gen (bspl _ prim.c), 
*     and torus_gen (bspl _ prim.c).
* 
* 9   Functions referenced by cylinder_gen() are:
*     add_vect1()
*     copyvector()
*     cross1()
*     fgeomalloc1()
*     mag()
*     midpoint1()
*     scale_vect1()
*     sub_vect()
*     unitvector1()
* 
* 10  Functions that reference cylinder_gen() are:
*     PrimitiveSurfCB()
* 
******************************************************************************/

ParSurf *cylinder_gen(vector **v, double r)
{
  ParSurf *fgm;   /* the biquadratic B-spline surface representing the
		   * cylinder */
  int i,j;        /* loop indices */
  vector *n;      /* vector connecting the two points v0 and v1 */

  /* allocate memory for the biquadratic B-spline surface. the number of
   * control points is 3 and 9 in u and v direction, respectively. */
  fgm = fgeomalloc1(2,3,2,9);

  /* set the knot vector in u as [0.0, 0.0, 0.0, 1.0, 1.0, 1.0]. */
  for (i=0;i<2;i++) {
    fgm->uknots[i] = 0.0; fgm->uknots[i+2] = 1.0;
  }

  /* set the knot vector in v as the vector at the beginning of the file. */
  for(i=0;i<12;i++) fgm->vknots[i] = U[i];

  /* get the vector connecting v[0] and v[1]. */
  n = sub_vect(v[1],v[0]);

  /* construct the quadratic B-spline form of a unit circle,
   * see circle_gen() */
  fgm->contpts[0][0]->x = 1.0; fgm->contpts[0][0]->y = 0.0; 
  fgm->contpts[0][0]->z = 0.0; fgm->contpts[0][0]->w = 1.0;
  cross1(fgm->contpts[0][0],n,fgm->contpts[0][2]);

  if(mag(fgm->contpts[0][2]) < ZERO) {
    fgm->contpts[0][0]->x = 0.0; fgm->contpts[0][0]->y = 1.0;
    fgm->contpts[0][0]->z = 0.0; fgm->contpts[0][0]->w = 1.0;
    cross1(fgm->contpts[0][0],n,fgm->contpts[0][2]);
  }

  if(mag(fgm->contpts[0][2]) < ZERO) {
    fgm->contpts[0][0]->x = 0.0; fgm->contpts[0][0]->y = 0.0; 
    fgm->contpts[0][0]->z = 1.0; fgm->contpts[0][0]->w = 1.0;
    cross1(fgm->contpts[0][0],n,fgm->contpts[0][2]);
  }

  unitvector1(fgm->contpts[0][2],fgm->contpts[0][2]);
  scale_vect1(-1.0,fgm->contpts[0][2],fgm->contpts[0][6]);

  cross1(fgm->contpts[0][2],n,fgm->contpts[0][0]);
  unitvector1(fgm->contpts[0][0],fgm->contpts[0][0]);
  copyvector(fgm->contpts[0][0],fgm->contpts[0][8]);
  scale_vect1(-1.0,fgm->contpts[0][0],fgm->contpts[0][4]);
  
  midpoint1(fgm->contpts[0][0],fgm->contpts[0][2],fgm->contpts[0][1]);
  unitvector1(fgm->contpts[0][1],fgm->contpts[0][1]);
  scale_vect1(M_SQRT2,fgm->contpts[0][1],fgm->contpts[0][1]);
  scale_vect1(-1.0,fgm->contpts[0][1],fgm->contpts[0][5]);

  midpoint1(fgm->contpts[0][4],fgm->contpts[0][2],fgm->contpts[0][3]);
  unitvector1(fgm->contpts[0][3],fgm->contpts[0][3]);
  scale_vect1(M_SQRT2,fgm->contpts[0][3],fgm->contpts[0][3]);
  scale_vect1(-1.0,fgm->contpts[0][3],fgm->contpts[0][7]);
  /* end of circle generation */

  /* set the control points of the cylinder. */
  for (i=0;i<9;i++) {
      /* temporarily set the control points at the other end as same as the
       * one obtained already. */
      copyvector(fgm->contpts[0][i],fgm->contpts[1][i]);
      /* scale all the control points such that the cylinder has its
       * radius = r. */
      for (j=0;j<2;j++) {
          scale_vect1(r,fgm->contpts[j][i],fgm->contpts[j][i]);
          /* translate the control points such that the two given points are
	   * on the edges of the cylinder. */
          add_vect1(v[j],fgm->contpts[j][i],fgm->contpts[j][i]);
          }
      if (i%2) 
         for (j=0;j<2;j++) {
	   scale_vect1(M_SQRT1_2,fgm->contpts[j][i],fgm->contpts[j][i]);
	   fgm->contpts[j][i]->w = M_SQRT1_2;
	 }
      }

  free((char *)n);
  return(fgm);
}

/****************************************************************************
*                               cone_gen()
*****************************************************************************
* 
* 1   Purpose
*     Given two vectors and a radius, this function generates the geometry for 
*     a cone in the form of a NURBS surface.
* 
* 2   Specification
*     #include "bspl.h"
*     ParSurf *cone_gen(vector **v, double r)
* 
* 3   Description
*     This function returns a NURBS data structure that contains the geometry 
*     for a linear by quadratic NURBS surface that represents a cone that is 
*     defined by the given two points V0, and V1, and a maximum radius r. The 
*     cone is oriented such that the axis connects V0 and V1.  The axis is 
*     expanded by a linear interpolation of a zero radius at V0 to a maximum 
*     radius of r at V1. The u parametric direction corresponds to the axial 
*     direction.
*     The v parametric direction corresponds to the radial direction. The
*     radius of the cone is defined as = ur, where:
*                        u = 0 at V0
* 		       u = 1 at V1
* 
* 4   References
*     [1] S. T. Tuohy. Sculptured shape creation, approximation, and 
*         interrogation. Engineer's Thesis, Massachusetts Institute of, 
* 	  Technology Department of Ocean Engineering Cambridge, Massachusetts,
*         1992.
* 
* 5   Parameters
*        1.vector ** v
*          On entry: an array of length two of vector data structures that 
* 	           contain the points: v[0] = V0, and v[1] = V1.
*        2.double r
*          On entry: the maximum radius of the cone
* 
* 6   Return Values, Error Indicators and Warnings
*     The address of the NURBS surface structure containing the conical 
*     surface is returned.
* 
* 7   Accuracy
*     Not applicable.
* 
* 8   Further Comments
*     Deallocate the surface with free_fgeom().
*     For the generation of other primitive surfaces, see functions plane_gen 
*     (bspl-prim.c), cylinder_gen (bspl-prim.c), sphere_gen (bspl-prim.c), 
*     and torus_gen (bspl-prim.c).
* 
* 9   Functions referenced by cone_gen() are:
*     add_vect1()
*     copyvector()
*     cross1()
*     fgeomalloc1()
*     mag()
*     midpoint1()
*     scale_vect1()
*     sub_vect()
*     unitvector1()
* 
* 10  Functions that reference cone_gen() are:
*     PrimitiveSurfCB()
* 
****************************************************************************/

ParSurf *cone_gen(vector **v, double r)
{
  ParSurf *fgm;     /* the biquadratic B-spline surface representing the
		     * cone. */
  int i,j;          /* loop index */
  vector *n;        /* vector connecting v0 and v1. */

  /* function declaration */
  void midpoint1(vector*,vector*,vector*); 

  /* allocate memory for the biquadratic B-spline */
  fgm = fgeomalloc1(2,3,2,9);

  /* set the knot vector in u as [0.0,0.0,1.0,1.0] */
  for(i=0;i<2;i++) {
    fgm->uknots[i] = 0.0; fgm->uknots[i+2] = 1.0;
  }

  /* set the knot vector in v as the vector at the beginning of the file. */
  for(i=0;i<12;i++) fgm->vknots[i] = U[i];

  /* get the vector connecting v0 and v1. */
  n = sub_vect(v[1],v[0]);

  /* generate a unit circle on the plane normal to n. */
  fgm->contpts[0][0]->x = 1.0; fgm->contpts[0][0]->y = 0.0; 
  fgm->contpts[0][0]->z = 0.0; fgm->contpts[0][0]->w = 1.0;
  cross1(fgm->contpts[0][0],n,fgm->contpts[0][2]);

  if(mag(fgm->contpts[0][2]) < ZERO) {
    fgm->contpts[0][0]->x = 0.0; fgm->contpts[0][0]->y = 1.0;
    fgm->contpts[0][0]->z = 0.0; fgm->contpts[0][0]->w = 1.0;
    cross1(fgm->contpts[0][0],n,fgm->contpts[0][2]);
  }

  if(mag(fgm->contpts[0][2]) < ZERO) {
    fgm->contpts[0][0]->x = 0.0; fgm->contpts[0][0]->y = 0.0; 
    fgm->contpts[0][0]->z = 1.0; fgm->contpts[0][0]->w = 1.0;
    cross1(fgm->contpts[0][0],n,fgm->contpts[0][2]);
  }

  unitvector1(fgm->contpts[0][2],fgm->contpts[0][2]);
  scale_vect1(-1.0,fgm->contpts[0][2],fgm->contpts[0][6]);

  cross1(fgm->contpts[0][2],n,fgm->contpts[0][0]);
  unitvector1(fgm->contpts[0][0],fgm->contpts[0][0]);
  copyvector(fgm->contpts[0][0],fgm->contpts[0][8]);
  scale_vect1(-1.0,fgm->contpts[0][0],fgm->contpts[0][4]);

  midpoint1(fgm->contpts[0][0],fgm->contpts[0][2],fgm->contpts[0][1]);
  unitvector1(fgm->contpts[0][1],fgm->contpts[0][1]);
  scale_vect1(M_SQRT2,fgm->contpts[0][1],fgm->contpts[0][1]);
  scale_vect1(-1.0,fgm->contpts[0][1],fgm->contpts[0][5]);

  midpoint1(fgm->contpts[0][4],fgm->contpts[0][2],fgm->contpts[0][3]);
  unitvector1(fgm->contpts[0][3],fgm->contpts[0][3]);
  scale_vect1(M_SQRT2,fgm->contpts[0][3],fgm->contpts[0][3]);
  scale_vect1(-1.0,fgm->contpts[0][3],fgm->contpts[0][7]);
  /* end of the generation of the unit circle. */

  /* scale the control points such that the radius is r at v1 and 0 at v0. */
  for(i=0;i<9;i++) {
/** the following line is wrong, it causes the apex to be at distance 2n
 ** from the base, where n is the diatance between v0 and v1.
 ** S. L. Abrams, 7/24/96 **/
/** copyvector(n,fgm->contpts[1][i]); **/                 /* radius=0 at v0 */
    scale_vect1(r,fgm->contpts[0][i],fgm->contpts[0][i]); /* radius=r at v1 */

    /* translate the control points such that v0 is the vertex of the cone
     *and v1
       is on the edge of the bottom. */
    for(j=0;j<2;j++) 
      add_vect1(v[j],fgm->contpts[j][i],fgm->contpts[j][i]);
    /* scale the four control points which are not on the circle. */
    if(i%2) 
      for(j=0;j<2;j++) {
	scale_vect1(M_SQRT1_2,fgm->contpts[j][i],fgm->contpts[j][i]);
	fgm->contpts[j][i]->w = M_SQRT1_2;
      }
  }

  free((char *)n);
  return(fgm);
}

/*****************************************************************************
*                                     sphere_gen()
******************************************************************************
* 
* 1   Purpose
*     Given a vector and a radius, this function generates the geometry for a 
*     sphere in the form of a NURBS surface.
* 
* 2   Specification
*     #include "bspl.h"
*     ParSurf *sphere_gen(vector *v, double r)
* 
* 3   Description
*     This function returns a NURBS data structure that contains the geometry 
*     for a biquadratic NURBS surface that represents a sphere that is defined by 
*     the center point V and the radius r.
* 
* 4   References
*     [1] S. T. Tuohy. Sculptured shape creation, approximation, and 
*         interrogation. Engineer's Thesis, Massachusetts Institute of, 
* 	  Technology Department of Ocean Engineering, Cambridge, Massachusetts,
*         1992.
* 
* 5   Parameters
*        1.vector * v
*          On entry: a vector data structures that contain the center point of 
* 	           the sphere.
*        2.double r
*          On entry: the radius of the sphere
* 
* 6   Return Values, Error Indicators and Warnings
*     The address of the NURBS surface structure containing the spherical
*     surface is returned.
* 
* 7   Accuracy
*     Not applicable.
* 
* 8   Further Comments
*     Deallocate the surface with free_fgeom().
*     For the generation of other primitive surfaces, see functions plane_gen 
*     (bspl _ prim.c), cone_gen (bspl _ prim.c), cylinder_gen (bspl _ prim.c), 
*     and torus_gen (bspl _ prim.c).
* 
* 9   Functions referenced by sphere_gen() are:
*     add_vect1()
*     copyvector()
*     fgeomalloc1()
*     midpoint1()
*     scale_vect1()
* 
* 10  Functions that reference sphere_gen() are:
*     PrimitiveSurfCB()
* 
****************************************************************************/

ParSurf *sphere_gen(vector *v, double r)
{
  ParSurf *fgm;   /* B-spline surface representing the sphere */
  int i,j;        /* loop indices */
  vector *n;      /* normal */

  /* function declaration */
  void midpoint1(vector*,vector*,vector*);

  /* allocate memory for the B-spline surface with order 3 and 5, number of
   * points 3 and 9 in u,v direction, respectively. */
  fgm = fgeomalloc1(3,3,5,9);

  /* set the knot vector in u as [0.0,0.0,0.0,0.5,0.5,1.0,1.0,1.0]. */
  for (i=0;i<3;i++) {
      fgm->uknots[i] = 0.0; fgm->uknots[i+5] = 1.0;
    }
  fgm->uknots[3] = fgm->uknots[4] = 0.5;

  /* set the knot vector in v as the vector at the beginning of the file. */
  for(i=0;i<12;i++) fgm->vknots[i] = U[i];

  /* set the control points */
  for(i=0;i<9;i++) {
    fgm->contpts[0][i]->x = 0.0;    fgm->contpts[0][i]->y = r;
    fgm->contpts[0][i]->z = 0.0;    fgm->contpts[0][i]->w = 1.0;

    fgm->contpts[4][i]->x = 0.0;    fgm->contpts[4][i]->y = -r;
    fgm->contpts[4][i]->z = 0.0;    fgm->contpts[4][i]->w = 1.0;
  }

  for(i=1;i<4;i++) {
    fgm->contpts[i][0]->x = r;    fgm->contpts[i][0]->y = 0.0;
    fgm->contpts[i][0]->z = 0.0;    fgm->contpts[i][0]->w = 1.0;
    copyvector(fgm->contpts[i][0],fgm->contpts[i][8]);
    scale_vect1(-1.0,fgm->contpts[i][0],fgm->contpts[i][4]);
    
    fgm->contpts[i][1]->x = r;    fgm->contpts[i][1]->y = 0.0;
    fgm->contpts[i][1]->z = r;    fgm->contpts[i][1]->w = 1.0;
    scale_vect1(-1.0,fgm->contpts[i][1],fgm->contpts[i][5]);
    
    fgm->contpts[i][2]->x = 0.0;    fgm->contpts[i][2]->y = 0.0;
    fgm->contpts[i][2]->z = r;      fgm->contpts[i][2]->w = 1.0;
    scale_vect1(-1.0,fgm->contpts[i][2],fgm->contpts[i][6]);
    
    fgm->contpts[i][3]->x = -r;     fgm->contpts[i][3]->y = 0.0;
    fgm->contpts[i][3]->z = r;    fgm->contpts[i][3]->w = 1.0;
    scale_vect1(-1.0,fgm->contpts[i][3],fgm->contpts[i][7]);
  }

  for (i=0;i<9;i++) {
      fgm->contpts[1][i]->y = r;    
      fgm->contpts[3][i]->y = -r;
      }

  /* translate, scale the control points such that the B-spline surface
   * represents the sphere centered at v and with its radius = r. */
  for (i=0;i<9;i++) {    
      for (j=0;j<5;j++) 
          add_vect1(v,fgm->contpts[j][i],fgm->contpts[j][i]); 
      if (i%2) 
         for (j=0;j<5;j++) {
	     scale_vect1(M_SQRT1_2,fgm->contpts[j][i],fgm->contpts[j][i]);
	     fgm->contpts[j][i]->w = M_SQRT1_2;
             }
      }
  for (j=0;j<5;j++) 
    if (j%2) 
      for (i=0;i<9;i++) {
	scale_vect1(M_SQRT1_2,fgm->contpts[j][i],fgm->contpts[j][i]);
	fgm->contpts[j][i]->w = fgm->contpts[j][i]->w * M_SQRT1_2;
      }
  
  return(fgm);
}

/****************************************************************************
*                              torus_gen()
*****************************************************************************
* 
* 1   Purpose
*     Given two vectors and two radii, this function generates the geometry 
*     for a torus in the form of a NURBS surface.
* 
* 2   Specification
*     #include "bspl.h"
*     ParSurf *torus_gen(vector **v, double *r)
* 
* 3   Description
*     This function returns a NURBS data structure that contains the geometry 
*     for a biquadratic NURBS surface that represents a torus that is defined 
*     by the given two points V0, and V1, and the radui r0 and r1. The torus 
*     is oriented such that a circle of radius r1 is swept around the axis 
*     formed by a segment connecting V0 and V1. The center of the torus is at
*     V0, and the radius of the sweep is r0. The u parametric direction 
*     corresponds to the sweep direction.
* 
* 4   References
*     [1]S. T. Tuohy. Sculptured shape creation, approximation, and 
*        interrogation. Engineer's Thesis, Massachusetts Institute of, 
*        Technology Department of Ocean Engineering, Cambridge, Massachusetts,
*        1992.
* 
* 5   Parameters
*       1.vector ** v
*         On entry: an array of length two of vector data structures that 
* 	          contain the points:
* 		  v[0] = V0, and v[1] = V1.
*       2.double* r
*         On entry: an array of length two that contains the radii: r[0] = r0, 
* 	          and r[1] = r1.
* 
* 6   Return Values, Error Indicators and Warnings
*     The address of the NURBS surface structure containing the toroidal
*     surface is returned.
* 
* 7   Accuracy
*     Not applicable.
* 
* 8   Further Comments
*     Deallocate the surface with free_fgeom().
*     For the generation of other primitive surfaces, see functions plane_gen 
*     (bspl _ prim.c), cylinder_gen (bspl _ prim.c), cone_gen (bspl _ prim.c), 
*     and sphere_gen (bspl _ prim.c).
* 
* 9   Functions referenced by torus_gen() are:
*     add_vect1()
*     copyvector()
*     fgeomalloc1()
*     midpoint1()
*     ParSurf_rotrans()
*     scale_vect1()
*     sub_vect()
*     unitvector1()
*     vectalloc()
* 
* 10  Functions that reference torus_gen() are:
*     PrimitiveSurfCB()
*  
**************************************************************************/

ParSurf *torus_gen(vector *v[2], double r[2])
{
  ParSurf *fgm;      /* B-spline surface representing the torus */
  int i,j;           /* loop indices */
  vector *n,*inner,*rot,*trn;
                     /* n -- vector connecting v0 and v1, which is the axis
		      *      of rotation of the torus
                      * inner -- a translating vector
                      * rot -- a rotating vector
                      * trn -- a translating vector
                      */

  /* function declaration, see midpoint1(), ParSurf_rotrans() */
  void midpoint1(vector*,vector*,vector*);
  void ParSurf_rotrans (ParSurf*,vector*,vector*);

  /* allocate memory for the B-spline surface and translation vector and
     rotation vector. */
  fgm = fgeomalloc1(3,3,9,9); inner = vectalloc();
  trn = vectalloc(); rot = vectalloc();

  /* set the knot vectors in u,v as the vector at the beginning of the file */ 
  for(i=0;i<12;i++) {
    fgm->uknots[i] = U[i]; fgm->vknots[i] = U[i];
  }

  /* construct a torus centered at the origin and with z-axis as its axis of 
   * rotation. */

  /* construct a circle with radius r1 on x-z plane. */ 
  fgm->contpts[0][0]->x = r[1];   fgm->contpts[0][0]->y = 0.0;
  fgm->contpts[0][0]->z = 0.0;    fgm->contpts[0][0]->w = 1.0;
  copyvector(fgm->contpts[0][0],fgm->contpts[0][8]);
  scale_vect1(-1.0,fgm->contpts[0][0],fgm->contpts[0][4]);
  
  fgm->contpts[0][1]->x = r[1];    fgm->contpts[0][1]->y = 0.0;
  fgm->contpts[0][1]->z = r[1];    fgm->contpts[0][1]->w = 1.0;
  scale_vect1(-1.0,fgm->contpts[0][1],fgm->contpts[0][5]);
  
  fgm->contpts[0][2]->x = 0.0;      fgm->contpts[0][2]->y = 0.0;
  fgm->contpts[0][2]->z = r[1];     fgm->contpts[0][2]->w = 1.0;
  scale_vect1(-1.0,fgm->contpts[0][2],fgm->contpts[0][6]);
  
  fgm->contpts[0][3]->x = -r[1];   fgm->contpts[0][3]->y = 0.0;
  fgm->contpts[0][3]->z =  r[1];   fgm->contpts[0][3]->w = 1.0;
  scale_vect1(-1.0,fgm->contpts[0][3],fgm->contpts[0][7]);
  
  inner->x = r[0];   /* set the x-component of the translating vector as
		      * r[0] */ 

  /* translate the control points such that the circle is centered at x=r[0]
   * on x-axis, and set them as the 1st column of the control points in v
   * direction of the torus. */
  for (i=0;i<9;i++) {
      add_vect1(inner,fgm->contpts[0][i],fgm->contpts[0][i]);
      /* the 9th column of control points ( in v direction ) are set as the
       * same as the 1th one. */
      copyvector(fgm->contpts[0][i],fgm->contpts[8][i]);
      /* the 5th column of control points ( in v direction ) are oppsite to
       * the 1th one about the z-axis. */
      copyvector(fgm->contpts[0][i],fgm->contpts[4][i]);
      fgm->contpts[4][i]->x = -1.0 * fgm->contpts[0][i]->x;
    }

  /* construct a circle with radius r1 on the plane which is obtained rotating
   * x-z plane 45 degree about z-axis. */ 
  fgm->contpts[1][0]->x = r[1];   fgm->contpts[1][0]->y = r[1];
  fgm->contpts[1][0]->z = 0.0;    fgm->contpts[1][0]->w = 1.0;
  copyvector(fgm->contpts[1][0],fgm->contpts[1][8]);
  scale_vect1(-1.0,fgm->contpts[1][0],fgm->contpts[1][4]);
  
  fgm->contpts[1][1]->x = r[1];    fgm->contpts[1][1]->y = r[1];
  fgm->contpts[1][1]->z = r[1];    fgm->contpts[1][1]->w = 1.0;
  scale_vect1(-1.0,fgm->contpts[1][1],fgm->contpts[1][5]);
  
  fgm->contpts[1][2]->x = 0.0;      fgm->contpts[1][2]->y = 0.0;
  fgm->contpts[1][2]->z = r[1];     fgm->contpts[1][2]->w = 1.0;
  scale_vect1(-1.0,fgm->contpts[1][2],fgm->contpts[1][6]);
  
  fgm->contpts[1][3]->x = -r[1];   fgm->contpts[1][3]->y = -r[1];
  fgm->contpts[1][3]->z =  r[1];   fgm->contpts[1][3]->w = 1.0;
  scale_vect1(-1.0,fgm->contpts[1][3],fgm->contpts[1][7]);
  
  inner->y = r[0];   /* now, the translating vector is [r0,r0,0] */

  /* translate the control points by r0 both along x and y-axis such that the 
   * circle is centered at [r0,r0,0.0], and also set them as the 2nd column of 
   * control points in v direction of the torus. */ 
  for (i=0;i<9;i++) {
      add_vect1(inner,fgm->contpts[1][i],fgm->contpts[1][i]);
      /* the 6th column of control points in v direction are opposite to the 
       * the 2nd column. */
      copyvector(fgm->contpts[1][i],fgm->contpts[5][i]);
      fgm->contpts[5][i]->x = -1.0 * fgm->contpts[1][i]->x;
      fgm->contpts[5][i]->y = -1.0 * fgm->contpts[1][i]->y;
      }

  /* construct a circle with radius r1 on y-z plane. */
  fgm->contpts[2][0]->x = 0.0;   fgm->contpts[2][0]->y = r[1];
  fgm->contpts[2][0]->z = 0.0;    fgm->contpts[2][0]->w = 1.0;
  copyvector(fgm->contpts[2][0],fgm->contpts[2][8]);
  scale_vect1(-1.0,fgm->contpts[2][0],fgm->contpts[2][4]);
  
  fgm->contpts[2][1]->x = 0.0;    fgm->contpts[2][1]->y = r[1];
  fgm->contpts[2][1]->z = r[1];    fgm->contpts[2][1]->w = 1.0;
  scale_vect1(-1.0,fgm->contpts[2][1],fgm->contpts[2][5]);
  
  fgm->contpts[2][2]->x = 0.0;      fgm->contpts[2][2]->y = 0.0;
  fgm->contpts[2][2]->z = r[1];     fgm->contpts[2][2]->w = 1.0;
  scale_vect1(-1.0,fgm->contpts[2][2],fgm->contpts[2][6]);
  
  fgm->contpts[2][3]->x = 0.0;   fgm->contpts[2][3]->y = -r[1];
  fgm->contpts[2][3]->z =  r[1];   fgm->contpts[2][3]->w = 1.0;
  scale_vect1(-1.0,fgm->contpts[2][3],fgm->contpts[2][7]);
  
  inner->x = 0.0;          /* the translating vector is [0.0,r0,0.0] */

  /* translate the control points such that the circle is centered at 
   * y=r0 on y-axis, and set them as the 3rd column of the control points
   * of the torus. */
  for(i=0;i<9;i++) {
    add_vect1(inner,fgm->contpts[2][i],fgm->contpts[2][i]);
    /* the 7th column of the control points of the torus is oppsite to
     * the 3rd column about z-axis. */
    copyvector(fgm->contpts[2][i],fgm->contpts[6][i]);
    fgm->contpts[6][i]->y = -1.0 * fgm->contpts[2][i]->y;
  }

  /* construct a circle with radius r1 on the plane which is 45 degree from
   * y-z plane. */
  fgm->contpts[3][0]->x = -r[1];   fgm->contpts[3][0]->y = r[1];
  fgm->contpts[3][0]->z = 0.0;    fgm->contpts[3][0]->w = 1.0;
  copyvector(fgm->contpts[3][0],fgm->contpts[3][8]);
  scale_vect1(-1.0,fgm->contpts[3][0],fgm->contpts[3][4]);
  
  fgm->contpts[3][1]->x = -r[1];    fgm->contpts[3][1]->y = r[1];
  fgm->contpts[3][1]->z = r[1];    fgm->contpts[3][1]->w = 1.0;
  scale_vect1(-1.0,fgm->contpts[3][1],fgm->contpts[3][5]);
  
  fgm->contpts[3][2]->x = 0.0;      fgm->contpts[3][2]->y = 0.0;
  fgm->contpts[3][2]->z = r[1];     fgm->contpts[3][2]->w = 1.0;
  scale_vect1(-1.0,fgm->contpts[3][2],fgm->contpts[3][6]);
  
  fgm->contpts[3][3]->x = r[1];   fgm->contpts[3][3]->y = -r[1];
  fgm->contpts[3][3]->z = r[1];   fgm->contpts[3][3]->w = 1.0;
  scale_vect1(-1.0,fgm->contpts[3][3],fgm->contpts[3][7]);
  
  inner->x = -r[0];              /* translating vector = { -r0, r0, 0.0 } */

  /* translate the control points such that the circle is centered at
   * (-r0,r0,0.0), and set them as the 4th column of the control points in v
   * direction of the torus. */
  for(i=0;i<9;i++) {
    add_vect1(inner,fgm->contpts[3][i],fgm->contpts[3][i]);
    /* the 8th column of the control points are oppsite to the 4th column
       about z-axis. */
    copyvector(fgm->contpts[3][i],fgm->contpts[7][i]);
    fgm->contpts[7][i]->x = -1.0 * fgm->contpts[3][i]->x;
    fgm->contpts[7][i]->y = -1.0 * fgm->contpts[3][i]->y;
  }

  /* obtain the real axis of rotation of the torus and normalize it as a unit
   * vector. */
  n = sub_vect(v[1],v[0]);
  unitvector1(n,n);

  /* determine the rotation */
  /* if the axis is not parallel to z-axis */
  if(n->z != 0.0) { 
    rot->x = (180.0/M_PI) * atan(n->y/n->z);
    rot->y = (180.0/M_PI) * atan(n->x/n->z);
  }
  /* if it is parallel to z-axis */
  else {
    rot->y = 90.0;
    rot->x = (180.0/M_PI) * atan(n->y/n->x);
  }

  /* rotate the torus in B-spline form to make its axis of rotation coincide
   * with the vector connecting v0 and v1. */
  ParSurf_rotrans(fgm,rot,trn); 

  /* free memory */
  free((char *) n);  free((char *) rot);  free((char *) trn);

  /* translate the torus to make it centered at v0. */
  for (i=0;i<9;i++) 
      for (j=0;j<9;j++) 
	  add_vect1(v[0],fgm->contpts[j][i],fgm->contpts[j][i]);

  /* scale the control points */
  for (i=0;i<9;i++) 
      if (i%2) 
         for (j=0;j<9;j++) {
	     scale_vect1(M_SQRT1_2,fgm->contpts[j][i],fgm->contpts[j][i]);
	     fgm->contpts[j][i]->w = M_SQRT1_2;
             }

  for (j=0;j<9;j++) 
      if (j%2) 
         for (i=0;i<9;i++) {
	     scale_vect1(M_SQRT1_2,fgm->contpts[j][i],fgm->contpts[j][i]);
	     fgm->contpts[j][i]->w = fgm->contpts[j][i]->w * M_SQRT1_2;
             }
  
  return(fgm);
}
  
/*
 *vector *midpoint(v1,v2)
 *vector *v1,*v2;
 *{
 *  vector *n;
 *  void midpoint1(vector*,vector*,vector*);
 *
 *  n = vectalloc();
 *
 *  midpoint1(v1,v2,n);
 *
 *  return(n);
 *}
 */

/*****************************************************************************
*                                 midpoint1()
******************************************************************************
* 
* 1   Purpose
*     Calculate midpoint of two vector coordinates.
* 
* 2   Specification
*     #include "bspl.h"
*     void midpoint(vector *v1, vector *v2, vector *n)
* 
* 3   Description
*     This function calculates the midpoint of two vector coordinates.
* 
*                                         v1 + v2
*                                     n = _______
*                                            2
* 
* 4   References
*     Not applicable.
* 
* 5   Parameters
*       1.vector * v1
*         On entry: the address of the vector structure containing the first 
* 	          coordinate.
*       2.vector * v2
*         On entry: the address of the vector structure containing the second 
* 	          coordinate.
*       3.vector * v1
*         On exit: the address of a vector structure containing the midpoint 
* 	          of the two coordinates.
* 
* 6   Return Values, Error Indicators and Warnings
*     Not applicable.
* 
* 7   Accuracy
*     Not applicable.
* 
* 8   Further Comments
*     Not applicable.
* 
* 9   Functions referenced by midpoint1() are:
* 
* 10  Functions that reference midpoint1() are:
*     circle_gen()
*     cone_gen()
*     cylinder_gen()
*     sphere_gen()
*     torus_gen()
* 
******************************************************************************/

void midpoint1(vector *v1, vector *v2, vector *n)
{
  n->x = (v1->x + v2->x)/2.0;
  n->y = (v1->y + v2->y)/2.0;
  n->z = (v1->z + v2->z)/2.0;
  n->w = 1.0;
}

/******************************************************************************
*                             ParCurv_rotrans()
*******************************************************************************
* 
* 1   Purpose
*     This routine translates and rotates a B-spline curve.
* 
* 2   Specification
*     #include "bspl.h"
*     void ParCurv_rotrans (ParCurv *bgp, vector *rot, vector *trn)
* 
* 3   Description
*     This function first translates a given B-spline curve by a translating
*     vector, then rotates the B-spline about an axis defined by a vector.
* 
* 4   References
*     Not applicable
* 
* 5   Parameters
*       1.ParCurv *bgp
*         On entry: the B-spline curve to be transformed.
* 	On exit:  the B-spline curve after transformation.
*       2.vector *rot
*         On entry: the rotating vector
*       3.vector *trn
*         On entry: the translating vector
* 
* 6   Return Values, Error Indicators and Warnings
*     Not applicable.
* 
* 7   Accuracy
*     Not applicable.
* 
* 8   Further Comments
*     Not applicable.
* 
* 9   Functions referenced by nodes_bsp() are: 
*     rotx1()
*     roty1()
*     rotz1()
*     translate1()
* 
* 10  Functions that reference nodes_bsp() are: None
* 
***************************************************************************/

void ParCurv_rotrans (ParCurv *bgp, vector *rot, vector *trn)
{
  int iu;
  /* declaration of transformation functions */
  void rotx1(double,vector*,vector*);   /* rotate about x-axis */
  void roty1(double,vector*,vector*);   /* rotate about y-axis */
  void rotz1(double,vector*,vector*);   /* rotate about z-axis */
  void translate1(double,double,double,vector*,vector*);   /* translate */

  /* calculate the control points after translation and rotation */
  for (iu = 0; iu < bgp->ncontpts; ++iu) {
    translate1 (trn->x / trn->w, trn->y / trn->w, trn->z / trn->w,
		bgp->contpts[iu], bgp->contpts[iu]);
    rotx1 (rot->x / rot->w, bgp->contpts[iu], bgp->contpts[iu]);
    roty1 (rot->y / rot->w, bgp->contpts[iu], bgp->contpts[iu]);
    rotz1 (rot->z / rot->w, bgp->contpts[iu], bgp->contpts[iu]);
  }
}

/****************************************************************************
*                            ParSurf_rotrans()
*****************************************************************************
*
* 1   Purpose
*     This function transforms a B-spline surface.
*
* 2   Specification
*     void ParSurf_rotrans (ParSurf *bgp, vector *rot, vector *trn)
*
* 3   Description
*     This function first translates the surface by a translating vector, then
*     rotates it about an axis defined by a vector.
*
* 4   References
*     Not applicable
*
* 5   Parameters
*       1.ParSurf *bgp
*         On entry: the B-spline surface to be transformed.
* 	On exit:  the B-spline surface after transformation.
*       2.vector *rot
*         On entry: the vector defining the rotating axis
*       3.vector *trn
*         On entry: the translating vector
* 
* 6   Return Values, Error Indicators and Warnings
*     Not applicable.
* 
* 7   Accuracy
*     Not applicable.
* 
* 8   Further Comments
*     Not applicable.
* 
* 9   Functions referenced by nodes_bsp() are:
*     rotx1()
*     roty1()
*     rotz1()
*     translate1()
* 
* 10  Functions that reference nodes_bsp() are: None
* 
******************************************************************************/

void ParSurf_rotrans (ParSurf *bgp, vector *rot, vector *trn)
{
  int iu,iv;    /* loop indices */

  /* function declaration */
  /* rotation */
  void rotx1(double,vector*,vector*);
  void roty1(double,vector*,vector*);
  void rotz1(double,vector*,vector*);
  /* translation */
  void translate1(double,double,double,vector*,vector*);
  
  /* calculate the control points after translation and rotation */
  for (iu = 0; iu < bgp->ucontpts; ++iu)
      for (iv = 0; iv < bgp->vcontpts; ++iv) {
          translate1 (trn->x / trn->w, trn->y / trn->w, trn->z / trn->w,
		      bgp->contpts[iu][iv], bgp->contpts[iu][iv]);
          rotx1 (rot->x / rot->w, bgp->contpts[iu][iv], bgp->contpts[iu][iv]);
          roty1 (rot->y / rot->w, bgp->contpts[iu][iv], bgp->contpts[iu][iv]);
          rotz1 (rot->z / rot->w, bgp->contpts[iu][iv], bgp->contpts[iu][iv]);
          }
}

/*****************************************************************************
*                                  rotx1()
******************************************************************************
*
* 1   Purpose
*     Rotate a point about the x-axis.
* 
* 2   Specification
*     #include "bspl.h"
*     void rotx1(double angle, vector *in, vector *out);
* 
* 3   Description
*     This routine rotates a point about the x-axis.
* 
* 4   References
*     Not applicable.
* 
* 5   Parameters
*       1.double angle
*         On entry: the angle of rotation, given in degrees.
*       2.vector * in
*         On entry: the address of a vector structure that contains the point 
* 	          in its original position.
*       3.vector * out
*         On exit: the address of a vector structure that contains the rotated 
* 	         point.
* 
* 6   Return Values, Error Indicators and Warnings
*     Not applicable.
* 
* 7   Accuracy
*     Not applicable.
* 
* 8   Further Comments
*     Not applicable.
* 
* 9   Functions referenced by rotx1() are:
* 
* 10  Functions that reference rotx1() are:
*     ParCurv_rotrans()
*     ParSurf_rotrans()
* 
******************************************************************************/

void rotx1(double angle, vector *in, vector *out)
{
  double y,z,rad,sine,cosine;

  rad = (angle * M_PI) / 180;   /* convert from degree to radian. */
  sine = sin(rad), cosine = cos(rad);

  /* recalculate y, z coordinates, x doesn't change. */
  y = in->y*cosine - in->z*sine, z = in->y*sine + in->z*cosine;

  out->x = in->x;
  out->y = y;
  out->z = z;
  out->w = in->w;
}

/*****************************************************************************
*                                    roty1()
******************************************************************************
* 
* 1   Purpose
*     Rotate a point about the y-axis.
* 
* 2   Specification
*     #include "bspl.h"
*     void roty1(double angle, vector *in, vector *out);
* 
* 3   Description
*     This routine rotates a point about the y-axis.
* 
* 4   References
*     Not applicable.
* 
* 5   Parameters
*        1.double angle
*          On entry: the angle of rotation, given in degrees.
*        2.vector * in
*          On entry: the address of a vector structure that contains the point 
* 	           in its original position.
*        3.vector * out
*          On exit: the address of a vector structure that contains the 
* 	           rotated point.
* 
* 6   Return Values, Error Indicators and Warnings
*     Not applicable.
* 
* 7   Accuracy
*     Not applicable.
* 
* 8   Further Comments
*     Not applicable.
* 
* 9   Functions referenced by roty1() are:
* 
* 10  Functions that reference roty1() are:
*     ParCurv_rotrans()
*     ParSurf_rotrans()
* 
****************************************************************************/

void roty1(double angle, vector *in, vector *out)
{
  double x,z,rad,sine,cosine;

  rad = (angle * M_PI) / 180;  /* convert degree to radian */
  sine = sin(rad), cosine = cos(rad);  

  /* recalculate the x, z coordinates, y doesn't changes. */
  x = in->x*cosine + in->z*sine, z = -in->x*sine + in->z*cosine;

  out->x = x;
  out->y = in->y;
  out->z = z;
  out->w = in->w;
}

/******************************************************************************
*                                   rotz1()
*******************************************************************************
* 
* 1   Purpose
*     Rotate a point about the z-axis.
* 
* 2   Specification
*     #include "bspl.h"
*     void rotx1(double angle, vector *in, vector *out);
* 
* 3   Description
*     This routine rotates a point about the z-axis.
* 
* 4   References
*     Not applicable.
* 
* 5   Parameters
*       1.double angle
*         On entry: the angle of rotation, given in degrees.
*       2.vector * in
*         On entry: the address of a vector structure that contains the point 
* 	          in its original position.
*       3.vector * out
*         On exit: the address of a vector structure that contains the rotated 
* 	          point.
* 
* 6   Return Values, Error Indicators and Warnings
*     Not applicable.
* 
* 7   Accuracy
*     Not applicable.
* 
* 8   Further Comments
*     Not applicable.
* 
* 9   Functions referenced by rotz1() are:
*  
* 10  Functions that reference rotz1() are:
*     ParCurv_rotrans()
*     ParSurf_rotrans()
* 
****************************************************************************/

void rotz1(double angle, vector *in, vector *out)
{
  double x,y,rad,sine,cosine;

  rad = (angle * M_PI) / 180;  /* convert from degree to radian */
  sine = sin(rad), cosine = cos(rad);

  /* recalculate x, y coordinates, z doesn't change. */
  x = in->x*cosine - in->y*sine, y = in->x*sine + in->y*cosine;

  out->x = x;
  out->y = y;
  out->z = in->z;
  out->w = in->w;
}

/*****************************************************************************
*                              translate1()
******************************************************************************
* 
* 1   Purpose
*     This routine translates a vector.
* 
* 2   Specification
*     #include "bspl.h"
*     void translate1(double tx, double ty, double tz, vector *in, vector *out)
* 
* 3   Decription
*     Not applicable
* 
* 4   References
*     Not applicable
* 
* 5   Parameters
*       1.double tx
*         On entry: the value of translation in x direction.
*       2.double ty
*         On entry: the value of translation in y direction.
*       3.double tz
*         On entry: the value of translation in z direction.
*       4.vector *in
*         On entry: address of the vector containing the point to be
*                 translated.
*       5.vector *out
*         On entry: address of the vector containing the point after
*                 translation.
* 
* 6   Return Values, Error Indicators and Warnings
*     Not applicable.
* 
* 7   Accuracy
*     Not applicable.
* 
* 8   Further Comments
*     Not applicable.
* 
* 9   Functions referenced by nodes_bsp() are: None
* 
* 10  Functions that reference nodes_bsp() are:
*     ParCurv_rotrans()
*     ParSurf_rotrans()
* 
****************************************************************************/

void translate1(double tx, double ty, double tz, vector *in, vector *out)
{
  out->x = in->x + in->w*tx;
  out->y = in->y + in->w*ty;
  out->z = in->z + in->w*tz;
  out->w = in->w;
}

ParCurv *arc_gen(vector **v, double r)
{
  /* vector *v[3] : v[0] is the center point   */
  /*                v[1] is the starting point */
  /*                v[2] is the terminal point */
  /* double r :     radius                     */

  /* Adopted from Algorithm A7.1 in Piegl and Tiller, The NURBS Book, */
  /* Springer, 1995, pp. 308-9.                                       */

  ParCurv
    *egeom;    /* the B-spline arc */
  vector
    *P0,      /* temporary vector holding the initial control point */
    *P1,      /* temporary vector holding the middle control point */
    *P2,      /* temporary vector holding the final control point */
    *T0,      /* temporary vector holding initial tangent vector */
    *T2;      /* temporary vector holding final tangent vector */
  double
    angle,    /* starting angle of each segment */
    dtheta,   /* sweep angle per segment */
    s, t,     /* dummy arguments */
    the,      /* the terminal angle */
    theta,    /* the sweep angle = the - ths */
    ths,      /* the starting angle */
    w1;       /* the weight of the middle control points */
              /* the initial and final control points have weight = 1.0 */
  int
    i, j,
    narcs;    /* number of Bezier segments */

  ths = atan2(v[1]->y-v[0]->y, v[1]->x-v[0]->x);     /* starting angle */
  if (ths < 0.0)
    ths += M_2_PI;

  the = atan2(v[2]->y-v[0]->y, v[2]->x-v[0]->x);     /* terminal angle */
  if (the < 0.0)
    the += M_2_PI;

  if (the < ths)
    the += 2.0*M_PI;
  theta = the - ths;                     /* sweep angle of the arc */

  if (theta <= M_PI_2)                   /*   0 to  90 deg, 1 Bezier segment */
    narcs = 1;
  else if (theta <= M_PI)                /*  90 to 180 deg, 2 segments */
    narcs = 2;
  else if (theta <= 30.*M_PI_2)          /* 180 to 270 deg, 3 segments */
    narcs = 3;
  else                                   /* 270 to 360 deg, 4 segments */
    narcs = 4;

  dtheta = theta/narcs;                  /* sweep angle per Bezier segment */

  /* number of control points is 1 + 2*narcs */
  /* 3 for the first segment and 2 more for each additional segment */

  egeom = egeomalloc1(3, 1 + 2*narcs);

  angle = ths;                           /* starting angle */
  w1 = cos(dtheta/2.0);                  /* weight for middle control points */

  P0 = vectalloc();                      /* initial control point of segment */
  P1 = vectalloc();                      /* middle control point of segment */
  P2 = vectalloc();                      /* final control point of segment */

  T0 = vectalloc();                      /* initial tangent of segment */
  T2 = vectalloc();                      /* final tangent of segment */

  P0->x = v[0]->x + r*cos(angle);        /* initial point for 1st segment */
  P0->y = v[0]->y + r*sin(angle);
  P0->z = v[0]->z;

  T0->x = -sin(angle);                   /* initial tangent for 1st segment */
  T0->y =  cos(angle);

  copyvector(P0, egeom->contpts[0]);     /* copy point to control polygon */

  for (i=1,j=0; i<=narcs; i++,j+=2) {    /* for each segment ... */
    angle += dtheta;                     /* terminal angle for this segment */

    P2->x = v[0]->x + r*cos(angle);      /* final point for this segment */
    P2->y = v[0]->y + r*sin(angle);
    P2->z = v[0]->z;

    T2->x = -sin(angle);                 /* final tangent for this segment */
    T2->y =  cos(angle);

    copyvector(P2, egeom->contpts[j+2]); /* copy point to control polygon */

    arc_gen_intersect(P0, T0, P2, T2, P1, &s, &t); /* middle point for this */
    P1->x = P1->x*w1;                    /* segment multiplied by the */
    P1->y = P1->y*w1;                    /* homogeneous weight */
    P1->z = P1->z*w1;
    P1->w = w1;

    copyvector(P1, egeom->contpts[j+1]); /* copy point to control polygon */
    
    if (i < narcs) {
      copyvector(P2, P0);                /* copy P2 to P0 for next segment */
      copyvector(T2, T0);                /* copy T2 to T0 for next segment */
    }
  }

  vectfree(P0);                          /* free vector storage */
  vectfree(P1);
  vectfree(P2);

  vectfree(T0);
  vectfree(T2);

  for (i=0; i<3; i++) {
    egeom->knots[i]                 = 0.0;           /* starting knots */
    egeom->knots[i+egeom->ncontpts] = 1.0;           /* terminal knots */
  }
  for (i=1; i<narcs; i++) {
    egeom->knots[i*2+1] = (double)i/(double)narcs;   /* internal knots */
    egeom->knots[i*2+2] = (double)i/(double)narcs;
  }

  return egeom;
}

int arc_gen_intersect(vector *P0, vector *T0, vector *P2, vector *T2,
		      vector *P1, double *s, double *t)
{
  double denom, P0T0x, P0T0y, P2T2x, P2T2y;

  /* Modified from: "Faster Line Segment Intersection," in David Kirk (ed.),
   * Graphics Gems III (Academic Press, 1992), pp. 199-202
   */

  /* P0 is a position vector and T0 is the tangent vector at P0 */
  /* P2 is a position vector and T2 is the tangent vector at P2 */

  /* let A = P0 and B = P0 + T0, and C = P2 and D = P2 + T2               */
  /* and s = [(A->y - C->y)(D->x - C->x) - (A->x - C->x)(D->y - C->y)] /  */
  /*         [(B->x - A->x)(D->y - C->y) - (B->y - A->y)(D->x - C->x)]    */
  /* and t = [(A->y - C->y)(B->x - A->x) - (A->x - C->x)(B->y - A->y)] /  */
  /*         [(B->x - A->x)(D->y - C->y) - (B->y - A->y)(D->x - C->x)]    */

  /* then the intersection P1 is: P1->x = A->x + s(B->x - A->x)           */
  /*                              P1->y = A->y + s(B->y - A->y)           */

  P0T0x = P0->x + T0->x;
  P0T0y = P0->y + T0->y;
  P2T2x = P2->x + T2->x;
  P2T2y = P2->y + T2->y;

  denom = (P0T0x - P0->x)*(P2T2y - P2->y) - (P0T0y - P0->y)*(P2T2x - P2->x);
  if (fabs(denom) < 1.0e-10)
    return 1;

  *s = ((P0->y - P2->y)*(P2T2x - P2->x) - (P0->x - P2->x)*(P2T2y - P2->y)) /
       denom;

  *t = ((P0->y - P2->y)*(P0T0x - P0->x) - (P0->x - P2->x)*(P0T0y - P0->y)) /
       denom;
  
  P1->x = P0->x + (*s)*T0->x;
  P1->y = P0->y + (*s)*T0->y;
  P1->z = P0->z;

  return 0;
}

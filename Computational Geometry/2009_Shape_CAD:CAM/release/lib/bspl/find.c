/************************************************************************
 *									*
			Copyright (C) 1996 by
	Massachusetts Institute of Technology, Cambridge, MA
			 All rights reserved

 *									*
 ************************************************************************/

#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include "bspl.h"
#include "gen.h"

/*****************************************************************************
*                                 find_points()
******************************************************************************
* 
* 1   Purpose
*     This function evaluates a curve on a suface at uniform parametric values.
* 
* 2   Specification
*     #include "bspl.h"
*     void find_points(ParCurv *eg, ParSurf *fgeom, double *tf, double **cq,, 
*                      int n, int cntr[])
* 
* 3   Description
*     This function calculates a set of points on a curve on a B-spline
*     surface.  The curve is obtained by mapping a B-spline curve in the
*     parametric domain onto the surface. The points are uniformly
*     distributed along the curve in the parametric domain. 
* 
* 4   References
*     Not applicable
* 
* 5   Parameters
*     1.ParCurv *eg
*       On entry: a B-spline curve in the parametric domain
*     2.ParSurf *fgeom
*       On entry: a B-spline surface on which the curve is mapped.
*     3.double *tf
*       On entry: an array of length n containing parametric values where the
*                 curve on surface is to be evaluated.
*     4.double *cq
*       On exit:  2D array of nx4, containing the homogeneous coordinates 
*                 x,y,z,w of all the evaluated points.
*     5.int n
*       On entry: the length of the array *tf.
*     6.int cntr[]
*       On entry: if cntr[0] = 1, the curve is non-perodic; else it's perodic.
*                 if cntr[1] = 1, the surface is non-perodic; else it is.
* 
* 6   Return Values, Error Indicators and Warnings
*     Not applicable
* 
* 7   Accuracy
*     Not applicable
* 
* 8   Further Comments
*     Not applicable
* 
* 9   Functions referenced by find_points() are:
*     rbspeval()
*     rbspeval_per()
*     revalderivsurf()
*     revalderivsurf_per()
* 
* 10  Functions that reference find_points() are: None
* 
*****************************************************************************/

find_points(ParCurv *eg, ParSurf *fgeom, double *tf, double **cq, int n, int cntr[])
{
  int i,j;                /* i -- the order of the curve
			   * j -- the number of control points of the curve */
  double u,v,dt;          /* u,v -- parametric values of the points on the
			   *        curve on the surface
			   * dt -- the uniform parametric difference between 
			   *       the points to be evaluated on the curve  */
  vector *v0,*v1;         /* v0 -- point of the curve in the domain of the,
			   *       surface consisting of x(u), y(v) and w.
			   * v1 -- point on the curve on the surface */

  /* if the curve is non-perodic */
  if (cntr[0]) {
     i = eg->order ; j = eg->ncontpts ;
     /* calculate the parametric distance between two points to be evaluated */
     dt=(eg->knots[i+j-1]-eg->knots[0])/(double)(n-1) ;
     /* tf varies from 0.0 to 1.0 */
     tf[0] = 0.0 ; tf[n-1] = 1.0 ;
     /* in the parametric domain of the surface, evaluate the curve at tf[i] */
     for (i=0;i<n;i++) {
         /* v0 is the point on the curve in the parametric domain of the
	  * surface. */
         v0 = rbspeval(eg,tf[i],0) ; 
         /* obtain the parametric values u and v of the corresponding point
	  * of the  curve on the surface. */
	 u = v0->x ; v = v0->y ;
         free((char *)v0) ;
         /* evaluate the homogeneous coordinates of the point */
	 /* if the surface is periodic. */
         if (cntr[1]) v1 = revalderivsurf_per(fgeom,u,v,0,0) ;
         /* if the surface is non-periodic. */
         else v1 = revalderivsurf(fgeom,u,v,0,0) ;
         /* put the coordinates into array cq */
         cq[i][0] = v1->x; cq[i][1] = v1->y; cq[i][2] = v1->z;
	 cq[i][3] = v1->w;
         free((char *)v1) ;
         /* if not exceed the range, evaluate next point. */
         if (i<n-2) tf[i+1] = tf[i]+dt ;
        }
     }
  /* if the curve is perodic */
  else {
     dt = 1.0/(double)n ; tf[0] = 0.0 ;
     for (i=0;i<n;i++) {
          /* evaluate the point on the curve in the domain, and assign u,v
	   * parametric values of the corresponding point on the curve on
	   * the surface. */
          v0 = rbspeval_per(eg,tf[i],0) ; u = v0->x ; v = v0->y ;
          free((char *)v0) ;
	  /* if the surface is perodic */
          if (cntr[1]) v1 = revalderivsurf_per(fgeom,u,v,0,0) ;
	  /* if the surface is non-perodic */
          else v1 = revalderivsurf(fgeom,u,v,0,0) ;
	  /* put the coordinate into cq[] */
          cq[i][0] = v1->x; cq[i][1] = v1->y; cq[i][2] = v1->z;
	  cq[i][3] = v1->w;
          free((char *)v1) ;
	  /* if not exceed the range, go on to the next point. */
          if (i<n-1) tf[i+1] = tf[i]+dt ;
          }
     }
}

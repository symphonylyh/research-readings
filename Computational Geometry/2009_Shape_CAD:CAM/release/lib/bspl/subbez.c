/* ***************************************************************************
 Copyright (C) 1996 Massachusetts Institute of Technology all rights reserved 
	Programmer: George A. Kriezis
**************************************************************************** */
#include <stdio.h>
#include "gen.h"
#include "bspl.h"

/******************************************************************************
*                                subbezier()
*******************************************************************************
* 
* 1   Purpose
*     This routine subdivides a Bezier patch out of a B-spline surface.
* 
* 2   Specification
*     #include "bspl.h"
*     void subbezier(ParSurf *fgeom, double ua, double ub, double va,
*                    double vb, ParSurf *subfgeom)
* 
* 3   Description
*     This function subdivides a small patch of Bezier surface out of a B- 
*     spline surface. The parametric range of the Bezier patch is specified by 
*     the lower and upper bounds both in u and v directions.
* 
* 4   References
*     Not applicable
* 
* 5   Parameters
*       1.ParSurf *fgeom
*         On entry: NURBS data structure containing the B-spline surface to be
*                   subdivided.
*       2.double ua
*         On entry: the lower bound in u direction of the subdivided region.
*       3.double ub
*         On entry: the upper bound in u direciton of the subdivided region.
*       4.double va
*         On entry: the lower bound in v direciton of the subdivided region.
*       5.double vb
*         On entry: the upper bound in v direction of the subdivided region.
*       6.ParSurf *subfgeom
*         On entry: NURBS data structure
*         On exit:  NURBS data structure containing the small Bezier patch.
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
* 9   Functions referenced by subbezier() are:
*     copyvector()
*     fgeomalloc2()
*     find()
*     free_fgeom()
*     lf2comp()
*     surfoslo3()
*     vectalloc()
* 
* 10 Functions that reference subbezier() are: None
* 
****************************************************************************/

void subbezier(ParSurf *fgeom, double ua, double ub, double va, double vb,
	       ParSurf *subfgeom)
{
  int is,js, i,j,uflag, vflag;           /* indices */
  int uorder, vorder,ucontpts, vcontpts; /* order & number of control points */
  struct fgeom *largefgeom;      /* the B-spline surface after inserting those
                                  * splitting points ua, ub, va, vb. This
                                  * surface contains the smaller one. */ 
  double eps,u1, v1, u2,v2;      /* the lower and upper bounds of the B-spline
                                  * surface */
  struct vector *v;   /* to be used to contain a control point */

  eps = 1.0e-13;      /* tolerance for comparison */

  /* allocate memory for the whole surface whose knot vector is to be augmented
   * by knot insertion. */
  uorder = fgeom->uorder;
  vorder = fgeom->vorder;
  largefgeom = fgeomalloc2((short)uorder,(short)vorder,
	       		   (short)(fgeom->ucontpts + 2*uorder),
			   (short)(fgeom->vcontpts + 2*vorder));

  /* the lower and upper bounds in u direction of the surface */
  u1 = fgeom->uknots[0];
  u2 = fgeom->uknots[fgeom->ucontpts+fgeom->uorder-1];

  /* the lower and upper bounds in v direction of the surface */
  v1 = fgeom->vknots[0];
  v2 = fgeom->vknots[fgeom->vcontpts+fgeom->vorder-1];

  /* first, subdivide the surface in u direction. */
  j = 0;
  uflag = 0;    
  /* find where the subdivided Bezier surface starts in u direction. */
  is = find(fgeom->ucontpts,fgeom->uknots,ua);  /* here uknots[is] <= ua */

  /* set the new knot vector of the B-spline surface. */ 
  /* from i=0 to is, the knot vector of the B-spline surface doesn't change */ 
  for (i=0; i<=is; i++)
      largefgeom->uknots[j++] = fgeom->uknots[i];
  /* now, make the startpoint ua as a knot of the B-spline with its
     mulitiplicity  = order in u. */
  js = is;
  if (lf2comp(u1, ua, eps) != 0)
     {   /* if the Bezier starts inside the range of the B-spline surface in
	  * u direction. */
      i = 0;
      while (i<uorder){
	    /* if the startpoint ua of the Bezier is one knot of the B-spline
	     * surface, count its multiplicity and move backward to the last
	     * knot which is smaller than ua. */ 
            if (lf2comp(ua,fgeom->uknots[js],eps) == 0){
	       i++; /* account for b-splines */
	       js--;
               }
	    /* make ua as a knot with multiplicity = order in u,  of the 
	     * B-spline surface. */ 
            else{
	       largefgeom->uknots[j++] = ua;
	       i++;
	     }
	  }
      uflag = js+1;  /* index for the start of the Bezier surface, which will 
		      * to be used to set the control points of the Bezier. */
     }
  else{
     /* if the Bezier surface starts at the beginning of the B-spline surface,
      * the index for the start is simply 0. */
     uflag = 0;
     }

  /* find where the subdivision ends. */
  js = find(fgeom->ucontpts,fgeom->uknots,ub);   /* here uknots[js] <= ub */
  /* copy the knot vector after the startpoint and before the endpoint */
  for (i=is+1; i<=js; i++)
      largefgeom->uknots[j++] = fgeom->uknots[i];
  /* make the endpoints ub as a knot of the B-spline with its multiplicity  = 
     order in u. */
  is = js;  
  if (lf2comp(ub, u2, eps) != 0)
     {    /* if the Bezier ends inside the range of the B-spline surface in u 
	   * direction. */
      i = 0;
      while (i<uorder){
            /* if the endpoint ub is a knot of the B-spline surface, count its
	     *  multiplicity, and move backward to the last knot which is
	     * smaller than ub. */
            if (lf2comp(ub,fgeom->uknots[is],eps) == 0){
	       i++; /* account for b-splines */
	       is--;
               }
	    /* make ub as a knot with mulitiplicity = order in u */
            else{
	       largefgeom->uknots[j++] = ub;
	       i++;
               }
            }
    }
  /* set the knot vector after the end point. */
  for (i=js+1; i<uorder+fgeom->ucontpts; i++)
      largefgeom->uknots[j++] = fgeom->uknots[i];
  /* set the number of control points in u direction. */
  largefgeom->ucontpts = j - fgeom->uorder;

  /************the end of subdivision in u direction*************/

  /* then, subdivide the surface in v direciton. the process is very similar
   * to the subdivision in u direction. */
  j = 0;
  vflag = 0;
  /* find the startpoint. */
  is = find(fgeom->vcontpts,fgeom->vknots,va);

  /* set the knot vector in v direction. both the startpoint and the endpoint
   * will be knots in the knot vector of the B-spline surface. their
   * mulitiplicities =  order in v. vflag tells where the Bezier starts. */
  /* set the knot vector before the startpoint. */
  for (i=0; i<=is; i++)
      largefgeom->vknots[j++] = fgeom->vknots[i];
  /* deal with the startpoint va. */
  js = is;
  if (lf2comp(v1, va, eps) != 0)
     {
      i = 0;
      while (i<vorder){
            if (lf2comp(va,fgeom->vknots[js],eps) == 0){
	       i++;  /* account for b-splines */
	       js--;
	     }
            else{
	       largefgeom->vknots[j++] = va;
	       i++;
	     }
	  }
      vflag = js+1;  /* index for point start in bezier piece */
     }
  else{
     vflag = 0;
     }
  /* deal with the endpoint vb. */
  js = find(fgeom->vcontpts,fgeom->vknots,vb);
  for (i=is+1; i<=js; i++)
      largefgeom->vknots[j++] = fgeom->vknots[i];
  is = js;  
  if (lf2comp(vb, v2, eps) != 0)
     {
      i = 0;
      while (i<vorder){
            if (lf2comp(vb,fgeom->vknots[is],eps) == 0){
	       i++; /* account for b-splines */
	       is--;
	     }
            else{
	       largefgeom->vknots[j++] = vb;
	       i++;
	     }
	  }
    }
  /* set the knot vector after the endpoint. */
  for (i=js+1; i<vorder+fgeom->vcontpts; i++)
      largefgeom->vknots[j++] = fgeom->vknots[i];
  /* set the number of control points in v direction */
  largefgeom->vcontpts = j - fgeom->vorder;

  /*************the end of subdivision in v direction*************/

  /* set the orders of the new B-spline. */    
  largefgeom->uorder = uorder;
  largefgeom->vorder = vorder;

  /* calculate the control points due to the change of the knot vectors. */
  surfoslo3(fgeom, largefgeom, 4);
 
  /* set the orders of the Bezier */
  subfgeom->uorder = fgeom->uorder;
  subfgeom->vorder = fgeom->vorder;
 
  /* set the number of control points of the Bezier */
  ucontpts = subfgeom->ucontpts = uorder;
  vcontpts = subfgeom->vcontpts = vorder;

  /* set the knot vector of the Bezier */
  for (i=0; i<uorder; i++)
      {
       subfgeom->uknots[i] = ua;
       subfgeom->uknots[i+uorder] = ub;
      }

  for (i=0; i<vorder; i++)
      {
       subfgeom->vknots[i] = va;
       subfgeom->vknots[i+vorder] = vb;
      }

  /* set the control points of the Bezier */
  for (i=0; i<ucontpts; i++)
      for (j=0; j<vcontpts; j++)
          {
           if ((v = subfgeom->contpts[i][j]) == NULL)
	      v = subfgeom->contpts[i][j] = vectalloc();
           copyvector(largefgeom->contpts[uflag+i][vflag+j], v);
          }
  /* free memory */
  free_fgeom(largefgeom);
}

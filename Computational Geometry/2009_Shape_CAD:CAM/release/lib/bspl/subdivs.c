/* ***************************************************************************
 Copyright (C) 1996 Massachusetts Institute of Technology all rights reserved 
	Programmer: George A. Kriezis
**************************************************************************** */
#include <stdio.h>
#include <math.h>
#include "bspl.h"

/*****************************************************************************
*                                 subdivids()
******************************************************************************
* 
* 1   Purpose
*     This function splits a NURBS surface to two NURBS surfaces at a given 
*     parameter value.
* 
* 2   Specification
*     #include "bspl.h"
*     void subdivids(ParSurf *fgeom1, ParSurf *fgeom2, double tnew, char dir, 
*                    int iadd)
* 
* 3   Description
*     This routine applies the Boehm knot insertion algorithm to a surface to 
*     split the surface in two pieces. It adds a number of knots equal to the 
*     order of the surface in the required parametric direction. It returns 
*     the two pieces of the surface. It keeps the same parametrization in the 
*     two surfaces to distinguish the pieces.
* 
* 4   References
*     [1] W. Boehm, Inserting New Knots Into B-Spline Curves, CAD, (12)4:199-
*       201, 1980.
* 
* 5   Parameters
*       1. ParSurf * fgeom1
*         On entry: contains the original NURBS surface data structure
*         On exit:  contains the NURBS surface data structure representing the 
* 	          lower portion, parametrically, of the split surface.
*       2. ParSurf * fgeom1
*         On entry: contains an allocated NURBS surface data structure
*         On exit:  contains the NURBS surface data structure representing the 
* 	          upper portion, parametrically, of the split surface.
*       3. double tnew
*         On entry: the parameter value at which the surface is to be split
*       4. int iadd
*         On entry: the number of knots to be added to split the surface.
*       5. char dir
*         On entry: parametric direction to insert knot at and split the,
* 	          surface 'u' and 'v'.
* 
* 6   Return Values, Error Indicators and Warnings
*     Not applicable.
* 
* 7   Accuracy
*     Not applicable.
* 
* 8   Further Comments
*     Adequate allocation of memory should be available on entry for both 
*     NURBS surface data structures. Unless it is known that there is one or
*     more knots already in the knot location at which the surface is to be
*     more knots  split, the value of iadd should equal the order of the
*     surface in the required parametric direction.  Otherwise, iadd should
*     equal the order minus the number of knots already present.
* 
* 9   Functions referenced by subdivids() are:
*     boehm_surface()
* 
* 10  Functions that reference subdivids() are:
* 
***************************************************************************/

void subdivids(ParSurf *fgeom1, ParSurf *fgeom2, double tnew, char dir,
	       int iadd)
/* split a b-spline surface to two b-spline surfaces */
/* keep the same parametrization to distinguish pieces  */
/* boehm's algorithm applied at a knot value of surface */
/* struct fgeom *fgeom1,*fgeom2;  original surface is changed and a
				 second surface is also retrurned */
/* double tnew;  knot value subdivide */
/* int iadd;  number of knots to add at a point = uorder if splitting */
/* char dir;  direction of knot addition */
{
  int jj,i,j,k,ihelp;  /* indices */
  ihelp = 0;

  /* insert tnew in fgeom using boehm's knot insertion algorithm */
  if (dir == 'u')
     /* if it is split in u direction, insert tnew to the u knot vector. */ 
     boehm_surface(fgeom1,tnew,iadd,dir);
  else 
     /* if it is split in v direction, insert tnew to the v knot vector. */
     if (dir == 'v')
       boehm_surface(fgeom1,tnew,iadd,dir);

  /* now, break fgeom to two different surfaces */

  if (dir == 'u')
     {     /* split in u direction */
      /* find the first knot with value tnew. its index is ihelp. */
      for (i=fgeom1->uorder; i<fgeom1->ucontpts; i++)
          if (FABS(fgeom1->uknots[i] - tnew) <1.e-10){
	     ihelp = i+fgeom1->uorder;
	     break;
             }
      /* set the values for the 2nd of the split surfaces, which has larger
       * parametric values. */
/* k is never used!!!! */
      k = fgeom1->ucontpts-ihelp+fgeom1->uorder;
      fgeom2->uorder = fgeom1->uorder;             /* order in u direction */
      fgeom2->ucontpts = fgeom1->ucontpts-ihelp+fgeom1->uorder;
                                /* number of control points in u direction */
      fgeom2->vorder = fgeom1->vorder;             /* order in v direction */
      fgeom2->vcontpts = fgeom1->vcontpts;   
                                /* number of control points in v direction */
      /* the v knot vector, unchanged. */
      for (i=0; i<fgeom1->vorder+fgeom1->vcontpts; i++)
          fgeom2->vknots[i]=fgeom1->vknots[i];
      /* set the knot vector in u of the second surface. */
      for (i=0; i<fgeom1->uorder; i++)
          fgeom2->uknots[i]=tnew;
      for (i=ihelp; i<fgeom1->ucontpts+fgeom1->uorder; i++){
          fgeom2->uknots[i-ihelp+fgeom1->uorder]=fgeom1->uknots[i];
          fgeom1->uknots[i]=0.0;  /* by the way, modify the knot vector in 
	 		      	   * u direction of the first surface. */
          }
     j= fgeom1->ucontpts;
     /* set the number of control points in u direction for 1st surface. */
     fgeom1->ucontpts = ihelp-fgeom1->uorder;
     /* set the control points of the 2nd surface. */
     for (jj=0; jj<fgeom1->vcontpts; jj++){
         for (i=fgeom1->ucontpts; i<j; i++)
	     fgeom2->contpts[i-fgeom1->ucontpts][jj] = fgeom1->contpts[i][jj];
         }
     }
  else if (dir == 'v')
     {   /* split in v direction */
         /* find the first knot with value tnew in the v knot vector. its
	  * index is ihelp. */
      for (i=fgeom1->vorder; i<fgeom1->vcontpts; i++)
          if (FABS(fgeom1->vknots[i] - tnew) <1.e-10){
	     ihelp = i+fgeom1->vorder;
	     break;
             }
/* again, k is never used!!!!!! */
      k = fgeom1->vcontpts-ihelp+fgeom1->vorder;
      j= fgeom1->vorder;
      /* set the values for the 2nd of the split surfaces, which has larger
	 parametric values. */
      fgeom2->vorder = fgeom1->vorder;           /* order in v direction */
      fgeom2->vcontpts = fgeom1->vcontpts-ihelp+fgeom1->vorder;
                              /* number of control points in v direction */
      fgeom2->uorder = fgeom1->uorder;           /* order in u direction */
      fgeom2->ucontpts = fgeom1->ucontpts;
                              /* number of control points in u direction */
      /* the u knot vector, unchanged. */
      for (i=0; i<fgeom1->uorder+fgeom1->ucontpts; i++)
          fgeom2->uknots[i]=fgeom1->uknots[i];
      /* set the v knot vector in v for the 2nd surface */
      for (i=0; i<fgeom1->vorder; i++)
          fgeom2->vknots[i]=tnew;
      for (i=ihelp; i<fgeom1->vcontpts+fgeom1->vorder; i++){
          fgeom2->vknots[i-ihelp+fgeom1->vorder]=fgeom1->vknots[i];
          fgeom1->vknots[i]=0.0;  /* modify the knot vector in v */
          }
      /* set the control points for the 2nd surface */
      j= fgeom1->vcontpts;
      fgeom1->vcontpts = ihelp-fgeom1->vorder;
      for (jj=0; jj<fgeom1->ucontpts; jj++){
          for (i=fgeom1->vcontpts; i<j; i++)
	      fgeom2->contpts[jj][i-fgeom1->vcontpts] = fgeom1->contpts[jj][i];
          }
      }
}

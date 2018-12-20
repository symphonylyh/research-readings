/* ***************************************************************************
 Copyright (C) 1996 Massachusetts Institute of Technology all rights reserved 
	Programmer: George A. Kriezis
**************************************************************************** */
#include <stdio.h>
#include <math.h>
#include "gen.h"
#include "bspl.h"

/******************************************************************************
*                                  subdivbs()
*******************************************************************************
* 
* 1   Purpose
*     This routine subdivides a B-spline curve to two B-spline curves
* 
* 2   Specification
*     #include "bspl.h"
*     void subdivbs(ParCurv *egeom1, ParCurv *egeom2, double tnew)
* 
* 3   Description
*     This routine subdivides a B-spline curve by using boehm's knot insertion
*     algorithm to insert a new knot tnew r times, where r is the order of the 
*     B-spline curve, and then breaking it to two B-spline curves with same 
*     order as the previous one.
* 
* 4   References
*     Not applicable
* 
* 5   Parameters
*       1.ParCurv * egeom1
*         On entry: NURBS curve data structure containing the B-spline curve to
* 	          be subdivided.
* 	On exit:  NURBS curve data structure containing the first curve of
* 	          the subdivision.
*       2.ParCurv * egeom2
*         On entry: NURBS curve data structure.
* 	On exit:  NURBS curve data structure containing the second curve of
* 	          the subdivision.
*       3.double tnew
*         On entry: the parametric value of the point where the subdivision is
* 	          to be proceeded.
* 
* 6   Return Values, Error Indicators and Warnings
* 
* 7   Accuracy
* 
* 8   Further Comments
* 
* 9   Functions referenced by subdivbs() are:
*     boehm_curve()
*     copyvector()
*     find()
*     vectalloc()
* 
* 10  Functions that reference subdivbs() are: None
****************************************************************************/

void subdivbs(ParCurv *egeom1, ParCurv *egeom2, double tnew)
/* subdivide a b-spline curve to two b-spline curves */
{
  int i,j,k,r;

  /* find the position to insert the new knot tnew.
   * here i is the index such that knots[i]<=tnew */
  i= find(egeom1->ncontpts,egeom1->knots,tnew);
  r=egeom1->order;
  j=0;
  /* if there are knots equal to tnew in the knot vector, count its
   * multiplicity. */
  while (fabs(egeom1->knots[i+j]-tnew) < 1.e-10){
      r--; j++;
  }   
  /* after the loop, r represents the number of times that tnew is to be
   * inserted. */

  /* insert tnew, r times, in egeom using boehm's knot insertion algorithm */
  for (k=0; k<r; k++)
    boehm_curve(egeom1,tnew,1);

  /* break egeom to two different curves */
  /* set the order of curve 2 */
  egeom2->order = egeom1->order;
  /* if insert tnew r times and r equals to the order of the curve,
     then the knot vector of curve 2 starts from the one next to the i-th
     * index; otherwise, it starts from the i-th knot. */
  if (r == egeom1->order)
     i++;
  /* set the knot vector of curve 2 */
  k=0;
  for (j=i; j<egeom1->order+egeom1->ncontpts; j++){
    egeom2->knots[k]=egeom1->knots[j];
    k++;
  }
  /* set the number of control points of curve 2 */
  egeom2->ncontpts = k-egeom2->order;       /* k is the number of knots */
  /* set the control points of curve 2 */
  k=1;
  for (j=egeom2->ncontpts-1; j>=0; j--){
    if (egeom2->contpts[j] == NULL)
      egeom2->contpts[j]=vectalloc();
    copyvector(egeom1->contpts[egeom1->ncontpts-k],egeom2->contpts[j]);
    k++;
  }
  /* set the number of control points of curve 1 */
  egeom1->ncontpts = i;
}

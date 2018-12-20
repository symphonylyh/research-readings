/* ***************************************************************************
 Copyright (C) 1996 Massachusetts Institute of Technology all rights reserved 
	Programmer: George A. Kriezis
**************************************************************************** */
#include <stdio.h>
#include <malloc.h>
#include "gen.h"
#include "bspl.h"

/*****************************************************************************
*                                boehm_curve()
******************************************************************************
* 1   Purpose
*     This function performs Boehm's knot insertion algorithm on a curve. It 
*     works for both the integral and rational curves.
* 
* 2   Specification
*     #include "bspl.h"
*     void boehm_curve(ParCurv *egeom, double tnew, int r)
* 
* 3   Description
*     This routine applies the Boehm knot insertion algorithm to a NURBS curve.
*     It adds a knot of value tnew, r times to the knot vector. It works for both
*     integral and rational curves.
* 
* 4   References
*     [1] W. Boehm, Inserting New Knots Into B-Spline Curves, CAD, 
*         (12)4:199-201, 1980.
* 
* 5   Parameters
*        1.ParCurv * egeom
*          On entry: Original curve in NURBS curve data format.
*          On exit: New curve with knots added in NURBS curve data format.
*        2.double tnew
*          On entry: The value of the knot to insert.
*        3.int r
*          On entry: number of times to insert knot ( spline order).
* 
* 6   Return Values, Error Indicators and Warnings
*     Not applicable.
* 
* 7   Accuracy
*     Not applicable.
* 
* 8   Further Comments
*     The curve structure contains on input the curve to be refined and on 
*     output the new refined curve.  Adequate allocation of memory should be 
*     accounted in advance for this structure before entry in the routine.
* 
* 9   Functions referenced by boehm_curv() are:
*     copyvector()
*     dbl_array1()
*     free_darray1()
*     free_varray1()
*     vectalloc()
*     vec_array1()
* 
* 10  Functions that reference boehm_curve() are:
*     subdivbs()
* 
******************************************************************************/

boehm_curve(ParCurv *egeom, double tnew, int r)
{
  int j,i,l,nknots,li,k;   /* i,j -- loop indices
			    * l -- the index of the position to insert the
			    *      new knot
			    * li -- index shift
			    * k -- order of the curve
			    * nknots -- number of knots  */
  double *knots,ai,ai1;    /* knots -- the knot vector
                            * ai,ai1 -- factors in the knot insertion
			    *           formula */
  vector *v,*v2,**a;       /* vectros and array of vectors */

  /* allocate memory */
  knots = dbl_array1((unsigned)(egeom->order+egeom->ncontpts+r)); 
  a = vec_array1((unsigned)(egeom->ncontpts+r));

  /* insert tnew r times into the knot vector. knots[] temporarily holds the 
     knot vector */
  nknots = egeom->order + egeom->ncontpts;
  li=0;      /* li indicates the distance that a knot is to be shifted after
	      * knot insertion. */
  for (i=0; i<nknots; i++){
      /* before the place of insertion, no need of shift, li = 0;
       * after the place of insertion, li = r. */
      knots[i+li]=egeom->knots[i];
      /* find the position of insertion and insert tnew r times. */
      if (egeom->knots[i]<= tnew && tnew < egeom->knots[i+1]){
         l = i;   /* remember the position of insertion for calculating 
		   * control points later. */
         for (j=1; j<r+1; j++){
	     knots[i+j]=tnew;
	     li=r;   /* all the following knots are to be shifted r
		      * positions */
             }
         }
      }

  /* allocate memory */
  v = vectalloc();
  v2 = vectalloc();

  k=egeom->order;
  /* the control points with indices l-k+1, ..., l, are temporarily copyed to
   * array a[i]. */     
  for (i=l-k+1; i<=l; i++)
      copyvector(egeom->contpts[i],a[i]);  

  /* calculate the control points after
   * inserting the new knot r times. */
  for (j=1; j<=r; j++){  /* inserting the new knot r times. */
      copyvector(a[l-k+j],v2);
      for (i=l-k+j+1; i<l+1; i++){
          ai = (tnew-egeom->knots[i])/(egeom->knots[i+k-1]-egeom->knots[i]);
          ai1 = 1.0-ai;
          /* evaluate control points */
          v->x = ai1*a[i-1]->x;
          v->y = ai1*a[i-1]->y;
          v->z = ai1*a[i-1]->z;
          v->w = ai1*a[i-1]->w;
          copyvector(v2,a[i-1]);
          v2->x = ai*a[i]->x + v->x;
          v2->y = ai*a[i]->y + v->y;
          v2->z = ai*a[i]->z + v->z;
          v2->w = ai*a[i]->w + v->w;
          }
      for (i=j; i>=1; i--)
          copyvector(a[i+l-1],a[i+l]);
      copyvector(v2,a[l]);
  }

  /* set the new knot vector */
  for (i=0; i<nknots+r; i++)
      egeom->knots[i]=knots[i];
  /* allocate memory for r new control points due to knot insertion if
   * necessary. */
  for (i=0; i<r; i++)
      if (egeom->contpts[egeom->ncontpts+i] == NULL)
         egeom->contpts[egeom->ncontpts+i]=vectalloc();
  
  /* set the new control points */
  for (i=egeom->ncontpts-1; i>=l; i--)
      copyvector(egeom->contpts[i],egeom->contpts[i+r]);
  for (i=l-k+2; i<l+r; i++)
      copyvector(a[i],egeom->contpts[i]);
  egeom->ncontpts += r;

/* free memory */
  free((char *)v);  free((char *)v2);
  free_varray1(a,(unsigned)egeom->ncontpts);
  free_darray1(knots);
}

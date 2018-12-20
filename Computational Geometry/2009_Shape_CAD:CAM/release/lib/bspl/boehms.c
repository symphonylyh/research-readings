/* ***************************************************************************
 Copyright (C) 1996 Massachusetts Institute of Technology all rights reserved 
	Programmer: George A. Kriezis
**************************************************************************** */
#include <stdio.h>
#include <malloc.h>
#include "bspl.h"

/****************************************************************************
*                          boehm_surface()
*****************************************************************************
* 
* 1   Purpose
*     This function performs Boehm's knot insertion algorithm on a NURBS
*     surface.
* 
* 2   Specification
*     #include "bspl.h"
*     void boehm_surface(ParSurf *fgeom, double tnew, int r, char dir)
* 
* 3   Description
*     This routine applies the Boehm knot insertion algorithm to a NURBS
*     surface. It adds a knot of value tnew, r times in one of the parametric 
*     directions of the surface. It works fo rboth integral and rational
*     B-spline surfaces.
* 
* 4   References
*     [1] W. Boehm, Inserting New Knots Into B-Spline Curves, CAD, 
*     (12)4:199-201, 1980.
* 
* 5   Parameters
*        1.ParSurf * fgeom
*          On entry: Original surface in NURBS surface data format.
*          On exit: New surface with knots added in NURBS surface data format.
*        2.double tnew
*          On entry: The value of the knot to insert.
*        3.int r
*          On entry: number of times to insert knot ( spline order).
*        4.char dir
*          On entry: parametric direction to insert knot, either 'u' or 'v'.
* 
* 6   Return Values, Error Indicators and Warnings
*     Not applicable.
* 
* 7   Accuracy
*     Not applicable.
* 
* 8   Further Comments
*     The surface structure contains on input the surface to be refined and on 
*     output the new refined surface.  Adequate allocation of memory should be 
*     accounted in advance for this structure before entry in the routine.
* 
* 9   Functions referenced by boehm_surface() are:
*     copyvector()
*     dbl_array1()
*     free_darray1()
*     free_varray1()
*     vectalloc()
*     vec_array1()
* 
* 10  Functions that reference boehm_surface() are:
*     subdivids()
* 
*****************************************************************************/

void boehm_surface(ParSurf *fgeom, double tnew, int r, char dir)
{
  int jj,j,i;          /* indices */
  int l,nknots,li,k;   /* l -- the position of insertion 
			* nknots -- number of knots
			* li -- the index shift of knots in knot vector
			*       due to knot insertion
			* k -- order           */
  double *knots,ai,ai1;   /* knot -- knot vector
			   * ai,ai1 -- the factors in the insertion formula */
  struct vector *v,*v2,**a;  /* vectors and array of vectors, to be used to
			      * temporarily hold control points   */

  /* initialize the index shift as 0 */
  li=0;

  /* knot insertion */
  /* if the insertion is in u direction */
  if (dir == 'u'){
     nknots = fgeom->uorder + fgeom->ucontpts;
     knots = dbl_array1((unsigned)(nknots+r));
     /* insert the new knot r times into the knot vector in u */
     for (i=0; i<nknots; i++){
         knots[i+li]=fgeom->uknots[i];
         if (fgeom->uknots[i]< tnew && tnew <= fgeom->uknots[i+1]){
	    l = i;
	    for (j=1; j<r+1; j++){
	        knots[i+j]=tnew;
	        li=r;
	        }
            }
         }
    a = vec_array1((unsigned)(fgeom->ucontpts+r));
    k=fgeom->uorder;
    }
  /* if the insertion is in v direction */
  else if (dir == 'v'){
    nknots = fgeom->vorder + fgeom->vcontpts;
    knots = dbl_array1((unsigned)(nknots+r));
    /* insert the new knot r times into the knot vector in v */
    for (i=0; i<nknots; i++){
        knots[i+li]=fgeom->vknots[i];
        if (fgeom->vknots[i]< tnew && tnew <= fgeom->vknots[i+1]){
	   l = i;
	   for (j=1; j<r+1; j++){
	       knots[i+j]=tnew;
	       li=r;
	       }
           }
        }
    a = vec_array1((unsigned)(fgeom->vcontpts+r));
    k=fgeom->vorder;
  }

  /* allocate memory */
  v = vectalloc();
  v2 = vectalloc();

  /* calculate the new control points after knot insertion */
  /* if the insertion is in u direction */
  if (dir == 'u'){
     for (jj=0; jj<fgeom->vcontpts; jj++){
         for (i=l-k+1; i<=l; i++)
	     /* copy those control points to be changed into array a */
	     copyvector(fgeom->contpts[i][jj],a[i]);
         for (j=1; j<=r; j++){
	     copyvector(a[l-k+j],v2);
	     for (i=l-k+j+1; i<l+1; i++){
	         /* calculate the factors */
	         ai = (tnew-fgeom->uknots[i]) /
		   (fgeom->uknots[i+k-1]-fgeom->uknots[i]);
	         ai1 = 1.0-ai;
                 /*  evaluate control points */
	         v->x = a[i-1]->x*ai1;
		 v->y = a[i-1]->y*ai1;
		 v->z = a[i-1]->z*ai1;
		 v->w = a[i-1]->w*ai1;
		 copyvector(v2,a[i-1]);
		 v2->x = ai*a[i]->x +v->x;
		 v2->y = ai*a[i]->y +v->y;
		 v2->z = ai*a[i]->z +v->z;
		 v2->w = ai*a[i]->w +v->w;
	        }
             for (i=j; i>=1; i--)
	         copyvector(a[i+l-1],a[i+l]);
	     copyvector(v2,a[l]);
            }
            /* set the new control points */
	    for (i=0; i<r; i++)
	        fgeom->contpts[fgeom->ucontpts+i][jj]=vectalloc();
            for (i=fgeom->ucontpts-1; i>=l; i--)
        	copyvector(fgeom->contpts[i][jj],fgeom->contpts[i+r][jj]);
            for (i=l-k+2; i<l+r; i++)
	        copyvector(a[i],fgeom->contpts[i][jj]);
        }
        /* set the new u knot vector and number of control points in u */
        for (i=0; i<nknots+r; i++)
            fgeom->uknots[i]=knots[i];
        fgeom->ucontpts += r;
  }
  /* if the insertion is in v direction */
  else if (dir == 'v'){
     for (jj=0; jj<fgeom->ucontpts; jj++){
         for (i=l-k+1; i<=l; i++)
	     /* copy the control points to be changed into array a */
	     copyvector(fgeom->contpts[jj][i],a[i]);
         for (j=1; j<=r; j++){
	     copyvector(a[l-k+j],v2);
	     for (i=l-k+j+1; i<l+1; i++){
	         /* calculate the factors */
	         ai = (tnew-fgeom->vknots[i]) /
		   (fgeom->vknots[i+k-1]-fgeom->vknots[i]);
	         ai1 = 1.0-ai;
		 /* calculate the new control points */
                 v->x = a[i-1]->x*ai1;
	         v->y = a[i-1]->y*ai1;
		 v->z = a[i-1]->z*ai1;
		 v->w = a[i-1]->w*ai1;
		 copyvector(v2,a[i-1]);
		 v2->x = ai*a[i]->x +v->x;
		 v2->y = ai*a[i]->y +v->y;
		 v2->z = ai*a[i]->z +v->z;
		 v2->w = ai*a[i]->w +v->w;
	        }
	    for (i=j; i>=1; i--)
	        copyvector(a[i+l-1],a[i+l]);
	    copyvector(v2,a[l]);
           }
        /* set the new control points */
        for (i=0; i<r; i++)
	    fgeom->contpts[jj][fgeom->vcontpts+i]=vectalloc();
        for (i=fgeom->vcontpts-1; i>=l; i--)
	    copyvector(fgeom->contpts[jj][i],fgeom->contpts[jj][i+r]);
        for (i=l-k+2; i<l+r; i++)
	    copyvector(a[i],fgeom->contpts[jj][i]);
       }
    /* set the new v knot vector */
    for (i=0; i<nknots+r; i++)
        fgeom->vknots[i]=knots[i];
    fgeom->vcontpts += r;
  }

  /* free memory */
  free((char *)v);  free((char *)v2);
  free_darray1(knots);
  if (dir == 'u'){
    free_varray1(a,(unsigned)fgeom->ucontpts);
  }
  else if (dir == 'v'){
    free_varray1(a,(unsigned)fgeom->vcontpts);
  }
}

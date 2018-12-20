/* ***************************************************************************
 Copyright (C) 1996 Massachusetts Institute of Technology all rights reserved 
	Programmer: George A. Kriezis
**************************************************************************** */
#include "gen.h"
#include "bspl.h"

/* NAG routine, which rearranges a vector of real numbers into ascending
 * or desecending order, see NAG manual, Vol. 8. */
void m01caf_(double *, int *, int *, char *, int *);   

/*****************************************************************************
*                               convexbox()
******************************************************************************
* 
* 1   Purpose
*     Calculate the convexbox of a NURBS surface.
* 
* 2   Specification
*     #include "bspl.h"
*     void convexbos(ParSurf *fgeom, double *xa, double *xb, double *ya, 
*                    double *yb, double *za, double *zb)
* 
* 3   Description
*     The convexbox of a NURBS surface is a rectangular parallelipiped
*     defined by the minimum and the maximum of x,y,z coordinates of the
*     control points of the surface.
* 
* 4   References
*     Not applicable.
* 
* 5   Parameters
*       1.ParCurv * fgeom
*         On entry: the address of a NURBS surface structure.
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
* 9   Functions referenced by convexbox() are:
*     dbl_array1()
*     free_darray1()
*     m01caf_()
* 
* 10  Functions that reference convexbox() are:
*  
******************************************************************************/

void convexbox(ParSurf *fgeom, double *xa, double *xb, double *ya, double *yb,
	       double *za, double *zb)
{
  int i,j;          /* loop indices */
  vector *vect;     /* vector holding x,y,z coordinates of each control
		     * point */
  double *x,*y,*z;  /* arrays holding x,y,z coordinates of control points 
		     * respectively */
  int m1,m2,n,ifail;  /* m1,m2 -- indicates the beginning and the end of the
		       *          array 
		       * n -- length of x, y, z arrays 
		       * ifail -- indicates the failure of calling NAG
		       *          routine */
  char order;         /* a -- indicate if the array is ascending */

  /* allocate memory for x,y,z arrays */
  x = dbl_array1((unsigned)(fgeom->ucontpts*fgeom->vcontpts));
  y = dbl_array1((unsigned)(fgeom->ucontpts*fgeom->vcontpts));
  z = dbl_array1((unsigned)(fgeom->ucontpts*fgeom->vcontpts));

  n = 0;
  /* put the coordinates of control points into the corresponding array */
  for(i=0; i<fgeom->ucontpts; i++)
    for(j=0; j<fgeom->vcontpts; j++)
      {
	vect = fgeom->contpts[i][j];
	
	x[n  ] = vect->x/vect->w;
	y[n  ] = vect->y/vect->w;
	z[n++] = vect->z/vect->w;
      }
  
  /* set up for calling NAG routine */
  m1 = 1;
  m2 = n;
  order = 'a';
  ifail = 0;
  
  /* sort x coordinates and get Xmin and Xmax */
  m01caf_(x, &m1, &m2, &order, &ifail);  
  *xa = x[0];
  *xb = x[n-1];

  /* sort y coordinates and get Ymin and Ymax */  
  m01caf_(y, &m1, &m2, &order, &ifail);  
  *ya = y[0];
  *yb = y[n-1];
  
  /* sort z coordinates and get Zmin and Zmax */
  m01caf_(z, &m1, &m2, &order, &ifail);  
  *za = z[0];
  *zb = z[n-1];

  /* free memory */
  free_darray1(x);
  free_darray1(y);
  free_darray1(z);

}

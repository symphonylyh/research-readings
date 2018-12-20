/* ***************************************************************************
 Copyright (C) 1996 Massachusetts Institute of Technology all rights reserved 
	Programmer: George A. Kriezis
**************************************************************************** */
#include <stdio.h>
#include <malloc.h>
#include "bspl.h"

/****************************************************************************
*                                 hodograph_surf()
*****************************************************************************
* 
* 1   Purpose
*     This function calculates the hodograph of a NURBS surface.
* 
* 2   Specification
*     #include "bspl.h"
*     void hodograph_surf(ParSurf *fgeom, ParSurf *hodog, char dir)
* 
* 3   Description
*     This routine calculates the hodograph (1st parametric derivative 
*     surface) of a NURBS surface. The NURBS hodograph is approximated using 
*     the integral hodograph as suggested by Sederberg.
* 
* 4   References
*     [1] T. W. Sederberg and X. Wang.  Rational Hodographs, Computer Aided 
*         Geometric Design, 5(4):333-335, 1988
* 
* 5   Parameters
*        1.ParSurf * fgeom
*          On entry: NURBS data structure containing geometry of the original 
* 	           surface.
*        2.ParSurf * hodog
*          On entry: NURBS data structure allocated for a surface of the same 
* 	           order and number of control points as the original surface.
*          On exit: NURBS data structure containing geometry of the surface 
* 	          hodograph.
*        3.char dir
*          On entry: direction of hodograph computation 'u' or 'v'
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
* 9   Functions referenced by hodograph_surf() are:
*     scale_vect1()
*     sub_vect1()
*     vectalloc()
* 
* 10 Functions that reference hodograph_surf() are: None
* 
*******************************************************************************/

void hodograph_surf(ParSurf *fgeom, ParSurf *hodog, char dir)
/* this routine calculates the hodograph of a b-spline surface */
/* calculates the integral hodograph of a rational b-spline using
 * sederberg's approximation - w = 1 */
/* dir is direction take hodograph */
{
  int i,k;            /* loop indices */
  int n,nm;           /* length of the u and the v knot vector,
		       * respectively */ 
  double prod;        /* the factor which is an inverse of the parametric
		       * difference between two knots */
  struct vector *v;   /* a vector */
  v = vectalloc();

  if (dir == 'u'){  /* hodograph w.r.t u */
    n = fgeom->uorder+fgeom->ucontpts;
    hodog->uorder = fgeom->uorder - 1;  /* order decreases by 1 in u */
    hodog->ucontpts = fgeom->ucontpts - 1;  /* number of control points 
					     * decreases by 1 in u */
    hodog->vorder = fgeom->vorder;      /* same order in v */
    hodog->vcontpts = fgeom->vcontpts;  /* same number of control points in
					 * v */
    nm = fgeom->vorder+fgeom->vcontpts;
    for (i=0; i<nm; i++)  /* same knot vector in v */
      hodog->vknots[i]=fgeom->vknots[i];
    for (i=1; i<n-1; i++) /* new knot vector in u */
      hodog->uknots[i-1] = fgeom->uknots[i];
    /* allocate memory for the new control points */
    for (i=0; i<hodog->ucontpts; i++)
      for (k=0; k<hodog->vcontpts; k++)
	if (hodog->contpts[i][k] == NULL)
	  hodog->contpts[i][k] = vectalloc();
    /* calculate the new control points */
    n = fgeom->uorder -1;
    for (i=1; i<fgeom->ucontpts; i++){
      prod = 1./(fgeom->uknots[i+n] - fgeom->uknots[i]);
      for (k=0; k<fgeom->vcontpts; k++){
	sub_vect1(fgeom->contpts[i][k],fgeom->contpts[i-1][k],v);
	scale_vect1(prod,v,hodog->contpts[i-1][k]);
      }
    }
  }
  else if (dir == 'v'){  /* hodograph w.r.t. v */
    n = fgeom->vorder+fgeom->vcontpts;
    hodog->vorder = fgeom->vorder - 1;      /* order decreases by 1 in v */
    hodog->vcontpts = fgeom->vcontpts - 1;  /* number of control points
					     * decreases by 1 in v */
    hodog->uorder = fgeom->uorder;      /* same order in u */
    hodog->ucontpts = fgeom->ucontpts;  /* same number of control points in
					 * u */
    nm = fgeom->uorder+fgeom->ucontpts;
    for (i=0; i<nm; i++)   /* same knot vector in u */
      hodog->uknots[i]=fgeom->uknots[i];
    for (i=1; i<n-1; i++)  /* new knot vector in v */
      hodog->vknots[i-1] = fgeom->vknots[i];
    /* allocate memory for the control points of hodograph */
    for (i=0; i<hodog->ucontpts; i++)
      for (k=0; k<hodog->vcontpts; k++)
	if (hodog->contpts[i][k] == NULL)
	  hodog->contpts[i][k] = vectalloc();
    /* calculate the control points of hodograph */
    n = fgeom->vorder -1;
    for (i=1; i<fgeom->vcontpts; i++){
      prod = 1./(fgeom->vknots[i+n] - fgeom->vknots[i]);
      for (k=0; k<fgeom->ucontpts; k++){
	sub_vect1(fgeom->contpts[k][i],fgeom->contpts[k][i-1],v);
	scale_vect1(prod,v,hodog->contpts[k][i-1]);
      }
    }
  }

  free((char *)v);
}

/****************************************************************************
*                               ParCurv_hodograph()
*****************************************************************************
* 
* 1   Purpose
*     This function calculates the hodograph of an integral B-spline curve.
* 
* 2   Specification
*     #include "bspl.h"
*     ParCurv *ParCurv_hodograph(ParCurv *egm)
* 
* 3   Description
*     This function returns the hodograph, 1st parametric derivative curve,
*     of an integral B-spline curve. 
* 
* 4   References
*     Not applicable
* 
* 5   Parameters
*        1.ParCurv * egm
*          On entry: data structure containing an integral B-spline curve.
*        
* 6   Return Values, Error Indicators and Warnings
*     The NURBS data structure containing a B-spline curve as the hodograph of
*     the input B-spline curve
* 
* 7   Accuracy
*     Not applicable.
* 
* 8   Further Comments
*     Not applicable.
* 
* 9   Functions referenced by ParCurv_hodograph() are:
*     copyvector()
*     egeomalloc1()
*     scale_vect1()
*     sub_vect1()
*     vectalloc()
* 
* 10 Functions that reference ParCurv_hodograph() are: None
* 
*******************************************************************************/

ParCurv *ParCurv_hodograph(ParCurv *egm)
{
  int j,i;             /* i -- loop index
			* j -- never used (glshen) */
  vector *vect;   /* working variable */
  ParCurv *hodo;  /* a NURBS data structure which will contain the hodograph
		   * of the B-spline curve. */

  /* allocate memory for hodo */
  hodo = egeomalloc1(egm->order-1,egm->ncontpts-1);

  /* set the knot vector of the hodograph as the knot vector of the original
   * curve but with multiplicity egm->order-1 at two ends. */
  for (i=0;i<hodo->order+hodo->ncontpts;i++)
      hodo->knots[i] = egm->knots[i+1];

  vect = vectalloc();

  /* calculate the control points of the hodograph */
  for(i=1;i<egm->ncontpts;i++) {
    sub_vect1(egm->contpts[i],egm->contpts[i-1],vect);
    scale_vect1((egm->order-1)/(egm->knots[i+egm->order-1] - egm->knots[i]),
		vect,vect);
    copyvector(vect,hodo->contpts[i-1]);
  }

  free((char *)vect);

  return(hodo);
}
















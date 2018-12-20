/************************************************************************
 *									*
			Copyright (C) 1996 by
	Massachusetts Institute of Technology, Cambridge, MA
			 All rights reserved

 *									*
 ************************************************************************/

#include <stdio.h>
#include <math.h>
#include "bspl.h"

static double ***derivpts;   /* global variable containing values of
			      * derivatives of a periodic B-spline surface */

static int NUMBERPOINTS = 0; /* this is so that memory is only allocated
                              * a single time. That is, when 0, memory is
                              * allocted and the variable set to the length 
                              * of the array derivpts. Subsequent calls will
                              * check this against the current required
                              * length, if greater, will  reallocate derivpts
                              * if less, will not. */

/****************************************************************************
*                                 evalrsurf_per()
*****************************************************************************
* 
* 1   Purpose
*     This function places the point and derivatives of a periodic NURBS
*     surface into an array of vectors.
* 
* 2   Specification
*     #include "bspl.h"
*     void evalrsurf_per(ParSurf *fgeom, double u, double v, int ider, 
*                        vector **vec)
* 
* 3   Description
*     This routine places into a supplied array of vector data structures the 
*     point and derivatives of the point on a B-spline surface corresponding 
*     to the parameters (u,v).
* 
* 4   References
*     [1]C. De Boor. A Practical Guide to Splines, Springer, New York, 1978.
* 
* 5   Parameters
*       1.ParSurf * fgeom
*         On entry: NURBS surface data structure containing geometry of 
* 	          surface to be evaluated.
*       2.double u
*         On entry: the parameter value at which the non-uniform integral 
* 	          B-spline surface is to be evaluated.
*       3.double v
*         On entry: the parameter value at which the non-uniform integral 
* 	          B-spline surface is to be evaluated.
*       4.int ider
*         On entry: an index referring to the required derivative with respect 
* 	          to u to be calculated such that if the surface is defined as 
* 		  R(u,v) then the routine will evaluate all derivatives up to 
* 		  and including:
* 		                      ider
*                                      @    R(u,v)
* 				  ------------------
*                                      uderiv  vderiv
* 				   @u      @v
* 		  ( @ represents partial derivative. )
* 		  where uderiv and vderiv are defined as:
* 		     uderiv = 0 and vderiv = 0 when ider = 0.
* 		     uderiv = 1 and vderiv = 0 when ider = 1.
* 		     uderiv = 0 and vderiv = 1 when ider = 2.
* 		     uderiv = 1 and vderiv = 1 when ider = 3.
* 		     uderiv = 2 and vderiv = 0 when ider = 4.
* 		     uderiv = 0 and vderiv = 2 when ider = 5.
* 		     uderiv = 3 and vderiv = 0 when ider = 6.
* 		     uderiv = 2 and vderiv = 1 when ider = 7.
* 		     uderiv = 1 and vderiv = 2 when ider = 8.
* 		     uderiv = 0 and vderiv = 3 when ider = 9.
*        5.vect ** vec
*          On entry: an array of pointers of length  ider.
*          On exit: the vector values of the surface evaluated at the given 
* 	          data points. The ith element of the array corresponds to:
* 		                       i
*                                       @ R(u,v)
* 				  ------------------
*                                      uderiv  vderiv
* 				   @u      @v
* 		  ( @ represents partial derivative. )
* 		  where uderiv and vderiv are defined as:
* 		     uderiv = 0 and vderiv = 0 when i = 0.
* 		     uderiv = 1 and vderiv = 0 when i = 1.
* 		     uderiv = 0 and vderiv = 1 when i = 2.
* 		     uderiv = 1 and vderiv = 1 when i = 3.
* 		     uderiv = 2 and vderiv = 0 when i = 4.
* 		     uderiv = 0 and vderiv = 2 when i = 5.
* 		     uderiv = 3 and vderiv = 0 when i = 6.
* 		     uderiv = 2 and vderiv = 1 when i = 7.
* 		     uderiv = 1 and vderiv = 2 when i = 8.
* 		     uderiv = 0 and vderiv = 3 when i = 9.
* 
* 6   Return Values, Error Indicators and Warnings
*     The index uderiv <= the order with respect to u of the non-uniform 
*     integral B-spline surface.
*     The index vderiv <= the order with respect to u of the non-uniform 
*     integral B-spline surface.
*     The value of u should lie between fgeom -> uknots[0] and fgeom ->
*     uknots[k] where k = (fgeom -> ucontpts) + (fgeom -> uorder) - 1.
*     The value of v should lie between fgeom -> vknots[0] and fgeom ->
*     vknots[j]  where j = (fgeom -> vcontpts) + (fgeom -> vorder) - 1.
* 
* 7   Accuracy
*     Not applicable.
* 
* 8   Further Comments
*     This basis function evaluation is performed using routine nbasisd (bspl).
*     This is a modified and more efficient version of subroutine 
*     revalderivsurf_per (bspl _ eval-surf.c).
*     See also eval_surface_bounded (bspl _ evalsurf.c).
* 
* 9   Functions referenced by evalrsurf_per() are:
*     dbl_array2()
*     evaldersurf_per()
*     evalsurf1_per()
*     find()
*     free_darray2()
*     free_iarray1()
*     int_array1()
*     nbasisd()
*     nbasisd_per()
* 
* 10  Functions that reference evalrsurf_per() are:
*     offset_geod()
*     offset_geod_par()
*     offset_normal()
* 
***************************************************************************/

void evalrsurf_per(ParSurf *fgeom, double u, double v, int fin, vector *vec[])
{
  vector *r0,*r1u,*r1v,*r2u,*r2v,*r2uv,*r3u,*r3v,*ruuv,*ruvv;
         /* the partial derivatives. here the numbers represent the orders, 
	  * u,v mean taking derivatives w.r.t. u and v */
  double q,r;     /* periodic B-spline basis */
  double **N1,**M1;
         /* matrices containing periodic B-spline basis in u,v directions */
  double w,w1u,w1v,w2u,w2v,w2uv,w3u,w3v,wuuv,wuvv;
         /* scaling factors in the transformation from homogeneous to cartesian
	  * coordinates */
  int i,j;  /* i -- loop index
	       j -- never used (glshen) */
  int ucontpts,vcontpts,uorder,vorder;   /* order and number of control
					  * points */
  int ilv,*P1;   /* ilv -- index
                  * P1 -- integer array containing indices of the control
		  * points whose associated B-spline basis are nonzero   */

  /* assignment and memory allocation */
  ucontpts = fgeom->ucontpts;  vcontpts = fgeom->vcontpts;
  uorder = fgeom->uorder;  vorder = fgeom->vorder;
  N1 = dbl_array2((unsigned)(uorder+1),(unsigned)(uorder+1));
  M1 = dbl_array2((unsigned)(vorder+1),(unsigned)(vorder+1));
  P1 = int_array1((unsigned)(uorder+1));

  /* calculate the periodic B-spline basis matrix in u */
  nbasisd_per(uorder,ucontpts,u,fgeom->uknots,N1,P1);
  /* find the index ilv, such that vknots[ilv] <= v <= vknots[ilv+1] */
  ilv = find(vcontpts,fgeom->vknots,v);
  /* calculate the periodic B-spline basis matrix in v */
  nbasisd(ilv, vorder, v, fgeom->vknots, M1); 
  
  for (i=0; i<=fin; i++){
    if(i==0) {    /* evaluate the position at u,v */
      r0 = evalsurf1_per(fgeom,ilv,N1,M1,P1);
      w = r0->w;
      /* transform from homogeneous to cartesian coordinates */
      r0->x = r0->x/w;
      r0->y = r0->y/w;
      r0->z = r0->z/w;
      r0->w = 1.0;
      vec[i]=r0;
    }
    else if(i == 1) { /* evaluate the 1st partial derivative w.r.t. u */
      r1u = evaldersurf_per(fgeom,1,0,ilv,N1,M1,P1);
      w1u = r1u->w;
      /* transform from homogeneous to cartesian coordinates */ 
      r1u->x = (r1u->x - r0->x*w1u)/w;
      r1u->y = (r1u->y - r0->y*w1u)/w;
      r1u->z = (r1u->z - r0->z*w1u)/w;
      r1u->w = 1.0;
      vec[i]=r1u;
    }
    else if(i==2) {  /* evaluate the 1st partial derivative w.r.t. v */	
      r1v = evaldersurf_per(fgeom,0,1,ilv,N1,M1,P1);
      w1v = r1v->w;
      /* transform from homogeneous to cartesian coordinates */
      r1v->x = (r1v->x - r0->x*w1v)/w;
      r1v->y = (r1v->y - r0->y*w1v)/w;
      r1v->z = (r1v->z - r0->z*w1v)/w;
      r1v->w = 1.0;
      vec[i]=r1v;
    }
    else if(i == 3) { /* evaluate the 2nd mixed partial derivative */
      r2uv = evaldersurf_per(fgeom,1,1,ilv,N1,M1,P1);
      w2uv = r2uv->w;
      /* transform from homogeneous to cartesian coordinates */
      r2uv->x = (r2uv->x - w1u*r1v->x - w1v*r1u->x - w2uv*r0->x)/w;
      r2uv->y = (r2uv->y - w1u*r1v->y - w1v*r1u->y - w2uv*r0->y)/w;
      r2uv->z = (r2uv->z - w1u*r1v->z - w1v*r1u->z - w2uv*r0->z)/w;
      r2uv->w = 1.0;
      vec[i]=r2uv;
    }
    else if(i == 4) { /* evaluate the 2nd partial derivative w.r.t. u */
      r2u = evaldersurf_per(fgeom,2,0,ilv,N1,M1,P1);
      w2u = r2u->w;
      /* transform from homogeneous to cartesian coordinates */
      r2u->x = (r2u->x - 2.0*w1u*r1u->x - w2u*r0->x)/w;
      r2u->y = (r2u->y - 2.0*w1u*r1u->y - w2u*r0->y)/w;
      r2u->z = (r2u->z - 2.0*w1u*r1u->z - w2u*r0->z)/w;
      r2u->w = 1.0;
      vec[i]=r2u;
    }
    else if(i == 5) { /* evaluate the 2nd partial derivative w.r.t. v */
      r2v = evaldersurf_per(fgeom,0,2,ilv,N1,M1,P1);
      w2v = r2v->w;
      /* transform from homogeneous to cartesian coordinates */
      r2v->x = (r2v->x - 2.0*w1v*r1v->x - w2v*r0->x)/w;
      r2v->y = (r2v->y - 2.0*w1v*r1v->y - w2v*r0->y)/w;
      r2v->z = (r2v->z - 2.0*w1v*r1v->z - w2v*r0->z)/w;
      r2v->w = 1.0;
      vec[i]=r2v;
    }
  }

  /* free memory */
  free_darray2(N1);  free_darray2(M1); free_iarray1(P1);
}

/****************************************************************************
*                               evalsurf1_per()
*****************************************************************************
* 
* 1   Purpose
*     Given the periodic B-spline basis functions, this function evaluates a 
*     point on a periodic B-spline surface.
* 
* 2   Specification
*     #include "bspl.h"
*     vector *evalsurf1_per(ParSurf *fgeom, int ilv, double **N1, double **M1, 
*                           int *P1)
* 
* 3   Description
*     Here only in v direction the index ilv such that vknot[ilv] <= v <= 
*     vknot[ilv+1], and only in u direction the P array containing the 
*     indices of the control points with nonzero associated basis functions,
*     are given.
* 
* 4   References
*     Not applicable
* 
* 5   Parameters
*     1.ParSurf *fgeom
*       On entry: the periodic B-spline surface to be evaluated.
*     2.int ilv
*       On entry: the index such that vknot[ilv] <= v <= vknot[ilv+1]
*     3.double **N1
*       On entry: the B-spline basis functions in u direction
*     4.double **M1
*       On entry: the B-spline basis functions in v direction
*     5.int *P1
*       On entry: the array containing the indices of the control points whose
*                 associated basis functions are nonzero.
* 
* 6   Return Values, Error Indicators and Warnings
*     This function returns the evaluated vector.
* 
* 7   Accuracy
*     Not applicable.
* 
* 8   Further comments
*     Not appicable
* 
* 9   Functions referenced by evalsurf1_per() are:
*     vectalloc()
* 
* 10  Functions that reference evalsurf1_per() are:
*     evalrsurf_per()
* 
****************************************************************************/

vector *evalsurf1_per(ParSurf *fgeom, int ilv, double **N1, double **M1,
		      int *P1)
{
  int ucontpts,uorder,vcontpts,vorder;   /* number of control points and order
					  * in u,v directions, respectively. */
  int i,j;                               /* loop indices */
  vector *p,*eval;               /* p -- to be used to contain a control point
				  * eval -- the evaluated point  */
  double q,r;                    /* values of basis functions */

  ucontpts = fgeom->ucontpts;vcontpts = fgeom->vcontpts;
  uorder = fgeom->uorder;vorder = fgeom->vorder;

  /* allocate memory and initialize */
  eval = vectalloc();
  eval->x = eval->y = eval->z = eval->w = 0.0;

  /* evaluate the point on the periodic B-spline surface */
  for(i=0; i<uorder; i++)
    for(j=ilv; j>ilv-vorder; j--) {
      p = fgeom->contpts[P1[i]][j];		
      q = N1[i+1][uorder];
      r = M1[j+vorder-ilv][vorder];
      eval->x += q*r*p->x;	
      eval->y += q*r*p->y;	
      eval->z += q*r*p->z;	
      eval->w += q*r*p->w;	
    }
  return(eval);
}

/******************************************************************************
*                                 evaldersurf_per()
*******************************************************************************
* 
* 1   Purpose
*     This function places the point and derivatives of a periodic integral 
*     NURBS surface into an array of vectors.
* 
* 2   Specification
*     #include "bspl.h"
*     void evaldersurf_per(ParSurf *fgeom, double u, double v, int ider, 
*                          vector **vec)
* 
* 3   Description
*     This routine places into a supplied array of vector data structures the 
*     point and derivatives of the point on a B-spline surface corresponding 
*     to the parameters (u,v).
* 
* 4   References
*     [1]C. De Boor. A Practical Guide to Splines, Springer, New York, 1978.
* 
* 5   Parameters
*       1.ParSurf * fgeom
*         On entry: NURBS surface data structure containing geometry of 
* 	          surface to be evaluated.
*       2.double u
*         On entry: the parameter value at which the non-uniform integral 
* 	          B-spline surface is to be evaluated.
*       3.double v
*         On entry: the parameter value at which the non-uniform integral 
* 	          B-spline surface is to be evaluated.
*       4.int ider
*         On entry: an index referring to the required derivative with respect 
* 	          to u to be calculated such that if the surface is defined as 
* 		  R(u,v) then the routine will evaluate all derivatives up to 
* 		  and including:
*                                       ider
*                                      @    R(u,v)
*                                   ------------------
*                                      uderiv  vderiv
* 				   @u      @v
* 		  ( @ represents partial derivative. )
* 		  where uderiv and vderiv are defined as:
* 		     uderiv = 0 and vderiv = 0 when ider = 0.
* 		     uderiv = 1 and vderiv = 0 when ider = 1.
* 		     uderiv = 0 and vderiv = 1 when ider = 2.
* 		     uderiv = 1 and vderiv = 1 when ider = 3.
* 		     uderiv = 2 and vderiv = 0 when ider = 4.
* 		     uderiv = 0 and vderiv = 2 when ider = 5.
* 		     uderiv = 3 and vderiv = 0 when ider = 6.
* 		     uderiv = 2 and vderiv = 1 when ider = 7.
* 		     uderiv = 1 and vderiv = 2 when ider = 8.
* 		     uderiv = 0 and vderiv = 3 when ider = 9.
*        5.vect ** vec
*          On entry: an array of pointers of length  ider.
*          On exit: the vector values of the surface evaluated at the given
* 	          data points. The ith element of the array corresponds to:
*                                        i
*                                       @ R(u,v)
*                                  ------------------
*                                     uderiv  vderiv
* 				  @u      @v
* 		  ( @ represents partial derivative. )
* 		  where uderiv and vderiv are defined as:
* 		     uderiv = 0 and vderiv = 0 when i = 0.
* 		     uderiv = 1 and vderiv = 0 when i = 1.
* 		     uderiv = 0 and vderiv = 1 when i = 2.
* 		     uderiv = 1 and vderiv = 1 when i = 3.
* 		     uderiv = 2 and vderiv = 0 when i = 4.
* 		     uderiv = 0 and vderiv = 2 when i = 5.
* 		     uderiv = 3 and vderiv = 0 when i = 6.
* 		     uderiv = 2 and vderiv = 1 when i = 7.
* 		     uderiv = 1 and vderiv = 2 when i = 8.
* 		     uderiv = 0 and vderiv = 3 when i = 9.
* 
* 6   Return Values, Error Indicators and Warnings
*     The index uderiv <= the order with respect to u of the non-uniform 
*     integral B-spline surface.
*     The index vderiv <= the order with respect to u of the non-uniform 
*     integral B-spline surface.
*     The value of u should lie between fgeom -> uknots[0] and fgeom ->
*     uknots[k]  where k = (fgeom -> ucontpts) + (fgeom -> uorder) - 1.
*     The value of v should lie between fgeom -> vknots[0] and fgeom ->
*     vknots[j]  where j = (fgeom -> vcontpts) + (fgeom -> vorder) - 1.
* 
* 7   Accuracy
*     Not applicable.
* 
* 8   Further Comments
*     This basis function evaluation is performed using routine nbasisd (bspl).
*     This is a modified and more efficient version of subroutine 
*     revalderivsurf_per (bspl _ evalsurf.c).
*     See also eval_surface_bounded (bspl _ evalsurf.c).
* 
* 9   Functions referenced by evaldersurf_per() are:
*     Aijqs_per()
*     dbl_array3()
*     free_darray3()
*     vectalloc()
* 
* 10  Functions that reference evaldersurf_per() are:
*     evalrsurf_per()
* 
***************************************************************************/

vector *evaldersurf_per(ParSurf *fgeom, int uderiv, int vderiv, int ilv,
			double **N1, double **M1, int *P1)
{
  int ucontpts,uorder,vcontpts,vorder; /* order and number of control points */
  int i,j,ii,jj;                       /* indices */
  int prod;                            /* a multiplying factor */
  vector *p,*eval;                     /* p -- to be used to contain a
                                        *      control point
					* eval -- the evaluated value */
  double q,r;                          /* periodic B-spline basis functions */
  static int flag = FALSE;             /* never used (glshen) */

  /* assignment */
  ucontpts = fgeom->ucontpts;vcontpts = fgeom->vcontpts;
  uorder = fgeom->uorder;vorder = fgeom->vorder;

  /* allocate memory, see the explaination at the head of the file */
  if (NUMBERPOINTS == 0){
    derivpts = dbl_array3((unsigned)(uorder*vorder),(unsigned)(uorder*vorder),
			  4);
    NUMBERPOINTS = uorder*vorder;
  }
  else if (NUMBERPOINTS < uorder*vorder){
    free_darray3(derivpts);
    derivpts = dbl_array3((unsigned)(uorder*vorder),(unsigned)(uorder*vorder),
			  4);
    NUMBERPOINTS = uorder*vorder;
  }

  /* calculate the control points of the hodograph. here some factors are not
     considered yet. */
  Aijqs_per(uderiv,uorder,ucontpts,P1,fgeom->uknots,ilv,vderiv,vorder,
	    fgeom->vknots,fgeom->contpts, derivpts);

  /* allocate memory for the evaluated value and initialize it. */
  eval = vectalloc();
  eval->x = eval->y = eval->z = eval->w = 0.0;
  
  /* calculate the factor (uorder-1)...(uorder-uderiv)(vorder-1)...
   * (vorder-vderiv) */
  prod = 1.0;
  for(i=1; i<=uderiv; i++)
    prod *= uorder-i;
  for(i=1; i<=vderiv; i++)
    prod *= vorder-i;

  /* calculate the derivative at u,v */
  for (i=0; i<uorder-uderiv; i++)
      for (j=ilv; j>ilv-vorder+vderiv; j--) {
          q = N1[i+1][uorder-uderiv];
          r = M1[j+vorder-vderiv-ilv][vorder-vderiv];
          /* index shift */
          ii = i+uderiv;
          jj = j-ilv+(vorder-1);
          ii = ii + uorder*jj;
          jj = uderiv + uorder*vderiv;
          eval->x  += q*r*derivpts[ii][jj][0];
          eval->y  += q*r*derivpts[ii][jj][1];
          eval->z  += q*r*derivpts[ii][jj][2];
          eval->w  += q*r*derivpts[ii][jj][3];
         }
  /* multiply the factor */  
  eval->x = prod*eval->x;
  eval->y = prod*eval->y;
  eval->z = prod*eval->z;
  eval->w = prod*eval->w;

  return(eval);
}

/****************************************************************************
*                                  Aijqs_per()
*****************************************************************************
* 
* 1   Purpose
*     This is a subroutine of function evalderivsurf_per(), which evaluates 
*     the derivatives of a periodic B-spline surface, see evalderivsurf_per().
*     It calculate the control points of the hodograph but without the scaling 
*     factor.
* 
* 2   Specification
*     #include "bspl.h"
*     void Aijqs_per(int q, int k, int o, int *P, double u[], int r, int s,
* 	             int l,double v[], vector ***contpts, double ***derivpts2)
* 
* 3   Description
*     This function calculates the control points of the derivative surface
*     of a periodic B-spline surface. The derivative surface itself is also a 
*     periodic B-spline surface with lower order. The control points of the 
*     derivative surface are recursively determined from the control points of 
*     the surface.
* 
* 4   References
*     [1]C. De Boor. A Practical Guide to Splines, Springer, New York, 1978.
* 
* 5   Parameters
*       1.int q
*         On entry: the order of the derivative w.r.t. u
*       2.int k
*         On entry: the order of the surface in u direciton
*       3.int o
*         On entry: the number of control points in u direction
*       4.int *P
*         On entry: an array containing the indices of the control points whose
* 	          associated periodic B-spline basis functions in u direction
* 		  are nonzero.
*       5.double u[]
*         On entry: the knot vector in u direction of the surface to be
*                 evaluated.
*       6.int r
*         On entry: the index of the knot such that vknots[r] <= v <=
*                   vknots[r+1], where vknots[] is the knot vector in v
*                   direction, v is the  parametric value to be evaluated.
*       7.int s
*         On entry: the order of the derivative w.r.t. v.
*       8.int l
*         On entry: the order in v direction of the surface to be evaluated.
*       9.double v[]
*         On entry: the knot vector in v direction of the surface to be
*                   evaluated.
*       10.vector ***contpts
*         On entry: the control points of the surface to be evaluated.
*       11.double ***derivpts2
*         On entry: a 3D array
*         On exit:  the 3D array containing the control points of the
*                   derivative surface.
* 
* 6   Return Values, Error Indicators and Warnings
*     Not applicable
* 
* 7   Accuracy
*     Not applicable.
* 
* 8   Further comments
*     Not appicable
* 
* 9   Functions referenced by Aijqs_per() are:
*     errormsg()
* 
* 10  Functions that reference Aijqs_per() are:
*     evalderivsurf_per()
*     evaldersurf_per()
* 
****************************************************************************/

void Aijqs_per(int q, int k, int o, int *P, double u[], int r, int s, int l,
	       double v[], vector ***contpts, double ***derivpts2)
/************************************************************************ 
 * The Aij are such that A sub(i) sup(j) = A[i-il+(k-1)][j]
 * Although the variable il is not available here, this value is useful
 * in indexing A[][] compactly and to calculate only Aij required for a 
 * a given il, see program evalderivbsp() where Aij is used as such. 
 **********************************************************************/
{
  vector *v1, *v2, *v3;  /* v2 -- to be used to contain a control point
			  * v1,v3 -- never used (glshen) */
  double du,dv,temp;     /* du,dv -- parametric differences in u,v, */
  int i,j,ireal,jreal,m,n,iv;  /* indices */

  /* for 0-th derivative, the control points are just the given ones. */
  for (i=0; i<k; i++) 
      for (j=0; j<l; j++) {
          v2 = contpts[P[i]][j+r-(l-1)];
          derivpts2[i+k*j][0][0] = v2->x; 
	  derivpts2[i+k*j][0][1] = v2->y;
	  derivpts2[i+k*j][0][2] = v2->z; 
	  derivpts2[i+k*j][0][3] = v2->w;
         }

  /* calculate the control points of the n-th derivative w.r.t. v, where n
   * varying from 1 to s. */
  for (n=1; n<=s; n++) 
      for (i=0; i<k; i++) 
	  for (j=n; j<l; j++) {
              jreal = j+r-(l-1); 
	      dv = v[jreal+l-n] - v[jreal]; 
	      if (dv < 0.0) 
		 errormsg(1,"dv is zero in routine Aij()\n");	
	      for (iv = 0; iv <=3; iv++){
		  temp = derivpts2[i+k*j][k*(n-1)][iv] -
		         derivpts2[i+k*(j-1)][k*(n-1)][iv];
		  derivpts2[i+k*j][k*n][iv] = temp/dv;
		 }
	    }

  /* calculate the control points of the m-th derivative w.r.t. u of the n-th 
   * derivative w.r.t. v, i.e., the result is the mixed (m+n)-th partial 
   * derivative with order m in u and n in v. */
  for (m=1; m<=q; m++) { 
      for (i=m; i<k; i++) { 
          for (j=0; j<l; j++) {
	      ireal = P[i]+k-m ;
              /* if ireal > o, du<0 and du = 1.0-du. */
	      if (ireal>o) { 
	         ireal = ireal%o ;
	         du = 1.0 ;
	        } 
	      else du = 0.0 ;
	      du += u[ireal] - u[P[i]];
	      if (du < 0.0) 
		 errormsg(1,"dv is zero in routine Aij()\n");	
	      for (iv = 0; iv <=3; iv++){
	          temp = derivpts2[i+k*j][m-1+k*s][iv] -
	                 derivpts2[i-1+k*j][m-1+k*s][iv];
	          derivpts2[i+k*j][m+k*s][iv] = temp/du;
	         }
             }
         }
     }
}

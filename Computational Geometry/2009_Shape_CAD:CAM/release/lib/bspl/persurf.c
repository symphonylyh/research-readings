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

/******************************************************************************
*                              revalderivsurf_per()
*******************************************************************************
* 
* 1   Purpose
*     This function returns the vector valued periodic NURBS surface point or 
*     derivative.
* 
* 2   Specification
*     #include "bspl.h"
*     vector *revalderivsurf_per(ParSurf *fgeom, double u, double v,
*                                int uderiv, int vderiv)
* 
* 3   Description
*     This routine returns a vector data structure which contains the point or 
*     the derivatives of the point on a periodic NURBS surface corresponding 
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
*         On entry: the parameter value at which the periodic non-uniform 
* 	          integral B-spline surface is to be evaluated.
*       3.double v
*         On entry: the parameter value at which the periodic non-uniform 
* 	          integral B-spline surface is to be evaluated.
*       4.int uderiv
*         On entry: an index referring to the required derivative with respect 
* 	          to u to be calculated such that if the surface is defined as 
* 		  R(u,v) then the routine will evaluate
*                                    uderiv+vderiv
* 				  @             R(u,v)
* 				-----------------------
*                                      uderiv  vderiv
* 				   @u      @v
* 		  ( @ represents partial derivative. )
* 		  where vderiv is defined below.
*       5.int vderiv
*         On entry: an index referring to the required derivative with respect 
* 	          to v to be calculated such that if the surface is defined as 
* 		  R(u,v) then the routine will evaluate
*                                   uderiv+vderiv
* 				  @             R(u,v)
* 				-----------------------
*                                      uderiv  vderiv
* 				   @u      @v
* 		  ( @ represents partial derivative. )
* 		  where uderiv is defined above.
* 
* 6   Return Values, Error Indicators and Warnings
*     The index uderiv <= the order with respect to u of the periodic non- 
*     uniform integral B-spline surface.
*     The index vderiv <= the order with respect to u of the periodic non- 
*     uniform integral B-spline surface.
*     The value of u should lie between egeom -> uknots[0] and fgeom ->
*     uknots[k]  where k = (fgeom -> ucontpts) + (fgeom -> uorder) - 1.
*     The value of v should lie between egeom -> vknots[0] and fgeom ->
*     vknots[j]  where j = (fgeom -> vcontpts) + (fgeom -> vorder) - 1.
* 
* 7   Accuracy
*     Not applicable.
* 
* 8   Further Comments
*     This basis function evaluation is performed using routine nbasisd_per
*     (bspl).
* 
* 9   Functions referenced by revalderivsurf_per() are:
*     evalderivsurf_per()
*     evalsurf_per()
*     vectfree()
* 
* 10  Functions that reference revalderivsurf_per() are:
*     cos_conv_nurbs()
*     find_points()
*     funcsurfp()
*     guncsurfp()
* 
***************************************************************************/

vector *revalderivsurf_per(ParSurf *fgeom, double u, double v, int uderiv, 
			   int vderiv)

/* evaluate a rational B-spline surface which is periodic in u and open
 * in v */
{
  vector *r0,*r1u,*r1v,*r2u,*r2v,*r2uv;   /* value and deirvatives at u,v */
  double w,w1u,w1v,w2u,w2v,w2uv;          /* scaling factors */

  /* if the orders of partial derivatives are no less than 0 */
  if(uderiv >= 0 && vderiv >= 0) {
    /* evaluate the point (u,v) */
    r0 = evalsurf_per(fgeom,u,v);
    w = r0->w;
    r0->x = r0->x/w;
    r0->y = r0->y/w;
    r0->z = r0->z/w;
    r0->w = 1.0;
    /* if only the value at u,v is wanted, return r0 */
    if(uderiv == 0 && vderiv == 0)
      return(r0);
  }

  /* if the order of partial derivative w.r.t. u is no less than 1 */
  if(uderiv >= 1) {
    /* evaluate the 1st partial derivative w.r.t. u by calling
     * evalderivsurf_per() */
    r1u = evalderivsurf_per(fgeom,u,v,1,0);
    w1u = r1u->w;
    r1u->x = (r1u->x - r0->x*w1u)/w;
    r1u->y = (r1u->y - r0->y*w1u)/w;
    r1u->z = (r1u->z - r0->z*w1u)/w;
    r1u->w = 1.0;
    /* if only the 1st partial derivative w.r.t. u is wanted, return r1u */
    if(uderiv == 1 && vderiv == 0) {
      vectfree(r0);
      return(r1u);
    }
  }

  /* if the order of partial derivative w.r.t. v is no less than 1 */
  if(vderiv >= 1) {	
    r1v = evalderivsurf_per(fgeom,u,v,0,1);
    w1v = r1v->w;
    r1v->x = (r1v->x - r0->x*w1v)/w;
    r1v->y = (r1v->y - r0->y*w1v)/w;
    r1v->z = (r1v->z - r0->z*w1v)/w;
    r1v->w = 1.0;
    /* if only the 1st partial derivative w.r.t. v is wanted, return r1v */
    if(uderiv == 0 && vderiv == 1) {
      vectfree(r0);
      return(r1v);
    }
  }

  /* if the orders of partial derivatives w.r.t. both u & v are no less than
   * 1 */
  if(uderiv >= 1 && vderiv >= 1) {
    r2uv = evalderivsurf_per(fgeom,u,v,1,1);
    w2uv = r2uv->w;
    r2uv->x = (r2uv->x - w1u*r1v->x - w1v*r1u->x - w2uv*r0->x)/w;
    r2uv->y = (r2uv->y - w1u*r1v->y - w1v*r1u->y - w2uv*r0->y)/w;
    r2uv->z = (r2uv->z - w1u*r1v->z - w1v*r1u->z - w2uv*r0->z)/w;
    r2uv->w = 1.0;
    /* if only the 2nd mixed partial derivative w.r.t. u & v is wanted,
     * return r2uv */
    if(uderiv == 1 && vderiv == 1) {
      vectfree(r0);
      vectfree(r1u);
      vectfree(r1v);
      return(r2uv);
    }
  }

  /* if the order of partial derivative w.r.t. u is no less than 2 */
  if(uderiv >= 2)	{
    r2u = evalderivsurf_per(fgeom,u,v,2,0);
    w2u = r2u->w;
    r2u->x = (r2u->x - 2.0*w1u*r1u->x - w2u*r0->x)/w;
    r2u->y = (r2u->y - 2.0*w1u*r1u->y - w2u*r0->y)/w;
    r2u->z = (r2u->z - 2.0*w1u*r1u->z - w2u*r0->z)/w;
    r2u->w = 1.0;
    /* if only the 2nd partial derivative w.r.t. u is wanted, return r2u */
    if(uderiv == 2 && vderiv == 0) {
      vectfree(r0);
      vectfree(r1u);
      return(r2u);
    }
  }

  /* if the order of partial derivative w.r.t. v is no less than 2 */
  if(vderiv >= 2) {
    r2v = evalderivsurf_per(fgeom,u,v,0,2);
    w2v = r2v->w;
    r2v->x = (r2v->x - 2.0*w1v*r1v->x - w2v*r0->x)/w;
    r2v->y = (r2v->y - 2.0*w1v*r1v->y - w2v*r0->y)/w;
    r2v->z = (r2v->z - 2.0*w1v*r1v->z - w2v*r0->z)/w;
    r2v->w = 1.0;
    /* if only the 2nd partial derivative w.r.t. v is wanted, return r2v */
    if(uderiv == 0 && vderiv == 2) {
      vectfree(r0);
      vectfree(r1v);
      return(r2v);
    }
  }
}

/****************************************************************************
*                                evalsurf_per()
*****************************************************************************
*
* 1   Purpose
*     This funtion returns the coordinates of a periodic non-uniform integral 
*     B-spline surface at a given set of parameter values.
* 
* 2   Specification
*     #include "bspl.h"
*     vector *evalsurf_per(ParSurf *fgeom, double u, double v)
* 
* 3   Description
*     This routine returns a vector data structure which contains the 
*     coordinates of the point on a periodic non-uniform integral B-spline 
*     surface corresponding to the supplied parameters.
* 
* 4   References
*     [1]C. De Boor. A Practical Guide to Splines, Springer, New York, 1978.
* 
* 5   Parameters
*       1.ParSurf * fgeom
*         On entry: NURBS surface data structure containing geometry of 
* 	          surface that is to be evaluated.
*       2.double u
*         On entry: parameter value at which the surface is evaluated.
*       3.double v
*         On entry: parameter value at which the surface is evaluated.
* 
* 6   Return Values, Error Indicators and Warnings
*     The value of u should lie between fgeom -> uknots[0] and fgeom ->
*     uknots[k]  where k = (fgeom -> ucontpts) + (fgeom -> uorder) - 1.
*     The value of v should lie between fgeom -> vknots[0] and fgeom ->
*     vknots[j]  where j = (fgeom -> vcontpts) + (fgeom -> vorder) - 1.
* 
* 7   Accuracy
*     Not applicable.
* 
* 8   Further Comments
*     This basis function evaluation is performed using routine nbasisd_per() 
*     This routine should not be used to evaluate derivatives of periodic
*     NURBS curves. 
*     For the evaluation of periodic NURBS curves, see evalrsurf_per().
* 
* 9   Functions referenced by evalsurf_per() are:
*     dbl_array2()
*     find()
*     free_darray2()
*     free_iarray1()
*     int_array1()
*     nbasisd()
*     nbasisd_per()
*     vectalloc()
* 
* 10  Functions that reference evalsurf_per() are:
*     revalderivsurf_per()
* 
***************************************************************************/

vector *evalsurf_per(ParSurf *fgeom, double u, double v)

/* evaluate an integral B-spline surface patch, which is periodic in u
 * and open in v */
{
  int ucontpts,uorder,vcontpts,vorder;  /* order and number of control
					 * points */
  int i,j,ilv,*P;    /* i,j -- loop indices 
		      *	ilv -- the index such that
		      *        vknot[ilv] <= v <= vknot[ilv+1]
		      * P -- the indices of the control points whose associated
		      *      basis functions are nonzero   */
  vector *p,*eval;   /* p -- to be used to contain a control point 
		      * eval -- the evaluated value */
  double q,r,**N,**M;  /* q,r -- values of B-spline basis
			* N -- the basis matrix in u direction
			* M -- the basis matrix in v direction  */

  /* assignment and memory allocation */
  ucontpts = fgeom->ucontpts; vcontpts = fgeom->vcontpts;
  uorder = fgeom->uorder; vorder = fgeom->vorder;
  N = dbl_array2((unsigned) uorder+1,(unsigned) uorder+1);
  M = dbl_array2((unsigned) vorder+1,(unsigned) vorder+1);
  P = int_array1((unsigned) uorder+1);

  /* evaluate the periodic B-spline basis in u */
  nbasisd_per(uorder,ucontpts,u,fgeom->uknots,N,P) ;

  /* evaluate the B-spline basis in v */	
  ilv = find(vcontpts,fgeom->vknots,v);
  nbasisd(ilv, vorder, v, fgeom->vknots, M); 

  /* allocate memory for eval and initialize it */
  eval = vectalloc();
  eval->x = eval->y = eval->z = eval->w = 0.0;

  /* evaluate the point */
  for (i=0;i<uorder;i++) 
      for (j=ilv; j>ilv-vorder; j--) {
          p=fgeom->contpts[P[i]][j] ; 
	  q=N[i+1][uorder] ; 
	  r=M[j+vorder-ilv][vorder] ;
	  eval->x += q*r*p->x ; eval->y += q*r*p->y;
	  eval->z += q*r*p->z ; eval->w += q*r*p->w;	
	  }
  
  /* free memory */
  free_darray2(N);
  free_darray2(M);
  free_iarray1(P);

  return(eval);
}

double ***derivpts;    /* global variable for the derivative of periodic
			* B-spline surface */

/*****************************************************************************
*                               evalderivsurf_per()
******************************************************************************
* 
* 1   Purpose
*     This function returns the vector valued periodic non-uniform integral 
*     B-spline surface point or derivative.
* 
* 2   Specification
*     #include "bspl.h"
*     vector *evalderivsurf_per(ParSurf *fgeom, double u, double v,
*                               int uderiv, int vderiv)
* 
* 3   Description
*     This routine returns a vector data structure which contains the point or 
*     the derivatives of the point on a periodic B-spline surface 
*     corresponding to the parameters (u,v).
* 
* 4   References
*     [1]C. De Boor. A Practical Guide to Splines, Springer, New York, 1978.
* 
* 5   Parameters
*       1.ParSurf * fgeom
*         On entry: NURBS surface data structure containing geometry of 
*                   surface to be evaluated.
*       2.double u
*         On entry: the parameter value at which the non-uniform integral 
*                   B-spline surface is to be evaluated.
*       3.double v
*         On entry: the parameter value at which the non-uniform integral 
*                   B-spline surface is to be evaluated.
*       4.int uderiv
*         On entry: an index referring to the required derivative with respect 
*                   to u to be calculated such that if the surface is defined 
*                   as R(u,v) then the routine will evaluate
*                                  uderiv+vderiv
*                                 @             R(u,v)
*                                ----------------------
*                                    uderiv  vderiv
*                                  @u      @v
*                   ( @ represents partial derivative. )
*                   where vderiv is defined below.
*       5.int vderiv
*         On entry: an index referring to the required derivative with respect 
*                   to v to be calculated such that if the surface is defined 
*                   as R(u,v) then the routine will evaluate
*                                  uderiv+vderiv
*                                 @             R(u,v)
*                                ----------------------
*                                    uderiv  vderiv
*                                  @u      @v
*                   ( @ represents partial derivative. )
*                   where uderiv is defined above.
* 
* 6   Return Values, Error Indicators and Warnings
*     The index <= uderiv  the order with respect to u of the periodic non- 
*     uniform integral B-spline surface.
*     The index <= vderiv  the order with respect to u of the periodic non- 
*     uniformintegral  B-spline surface.
*     The value of u should lie between egeom -> uknots[0] and fgeom ->
*     uknots[k], where k = (fgeom -> ucontpts) + (fgeom -> uorder) - 1.
*     The value of v should lie between egeom -> vknots[0] and fgeom ->
*     vknots[j], where j = (fgeom -> vcontpts) + (fgeom -> vorder) - 1.
* 
* 7   Accuracy
*     Not applicable.
* 
* 8   Further Comments
*     This basis function evaluation is performed using routine nbasisd_per
*     (bspl).
*     This routine should not be used to evaluate derivatives of periodic
*     NURBS curves. 
*     For the evaluation of NURBS curves, see evalrsurf_per (bspl) for details.
* 
* 9   Functions referenced by evalderivsurf_per() are:
*     Aijqs_per()
*     dbl_array2()
*     dbl_array3()
*     find()
*     free_darray2()
*     free_iarray1()
*     int_array1()
*     nbasisd()
*     nbasisd_per()
*     vectalloc()
* 
* 10  Functions that reference evalderivsurf_per() are:
*     revalderivsurf_per()
* 
***************************************************************************/

vector *evalderivsurf_per(ParSurf *fgeom, double u, double v, int uderiv,
			  int vderiv)
{
  int i,j;        /* loop indices */
  int ilv,ii,jj;  /* indices */
  int prod;       /* a scalor factor */
  int ucontpts,uorder,vcontpts,vorder;  /* order and number of control
					 * points */
  int *P;   /* the array containing the indices of the control points whose
	     * associated basis functions are nonzero */
  double q,r;       /* values of B-spline basis */
  double **N,**M;   /* the matrices of periodic B-spline basis */ 
  vector *p,*eval;  /* vectors containing a control point and the evaluated
		     * value */
  static int flag = FALSE; /* flag indicating whether to allocate memory for 
			    * the derivatives */

  /* assignment and memory allocation */
  ucontpts = fgeom->ucontpts; vcontpts = fgeom->vcontpts;
  uorder = fgeom->uorder; vorder = fgeom->vorder;
  if (flag==FALSE){
     derivpts = dbl_array3((unsigned) uorder*vorder,(unsigned) uorder*vorder,
			   4);
     flag = TRUE;
     }
  N = dbl_array2((unsigned) uorder+1,(unsigned) uorder+1);
  M = dbl_array2((unsigned) vorder+1,(unsigned) vorder+1);
  P = int_array1((unsigned) uorder+1);

  /* evaluate the periodic B-spline basis at u */
  nbasisd_per(uorder,ucontpts,u,fgeom->uknots,N,P) ;
	
  /* evaluate the B-spline basis at v */
  ilv = find(vcontpts,fgeom->vknots,v);
  nbasisd(ilv, vorder, v, fgeom->vknots, M); 

  /* calculate the control points of the derivative without the scalor
   * factor */
  Aijqs_per(uderiv,uorder,ucontpts,P,fgeom->uknots,ilv,vderiv,vorder,
            fgeom->vknots,fgeom->contpts, derivpts);

  /* initialize the evaluate value */
  eval = vectalloc(); eval->x = eval->y = eval->z = eval->w = 0.0;

  /* calculate the scalor factor (uorder-1)...(uorder-uderiv)
     (vorder-1)...(vorder-vderiv) */
  prod = 1; 
  for(i=1; i<=uderiv; i++) prod *= uorder-i;
  for(i=1; i<=vderiv; i++) prod *= vorder-i;

  /* evaluate the value at u,v */
  for (i=0;i<uorder-uderiv;i++) 
      for (j=ilv; j>ilv-vorder+vderiv; j--) {
          q = N[i+1][uorder-uderiv]; 
	  r = M[j+vorder-vderiv-ilv][vorder-vderiv];
	  ii = i+uderiv; 
	  jj = j-ilv+(vorder-1);
	  ii = ii + uorder*jj;
	  jj = uderiv+uorder*vderiv;
	  eval->x  += q*r*derivpts[ii][jj][0];
	  eval->y  += q*r*derivpts[ii][jj][1];
	  eval->z  += q*r*derivpts[ii][jj][2];
	  eval->w  += q*r*derivpts[ii][jj][3];
         }

  /* multiply the scalor factor */
  eval->x = prod*eval->x; eval->y = prod*eval->y;
  eval->z = prod*eval->z; eval->w = prod*eval->w;

  /* free memory */
  free_darray2(N);
  free_darray2(M);
  free_iarray1(P);

  return(eval);
}

/************************************************************************/
/***** >>>>>> Aijqs_per() already exists in evalsurf_per.c <<<<<<< ******/
/************************************************************************/
/*
void Aijqs_per(int q, int k, int o, int *P, double u[],
		 int r, int s, int l, double v[], vector ***contpts,
		 double ***derivpts2)
*/
 /************************************************************************ 
 * The Aij are such that A sub(i) sup(j) = A[i-il+(k-1)][j]
 * Although the variable il is not available here, this value is useful
 * in indexing A[][] compactly and to calculate only Aij required for a 
 * a given il, see program evalderivbsp() where Aij is used as such. 
 **********************************************************************/
/*
 * {
 *  vector *v1, *v2, *v3; 
 *  double du,dv,temp; 
 *  int i,j,ireal,jreal,m,n,iv;
 *
 *  for(i=0; i<k; i++) for(j=0; j<l; j++) {
 *     v2 = contpts[P[i]][j+r-(l-1)];
 *     derivpts2[i+k*j][0][0] = v2->x; derivpts2[i+k*j][0][1] = v2->y;
 *     derivpts2[i+k*j][0][2] = v2->z; derivpts2[i+k*j][0][3] = v2->w;
 *   }
 *   for(n=1; n<=s; n++) for(i=0; i<k; i++) for(j=n; j<l; j++) {
 *     jreal = j+r-(l-1); dv = v[jreal+l-n] - v[jreal]; 
 *     if(dv < 0.0) errormsg(1,"dv is zero in routine Aij()\n");	
 *     for(iv = 0; iv <=3; iv++){
 *       temp = derivpts2[i+k*j][k*(n-1)][iv] -
 * 	derivpts2[i+k*(j-1)][k*(n-1)][iv];
 *       derivpts2[i+k*j][k*n][iv] = temp/dv;
 *     }
 *   }
 *   for(m=1; m<=q; m++) { 
 *     for(i=m; i<k; i++) { 
 *       for(j=0; j<l; j++) {
 * 	ireal = P[i]+k-m ;
 * 	if (ireal>o) { 
 * 	  ireal = ireal%o ;
 * 	  du = 1.0 ;
 * 	} 
 * 	else du = 0.0 ;
 * 	du += u[ireal] - u[P[i]];
 * 	if(du < 0.0) errormsg(1,"dv is zero in routine Aij()\n");	
 * 	for(iv = 0; iv <=3; iv++){
 * 	  temp = derivpts2[i+k*j][m-1+k*s][iv] -
 * 	    derivpts2[i-1+k*j][m-1+k*s][iv];
 * 	  derivpts2[i+k*j][m+k*s][iv] = temp/du;
 * 	}
 *       }
 *     }
 *   }
 * }
*/

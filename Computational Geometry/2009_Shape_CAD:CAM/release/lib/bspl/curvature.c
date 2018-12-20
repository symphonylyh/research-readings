/* ***************************************************************************
 Copyright (C) 1996 Massachusetts Institute of Technology all rights reserved 
	Programmer: George A. Kriezis
**************************************************************************** */
# include <math.h>
# include "gen.h"
# include "bspl.h"

/*****************************************************************************
*                               curvature1()
******************************************************************************
* 
* 1   Purpose
*     This function calculates the maximum and minimum principal curvatures 
*     for a NURBS surface at a specified point on the surface.
* 
* 2   Specification
*     #include "bspl.h"
*     int curvature1(ParSurf *fgeom, double u, double v, double *kmax, 
*                    double *kmin)
* 
* 3   Description
*     This routine computes the maximum and minimum principal curvatures of a 
*     NURBS surface point by solving for the roots of a quadratic equation. 
*     The coefficients of the quadratic equation are obtained from a 
*     minimization of all the normal curvatures of a surface at a given point. 
*     These coefficients are functions of the elements of the first and
*     second fundamental tensors of the surface at the given point.
* 
* 4   References
*     [1] P. G. Alourdas. Shape Creation, Interrogation and Fairing Using 
*         B-Splines, Engineer's Thesis, Massachusetts Institute of Technology, 
* 	Department of Ocean Engineering, Cambridge, Massachusetts, 1989.
* 
* 5   Parameters
*        1.ParSurf * fgeom
*          On entry: NURBS surface data structure containing the geometry of 
* 	           the surface to be evaluated.
*        2.double u
*          On entry: u parametric value of the surface point at which the 
* 	           curvatures are to be computed.
*        3.double v
*          On entry: v parametric value of the surface point at which the 
* 	           curvatures are to be computed.
*        4.double * kmax
*          On exit: the evaluated maximum principal curvature of the NURBS 
* 	          surface.
*        5.double * kmin
*          On exit: the evaluated minimum principal curvature of the NURBS 
* 	          surface.
* 
* 6   Return Values, Error Indicators and Warnings
*     Not applicable.
* 
* 7   Accuracy
*     Not applicable.
* 
* 8   Further Comments
*     This is a modified and more efficient version of curvature 
*     (bspl _ curvature.c). This routine does not return a value.
* 
* 9   Functions referenced by curvature1() are:
*     dot()
*     errormsg()
*     evalrsurf()
*     mag()
*     normalsurf1()
* 
* 10  Functions that reference curvature1() are:
*     lefunct1()
* 
******************************************************************************/

int curvature1(ParSurf *fgeom, double u, double v, double *kmax, double *kmin)
/* more efficient subroutine for calculating curvatures of parametric
   patches */
{
  vector *temp,*r[6],*ru,*rv,*ruu,*rvv,*ruv,*n;
         /* temp -- never used (glshen) 
	  * r[6] -- array containing derivatives
	  * ru, rv, ruu, rvv, ruv -- derivatives
	  * n -- normal    */
  double H,K,Gdet,Ddet,d11,d12,d21,d22,g11,g12,g21,g22,tempmod;
         /* H -- gaussian curvature
	  * K -- mean curvature
	  * others -- intermediate results */

  /* calculate the derivatives */
  evalrsurf(fgeom,u,v,5,r);
  ru = r[1];                /* 1st derivative with respect to u */
  rv = r[2];                /* 1st derivative with respect to v */
  /* if both the derivatives are 0 */
  if ( mag(ru)*mag(rv) < ZERO)
     {
      errormsg(1,"degenerate point, curvature not being calculated");
      return(0);
     }
  ruu = r[4];               /* 2nd derivative with respect to u */
  rvv = r[5];               /* 2nd derivative with respect to v */
  ruv = r[3];               /* 2nd comb. derivative with respect to u & v */
  
  /* normal at the point */
  n = normalsurf1(fgeom,u,v);

  /* intermediate results */
  g11 = dot(ru,ru);
  g12 = g21 = dot(ru,rv);
  g22 = dot(rv,rv);

  d11 = dot(n,ruu);
  d12 = d21 = dot(n,ruv);
  d22 = dot(n,rvv);

  Gdet = g11*g22 - g12*g21;
  Ddet = d11*d22 - d12*d21;

  /* gaussian curvature */
  H = ((g12*d21 + g21*d12) - (g11*d22 + g22*d11))/(2.0*Gdet);

  /* mean curvature */
  K = Ddet/Gdet;

  tempmod = H*H - K;
  tempmod = FABS(tempmod);

  /* maximum and minimum principal curvatures */
  *kmax = H + SQRT(tempmod);
  *kmin = H - SQRT(tempmod);
}

/****************************************************************************
*                              curvature()
*****************************************************************************
*
* 1   Purpose
*     This function calculates the maximum and minimum principal curvatures 
*     for a NURBS surface at a specified point on the surface.
* 
* 2   Specification
*     #include "bspl.h"
*     int curvature(ParSurf *fgeom, double u, double v, double *kmax, 
*                   double *kmin)
* 
* 3   Description
*     This routine computes the maximum and minimum principal curvatures of a 
*     NURBS surface point by solving for the roots of a quadratic equation. 
*     The coefficients of the quadratic equation are obtained from a 
*     minimization of all the normal curvatures of a surface at a given point. 
*     These coefficients are functions of the elements of the first and 
*     second fundamental tensors of the surface at the given point.
* 
* 4   References
*     [1] P. G. Alourdas. Shape Creation, Interrogation and Fairing Using 
*         B-Splines, Engineer's Thesis, Massachusetts Institute of Technology, 
* 	Department of Ocean Engineering, Cambridge, Massachusetts, 1989.
* 
* 5   Parameters
*        1.ParSurf * fgeom
*          On entry: NURBS surface data structure containing the geometry of
* 	           the surface to be evaluated.
*        2.double u
*          On entry: u parametric value of the surface point at which the 
* 	           curvatures are to be computed.
*        3.double v
*          On entry: v parametric value of the surface point at which the 
* 	           curvatures are to be computed.
*        4.double * kmax
*          On exit: the evaluated maximum principal curvature of the NURBS 
* 	          surface.
*        5.double * kmin
*          On exit: the evaluated minimum principal curvature of the NURBS 
* 	          surface.
* 
* 6   Return Values, Error Indicators and Warnings
*     Not applicable.
* 
* 7   Accuracy
*     Not applicable.
* 
* 8   Further Comments
*     This is not a very efficient routine. For a more efficient routine for 
*     computing the surface principal curvatures, see curvature1 
*     (bspl _ curvature.c). This routine does not return a value.
* 
* 9   Functions referenced by curvature() are:
*     dot()
*     errormsg()
*     mag()
*     normalsurf()
*     revalderivsurf()
* 
* 10  Functions that reference curvature() are: None
* 
*****************************************************************************/

int curvature(ParSurf *fgeom, double u, double v, double *kmax, double *kmin)
{
  /* see curvature1() about these variables */
  vector *temp, *ru,*rv,*ruu,*rvv,*ruv,*n;
  double H,K,Gdet,Ddet,d11,d12,d21,d22,g11,g12,g21,g22,tempmod;

  /* calculate the derivatives */
  ru = revalderivsurf(fgeom,u,v,1,0);   /* 1st derivative with respect to u */
  rv = revalderivsurf(fgeom,u,v,0,1);   /* 1st derivative with respect to v */

  /* if both the derivatives are 0 */
  if ( mag(ru)*mag(rv) < ZERO)
     {
      errormsg(1,"degenerate point, curvature not being calculated");
      return(0);
     }

  ruu = revalderivsurf(fgeom,u,v,2,0);   /* 2nd derivative with respect to u */
  rvv = revalderivsurf(fgeom,u,v,0,2);   /* 2nd derivative with respect to v */
  ruv = revalderivsurf(fgeom,u,v,1,1);   /* 2nd comb. derivative with respect
					  * to u & v */

  n = normalsurf(fgeom,u,v);       /* normal at the point */

  /* intermediate results */
  g11 = dot(ru,ru);
  g12 = g21 = dot(ru,rv);
  g22 = dot(rv,rv);

  d11 = dot(n,ruu);
  d12 = d21 = dot(n,ruv);
  d22 = dot(n,rvv);

  Gdet = g11*g22 - g12*g21;
  Ddet = d11*d22 - d12*d21;

  /* gaussian curvature */
  H = ((g12*d21 + g21*d12) - (g11*d22 + g22*d11))/(2.0*Gdet);

  /* mean curvature */
  K = Ddet/Gdet;

  tempmod = H*H - K;
  tempmod = FABS(tempmod);

  /* maximum and minimum principal curvatures */
  *kmax = H + SQRT(tempmod);
  *kmin = H - SQRT(tempmod);
}






















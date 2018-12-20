/* ***************************************************************************
 Copyright (C) 1996 Massachusetts Institute of Technology all rights reserved 
	Programmer: George A. Kriezis
**************************************************************************** */
#include <malloc.h>
#include "gen.h"
#include "bspl.h"

static double ***derivpts;

static int NUMBERPOINTS = 0; /* this is so that memory is only allocated a
			      * single time. That is, when 0, memory is
			      * allocted and the variable set to the length 
			      * of the array derivpts. Subsequent calls will
                              * check this against the current required
                              * length, if greater, will reallocate derivpts
                              * if less, will not. */
#define ZERO1 1.e-10

/*****************************************************************************
*                                 normalsurf()
******************************************************************************
* 
* 1   Purpose
*     This function returns the computed unit vector normal of a NURBS surface.
* 
* 2   Specification
*     #include "bspl.h"
*     vector *normalsurf(ParSurf *fgeom, double u, double v)
* 
* 3   Description
*     This routine computes the normal to a surface at a given parametric 
*     point by finding the unit vector in the direction of the cross product 
*     of the first order derivatives in the u and v directions.
*     If one of the derivatives is zero at a point due to the parametrization, 
*     then the above calculation becomes undefined.  In such cases, the second 
*     order mixed derivative in the corresponding direction is used in the 
*     cross product.  Such cases normally occur when quadrilateral patches 
*     degenerate to triangular patches (see section on Triangular Patches in
*     Faux and Pratt (1979)). The computed normal is returned as a vector
*     whose w = 1.
* 
* 4   References
*     [1] P. G. Alourdas. Shape Creation, Interrogation and Fairing Using 
*         B-Splines, Engineer's Thesis, Massachusetts Institute of Technology, 
* 	Department of Ocean Engineering, Cambridge, Massachusetts, 1989.
*     [2] I. D. Faux and M. J. Pratt.  Computational Geometry for Design and 
*         Manufacture, Ellis Horwood: Chichester, England, 1979.
* 
* 5   Parameters
*        1.ParSurf * fgeom
*          On entry: NURBS surface data structure containing geometry of
* 	           surface that is to be evaluated.
*        2.double u
*          On entry: parameter value at which to compute the surface normal.
*        3.double v
*          On entry: parameter value at which to compute the surface normal.
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
* 9   Functions referenced by normalsurf() are:
*     cross()
*     mag()
*     revalderivsurf()
*     unitvector1()
* 
* 10  Functions that reference normalsurf() are:
*     build_offsets_int()
*     check_patch_u()
*     compute_diffeq()
*     compute_ortho_proj()
*     con_confun()
*     curvature()
*     DrawFacetEvaluateMap()
*     DrawSurfEvaluate()
*     DrawSurfPivot()
*     eval_cylsect_param()
*     eval_cylsect_param3()
*     frenet_tr()
*     funcsurf()
*     FuncSurfAndCyl()
*     FuncSurfAndCyl3d()
*     guncsurf()
*     init_ortho_proj()
*     init_proj()
*     integral_build_offset()
*     integral_sample_offset()
*     localize_diagnostic()
*     localize_sumsq()
*     localize_sumsq_opt()
*     main()
*     MinDistance()
*     PostScriptFacetEvaluateMap()
*     PostScriptSurfEvaluate()
*     sample_mod_surf()
*     sample_offset()
*     signed_distance()
*     SubDistance()
*     tol_sample_offsets()
* 
******************************************************************************/

vector *normalsurf(ParSurf *fgeom, double u, double v)
{
  vector *ru,*rv,*ruv;   /* derivatives */
  vector *n;             /* the normal */

  ru = revalderivsurf(fgeom,u,v,1,0); /* 1st partial derivative w.r.t. u */
  rv = revalderivsurf(fgeom,u,v,0,1); /* 1st partial derivative w.r.t. v */

  if(mag(ru) < ZERO1)                /* if ru = 0 within tolerence ZERO1 */
	{
	ruv = revalderivsurf(fgeom,u,v,1,1);   /* 2nd mixed derivative is used 
                                                * instead of ru */
	n = cross(ruv,rv);	/* normal is the cross product of ruv and rv */
	free((char *) ruv);
	}
  else if(mag(rv) < ZERO1)      /* if rv = 0 within tolerence ZERO1 */
	{
	ruv = revalderivsurf(fgeom,u,v,1,1);   /* 2nd mixed derivative is used
                                                * instead of rv */
	n = cross(ru,ruv);      /* normal is the cross product of ru and ruv */
	free((char *) ruv);
	}
  else
	n = cross(ru,rv);       /* if both 1st partial derivatives are nonzero,
                                 * the normal is the cross product of ru and
				 * uv */

  unitvector1(n,n);   /* normalize to unit vector */

  free((char *) ru);
  free((char *) rv);

  return(n);
}

/*****************************************************************************
*                                 revalderivsurf()
******************************************************************************
* 
* 1   Purpose
*     This function returns the vector valued NURBS surface point or
*     derivative.
* 
* 2   Specification
*     #include "bspl.h"
*     vector *revalderivsurf(ParSurf *fgeom, double u, double v, int uderiv, 
*                            int vderiv)
* 
* 3   Description
*     This routine returns a vector data structure which contains the point or 
*     the derivatives of the point on a NURBS surface corresponding to the 
*     parameters (u,v).
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
*                                    uderiv+vderiv
* 				  @             R(u,v)
* 				-----------------------
*                                      uderiv  vderiv
* 				   @u      @v
* 		  ( @ represents partial derivative. )
* 		  where uderiv is defined above.
* 
* 6   Return Values, Error Indicators and Warnings
*     The index uderiv <= the order with respect to u of the non-uniform 
*     integral B-spline surface.
*     The index vderiv <= the order with respect to u of the non-uniform 
*     integral B-spline surface.
*     The value of u should lie between fgeom->uknots[fgeom->uorder-1] and
*     fgeom->uknots[ucontpts].
*     The value of v should lie between fgeom->vknots[fgeom->vorder-1] and
*     fgeom->vknots[vcontpts].
* 
* 7   Accuracy
*     Not applicable.
* 
* 8   Further Comments
*     This basis function evaluation is performed using routine nbasisd (bspl).
* 
* 9   Functions referenced by revalderivsurf() are:
*     evalderivsurf()
*     evalsurf()
*     vectfree()
*  
* 10  Functions that reference revalderivsurf() are:
*     build_offsets_int()
*     CalcHubTip()
*     check_patch_u()
*     cos_conv_nurbs()
*     curvature()
*     cyl_sample_gencyl()
*     DrawGeodesics()
*     eval_LE3D()
*     find_points()
*     frenet_tr()
*     funcsurf()
*     generatrix()
*     guncsurf()
*     integral_build_offset()
*     integral_sample_offset()
*     int_sample_gencyl()
*     normalsurf()
*     position_err()
*     PowParSurf_compare()
*     sample_gencyl()
*     sample_mod_surf()
*     sample_offset()
*     SpanwiseU()
*     test_tangent()
*     tol_sample_offsets()
* 
******************************************************************************/

vector *revalderivsurf(ParSurf *fgeom, double u, double v, int uderiv,
		       int vderiv)
{
  vector *r0,*r1u,*r1v,*r2u,*r2v,*r2uv;  /* value and deirvatives at u,v */
  double w,w1u,w1v,w2u,w2v,w2uv;         /* scaling factors */

  /* if the orders of partial derivatives are no less than 0 */
  if(uderiv >= 0 && vderiv >= 0) {
    r0 = evalsurf(fgeom,u,v);            /* evaluate the point u,v */
    w = r0->w;
    /* make w-coordinate  = 1 */
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
    /* evaluate the 1st partial derivative w.r.t. u by calling evalderivsurf */
    r1u = evalderivsurf(fgeom,u,v,1,0);        
    w1u = r1u->w;
    /* make w-coordinate = 1 */
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
    /* evaluate the 1st partial derivative w.r.t. v by calling evalderivsurf */
    r1v = evalderivsurf(fgeom,u,v,0,1);
    w1v = r1v->w;
    /* make w-coordinate = 1 */
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
    /* evaluate the 2nd mixed partial derivative w.r.t. u & v */
    r2uv = evalderivsurf(fgeom,u,v,1,1);
    w2uv = r2uv->w;
    /* make w-coordinate = 1 */
    r2uv->x = (r2uv->x - w1u*r1v->x - w1v*r1u->x - w2uv*r0->x)/w;
    r2uv->y = (r2uv->y - w1u*r1v->y - w1v*r1u->y - w2uv*r0->y)/w;
    r2uv->z = (r2uv->z - w1u*r1v->z - w1v*r1u->z - w2uv*r0->z)/w;
    r2uv->w = 1.0;
    /* if only the 2nd mixed partial derivative w.r.t. u & v is wanted, return
     * r2uv */
    if(uderiv == 1 && vderiv == 1) {
      vectfree(r0);
      vectfree(r1u);
      vectfree(r1v);
      return(r2uv);
    }
  }

  /* if the order of partial derivative w.r.t. u is no less than 2 */
  if(uderiv >= 2) {
    /* evaluate the 2nd partial derivative w.r.t. u by calling evalderivsurf */
    r2u = evalderivsurf(fgeom,u,v,2,0);
    w2u = r2u->w;
    /* make w-coordinate = 1 */
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
    /* evaluate the 2nd partial derivative w.r.t. v by calling evalderivsurf */
    r2v = evalderivsurf(fgeom,u,v,0,2);
    w2v = r2v->w;
    /* make w-coordinate = 1 */
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

/******************************************************************************
*                                evalsurf()
*******************************************************************************
* 
* 1   Purpose
*     This funtion returns the coordinates of a non-uniform integral B-spline 
*     surface at a given set of parameter values.
* 
* 2   Specification
*     #include "bspl.h"
*     vector *evalsurf(ParSurf *fgeom, double u, double v)
* 
* 3   Description
*     This routine returns a vector data structure which contains the
*     coordinates of the point on a non-uniform integral B-spline surface 
*     corresponding to the supplied parameters.
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
*     The value of u should lie between fgeom->uknots[fgeom->uorder-1] and
*     fgeom->uknots[fgeom->ucontpts].
*     The value of v should lie between fgeom->vknots[fgeom->vorder-] and
*     fgeom->vknots[fgeom->vcontpts].
* 
* 7   Accuracy
*     Not applicable.
* 
* 8   Further Comments
*     This basis function evaluation is performed using routine nbasisd 
*     (bspl _ cxdeb.c).
*     This routine should not be used to evaluate derivatives of NURBS curves. 
*     For the evaluation of NURBS curves, see evalrsurf (bspl _ evalsurf1.c) 
*     for details.
* 
* 9   Functions referenced by evalsurf() are:
*     dbl_array2()
*     find()
*     free_darray2()
*     nbasisd()
*     vectalloc()
*  
* 10  Functions that reference evalsurf() are:
*     build_offset()
*     DrawCosKnots()
*     DrawCosPolygon()
*     DrawFacetHeightMap()
*     DrawFacetShadedMap()
*     DrawFacetSlopeMap()
*     DrawFacetWireframeMap()
*     DrawSurfKnots()
*     DrawSurfSplits()
*     DrawSurfWireframe()
*     DrawTrimWireframe()
*     DrawUv()
*     DrawUvPolygon()
*     EvalSurfAndCyl()
*     EvalSurfAndCyl3d()
*     find_error_surf()
*     find_error_surf1()
*     FuncSurfAndCyl()
*     FuncSurfAndCyl3d()
*     init_ortho_proj()
*     init_proj()
*     output_3D()
*     PostScriptCosKnots()
*     PostScriptCosPolygon()
*     PostScriptFacetMap()
*     PostScriptSurfKnots()
*     PostScriptSurfWireframe()
*     PostScriptTrimWireframe()
*     PostScriptUv()
*     PostScriptUvPolygon()
*     ps_ParSurf()
*     ps_ParSurf2()
*     revalderivsurf()
*     RobustDistCB()
*     SampleSurfGridCB()
*     SampleSurfListCB()
*     SaveIgesCos()
*     SaveIgesSurf()
*     scattered_fit()
*     signed_distance()
*     TrimSection()
*     TrimSection3d()
*     trim_section()
*     trim_section3()
*     UvInputCB()
*     UvPtsToListCB()
* 
*******************************************************************************/

vector *evalsurf(ParSurf *fgeom, double u, double v)
{
  int ucontpts,uorder,vcontpts,vorder;  /* orders and numbers of control 
                                         * points in u,v directions */
  int i,j,ilu,ilv;                      /* index */
  vector *p;        /* a control point */
  vector *eval;     /* vector containing the coordinates of evaluated point */
  double q,r;       /* B-spline basis */
  double **N,**M;   /* matrices containing B-spline basis in u,v directions */ 

  /* initialize */
  ucontpts = fgeom->ucontpts;
  vcontpts = fgeom->vcontpts;
  uorder = fgeom->uorder;
  vorder = fgeom->vorder;
  N = dbl_array2((unsigned)(uorder+1),(unsigned)(uorder+1));
  M = dbl_array2((unsigned)(vorder+1),(unsigned)(vorder+1));

  /* find the index ilu, such that uknots[ilu] <= u <= uknots[ilu+1] */
  ilu = find(ucontpts,fgeom->uknots,u);
  /* calculate the matrix of B-spline basis in u direction */
  nbasisd(ilu, uorder, u, fgeom->uknots, N); 
	 
  /* find the index ilv, such that vknots[ilv] <= v <= vknots[ilv+1] */
  ilv = find(vcontpts,fgeom->vknots,v);
  /* calculate the matrix of B-spline basis in v direction */
  nbasisd(ilv, vorder, v, fgeom->vknots, M); 

  /* allocate memory for vector eval, and initialize it */
  eval = vectalloc();
  eval->x = eval->y = eval->z = eval->w = 0.0;

  /* summation over the control points in u,v directions */
  for(i=ilu; i>ilu-uorder; i--)
	for(j=ilv; j>ilv-vorder; j--)
		{
		p = fgeom->contpts[i][j];		
		q = N[i+uorder-ilu][uorder];
		r = M[j+vorder-ilv][vorder];
		eval->x += q*r*p->x;	
		eval->y += q*r*p->y;	
		eval->z += q*r*p->z;	
		eval->w += q*r*p->w;	
		}

  free_darray2(M);
  free_darray2(N);	

  return(eval);
}

/*****************************************************************************
*                                evalderivsurf()
******************************************************************************
* 
* 1   Purpose
*     This function returns the vector valued non-uniform integral B-spline 
*     surface point or derivative.
* 
* 2   Specification
*     #include "bspl.h"
*     vector *evalderivsurf(ParSurf *fgeom, double u, double v, int uderiv, 
*                           int vderiv)
* 
* 3   Description
*     This routine returns a vector data structure which contains the point or 
*     the derivatives of the point on a B-spline surface corresponding to the 
*     parameters (u,v).
* 
* 4   References
*     [1]C. De Boor. A Practical Guide to Splines, Springer, New York, 1978.
* 
* 5   Parameters
*       1.ParSurf * fgeom
*         On entry: NURBS surface data structure containing geometry of 
* 	          surface points to be evaluated.
*       2.double u
*         On entry: the parameter value at which the non-uniform integral 
* 	          B-spline surface is to be evaluated.
*       3.double v
*         On entry: the parameter value at which the non-uniform integral 
* 	          B-spline surface is to be evaluated.
*       4.int uderiv
*         On entry: an index referring to the required derivative with respect 
* 	          to u to be calculated such that if the surface is defined 
* 		  as R(u,v) then the routine will evaluate
*                                  uderiv+vderiv
*                                 @             R(u,v)
*                                ----------------------
*                                    uderiv  vderiv
*                                  @u      @v
* 		  ( @ represents partial derivative. )
* 		  where vderiv is defined below.
*       5.int vderiv
*         On entry: an index referring to the required derivative with respect 
* 	          to v to be calculated such that if the surface is defined as 
* 		  R(u,v) then the routine will evaluate
*                                  uderiv+vderiv
* 			        @             R(u,v)
*                                ----------------------
*                                    uderiv  vderiv
*                                  @u      @V
* 		  ( @ represents partial derivative. )
* 		  where uderiv is defined above.
* 
* 6   Return Values, Error Indicators and Warnings
*     The index uderiv <= the order with respect to u of the non-uniform 
*     integral B-spline surface.
*     The index vderiv <= the order with respect to u of the non-uniform 
*     integral B-spline surface.
*     The value of u should lie between fgeom->uknots[fgeom->uorder-1] and
*     fgeom->uknots[fgeom->ucontpts].
*     The value of v should lie between egeom->vknots[fgeom->vorder-1] and
*     fgeom->vknots[fgeom->vcontpts].
* 
* 7   Accuracy
*     Not applicable.
* 
* 8   Further Comments
*     This basis function evaluation is performed using routine nbasisd (bspl).
*     This routine should not be used to evaluate derivatives of NURBS curves. 
*     For the evaluation of NURBS curves, see evalrsurf (bspl _ ) for details.
* 
* 9   Functions referenced by evalderivsurf() are:
*     Aijqs()
*     dbl_array2()
*     dbl_array3()
*     find()
*     free_darray2()
*     free_darray3()
*     nbasisd()
*     vectalloc()
* 
* 10  Functions that reference evalderivsurf() are:
*     ParSurf_to_PowSurf()
*     revalderivsurf()
* 
*****************************************************************************/

vector *evalderivsurf(ParSurf *fgeom, double u, double v, int uderiv,
		      int vderiv)
{
  int ucontpts,uorder,vcontpts,vorder; /* orders and numbers of control points
                                        * in u,v directions */
  int i,j,ilu,ilv,ii,jj;               /* indices */
  int prod;                            /* factor of the derivative surface */
  vector *p;         /* a control point */
  vector *eval;      /* vector containing the evaluated derivative */
  double q,r;        /* B-spline basis */
  double **N,**M;    /* matrices containing B-spline basis in u,v directions */

  /* initialize */
  ucontpts = fgeom->ucontpts;
  vcontpts = fgeom->vcontpts;
  uorder = fgeom->uorder;
  vorder = fgeom->vorder;
  N = dbl_array2((unsigned)(uorder+1),(unsigned)(uorder+1));
  M = dbl_array2((unsigned)(vorder+1),(unsigned)(vorder+1));

  /* allocate memory for derivatives */
  if (NUMBERPOINTS == 0){
     derivpts = dbl_array3((unsigned)(uorder*vorder),
			   (unsigned)(uorder*vorder),4);
     NUMBERPOINTS = uorder*vorder;
     }
  else if (NUMBERPOINTS < uorder*vorder){
     free_darray3(derivpts);
     derivpts = dbl_array3((unsigned)(uorder*vorder),
			   (unsigned)(uorder*vorder),4);
     NUMBERPOINTS = uorder*vorder;
  }

  /* find the index ilu, such that uknots[ilu] <= u <= uknots[ilu+1] */
  ilu = find(ucontpts,fgeom->uknots,u);
  /* calculate the matrix of B-spline basis in u direction */
  nbasisd(ilu, uorder, u, fgeom->uknots, N); 
	
  /* find the index ilv, such that vknots[ilv] <= v <= vknots[ilv+1] */
  ilv = find(vcontpts,fgeom->vknots,v);
  /* calculate the matrix of B-spline basis in v direction */
  nbasisd(ilv, vorder, v, fgeom->vknots, M); 

  /* calculate the control points of the derivative surface (without factor) */
  Aijqs(ilu,uderiv,uorder,fgeom->uknots,ilv,vderiv,vorder,fgeom->vknots,
	fgeom->contpts, derivpts);

  /* allocate memory for vector eval, and initialize it */
  eval = vectalloc();
  eval->x = eval->y = eval->z = eval->w = 0.0;

  /* calculate the factor (uorder-1)...(uorder-uderiv)(vorder-1)...
   * (vorder-vderiv) of the derivative surface */
  prod = 1.0;
  for(i=1; i<=uderiv; i++)
	prod *= uorder-i;
  for(i=1; i<=vderiv; i++)
	prod *= vorder-i;

  /* evaluate the derivative surface at u,v (without factor) */
  for(i=ilu; i>ilu-uorder+uderiv; i--)
	for(j=ilv; j>ilv-vorder+vderiv; j--)
		{
		q = N[i+uorder-uderiv-ilu][uorder-uderiv];
		r = M[j+vorder-vderiv-ilv][vorder-vderiv];
		ii = i-ilu+(uorder-1);
		jj = j-ilv+(vorder-1);
		ii = ii+uorder*jj;
		jj = uderiv+uorder*vderiv;
		eval->x  += q*r*derivpts[ii][jj][0];
		eval->y  += q*r*derivpts[ii][jj][1];
		eval->z  += q*r*derivpts[ii][jj][2];
		eval->w  += q*r*derivpts[ii][jj][3];
		}
	
  /* multiply by factor */
  eval->x = prod*eval->x;
  eval->y = prod*eval->y;
  eval->z = prod*eval->z;
  eval->w = prod*eval->w;

  free_darray2(M);
  free_darray2(N);

  return(eval);
}

/****************************************************************************
*                                 zevalderivsurf()
*****************************************************************************
* 
* 1   Purpose
*     This function returns z-component of the vector valued non-uniform 
*     integral B-spline surface point or derivative.
* 
* 2   Specification
*     #include "bspl.h"
*     vector *zevalderivsurf(ParSurf *fgeom, double u, double v, int uderiv, 
*                           int vderiv)
* 
* 3   Description
*     This routine returns z-component of the point or the derivatives of the 
*     point on a B-spline surface corresponding to the parameters (u,v).
* 
* 4   References
*     [1]C. De Boor. A Practical Guide to Splines, Springer, New York, 1978.
* 
* 5   Parameters
*       1.ParSurf * fgeom
*         On entry: NURBS surface data structure containing geometry of 
*         surface to be evaluated.
*       2.double u
*         On entry: the parameter value at which the non-uniform integral 
*                   B-spline surface is to be evaluated.
*       3.double v
*         On entry: the parameter value at which the non-uniform integral 
*                   B-spline surface is to be evaluated.
*       4.int uderiv
*         On entry: an index referring to the required derivative with respect 
*                   to u to be calculated such that if the surface is defined 
*                   as R(u,v) then the routine will evaluate z-component of
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
*                   as R(u,v) then the routine will evaluate z-component of
*                                  uderiv+vderiv
*                                 @             R(u,v)
*                                ----------------------
*                                    uderiv  vderiv
*                                  @u      @V
*                   ( @ represents partial derivative. )
*                   where uderiv is defined above.
* 
* 6   Return Values, Error Indicators and Warnings
*     The index uderiv <= the order with respect to u of the non-uniform 
*     integral B-spline surface.
*     The index vderiv <= the order with respect to u of the non-uniform 
*     integral B-spline surface.
*     The value of u should lie between fgeom->uknots[fgeom->uorder-1] and
*     fgeom->uknots[fgeom->ucontpts].
*     The value of v should lie between fgeom->vknots[fgeom->vorder-1] and
*     fgeom->vknots[fgeom->vcontpts].
* 
* 7   Accuracy
*     Not applicable.
* 
* 8   Further Comments
*     This basis function evaluation is performed using routine nbasisd (bspl).
*     For evaluation of the point of the derivative of a vector-valued B-spline
*     surface, see evalderivsurf().
* 
* 9   Functions referenced by zevalderivsurf() are:
*     Aijqs_z()
*     dbl_array2()
*     dbl_array3()
*     find()
*     free_darray2()
*     free_darray3()
*     nbasisd()
* 
* 10  Functions that reference zevalderivsurf() are: None
* 
***************************************************************************/

double zevalderivsurf(ParSurf *fgeom, double u, double v, int uderiv,
		      int vderiv)
{
  int ucontpts,uorder,vcontpts,vorder; /* orders and numbers of control points
                                        * in u,v directions */
  int i,j,ilu,ilv,ii,jj;               /* indices */
  int prod;                            /* factor of the derivative surface */
  vector *p;                   
  double eval,q,r;   /* eval is to contain z-component of the derivative. */
  double **N,**M;    /* matrices containing B-spline basis in u,v directions */

  ucontpts = fgeom->ucontpts;
  vcontpts = fgeom->vcontpts;
  uorder = fgeom->uorder;
  vorder = fgeom->vorder;
  N = dbl_array2((unsigned)(uorder+1),(unsigned)(uorder+1));
  M = dbl_array2((unsigned)(vorder+1),(unsigned)(vorder+1));
  
  /* allocate memory */
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

  /* find the index ilu, such that uknots[ilu] <= u <= uknots[ilu+1] */
  ilu = find(ucontpts,fgeom->uknots,u);
  /* calculate the matrix of B-spline basis in u direction */
  nbasisd(ilu, uorder, u, fgeom->uknots, N); 
	
  /* find the index ilv, such that vknots[ilv] <= v <= vknots[ilv+1] */
  ilv = find(vcontpts,fgeom->vknots,v);
  /* calculate the matrix of B-spline basis in v direction */
  nbasisd(ilv, vorder, v, fgeom->vknots, M); 

  /* calculate z-components of the control points of the derivative surface 
     (without factor) */
  Aijqs_z(ilu,uderiv,uorder,fgeom->uknots,ilv,vderiv,vorder,fgeom->vknots,
	  fgeom->contpts, derivpts);

  /* calculate the factor (uorder-1)...(uorder-uderiv)(vorder-1)...
   *(vorder-vderiv) of the derivative surface */
  prod = 1.0;
  for(i=1; i<=uderiv; i++)
	prod *= uorder-i;
  for(i=1; i<=vderiv; i++)
	prod *= vorder-i;

  eval = 0.0;

  /* evaluate z-component of the derivative surface at u,v (without factor) */
  for(i=ilu; i>ilu-uorder+uderiv; i--)
	for(j=ilv; j>ilv-vorder+vderiv; j--)
		{
		q = N[i+uorder-uderiv-ilu][uorder-uderiv];
		r = M[j+vorder-vderiv-ilv][vorder-vderiv];
		ii = i-ilu+(uorder-1);
		jj = j-ilv+(vorder-1);
		ii = ii +uorder*jj;
		jj = uderiv + uorder*vderiv;
		eval += q*r*derivpts[ii][jj][2];
		}
	
  /* multiply by the factor */
  eval = prod*eval;

  free_darray2(M);
  free_darray2(N);

  return(eval);
}

/****************************************************************************
*                                  Aijqs()
*****************************************************************************
* 
* 1   Purpose
*     This is a subroutine of function evalderivsurf(), which evaluates the 
*     derivatives of a B-spline surface, see evalderivsurf(). 
* 
* 2   Specification
*     #include "bspl.h"
*     void Aijqs(int p, int q, int k, double u[], int r, int s, int l, 
*                double v[], vector ***contpts, double ***derivptsa)
* 
* 3   Description
*     This function calculates the control points of the derivative surface
*     of a B-spline surface. The derivative surface itself is also a B-spline 
*     surface with lower order. The control points of the derivative surface 
*     are recursively determined by the control points of the surface.
* 
* 4   References
*     [1]C. De Boor. A Practical Guide to Splines, Springer, New York, 1978.
* 
* 5   Parameters
*       1.int p
*         On entry: the index of the knot such that uknots[p] <= u <=
*                 uknots[p+1], where uknots[] is the knot vector in u
* 		  direction, u is the  parametric value to be evaluated.
*       2.int q
*         On entry: the order of the derivative w.r.t. u.
*       3.int k
*         On entry: the order in u direction of the surface to be evaluated.
*       4.double u[]
*         On entry: the knot vector in u direction of the surface to be
*                 evaluated.
*       5.int r
*         On entry: the index of the knot such that vknots[r] <= v <=
*                 vknots[r+1], where vknots[] is the knot vector in v
*                 direction, v is the parametric value to be evaluated.
*       6.int s
*         On entry: the order of the derivative w.r.t. v.
*       7.int l
*         On entry: the order in v direction of the surface to be evaluated.
*       8.double v[]
*         On entry: the knot vector in v direction of the surface to be
*                 evaluated.
*       9.vector ***contpts
*         On entry: the control points of the surface to be evaluated.
*       6.double ***derivptsa
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
* 9   Functions referenced by Aijqs() are:
*     errormsg()
* 
* 10  Functions that reference Aijqs() are:
*     evalderivsurf()
*     evaldersurf()
* 
*****************************************************************************/

void Aijqs(int p, int q, int k, double u[], int r, int s, int l, double v[],
	   vector ***contpts, double ***derivptsa)
/************************************************************************ 
 * The Aij are such that A sub(i) sup(j) = A[i-il+(k-1)][j]
 * Although the variable il is not available here, this value is useful
 * in indexing A[][] compactly and to calculate only Aij required for a 
 * a given il, see program evalderivbsp() where Aij is used as such. 
 **********************************************************************/
{
  vector *v1, *v2, *v3;          /* working vectors for temporarily carrying 
                                  * a control point
	        	       	  * v1,v3 are never used (glshen) */
  double du,dv,temp;             /* denominators in the derivative formula,
                                  * which are differences between two knots
                                  * in u,v directions respectively. */  
  int i,j,ireal,jreal,m,n,iv;    /* indices */

  /* for 0-th derivative, the control points are just the original ones. */
  for (i=0; i<k; i++)
      for (j=0; j<l; j++)
	  {
	   v2 = contpts[i+p-(k-1)][j+r-(l-1)];
	   derivptsa[i+k*j][0][0] = v2->x;
	   derivptsa[i+k*j][0][1] = v2->y;
	   derivptsa[i+k*j][0][2] = v2->z;
	   derivptsa[i+k*j][0][3] = v2->w;
	   }

  /* calculate the control points of the n-th derivative w.r.t. v, n varying 
   * from 1 to s. */
  for (n=1; n<=s; n++)
      for (i=0; i<k; i++)
	  for (j=n; j<l; j++)
	      {
	       jreal = j+r-(l-1);
               dv = v[jreal+l-n] - v[jreal]; 
	       if (dv < ZERO1)
		   errormsg(1,"dv is ZERO1 in routine Aij()");	
	       for (iv=0; iv<4; iv++){
		   temp = derivptsa[i+k*j][k*(n-1)][iv] - 
		          derivptsa[i+k*(j-1)][k*(n-1)][iv];
		   derivptsa[i+k*j][k*n][iv] = temp/dv;
		   }
	      }

  /* calculate the control points of the m-th derivative w.r.t. u of the n-th 
   * derivative with respect to v, i.e., the result is the mixed (m+n)-th
   * partial derivative with order m in u and n in v. */
  for (m=1; m<=q; m++)
      for (i=m; i<k; i++)
	  for (j=0; j<l; j++)
	      {
	       ireal = i+p-(k-1);
               du = u[ireal+k-m] - u[ireal]; 
 	       if (du < ZERO1)
		  errormsg(1,"du is ZERO1 in routine Aij()");	
	       for (iv=0; iv<4; iv++){
		   temp = derivptsa[i+k*j][m-1+k*s][iv] - 
		          derivptsa[i-1+k*j][m-1+k*s][iv];
		   derivptsa[i+k*j][m+k*s][iv] = temp/du;
		   }
	      }
}

/****************************************************************************
*                                  Aijqs_z()
*****************************************************************************
* 
* 1   Purpose
*     This is a subroutine of function zevalderivsurf(), which evaluates z-
*     component of the derivatives of a B-spline surface, see zevalderivsurf().
* 
* 2   Specification
*     #include "bspl.h"
*     void Aijqs_z(int p, int q, int k, double u[], int r, int s, int l, 
*                double v[], vector ***contpts, double ***derivptsa)
* 
* 3   Description
*     This function calculates z-components of the control points of the 
*     derivative surface of a B-spline surface. The derivative surface itself
*     is also a B-spline surface with lower order. The z-components of the
*     control points of the derivative surface are recursively determined
*     by z-components of the control points of the surface.
* 
* 4   References
*     [1]C. De Boor. A Practical Guide to Splines, Springer, New York, 1978.
* 
* 5   Parameters
*       1.int p
*         On entry: the index of the knot such that uknots[p] <= u <=
*                   uknots[p+1], where uknots[] is the knot vector in u
*                   direction, u is the  parametric value to be evaluated.
*       2.int q
*         On entry: the order of the derivative w.r.t. u.
*       3.int k
*         On entry: the order in u direction of the surface to be evaluated.
*       4.double u[]
*         On entry: the knot vector in u direction of the surface to be
*                   evaluated.
*       5.int r
*         On entry: the index of the knot such that vknots[r] <= v <=
*                   vknots[r+1], where vknots[] is the knot vector in v
*                   direction, v is the  parametric value to be evaluated.
*       6.int s
*         On entry: the order of the derivative w.r.t. v.
*       7.int l
*         On entry: the order in v direction of the surface to be evaluated.
*       8.double v[]
*         On entry: the knot vector in v direction of the surface to be
*                   evaluated.
*       9.vector ***contpts
*         On entry: the control points of the surface to be evaluated.
*       6.double ***derivptsa
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
* 9   Functions referenced by Aijqs_z() are:
*     errormsg()
* 
* 10  Functions that reference Aijqs_z() are:
*     evalderivsurf()
*     evaldersurf()
* 
*****************************************************************************/

void Aijqs_z(int p, int q, int k, double u[], int r, int s, int l,
	     double v[], vector *** contpts, double ***derivptsa)
/************************************************************************ 
 * This is a modification of routine Aijqs and computes only the z component.
 * The Aij are such that A sub(i) sup(j) = A[i-il+(k-1)][j]
 * Although the variable il is not available here, this value is useful
 * in indexing A[][] compactly and to calculate only Aij required for a 
 * a given il, see program evalderivbsp() where Aij is used as such. 
 **********************************************************************/
{
  vector *v2;          /* working vector temporarily containing a control
			* point */
  double du,dv,temp;   /*denominators in the derivative formula,
                        * which are differences between two knots in u,v
			* directions respectively. */   
  int i,j,ireal,jreal,m,n;  /* indices */

  /* for 0-th derivative, the control points are just original ones. */
  for(i=0; i<k; i++)
	for(j=0; j<l; j++)
		{
		v2 = contpts[i+p-(k-1)][j+r-(l-1)];
		derivptsa[i+k*j][0][2] = v2->z;
		}

  /* calculate z-components of the control points of n-th partial derivative
   * with respect to v */
  for(n=1; n<=s; n++)
	for(i=0; i<k; i++)
		for(j=n; j<l; j++)
		{
		jreal = j+r-(l-1);
           	dv = v[jreal+l-n] - v[jreal]; 
		if(dv < ZERO1)
		  errormsg(1,"dv is ZERO1 in routine Aij()");	
		temp = derivptsa[i+k*j][k*(n-1)][2] - 
		  derivptsa[i+k*(j-1)][k*(n-1)][2];
		derivptsa[i+k*j][k*n][2] = temp/dv;
		}

  /* calculate z-components of the control points of (m+n)-th mixed derivative
   * with order m in u and n in v. */
  for(m=1; m<=q; m++)
	for(i=m; i<k; i++)
		for(j=0; j<l; j++)
		{
		ireal = i+p-(k-1);
           	du = u[ireal+k-m] - u[ireal]; 

		if(du < ZERO1)
		  errormsg(1,"du is ZERO1 in routine Aij()");	
		temp = derivptsa[i+k*j][m-1+k*s][2] - 
		  derivptsa[i-1+k*j][m-1+k*s][2];
		derivptsa[i+k*j][m+k*s][2] = temp/du;
		}
}

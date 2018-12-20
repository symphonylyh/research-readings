/* ***************************************************************************
 Copyright (C) 1996 Massachusetts Institute of Technology all rights reserved 
	Programmer: George A. Kriezis
**************************************************************************** */
#include <malloc.h>
#include "gen.h"
#include "bspl.h"

static double ***derivpts;
static int NUMBERPOINTS = 0; /* this is so that memory is only allocated
                              * a single time. That is, when 0, memory is
                              * allocted and the variable set to the length
                              * of the array  derivpts. Subsequent calls will
                              * check this against the current required
                              * length, if greater, will reallocate derivpts
                              * if less, will not. */
#define ZERO1 1.e-10

/******************************************************************************
*                                 normalsurf1()
*******************************************************************************
* 
* 1   Purpose
*     This function returns the computed unit vector normal of a NURBS surface.
* 
* 2   Specification
*     #include "bspl.h"
*     vector *normalsurf1(ParSurf *fgeom, double u, double v)
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
*      whose w = 1.
* 
* 4   References
*     [1] P. G. Alourdas. Shape Creation, Interrogation and Fairing Using 
*         B-Splines, Engineer's Thesis, Massachusetts Institute of Technology, 
* 	Department of Ocean Engineering, Cambridge, Massachusetts, 1989.
*     [2] I. D. Faux and M. J. Pratt.  Computational Geometry for Design and 
*         Manufacture, Ellis Horwood: Chichester, England, 1979.
* 
* 5   Parameters
*       1.ParSurf * fgeom
*         On entry: NURBS surface data structure containing geometry of 
* 	          surface that is to be evaluated.
*       2.double u
*         On entry: parameter value at which to compute the surface normal.
*       3.double v
*         On entry: parameter value at which to compute the surface normal.
* 
* 6   Return Values, Error Indicators and Warnings
*     Not applicable.
* 
* 7   Accuracy
*     Not applicable.
* 
* 8   Further Comments
*     This is a modified and mor efficient version of normalsurf (bspl -
*     evalsurf.c).
* 
* 9   Functions referenced by normalsurf1() are:
*     cross()
*     evalrsurf()
*     mag()
*     unitvector1()
* 
* 10  Functions that reference normalsurf1() are:
*     CalcIsophotes()
*     curvature1()
*     fcn()
* 
*******************************************************************************/

vector *normalsurf1(ParSurf *fgeom, double u, double v)
/* evaluate normal to parametric surface */
{
  vector *ru,*rv,*ruv,*n,*vec[4];  /* derivatives and normal */
  int i;

  /* evaluate the position, 1st partial derivative w.r.t. u and v, 2nd mixed
     partial derivative at parametric values u and v. */
  evalrsurf(fgeom,u,v,3,vec);
  ru = vec[1];           /* the 1st partial derivative w.r.t. to u */
  rv = vec[2];           /* the 1st partial derivative w.r.t. to v */

  if(mag(ru) < ZERO1) { /* if ru is zero within tolerance, 
			 * use ruv to calculate the normal */
    ruv = vec[3];       /* the 2nd mixed partial derivative */
    n = cross(ruv,rv);	/* normal is the cross product of rv and ruv */
  }
  else if(mag(rv) < ZERO1) { /* if rv is zero within tolerance,
			      * use ruv to calculate the normal */
    ruv = vec[3];
    n = cross(ru,ruv);       /* normal is the cross product of ru and ruv */
  }
  else                       /* if both ru and rv are not zero */
    n = cross(ru,rv);        /* normal is the cross product of ru and rv */

  unitvector1(n,n);
  for (i=0; i<4; i++)
    free((char *)vec[i]);
  return(n);
}

/******************************************************************************
*                                 evalrsurf()
*******************************************************************************
* 
* 1   Purpose
*     This function places the point and derivatives of a NURBS surface into 
*     an array of vectors.
* 
* 2   Specification
*     #include "bspl.h"
*     void evalrsurf(ParSurf *fgeom, double u, double v, int ider,
*                    vector **vec)
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
*                                    ---------------
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
* 
*                                        i
*                                       @ R(u,v)
* 				   ----------------
*                                       uderiv  vderiv
* 				    @u      @v
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
*     This is a modified and more efficient version of subroutine 
*     revalderivsurf (bspl _ evalsurf.c).
*     See also eval_surface_bounded (bspl _ evalsurf.c).
* 
* 9   Functions referenced by evalrsurf() are:
*     dbl_array2()
*     evaldersurf()
*     evalsurf1()
*     find()
*     free_darray2()
*     nbasisd()
*  
* 10  Functions that reference evalrsurf() are:
*     axis_and_terminate()
*     build_offset()
*     CalcIsophotes()
*     CalcReflectionLines()
*     compute_diffeq()
*     compute_ortho_proj()
*     ConstrainedCB()
*     curvature1()
*     DrawFacetEvaluateMap()
*     DrawIsophotes()
*     DrawLinesOfCurvature()
*     DrawReflectionLines()
*     DrawSurfEvaluate()
*     DrawSurfPivot()
*     EvaluateSurf()
*     eval_cylsect_param()
*     eval_cylsect_param3()
*     eval_surface_bounded()
*     fcn()
*     fdist()
*     find_estim()
*     find_estim2()
*     find_init()
*     fn()
*     GetSurfValues()
*     main()
*     MinDistance()
*     normaldu()
*     normalsurf1()
*     offset_geod()
*     offset_geod_par()
*     offset_normal()
*     OrthoDistance()
*     PostScriptFacetEvaluateMap()
*     PostScriptSurfEvaluate()
*     ReadDeslabLocal()
*     ReadDeslabMinDist()
*     SubDistance()
*     SurfPivotUpdate()
*     surf_data()
*     Unconloc()
*
****************************************************************************/

void evalrsurf(ParSurf *fgeom, double u, double v, int fin, vector *vec[])
{
  double **N1,**M1;           /* matrices containing B-spline basis in u,v
                               * directions respectively. */
  double w,w1u,w1v,w2u,w2v,w2uv,w3u,w3v,wuuv,wuvv,w2u2v;
                              /* scaling factors in transformation from
			       * homogeneous to cartesian coordinates. */
  int ilu,ilv;                /* indices */
  
  /* allocate memories for N1 and M1. */
  N1 = dbl_array2((unsigned)(fgeom->uorder+1),(unsigned)(fgeom->uorder+1));
  M1 = dbl_array2((unsigned)(fgeom->vorder+1),(unsigned)(fgeom->vorder+1));

  /* find the index ilu, such that uknots[ilu] <= u <= uknots[ilu+1] */
  ilu = find(fgeom->ucontpts,fgeom->uknots,u);
  /* calculate the matrix of B-spline basis in u direction */
  nbasisd(ilu, fgeom->uorder, u, fgeom->uknots, N1); 
	
  /* find the index ilv, such that vknots[ilv] <= v <= vknots[ilu+1] */
  ilv = find(fgeom->vcontpts,fgeom->vknots,v);
  /* calculate the matrix of B-spline basis in v direction */
  nbasisd(ilv, fgeom->vorder, v, fgeom->vknots, M1); 

  /* if the index of derivative is no less than 0, calculate the position with 
     parametric values u & v. */
  if(fin>=0) {
    /* evaluate the point with parametric values u & v on the surface by
     * calling evalsurf1() */
    vec[0] = evalsurf1(fgeom,ilu,ilv,N1,M1);
    w = vec[0]->w;
    /* transform from homogeneous to cartesian coordinates. */
    vec[0]->x = vec[0]->x/w;      vec[0]->y = vec[0]->y/w;
    vec[0]->z = vec[0]->z/w;      vec[0]->w = 1.0;
  }

  /* if the index of derivative is no less than 1, calculate the 1st partial
   * derivative w.r.t. u */
  if(fin >= 1) {
    /* evaluate the 1st partial derivative w.r.t. u by calling evaldersurf() */
    vec[1] = evaldersurf(fgeom,1,0,ilu,ilv,N1,M1);
    w1u = vec[1]->w;
    /* transform from cartesian to homogeneous coordinates */
    vec[1]->x = (vec[1]->x - vec[0]->x*w1u)/w; 
    vec[1]->y = (vec[1]->y - vec[0]->y*w1u)/w;
    vec[1]->z = (vec[1]->z - vec[0]->z*w1u)/w;      vec[1]->w = 1.0;
  }

  /* if the index of derivative is no less than 2, calculate the 1st partial
   *derivative w.r.t. v */
  if(fin >= 2) {
    /* evaluate the 1st partial derivative w.r.t. v by calling evaldersurf() */
    vec[2] = evaldersurf(fgeom,0,1,ilu,ilv,N1,M1);
    w1v = vec[2]->w;
    /* transform from homogeneous to cartesian coordinates */
    vec[2]->x = (vec[2]->x-vec[0]->x*w1v)/w;
    vec[2]->y =(vec[2]->y-vec[0]->y*w1v)/w;
    vec[2]->z = (vec[2]->z - vec[0]->z*w1v)/w;	vec[2]->w = 1.0;
  }

  /* if the index of derivative is no less than 3, calculate the 2nd mixed
   * partial derivative w.r.t. u & v */
  if(fin >= 3) {
    /* evaluate the 2nd mixed partial derivative w.r.t. u & v by calling 
       evaldersurf() */
    vec[3] = evaldersurf(fgeom,1,1,ilu,ilv,N1,M1);
    w2uv = vec[3]->w;
    /* transform from homogeneous to cartesian coordinates */
    vec[3]->x = (vec[3]->x - w1u*vec[2]->x - w1v*vec[1]->x - w2uv*vec[0]->x)/w;
    vec[3]->y = (vec[3]->y - w1u*vec[2]->y - w1v*vec[1]->y - w2uv*vec[0]->y)/w;
    vec[3]->z = (vec[3]->z - w1u*vec[2]->z - w1v*vec[1]->z - w2uv*vec[0]->z)/w;
    vec[3]->w = 1.0;
  }

  /* if the index of derivative is no less than 4, calculate the 2nd partial 
   * derivative w.r.t. u */
  if(fin >= 4) {
    /* evaluate the 2nd partial derivative w.r.t. u by calling evaldersurf() */
    vec[4] = evaldersurf(fgeom,2,0,ilu,ilv,N1,M1);
    w2u = vec[4]->w;
    /* transform from homogeneous to cartesiancoordinates */
    vec[4]->x = (vec[4]->x - 2.0*w1u*vec[1]->x - w2u*vec[0]->x)/w;
    vec[4]->y = (vec[4]->y - 2.0*w1u*vec[1]->y - w2u*vec[0]->y)/w;
    vec[4]->z = (vec[4]->z - 2.0*w1u*vec[1]->z - w2u*vec[0]->z)/w;
    vec[4]->w = 1.0;
  }

  /* if the index of derivative is no less than 5, calculate the 2nd partial
     derivative w.r.t. v */
  if(fin >= 5) {
    /* evaluate the 2nd partial derivative w.r.t. v by calling evaldersurf() */
    vec[5] = evaldersurf(fgeom,0,2,ilu,ilv,N1,M1);
    w2v = vec[5]->w;
    /* transform from homogeneous to cartesian coordinates */
    vec[5]->x = (vec[5]->x - 2.0*w1v*vec[2]->x - w2v*vec[0]->x)/w;
    vec[5]->y = (vec[5]->y - 2.0*w1v*vec[2]->y - w2v*vec[0]->y)/w;
    vec[5]->z = (vec[5]->z - 2.0*w1v*vec[2]->z - w2v*vec[0]->z)/w;
    vec[5]->w = 1.0;
  }

  /* if the index of derivative is no less than 6, calculate the 3rd partial
     derivative w.r.t. u */
  if(fin >= 6) {
    /* evaluate the 3rd partial derivative w.r.t. u by calling evaldersurf() */
    vec[6] = evaldersurf(fgeom,3,0,ilu,ilv,N1,M1);
    w3u = vec[6]->w;
    /* transform from homogeneous to cartesiancoordinates */
    vec[6]->x = (vec[6]->x-3.0*w1u*vec[4]->x-3.0*w2u*vec[1]->x -
		 w3u*vec[0]->x)/w;
    vec[6]->y = (vec[6]->y-3.0*w1u*vec[4]->y-3.0*w2u*vec[1]->y -
		 w3u*vec[0]->y)/w;
    vec[6]->z = (vec[6]->z-3.0*w1u*vec[4]->z-3.0*w2u*vec[1]->z -
		 w3u*vec[0]->z)/w;
    vec[6]->w = 1.0;
  }

  /* if the index of derivative is no less than 7, calculate the 3rd mixed
     partial derivative with order 2 in u and order 1 in v. */
  if(fin >= 7) {
    /* evaluate the 3rd mixed partial derivative with order 2 in u and 1 in
     * v */
    vec[7] = evaldersurf(fgeom,2,1,ilu,ilv,N1,M1);
    wuuv = vec[7]->w;
    /* transform from homogeneous to cartesian coordinates */
    vec[7]->x = (vec[7]->x - 2.0*w1u*vec[3]->x - w1v*vec[4]->x -
		 2.0*w2uv*vec[1]->x - w2u*vec[2]->x - wuuv*vec[0]->x)/w;
    vec[7]->y = (vec[7]->y - 2.0*w1u*vec[3]->y - w1v*vec[4]->y -
		 2.0*w2uv*vec[1]->y - w2u*vec[2]->y - wuuv*vec[0]->y)/w;
    vec[7]->z = (vec[7]->z - 2.0*w1u*vec[3]->z - w1v*vec[4]->z -
		 2.0*w2uv*vec[1]->z - w2u*vec[2]->z - wuuv*vec[0]->z)/w;
    vec[7]->w = 1.0;
  }

  /* if the index of derivative is no less than 8, calculate the 3rd mixed
     partial derivative with order 1 in u and order 2 in v. */
  if(fin >= 8) {
    /* evaluate the 3rd mixed partial derivative with order 1 in u and 2 in
     * v */
    vec[8] = evaldersurf(fgeom,1,2,ilu,ilv,N1,M1);
    wuvv = vec[8]->w;
    /* transform from homogeneous to cartesian coordinates */
    vec[8]->x = (vec[8]->x - 2.0*w1v*vec[3]->x - w1u*vec[5]->x -
		 2.0*w2uv*vec[2]->x - w2v*vec[1]->x - wuvv*vec[0]->x)/w;
    vec[8]->y = (vec[8]->y - 2.0*w1v*vec[3]->y - w1u*vec[5]->y -
		 2.0*w2uv*vec[2]->y - w2v*vec[1]->y - wuvv*vec[0]->y)/w;
    vec[8]->z = (vec[8]->z - 2.0*w1v*vec[3]->z - w1u*vec[5]->z -
		 2.0*w2uv*vec[2]->z - w2v*vec[1]->z - wuvv*vec[0]->z)/w;
    vec[8]->w = 1.0;
  }

  /* if the index of derivative is no less than 9, calculate the 3rd partial
     derivative w.r.t. v */
  if(fin >= 9) {
    /* evaluate the 3rd partial derivative w.r.t. v */
    vec[9] = evaldersurf(fgeom,0,3,ilu,ilv,N1,M1);
    w3v = vec[9]->w;
    /* transform from homogeneous to cartesian coordinates */
    vec[9]->x = (vec[9]->x-3.0*w1v*vec[5]->x-3.0*w2v*vec[2]->x -
		 w3v*vec[0]->x)/w;
    vec[9]->y = (vec[9]->y-3.0*w1v*vec[5]->y-3.0*w2v*vec[2]->y -
		 w3v*vec[0]->y)/w;
    vec[9]->z = (vec[9]->z-3.0*w1v*vec[5]->z-3.0*w2v*vec[2]->z -
		 w3v*vec[0]->z)/w;
    vec[9]->w = 1.0;
  }
  
  /* if the index of derivative is no less than 10, calculate the 4th mixed 
     partial derivative with order 2 in u and order 2 in v. */
  if(fin == 10) {
    /* evaluate the 4th mixed partial derivative with order 2 in u and v */
    vec[10] = evaldersurf(fgeom,2,2,ilu,ilv,N1,M1);
    w2u2v = vec[10]->w;
    /* transform from homogeneous to cartesian coordinates. */
    vec[10]->x = (vec[10]->x-2.0*w1v*vec[7]->x-2.0*w1u*vec[8]->x-w2v*vec[4]->x
		  - 4.0*w2uv*vec[3]->x - w2u*vec[5]->x - 2.0*wuvv*vec[1]->x - 
		  2.0*wuuv*vec[2]->x - w2u2v*vec[0]->x)/w;
    vec[10]->y = (vec[10]->y-2.0*w1v*vec[7]->y-2.0*w1u*vec[8]->y-w2v*vec[4]->y
		  - 4.0*w2uv*vec[3]->y - w2u*vec[5]->y - 2.0*wuvv*vec[1]->y - 
		  2.0*wuuv*vec[2]->y - w2u2v*vec[0]->y)/w;
    vec[10]->z = (vec[10]->z-2.0*w1v*vec[7]->z - 2.0*w1u*vec[8]->z -
		  w2v*vec[4]->z - 4.0*w2uv*vec[3]->z - w2u*vec[5]->z -
		  2.0*wuvv*vec[1]->z - 2.0*wuuv*vec[2]->z - w2u2v*vec[0]->z)/w;
    vec[10]->w = 1.0;
  }

  free_darray2(N1);  free_darray2(M1);
}

/******************************************************************************
*                               evalsurf1()
*******************************************************************************
* 
* 1   Purpose
*     This function returns the vector valued non-uniform integral B-spline 
*     surface point.
* 
* 2   Specification
*     #include "bspl.h"
*     vector *evalsurf1(ParSurf *fgeom, int ilu, int ilv, double **N,
*                       double **M)
* 
* 3   Description
*     This routine returns a vector data structure which contains the point on 
*     a B-spline surface corresponding to the supplied values for ilu and ilv, 
*     and matrices N abd M.
* 
* 4   References
*     [1]C. De Boor. A Practical Guide to Splines, Springer, New York, 1978.
* 
* 5   Parameters
*       1.ParSurf * fgeom
*         On entry: NURBS surface data structure containing geometry of 
* 	          surface to be evaluated.
*       2.int ilu
*         On entry: index referring to the indicial location of the parameter 
* 	          value of u with respect to the array representing the knot 
* 		  vector in the u direction.
*       3.int ilv
*         On entry: index referring to the indicial location of the parameter 
* 	          value of v with respect to the array representing the knot 
* 		  vector in the v direction.
*       4.double ** N
*         On entry: matrix of length (order in the u direction x order in the 
* 	          u direction) which contains the computed B-spline basis 
* 		  values with repsect to the u parameter.
*       5.double ** N
*         On entry: matrix of length (order in the v direction x order in the 
* 	          v direction) which contains the computed B-spline basis 
* 		  values with repsect to the v parameter.
* 
* 6   Return Values, Error Indicators and Warnings
*     The value of u should lie between fgeom->uknots[fgeom->uorder-1] and
*     fgeom->uknots[fgeom->ucontpts].
*     The value of v should lie between fgeom->vknots[fgeom->vorder-1] and
*     fgeom->vknots[vcontpts].
* 
* 7   Accuracy
*     Not applicable.
* 
* 8   Further Comments
*     This modified version of evalsurf (bspl _ evalderivsurf.c) is more 
*     efficient when multiple evaluation calls are made such that the 
*     parameter values do not entail recomputing the basis matrices or that 
*     the basis matrices are computed elsewhere.
*     This basis function evaluation is performed using routine nbasisd (bspl).
*     This routine should not be used to evaluate a NURBS surface. For the 
*     evaluation of NURBS surfaces, see evalrsurf (bspl _ ) for details.
*     For evaluation of ilu and ilv, see routine find (bspl _ util.c).
*     For evaluation of N and N, see routine nbasisd (bspl _ ).
* 
* 9   Functions referenced by evalsurf1() are:
*     vectalloc()
* 
* 10  Functions that reference evalsurf1() are:
*     evalrsurf()
* 
****************************************************************************/

vector *evalsurf1(ParSurf *fgeom, int ilu, int ilv, double **N1, double **M1)
{
  int ucontpts,uorder,vcontpts,vorder;  /* orders and numbers of control points
                                         * in u,v directions */ 
  int i,j;
  vector *p,*eval;       /* working vector temporarily carrying a control point
                          * and the vector containing the evaluated point */
  double q,r;

  ucontpts = fgeom->ucontpts;
  vcontpts = fgeom->vcontpts;
  uorder = fgeom->uorder;
  vorder = fgeom->vorder;

  /* allocate memory for eval, and initialize it. */
  eval = vectalloc();
  eval->x = eval->y = eval->z = eval->w = 0.0;

  /* evaluate the point on the surface */
  for(i=ilu; i>ilu-uorder; i--)
	for(j=ilv; j>ilv-vorder; j--)
		{
		p = fgeom->contpts[i][j];		
		q = N1[i+uorder-ilu][uorder];
		r = M1[j+vorder-ilv][vorder];
		eval->x += q*r*p->x;	
		eval->y += q*r*p->y;	
		eval->z += q*r*p->z;	
		eval->w += q*r*p->w;	
		}
	
  return(eval);
}

/******************************************************************************
*                                evaldersurf()
*******************************************************************************
* 
* 1   Purpose
*     This function returns the vector valued non-uniform integral B-spline 
*     surface point or derivative.
* 
* 2   Specification
*     #include "bspl.h"
*     vector *evaldersurf(ParSurf *fgeom, int uderiv, int vderiv, int ilu, 
*                         int ilv, double **N, double **M)
* 
* 3   Description
*     This routine returns a vector data structure which contains the point or 
*     derivatives of the point on a B-spline surface corresponding to the 
*     supplied values for ilu and ilv, and matrices N abd M.
* 
* 4   References
*     [1]C. De Boor. A Practical Guide to Splines, Springer, New York, 1978.
* 
* 5   Parameters
*       1.ParSurf * fgeom
*         On entry: NURBS surface data structure containing geometry of 
* 	          surface to be evaluated.
*       2.int uderiv
*         On entry: an index referring to the required derivative with respect 
* 	          to u to be calculated such that if the surface is defined as 
* 		  R(u,v) then the routine will evaluate
*                                  uderiv+vderiv
*                                 @             R(u,v)
*                                ----------------------
*                                     uderiv  vderiv
* 				  @u      @v
* 		  ( @ represents partial derivative. )
* 		  where vderiv is defined below.
*       3.int vderiv
*         On entry: an index referring to the required derivative with respect 
* 	          to v to be calculated such that if the surface is defined as 
* 		  R(u,v) then the routine will evaluate
*                                  uderiv+vderiv
*                                 @             R(u,v)
*                                ----------------------
*                                     uderiv  vderiv
*                                   @u      @v
* 		  ( @ represents partial derivative. )
* 		  where uderiv is defined above.
*       4.int ilu
*         On entry: index referring to the indicial location of the parameter 
* 	          value of u with respect to the array representing the knot 
* 		  vector in the u direction.
*       5.int ilv
*         On entry: index referring to the indicial location of the parameter 
* 	          value of v with respect to the array representing the knot 
* 		  vector in the v direction.
*       6.double ** N
*         On entry: matrix of length (order in the u direction x order in the 
* 	          u direction) which contains the computed B-spline basis 
* 		  values with repsect to the u parameter.
*       7.double ** N
*         On entry: matrix of length (order in the v direction x order in the 
* 	          v direction) which contains the computed B-spline basis 
* 		  values with repsect to the v parameter.
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
*     This modified version of evalderivsurf (bspl _ evalderivsurf.c) is more 
*     efficient when multiple evaluation calls are made such that the 
*     parameter values do not entail recomputing the basis matrices or that
*     the basis matrices are computed elsewhere.
*     This basis function evaluation is performed using routine nbasisd (bspl).
*     This routine should not be used to evaluate derivatives of NURBS
*     surfaces. For the evaluation of NURBS surfaces, see evalrsurf (bspl)
*     for details.
*     For evaluation of ilu and ilv, see routine find (bspl _ util.c).
*     For evaluation of N and N, see routine nbasisd (bspl _ ).
* 
* 9   Functions referenced by evaldersurf() are:
*     Aijqs()
*     dbl_array3()
*     free_darray3()
*     vectalloc()
* 
* 10  Functions that reference evaldersurf() are:
*     evalrsurf()
* 
****************************************************************************/

vector *evaldersurf(ParSurf *fgeom, int uderiv, int vderiv, 
		    int ilu, int ilv, double **N1, double **M1)
/* more efficient evaluation routine */
/* evaluates the rational b-spline derivatives using Farouki's offset paper
 * it calculates up to third derivatives */
{
  int ucontpts,uorder,vcontpts,vorder;   /* orders and numbers of control 
					  * points in u,v directions */
  int i,j,ii,jj;                         /* indices */
  int prod;                              /* factor in the derivative formula */
  vector *eval;                          /* the evaluated result */
  double q,r;

  ucontpts = fgeom->ucontpts; vcontpts = fgeom->vcontpts;
  uorder = fgeom->uorder;     vorder = fgeom->vorder;

  /* if the order of derivative is higher than order, set the derivative = 0 */
  if (uorder <= uderiv || vorder <= vderiv){
    eval = vectalloc();
    eval->x = 0.0; eval->y = 0.0; eval->z = 0.0; eval->w = 0.0;
    return(eval);
  }
  
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

  /* calculate the control points of the derivative surface (without factor) */
  Aijqs(ilu,uderiv,uorder,fgeom->uknots,ilv,vderiv,vorder,fgeom->vknots,
	fgeom->contpts, derivpts);

  /* allocate memory for eval, and initialize it. */
  eval = vectalloc();  eval->x = eval->y = eval->z = eval->w = 0.0;
  
  /* calculate the factor (uorder-1)...(uorder-uderiv)(vorder-1)...
   * (vorder-vderive) */
  prod = 1.0;
  for(i=1; i<=uderiv; i++) prod *= uorder-i;
  for(i=1; i<=vderiv; i++) prod *= vorder-i;
  
  /* evaluate the derivative */
  for(i=ilu; i>ilu-uorder+uderiv; i--)
    for(j=ilv; j>ilv-vorder+vderiv; j--) {
      q = N1[i+uorder-uderiv-ilu][uorder-uderiv];
      r = M1[j+vorder-vderiv-ilv][vorder-vderiv];
      ii = i-ilu+(uorder-1);      jj = j-ilv+(vorder-1);
      ii = ii + uorder*jj;        jj = uderiv + uorder*vderiv;
      eval->x  += q*r*derivpts[ii][jj][0]; eval->y  += q*r*derivpts[ii][jj][1];
      eval->z  += q*r*derivpts[ii][jj][2]; eval->w  += q*r*derivpts[ii][jj][3];
    }
  
  /* multiply by the factor */
  eval->x = prod*eval->x;  eval->y = prod*eval->y;
  eval->z = prod*eval->z;  eval->w = prod*eval->w;

  return(eval);
}

/*****************************************************************************
*                                zevaldsurf()
******************************************************************************
* 
* 1   Purpose
*     This routine is a modified version of evaldsurf(). Instead of returning 
*     the derivative, it returns z-component of the derivative. 
* 
* 2   #include "bspl.h"
*     double zevaldsurf(ParSurf *fgeom, int uderiv, int vderiv, int ilu,
*                       int ilv, double **N1, double **M1)
* 
* 3   Description
*     see Description of evaldsurf().
* 
* 4   References
*     Not applicable
*  
* 5   Parameters
*     see evaldsurf().
* 
* 6   Return Values, Error Indicators and Warnings
*     This function a double value of z-component of the derivative.
* 
* 7   Accuracy
* 
* 8   Further Comments
*     See also evaldsurf().
* 
* 9   Functions referenced by bernrs() are:
*     Aijqs_z()
*     dbl_array3()
*     free_darray3()
* 
* 10  Functions that reference bernrs() are: None
* 
****************************************************************************/

double zevaldsurf(ParSurf *fgeom, int uderiv, int vderiv, int ilu, int ilv,
		  double **N1, double **M1)
{
  int ucontpts,uorder,vcontpts,vorder;  /* order and number of control points 
					 * in u,v directions */ 
  int i,j,ii,jj;                        /* indices */
  int prod;                             /* factor in the derivative formula 
					 * of a B-spline surface */
  double eval,q,r;                      /* eval is to contain z-component of 
					 * the derivative */

  ucontpts = fgeom->ucontpts;
  vcontpts = fgeom->vcontpts;
  uorder = fgeom->uorder;
  vorder = fgeom->vorder;

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

  /* calculate z-components of the control points of the derivative surface
     (without factor) */
  Aijqs_z(ilu,uderiv,uorder,fgeom->uknots,ilv,vderiv,vorder,fgeom->vknots,
	  fgeom->contpts, derivpts);

  /* calculate the factor (uorder-1)...(uorder-uderiv)(vorder-1)...
   * (vorder-vderiv) */
  prod = 1.0;
  for(i=1; i<=uderiv; i++)
	prod *= uorder-i;
  for(i=1; i<=vderiv; i++)
	prod *= vorder-i;

  eval = 0.0;

  /* evaluate z-component of the derivative */
  for(i=ilu; i>ilu-uorder+uderiv; i--)
	for(j=ilv; j>ilv-vorder+vderiv; j--)
		{
		q = N1[i+uorder-uderiv-ilu][uorder-uderiv];
		r = M1[j+vorder-vderiv-ilv][vorder-vderiv];
		ii = i-ilu+(uorder-1);
		jj = j-ilv+(vorder-1);
		ii = ii +uorder*jj;
		jj = uderiv + uorder*vderiv;
		eval += q*r*derivpts[ii][jj][2];
		}
	
  /* multiply by the factor */
  eval = prod*eval;

  return(eval);
}

/***************************************************************************
*                        eval_surface_bounded()
****************************************************************************
* 
1   Purpose
*     This function places the point and derivatives of a NURBS surface into 
*     an array of vectors. It is the recommended routine for use with the NAG 
*     library minimization routines.
* 
* 2   Specification
*     #include "bspl.h"
*     void eval_surface_bounded(vector **vec, ParSurf *fgeom, int ider, 
*                               double *xc)
* 
* 3   Description
*     This routine places into a supplied array of vector data structures the 
*     point and derivatives of the point of a NURBS surface corresponding to 
*     the parameters (u,v). It is very good if used with the NAG minimization 
*     routines, since it ensures the evaluation of the the function and deri-
*     vatives at a feasible point inside the domain.  If a point is outside the
*     domain, linear extrapolation is performed to evaluate the function and 
*     derivatives at this point (extrapolation from points within the domain).
* 
* 4   References
*     [1]C. De Boor. A Practical Guide to Splines, Springer, New York, 1978.
* 
* 5   Parameters
*       1.ParSurf * fgeom
*         On entry: NURBS surface data structure containing geometry of 
* 	          surface to be evaluated.
*       2.double * xc
*         On entry: the parameter values (xc[0] = u and xc[1] = v) at which 
* 	          the non-uniform integral B-spline surface is to be evaluated.
*       3.int ider
*         On entry: an index referring to the required derivative with respect 
* 	          to u to be calculated such that if the surface is defined 
* 		  as R(u, v) then the routine will evaluate all derivatives 
* 		  up to and including:
*                                  ider
*                                 @    R(u,v)
*                             -------------------
*                                 uderiv   vderiv
*                               @u       @v
* 		  ( @ represents partial derivatives )
* 
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
* 
*        4.vect ** vec
*          On entry: an array of pointers of length  ider.
*          On exit: the vector values of the surface evaluated at the given
* 	          data points. The ith element of the array corresponds to:
*                                    i
*                                   @ R(u,v)
*                               -----------------
*                                   uderiv  vderiv
*                                 @u      @v
* 		   ( @ represents partial derivative )
*                    where uderiv and vderiv are defined as:
* 		      uderiv = 0 and vderiv = 0 when i = 0.
* 		      uderiv = 1 and vderiv = 0 when i = 1.
* 		      uderiv = 0 and vderiv = 1 when i = 2.
* 		      uderiv = 1 and vderiv = 1 when i = 3.
* 		      uderiv = 2 and vderiv = 0 when i = 4.
* 		      uderiv = 0 and vderiv = 2 when i = 5.
* 		      uderiv = 3 and vderiv = 0 when i = 6.
* 		      uderiv = 2 and vderiv = 1 when i = 7.
* 		      uderiv = 1 and vderiv = 2 when i = 8.
* 		      uderiv = 0 and vderiv = 3 when i = 9.
* 
* 6   Return Values, Error Indicators and Warnings
*     The index uderiv <= the order with respect to u of the non-uniform 
*     integral B-spline surface.
*     The index vderiv <= the order with respect to u of the non-uniform 
*     integral B-spline surface.
* 
* 7   Accuracy
*     Not applicable.
* 
* 8   Further Comments
*     This basis function evaluation is performed using routine nbasisd 
*     (bspl _ cxdeb.c). This is a modified and more robust version of
*     subroutine evalrsurf (bspl _ evalsurf.c).
* 
* 9   Functions referenced by eval_surface_bounded() are:
*     evalrsurf()
*     scale_vect1()
*     sub_vect1()
* 
* 10  Functions that reference eval_surface_bounded() are:
*     fdist_surf()
*     fdist_surf1()
*     fdist_surf2()
*     footpoint_fcn()
*     localize_diagnostic()
*     localize_sumsq()
*     localize_sumsq_opt()
*     point_eval()
* 
*****************************************************************************/

/* TO IMPROVE PERFORMANCE WITH NAG MINIMIZATION ROUTINES */
void eval_surface_bounded(vector *r[], ParSurf *fgeom1, int index, double xc[])
/* This routine evaluates the function and derivatives of a b-Spline
 * surface in a bounded domain. It is used with the NAG minimization
 * routines and ensures that we evaluate the surface at a feasible
 * point in the domain. If a point is outside the domain linear
 * extrapolation is performed to evaluate the function and derivatives
 * at this point. (extrapolate from points within the domain)  
 * index = number use in evalrsurf*/
{
  int i;
  vector *r1[11];      /* will contain the evaluated derivatives. */
  double help1,help2;  /* minimum/maximum knot in u,v directions respectively,
			* depends on which side u,v are out of the range. */
  double eps,del;      /* the distance in u,v directions, which measures how 
		        * far the parametric values u,v are out of the range */

/* evaluate positions and derivatives up to order 3 */

  /* check exceeding bounds and correct evaluation*/
  if (xc[0] < fgeom1->uknots[fgeom1->uorder-1] ||
      xc[0] > fgeom1->uknots[fgeom1->ucontpts] ||
      xc[1] < fgeom1->vknots[fgeom1->vorder-1] ||
      xc[1] > fgeom1->vknots[fgeom1->vcontpts])
     {  /* if u and/or v do exceed the range */
        if (xc[0] < fgeom1->uknots[fgeom1->uorder-1])
	   { /* u is left to the range, linear extrapolation will be from
	      * a point inside the range uknots[uorder-1]+eps, via
	      * uknots[uorder-1] */
            eps = fgeom1->uknots[fgeom1->uorder-1] - xc[0];
            help1= fgeom1->uknots[fgeom1->uorder-1];
           }
        else if (xc[0] > fgeom1->uknots[fgeom1->ucontpts])
	  { /* u is right to the range, linear extrapolation will be from
	     * a point inside the range uknots[ucontpts]+eps, via
	     * uknots[ucontpts] */
	    eps = fgeom1->uknots[fgeom1->ucontpts] - xc[0];
	    help1 = fgeom1->uknots[fgeom1->ucontpts];
	  }
	else{ /* u is inside the range */
	  eps = 0.0;
	  help1 = xc[0];
	}
        if (xc[1] < fgeom1->vknots[fgeom1->vorder-1])
	   { /* v is left to the range, linear extrapolation will be from 
	      * a point inside the range vknots[vorder-1]+del, via
	      * vknots[vorder-1] */
            del = fgeom1->vknots[fgeom1->vorder-1] - xc[1];
            help2= fgeom1->vknots[fgeom1->vorder-1];
           }
        else if (xc[1] > fgeom1->vknots[fgeom1->vcontpts])
	  { /* v is right to the range, linear extrapolation will be from
	     * a point inside the range vknots[vcontpts]+del, via
	     * vknots[vcontpts] */
	    del = fgeom1->vknots[fgeom1->vcontpts] - xc[1];
	    help2 = fgeom1->vknots[fgeom1->vcontpts];
	  }
	else{ /* v is inside the range */
	  del = 0.0;
	  help2 = xc[1];
	}

       /* evaluate linearly if outside domain interpolation */
          evalrsurf(fgeom1,help1,help2,index,r);   /* evaluate the point on
						    * the boundary */
        evalrsurf(fgeom1,help1+eps,help2+del,index,r1);
	                             /* evaluate the point inside the range */
        /* do linear extrapolation */
        for (i=0; i<index+1; i++){
            scale_vect1(2.0,r[i],r[i]);
            sub_vect1(r[i],r1[i],r[i]);
            free((char *)r1[i]);
           }
     }
  else  /* if both u & v are inside the range, evaluate the derivatives by 
	   simply calling evalrsurf(). */
    evalrsurf(fgeom1,xc[0],xc[1],index,r);
}

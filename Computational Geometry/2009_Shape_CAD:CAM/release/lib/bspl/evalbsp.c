/* ***************************************************************************
 Copyright (C) 1996 Massachusetts Institute of Technology all rights reserved 
	Programmer: George A. Kriezis
**************************************************************************** */
#include <stdlib.h>
#include <malloc.h>
#include "gen.h"
#include "bspl.h"
#define ZERO1 1.e-10

/*****************************************************************************
*                                evalbsp()
******************************************************************************
* 
* 1   Purpose
*     This function evaluates a non-uniform integral B-spline curve.
* 
* 2   Specification
*     #include "bspl.h"
*     vector *evalbsp(ParCurv *egeom, double u)
* 
* 3   Description
*     This routine calculates and returns a vector data structure that contains 
*     the point value of the given non-uniform integral B-spline curve at the 
*     specific parameter value.
* 
* 4   References
*     [1]C. De Boor. A Practical Guide to Splines, Springer, New York, 1978.
* 
* 5   Parameters
*       1.ParCurv * egeom
*         On entry: NURBS data structure containing geometry of curve to be 
* 	          evaluated.
*         On exit:
*       2.double u
*         On entry: the parameter value at which the non-uniform integral 
* 	          B-spline curve is to be evaluated.
* 
* 6   Return Values, Error Indicators and Warnings
*     The value of u should lie between egeom -> knots[0] and egeom -> knots[k] 
*     where k = (egeom -> ncontpts) + (egeom -> order) - 1.  If not, it will 
*     print an error message and stop.
* 
* 7   Accuracy
*     Not applicable.
* 
* 8   Further Comments
*     This basis function evaluation is performed using routine nbasisd (bspl _ ).
*     This routine should not be used to evaluate derivatives of NURBS curves. 
*     For the evaluation of NURBS curves, see rbspeval (bspl _ ) for details.
* 
* 9   Functions referenced by evalbsp() are:
*     dbl_array2()
*     errormsg()
*     find()
*     free_darray2()
*     nbasisd()
*     vectalloc()
* 
* 10  Functions that reference evalbsp() are:
*     arc_length()
*     CamberBisect()
*     CamberBisect3d()
*     checkoff()
*     ContourFacets()
*     CurveOrientation()
*     DrawCosEvaluate()
*     DrawCosKnots()
*     DrawCurvKnots()
*     DrawTrimWireframe()
*     EvaluateCos()
*     EvaluateCurv()
*     funccurve()
*     GetCosValues()
*     GetCurvValues()
*     IsoLoopIntersect()
*     LinearCosNormal()
*     localize_sumsq_2d()
*     MinDistCurvCB()
*     normalcurv()
*     PostScriptCosEvaluate()
*     PostScriptCosKnots()
*     PostScriptCurvKnots()
*     PostScriptTrimWireframe()
*     ps_ParCurv()
*     rbspeval()
*     ReadDeslabLocal()
*     ReadDeslabMinDist()
*     RobustCurvCB()
*     sample()
*     SampleCurvCB()
*     sample_proj()
*     SaveIgesCos()
*     SaveIgesCurv()
*     TrimPoints()
*     Unconstrained2dCB()
*   
*****************************************************************************/

vector *evalbsp(ParCurv *egeom, double u)
{
  int ncontpts,order;  /* number of control points and order */
  int i,il;            /* indices */
  vector *p,*eval;     /* control point and evaluated point */
  double q;            /* B-spline basis */
  double **N;          /* matrix containing B-spline basis */

  /* allocate memory */
  N = dbl_array2((unsigned)(egeom->order+1),(unsigned)(egeom->order+1));
  ncontpts = egeom->ncontpts;
  order = egeom->order;

  /*printf("%f\n",u);*/

  if (u > egeom->knots[ncontpts] || u < egeom->knots[egeom->order-1])
  /* if u is out of the range of the knot vector. */
	{
	errormsg(3,"outside  the knot range in function evalbsp");
	}
  else  
	{
	il = find(ncontpts,egeom->knots,u); /* find index il such that
                                             * knots[il] <= u <= knots[il+1] */
	
	nbasisd(il, order, u, egeom->knots, N);  /* evaluate the B-spline
						  * basis */
	
	eval = vectalloc();
	eval->x = eval->y = eval->z = eval->w = 0.0;
 
        /* evaluate the point with parametric value u */
	for(i=il; i>il-order; i--)   /* summation over control points */
		{
		p = egeom->contpts[i];		
		q = N[i+order-il][order];
		eval->x += q*p->x;	
		eval->y += q*p->y;	
		eval->z += q*p->z;	
		eval->w += q*p->w;	
		}
	}

  free_darray2(N);

  return(eval);
}

/****************************************************************************
*                            evalderivbsp()
*****************************************************************************
* 
* 1   Purpose
*     This function evaluates a non-uniform integral B-spline curve or 
*     derivative.
* 
* 2   Specification
*     #include "bspl.h"
*     vector *evalderivbsp(ParCurv *egeom, double u, int deriv)
* 
* 3   Description
*     This routine calculates and returns a vector data structure that 
*     contains the point or derivative value of the given non-uniform integral
*     B-spline  curve at the specific parameter value.
* 
* 4   References
*     [1]C. De Boor. A Practical Guide to Splines, Springer, New York, 1978.
* 
* 5   Parameters
*       1.ParCurv * egeom
*         On entry: NURBS data structure containing geometry of curve to be 
* 	          evaluated.
*       2.double u
*         On entry: the parameter value at which the non-uniform integral 
* 	          B-spline curve is to be evaluated.
*       3.int deriv
*         On entry: an index referring to the required derivative to be 
* 	          calculated such that if the curve is defined as R(u) then 
* 		  the routine will evaluate
*                                    deriv
*                                   @     R(u)
*                                  ------------
*                                      deriv
*                                    @u
* 		  ( @ represents partial derivative. )
* 
* 6   Return Values, Error Indicators and Warnings
*     The index deriv <= the order of the non-uniform integral B-spline curve.
*     The value of u should lie between (egeom -> knots[0]) and (egeom ->
*     knots[k]) where k = (egeom -> ncontpts) + (egeom -> order) - 1.  If not,
*     it will  print an error message and stop.
* 
* 8   Further Comments
*     This basis function evaluation is performed using routine nbasisd.
*     This routine should not be used to evaluate derivatives of NURBS curves. 
*     For the evaluation of NURBS curves, see rbspeval (bspl _ ) for details.
* 
* 9   Functions referenced by evalderivbsp() are:
*     Aij()
*     dbl_array2()
*     dbl_array3()
*     find()
*     free_darray2()
*     free_darray3()
*     nbasisd()
*     vectalloc()
* 
* 10  Functions that reference evalderivbsp() are:
*     funccurve()
*     knot_holzle()
*     next_u()
*     point_bspl()
*     rbspeval()

******************************************************************************/

vector *evalderivbsp(ParCurv *egeom, double u, int deriv)
{
  int ncontpts,order,i,il;    /* number of control points, order, indices */
  vector *p,*eval;            /* a control point and the evaluated point */
  double q;                   /* B-spline basis */
  int prod;                   /* the factor of the derivative curve */
  double **N;                 /* matrix containing B-spline basis */
  double ***derivpts;         /* control points of derivative curves 
                               * without factors */

  /* allocate memory */
  N = dbl_array2((unsigned)(egeom->order+1),(unsigned)(egeom->order+1));
  derivpts = dbl_array3((unsigned)egeom->order,(unsigned)egeom->order,4);

  prod = 1.0;
  ncontpts = egeom->ncontpts;
  order = egeom->order;

  /* if u is out of the range of the knot vector. */
  if( u > egeom->knots[ncontpts+order-1] || u < egeom->knots[0] )
	{
	printf("u = %f outside the knot range in fn. evalderivbsp\n", u);
	exit(0);
	} 

  /* calculate the factor (order-1)(order-2)...(order-deriv) of the 
     derivative curve */
  for(i=1; i<=deriv; i++)
	prod *= order-i;

  /* find index il such that knots[il] <= u <= knots[il+1] */
  il = find(ncontpts,egeom->knots,u);

  /* evaluate the matrix of basis */	
  nbasisd(il, order, u, egeom->knots, N);

  /* calculate the control points of the derivative curves */ 
  Aij(il, deriv, order, egeom->knots, egeom->contpts, derivpts);

  eval = vectalloc(); 
  eval->x = eval->y = eval->z = eval->w = 0.0;

  /* evaluate the derivative curve without the factor */
  for(i=il; i>il-order+deriv; i--)
	{
	q = N[i+order-deriv-il][order-deriv];
	eval->x  += q*derivpts[i-il+(order-1)][deriv][0];
	eval->y  += q*derivpts[i-il+(order-1)][deriv][1];
	eval->z  += q*derivpts[i-il+(order-1)][deriv][2];
	eval->w  += q*derivpts[i-il+(order-1)][deriv][3];
	}

  /* multiply by the factor */
  eval->x = prod*eval->x;
  eval->y = prod*eval->y;
  eval->z = prod*eval->z;
  eval->w = prod*eval->w;

  /* free memories */
  free_darray2(N);
  free_darray3(derivpts);

  return(eval);
}

/*****************************************************************************
*                                  rbspeval()
******************************************************************************
* 
* 1   Purpose
*     This function evaluates a NURBS curve and computes derivatives.
* 
* 2   Specification
*     #include "bspl.h"
*     vector *rbspeval(ParCurv *egeom, double u, int deriv)
* 
* 3   Description
*     This routine calculates and returns a vector data structure that 
*     contains the point or derivative of a point value of the given NURBS
*     curve at the  specific parameter value. The returned vector data
*     structure contains  values for all four homogeneous coordinates.
* 
* 4   References
*     [1]C. De Boor. A Practical Guide to Splines, Springer, New York, 1978.
* 
* 5   Parameters
*       1.ParCurv * egeom
*         On entry: NURBS data structure containing geometry of curve to be 
* 	          evaluated.
*       2.double u
*         On entry: the parameter value at which the non-uniform integral 
* 	          B-spline curve is to be evaluated.
*       3.int deriv
*         On entry: index denoting required derivative where if the curve is 
* 	          defined as R(u) then the routine will return
*                                        deriv
*                                       d     R(u)
*                                      ------------
*                                           deriv
*                                         du
* 
* 6   Return Values, Error Indicators and Warnings
*     The value of u should lie between (egeom -> knots[0]) and (egeom ->
*     knots[k])  where k = (egeom -> ncontpts)+(egeom -> order)-1. If not, it
*     will print an  error message and stop.
*     The value for deriv <= 4.
* 
* 7   Accuracy
*     Not applicable.
* 
* 8   Further Comments
*     This basis function evaluation is performed using routine nbasisd 
*     (bspl _ cxdeb.c).
* 
* 9   Functions referenced by rbspeval() are:
*     evalbsp()
*     evalderivbsp()
* 
* 10  Functions that reference rbspeval() are:
*     bone_intersect()
*     bone_intersect3()
*     calc_disc()
*     calc_disc_curv()
*     camber_line_2d()
*     camber_line_3d()
*     camber_surface()
*     compute_diffeq()
*     compute_ortho_proj()
*     cos_conv_nurbs()
*     cross_link()
*     CurveOrientation()
*     curv_dev()
*     cyl_frenet_tr()
*     cyl_generatrix()
*     cyl_sample_gencyl()
*     deriv_dev()
*     DrawCurvFunction()
*     DrawCurvTorsion()
*     EvalMaxCamber()
*     EvalMaxCamber3d()
*     EvaluateCurv()
*     eval_curve_bounded()
*     eval_iso()
*     eval_ParCurv()
*     fcn_crit()
*     fcn_crit3()
*     fdist_curv()
*     FilterNegative()
*     find_error_curv()
*     find_error_curv1()
*     find_estim_curv()
*     find_points()
*     funccurve()
*     fun_camber()
*     fun_camber3()
*     GetCurvTang()
*     GetCurvValues()
*     g_camber()
*     g_camber3()
*     init_ortho_proj()
*     init_proj()
*     int_sample_gencyl()
*     lead_chord()
*     loft_integral()
*     max_chord()
*     max_cmbr()
*     max_cmbr3()
*     max_curve_curvature()
*     normalcurv()
*     offset_curv2d()
*     offset_geod()
*     offset_geod_par()
*     offset_normal()
*     ParCurv_to_PowCurv()
*     position_err()
*     PostScriptCurvFunction()
*     PostScriptCurvTorsion()
*     PowParCurv_compare()
*     rbspcur()
*     sample_blend()
*     sample_gencyl()
*     sample_mod_surf()
*     surf_data()
*     trim_start_par()
*     trim_start_par3()
* 
******************************************************************************/
vector *rbspeval(ParCurv *egeom, double u, int deriv)
{
  vector *eval0,*eval1,*eval2,*eval3,*eval4; /* derivatives */
  double w,wp,wpp,wppp,wpppp;                /* scaling factors */

  /* if the order of derivative is 0, just evaluate the curve at u by
   * calling function evalbsp(). */
  eval0 = evalbsp(egeom,u);
  w = eval0->w;

  /* if 0th derivative, return eval0. */
  if(deriv == 0)
    return(eval0);
  
  /* if the order of derivative is no less than 1 */
  if(deriv >= 1)
    {
      /* evaluate the 1st derivative by calling evalderivbsp() */
      eval1 = evalderivbsp(egeom,u,1);
      wp = eval1->w;
      /* transform to homogeneous coordinates */
      eval1->x = (eval1->x - (eval0->x/w)*wp)/w;
      eval1->y = (eval1->y - (eval0->y/w)*wp)/w;
      eval1->z = (eval1->z - (eval0->z/w)*wp)/w;
      eval1->w = 1.0;
      /* if only the 1st derivative is wanted, free eval0 and return eval1 */
      if(deriv == 1)
	{
	  free((char *) eval0);
	  return(eval1);
	}
    }

  /* if the order of derivative is no less than 2 */
  if(deriv >= 2)
    {
      /* evaluate the 2nd derivative by calling evalderivbsp() */
      eval2 = evalderivbsp(egeom,u,2);
      wpp = eval2->w;
      /* transform to homogeneous coordinates */ 
      eval2->x = (eval2->x - 2.0*eval1->x*wp - (wpp*eval0->x/w))/w;
      eval2->y = (eval2->y - 2.0*eval1->y*wp - (wpp*eval0->y/w))/w;
      eval2->z = (eval2->z - 2.0*eval1->z*wp - (wpp*eval0->z/w))/w;
      eval2->w = 1.0;
      /* if only the 2nd derivative is wanted, free eval0, eval1 and return
       * eval2 */ 
      if(deriv == 2)
	{
	  free((char *) eval0);
	  free((char *) eval1);
	  return(eval2);
	}
    }

  /* if the order of derivative is no less than 3 */
  if(deriv >= 3)
    {
      /* evaluate the 3rd derivative by calling evalderivbsp() */
      eval3 = evalderivbsp(egeom,u,3);
      wppp = eval3->w;
      /* transform to homogeneous coordinates */
      eval3->x = (eval3->x - 3.0*eval2->x*wp - 3.0*eval1->x*wpp
		  - wppp*eval0->x/w)/w;
      eval3->y = (eval3->y - 3.0*eval2->y*wp - 3.0*eval1->y*wpp
		  - wppp*eval0->y/w)/w;
      eval3->z = (eval3->z - 3.0*eval2->z*wp - 3.0*eval1->z*wpp
		  - wppp*eval0->z/w)/w;
      eval3->w = 1.0;
      /* if only the 3rd derivative is wanted, return eval3 and free others */
      if(deriv == 3)
	{
	  free((char *) eval0);
	  free((char *) eval1);
	  free((char *) eval2);
	  return(eval3);
	}
    }

  /* if the order of derivative is no less than 4 */
  if(deriv >= 4)
    {
      /* evaluate the 4th derivative by calling evalderivbsp() */
      eval4 = evalderivbsp(egeom,u,4);
      wpppp = eval4->w;
      /* transform to homogeneous coordinates */
      eval4->x = (eval4->x - 4.0*eval3->x*wp - 6.0*wpp*eval2->x
		  - 4.0*wppp*eval1->x - (wpppp*eval0->x/w))/w;
      eval4->y = (eval4->y - 4.0*eval3->y*wp - 6.0*wpp*eval2->y
		  - 4.0*wppp*eval1->y - (wpppp*eval0->y/w))/w;
      eval4->z = (eval4->z - 4.0*eval3->z*wp - 6.0*wpp*eval2->z
		  - 4.0*wppp*eval1->z - (wpppp*eval0->z/w))/w;
      eval4->w = 1.0;
      /* if only the 4th derivative is wanted, return eval4 and free others */
      if(deriv == 4)
	{
	  free((char *) eval0);
	  free((char *) eval1);
	  free((char *) eval2);
	  free((char *) eval3);
	  return(eval4);
        }
    }
}

/******************************************************************************
*                                 Aij()
*******************************************************************************
* 
* 1   Purpose
*     This is a subroutine of function evalderivbsp(), which evaluates the 
*     derivatives of a B-spline curve, see evalderivbsp(). 
* 
* 2   Specification
*     #include "bspl.h"
*     void Aij(int p, int q, int k, double t[], vector *contpts[], 
*              double ***derivpts)
* 
* 3   Description
*     This function calculates the control points of the q-th derivative curve
*     of a B-spline curve. The derivative curve itself is also a B-spline curve
*     with lower order (order of the curve - q ). The control points of the
*     derivative curve are recursively determined from the control points of
*     the curve.
* 
* 4   References
*     [1]C. De Boor. A Practical Guide to Splines, Springer, New York, 1978.
* 
* 5   Parameters
*       1.int p
*         On entry: the index of the knot such that knots[p]<=u<=knots[p+1],
* 	          where knots[] is the knot vector, u is the parametric value
* 		  to be evaluated.
*       2.int q
*         On entry: the order of the derivative.
*       3.int k
*         On entry: the order of the curve to be evaluated.
*       4.double t[]
*         On entry: the knot vector of the curve to be evaluated.
*       5.vector *contpts[]
*         On entry: the control points of the curve to be evaluated.
*       6.double ***derivpts
*         On entry: a 3D array
*         On exit:  the 3D array containing the control points of the
*                   derivative curves with derivative order varying from 0 to
*                   deriv.
* 
* 6   Return Values, Error Indicators and Warnings
*     Only derivpts[i][deriv][], deriv <= i <= n, are used by function
*     evalderivbsp(). Here n is # of control points of the original curve.
* 
* 7   Accuracy
*     Not applicable.
* 
* 8   Further comments
*     Not appicable
* 
* 9   Functions referenced by Aij() are:
*     errormsg()
* 
* 10  Functions that reference Aij() are:
*     evalderivbsp()
* 
*****************************************************************************/

void Aij(int p, int q, int k, double t[], vector *contpts[],
	 double ***derivpts)
/************************************************************************ 
 * The Aij are such that A sub(i) sup(j) = A[i-il+(k-1)][j]
 * Although the variable il is not available here, this value is useful
 * in indexing A[][] compactly and to calculate only Aij required for a 
 * a given il, see program evalderivbsp() where Aij is used as such. 
 **********************************************************************/
{
  vector *v1, *v2, *v3;   /* v1,v3 -- never used (glshen) 
			   * v2 -- control point */
  double dt;              /* the factor in hodograph of a B-spline curve */
  int i,j,ireal;          /* ireal -- index shift */

  /* set the control points of 0-th derivative curve 
     as the same as the originals */
  for(i=0; i<k; i++)
	{
	v2 = contpts[i+p-(k-1)];
	derivpts[i][0][0] = v2->x;
	derivpts[i][0][1] = v2->y;
	derivpts[i][0][2] = v2->z;
	derivpts[i][0][3] = v2->w;
	}

  /* calculate the control points of j-th derivative curves , 1 <= j <= q */
  for(j=1; j<=q; j++)
     for(i=j; i<k; i++) {
	ireal = i+p-(k-1);      /* index shift */
        dt = t[ireal+k-j] - t[ireal]; 
        if(dt < ZERO1)
	  errormsg(1,"dt is zero in routine Aij()");	
        derivpts[i][j][0] = (derivpts[i][j-1][0] - derivpts[i-1][j-1][0])/dt;
        derivpts[i][j][1] = (derivpts[i][j-1][1] - derivpts[i-1][j-1][1])/dt;
        derivpts[i][j][2] = (derivpts[i][j-1][2] - derivpts[i-1][j-1][2])/dt;
        derivpts[i][j][3] = (derivpts[i][j-1][3] - derivpts[i-1][j-1][3])/dt;
     }
}

/*****************************************************************************
*                               normalcurv()
******************************************************************************
* 
* 1   Purpose
*     This function returns the computed unit normal vector of a NURBS curve.
* 
* 2   Specification
*     #include "bspl.h"
*     vector *normalcurv(ParCurv *egeom, double t, int is_open)
* 
* 3   Description
*     This routine computes the normal to a curve at a given parametric point. 
*     The computed normal is returned as a vector whose w = 1.
* 
* 4   References
*     [1] P. G. Alourdas. Shape Creation, Interrogation and Fairing Using 
*         B-Splines, Engineer's Thesis, Massachusetts Institute of Technology, 
* 	Department of Ocean Engineering, Cambridge, Massachusetts, 1989.
*     [2] I. D. Faux and M. J. Pratt.  Computational Geometry for Design and 
*         Manufacture, Ellis Horwood: Chichester, England, 1979.
* 
* 5   Parameters
*       1.ParCurv * egeom
*         On entry: NURBS curve data structure containing geometry of curve
* 	          that is to be evaluated.
*       2.double t
*         On entry: parameter value at which to compute the curve normal.
*       3.int is_open
*         On entry: flag specifying whether the curve is periodic (is_open=0) 
* 	          or non-periodic (is_open = 1).
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
* 9   Functions referenced by normalcurv() are:
*     cross()
*     dot()
*     evalbsp()
*     evalbsp_per()
*     mag()
*     rbspeval()
*     rbspeval_per()
*     scale_vect1()
*     sub_vect()
*     unitvector1()
*     vectfree()
* 
* 10  Functions that reference normalcurv() are:
*     localize_sumsq_2d()
*     MinDistCurvCB()
* 
******************************************************************************/

vector *normalcurv(ParCurv *egeom, double t, int is_open)
{
  vector *c, *norm, *r, *r1, *r2, *u;
         /* c -- never used (glshen)
	    norm -- the normal
	    r -- vector value at t
            r1,r2 - derivatives at t
            u -- unit vector */
  double a, b, k;     /* a,b -- scalor
			 k -- never used (glshen) */

  if (is_open)     /* if the curve is non-periodic */
     {
      r = evalbsp(egeom, t);                /* value at t */
      r1 = rbspeval(egeom, t, 1);           /* 1st derivative at t */
      r2 = rbspeval(egeom, t, 2);           /* 2nd derivative at t */
     }
  else             /* if the curve is periodic */
     {
      r = evalbsp_per(egeom, t);            /* value at t */
      r1 = rbspeval_per(egeom, t, 1);       /* 1st derivative at t */
      r2 = rbspeval_per(egeom, t, 2);       /* 2nd derivative at t */
     }

/** these lines commented out 4/10/96 because they seem to have no purpose **/
/*  a = mag(r1);       */                     /* norm of 1st derivative */
/*  u = cross(r1, r2); */                     /* cross product of r1 and r2 */
/*  b = mag(u);        */                     /* norm of the cross product */
/** end of change **/

  /* calculate the normal */ 
  a = dot(r2, r1);                          /* dot product of r1 and r2 */
         
  b = (mag(r1)*mag(r1));                    /* square of the norm of r1 */
/** replace next line since u is no longer allocated, 4/10/96 **/
/** scale_vect1(a/b, r1, u); **/            /* scale r1 with factor a/b */
  u = scale_vect(a/b, r1);                  /* scale r1 with factor a/b */
  norm = sub_vect(r2, u);                    

  vectfree(r2);
  vectfree(u);

  unitvector1(norm, norm);                  /* normalize to unit vector */

  return norm;
}

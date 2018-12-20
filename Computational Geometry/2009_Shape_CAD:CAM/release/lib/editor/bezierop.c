/* Copyright (C) 1993 Massachusetts Institute of Technology, Cambridge, MA
   All rights reserved
*/

/* bezierop.c */

/* bezier_clip_decision_egeom()
   bezier_knots()
   calc_interval_egeom()
   checknot()
   clip_bezier_curve()
   eval_curve_bounded()
   normalize_egeom()
   prodknot()
*/

#include <math.h>
#include <malloc.h>
#include "gen.h"
#include "bspl.h"
#include "editor.h"

/*---------------------------------------------------------------------------
   Calculates the knot vector of a Bezier curve/surface of the order
   specified in argument list.
  -------------------------------------------------------------------------*/

void bezier_knots(int order, double first, double last, double *t)
{
  short i;
  
  for (i=0; i<order; i++)         
    t[i] = first;

  for (i=order; i<order+order; i++)
    t[i] = last;  
}

/************************************************************************
  NOTICE that egeom(u,v) is a Bezier function, egeom(u,v)=[u,v,egeom(u,v)]
  dir (u,or v) is the clipping direction.
  TOL is a tolerance for terminating clipping.
  MRVAL is the minimum reduction value for each clipping; 

     0 : If either fgeom1, or fgeom2 are never 0; 
     1 : If the required tolerance is reached;
     2 : Do subdivision;
     3 : Bezier clip;
 ************************************************************************/

int bezier_clip_decision_egeom(ParCurv *egeom, int dir, double *Tmin,
			       double *Tmax)
{
  int id, m=0;
  double2 tmin, tmax, dt, Dt, dtold;

     /*---------------------------------------------*/
     /* tmin < parameter range of egeom < tmax      */
     /*---------------------------------------------*/

  tmin = egeom->knots[0]; 
  tmax = egeom->knots[egeom->kmem-1];

  dt = tmax - tmin; 

  dtold = dt;

    /*------------------------------------------------------------*/
    /*       Calculate [umin,umax], [vmin,vmax] intervals         */
    /*------------------------------------------------------------*/
    
  if (dir == EDITOR_U)
    id = calc_interval_egeom(egeom, EDITOR_U, &tmin, &tmax); 
  else
    id = calc_interval_egeom(egeom, EDITOR_V, &tmin, &tmax); 

    /*------------------------------------------------------------*/
    /*  No intersection is found                                  */
    /*------------------------------------------------------------*/

  if (id == 0) {
    *Tmin = 0.0;
    *Tmax = 0.0;
    return(EDITOR_NO_ROOT);
  }

    /*------------------------------------------------------------*/
    /*  Check if tmax - tmin < TOL if TRUE return WITHIN_TOL      */
    /*------------------------------------------------------------*/

  if ((tmax - tmin)<EDITOR_TOL) {
    *Tmin = tmin;
    *Tmax = tmax;
    return(EDITOR_WITHIN_TOL);
  }

    /*------------------------------------------------------------*/
    /*             Calculate interval reduced value               */
    /*------------------------------------------------------------*/

  dt = (tmax - tmin);      
  Dt = dt/dtold;
  *Tmin = tmin; 
  *Tmax = tmax; 
  
  if (Dt < EDITOR_MRVAL) 
    return(EDITOR_BEZIER_CLIP);
  else
    return(EDITOR_SUBDIVISION);
}

/*--------------------------------------------------------------------------
   Calculates interval [Tmin,Tmax], where f(u)=0 or f(v)=0.
   NOTICE that f(u) = [u, f(u)] or f(v) = [v, f(v)] .
   0 : convex hull does not intersect.
   1 : convex hull intersects axis

  --------------------------------------------------------------------------*/

int calc_interval_egeom(ParCurv *egeom, int dir, double *Tmin, double *Tmax)
{
  int i, j, k, nk, low_cvx_hull, upp_cvx_hull, order, ncontpts;
  int cvx_intrsct;
  double knot_s, knot_e;
  double *t, *d, *t_cvxl, *t_cvxu;
  double  *d_cvxl, *d_cvxu, *t_param;

  order = egeom->order;
  ncontpts = egeom->ncontpts;
  knot_s = egeom->knots[0];
  knot_e = egeom->knots[egeom->kmem-1];

    /*------------------------------------------------------------*/
    /*        Allocation                                          */
    /*------------------------------------------------------------*/
     
  t = dbl_array1 (order);
  d = dbl_array1 (order);
  t_cvxl  = dbl_array1 (order);
  t_cvxu  = dbl_array1 (order);
  d_cvxl  = dbl_array1 (order);
  d_cvxu  = dbl_array1 (order);
  t_param = dbl_array1 (order);
     
     /*------------------------------------------------------------------*/
     /* Calculate  distance from axis at each i, i=0,...,contpts-1 */
     /*------------------------------------------------------------------*/

  for (i=0; i<ncontpts; i++) {
    d[i] = egeom->contpts[i]->z;

    if(dir == EDITOR_V)
      t[i] = egeom->contpts[i]->x;
    else
      t[i] = egeom->contpts[i]->y;
  }
     
     
     /*------------------------------------------------------------------*/
     /* Calculate convex hull                                            */
     /*------------------------------------------------------------------*/
     
  low_cvx_hull = calc_convex_hull(EDITOR_LOW, ncontpts, t, d, t_cvxl, d_cvxl);
  upp_cvx_hull = calc_convex_hull(EDITOR_UPP, ncontpts, t, d, t_cvxu, d_cvxu);
     
     /*------------------
       draw convex hull
       -----------------*/

     /*-----------------------------
       convex hull axis intersection 
       -----------------------------*/

  cvx_intrsct = cvx_hull_axis_intrsct(low_cvx_hull, t_cvxl, d_cvxl,
				      upp_cvx_hull, t_cvxu, d_cvxu,
				      t_param, knot_s, knot_e, Tmin, Tmax);
  
     /*------------------------------------------------------------------*/
     /* Free memory                                                      */
     /*------------------------------------------------------------------*/

  free_darray1(t);
  free_darray1(d);
  free_darray1(t_cvxl);
  free_darray1(t_cvxu);
  free_darray1(d_cvxl);
  free_darray1(d_cvxu);
  free_darray1(t_param);
  
  return(cvx_intrsct);
}

/*--------------------------------------------------------------------------*/
  
ParCurv *clip_bezier_curve(ParCurv *egeom, double min, double max)
{
  ParCurv *newegeom, *kid1, *kid2;
  int knot_s_comp, knot_e_comp;
  
  knot_s_comp = lf2comp(min, egeom->knots[0], EDITOR_ZEROD);
  knot_e_comp = lf2comp(egeom->knots[egeom->kmem-1], max,  EDITOR_ZEROD);
     
  if (knot_s_comp == 1) {
    decasteljau_curve(egeom, min, &kid1, &kid2);
    free_egeom(kid1);
  }
  else if (knot_s_comp == 0) {
    kid2 = egeomalloc1(egeom->order, egeom->ncontpts);
    copyegeom(egeom, kid2);
  }

  if (knot_e_comp == 1) {
    decasteljau_curve(kid2, max, &newegeom, &kid1);
    free_egeom(kid1);
    free_egeom(kid2);
  }
  else if (knot_e_comp == 0) {
    newegeom = egeomalloc1(egeom->order, egeom->ncontpts);
    copyegeom(kid2, newegeom);
    free_egeom(kid2);
  }

  return(newegeom);
}

/************************************************************************
   This routine will calculate the number of new knots.
 ************************************************************************/

int checknot(ParCurv *egeom)
{
  short i = 0, j, k = 0, kk = 1, n, order;

  order = egeom->order;
  n = egeom->ncontpts;

  k = 0;
  kk =1;
  i = 0;
  while (1) {
    for (j=0; j<order; j++) {
      if (i+j+1 > order+n-1) {
	k = k+order;
	return k;
      } 
      if (fabs(egeom->knots[i+j+1]- egeom->knots[i+j]) < EDITOR_TOLKNOT)
	kk++;
      else
	break;
    }
    k = k+order;
    i=i+kk;
    kk=1;
  }
}

/************************************************************************
   This routine will produce new knot vector kn from original B-spline
   curve egeom
 ************************************************************************/

void prodknot(ParCurv *egeom, double *kn)
{
  int i = 0, ii, j, k = 0, kk = 1, n, order;

  order = egeom->order;
  n = egeom->ncontpts;

  k=0;
  kk=1;
  i=0;
  while (1) {
    for (j=0; j<order; j++) {
      if (i+j+1> order+n-1)  {
	for (ii=k; ii<k+order; ii++) {
	  kn[ii] = egeom->knots[i];
	}
	k = k+order;
	return;
      }
      if (fabs(egeom->knots[i+j+1]- egeom->knots[i+j]) < EDITOR_TOLKNOT)
	kk++;
      else
	break;
    }
    for (ii=k; ii<k+order; ii++)
      kn[ii] = egeom->knots[i];

    k = k+order;
    i = i+kk;
    kk = 1;
  }
}

void normalize_egeom(ParCurv *egeom)
{
  double max_val;
  short i;
     
  max_val = 1.0e-10;
  for (i=0; i<egeom->ncontpts; i++)
    if(fabs(egeom->contpts[i]->z) > max_val)
      max_val = fabs(egeom->contpts[i]->z);
     
  if (max_val > 1.0e-10)
    for (i=0; i<egeom->ncontpts; i++)
      egeom->contpts[i]->z /= max_val;
}

vector *eval_curve_bounded(ParCurv *egeom, double t, int index)
{
  vector *vect, *vect1;
  double help, eps;

  if (t<egeom->knots[0] || t>egeom->knots[egeom->kmem-1]) {
    if (t<egeom->knots[0]) {
      eps = egeom->knots[0] - t;
      help = egeom->knots[0];
    }
    else if (t>egeom->knots[egeom->kmem-1]) {
      eps = egeom->knots[egeom->kmem-1] - t;
      help = egeom->knots[egeom->kmem-1];
    }
    else{
      eps = 0.0;
      help = t;
    }

    vect = rbspeval(egeom, help, index);
    vect1 = rbspeval(egeom, help+eps, index);

    scale_vect1(2.0, vect, vect);
    sub_vect1(vect, vect1, vect);
    vectfree(vect1);
  }
  else
    vect = rbspeval(egeom, t, index);

  return (vect);
}

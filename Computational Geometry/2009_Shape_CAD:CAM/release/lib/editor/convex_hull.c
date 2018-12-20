/* Copyright (C) 1993 Massachusetts Institute of Technology, Cambridge, MA
   All rights reserved
*/

/* convex_hull.c */

/* calc_convex_hull()
   calc_intsct()
   common_interval()
   cvx_hull_axis_intrsct()
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gen.h"
#include "bspl.h"
#include "editor.h"

int cvx_hull_axis_intrsct(int low_cvx_hull, double *t_cvxl, double *d_cvxl,
			  int upp_cvx_hull, double *t_cvxu, double *d_cvxu,
			  double *t_param, double knot_s, double knot_e,
			  double *Tmin, double *Tmax)
{
  int low_cvx_comp, upp_cvx_comp;
  int dcvxl_s_comp, dcvxl_e_comp, dcvxu_s_comp, dcvxu_e_comp;
  int ilow, low_cvx_zero, low_cvx_pos, low_cvx_neg;
  int iupp, upp_cvx_zero, upp_cvx_pos, upp_cvx_neg;
  int intsct, low_intsct, return_value;
  short i;

/*------------------------------------------------------------------
  count numbers of positive convex hull vertex
  count numbers of negative convex hull vertex
  count numbers of zero     convex hull vertex
  -----------------------------------------------------------------*/

/*---------------
  Initialization
  ---------------*/

  ilow = 0;
  iupp = 0;
  low_cvx_zero = 0;
  low_cvx_pos  = 0;
  low_cvx_neg  = 0;
  upp_cvx_zero = 0;
  upp_cvx_pos  = 0;
  upp_cvx_neg  = 0;
  dcvxl_s_comp  = lf2comp(d_cvxl[0], 0., EDITOR_ZERO_CVX);
  dcvxl_e_comp  = lf2comp(d_cvxl[low_cvx_hull-1], 0., EDITOR_ZERO_CVX);
  dcvxu_s_comp  = lf2comp(d_cvxu[0], 0., EDITOR_ZERO_CVX);
  dcvxu_e_comp  = lf2comp(d_cvxu[upp_cvx_hull-1], 0., EDITOR_ZERO_CVX);
  
     /*----------------
       low convex hull
       ---------------*/
  for(i = 0; i < low_cvx_hull; i++) {
    low_cvx_comp = lf2comp(d_cvxl[i], 0., EDITOR_ZERO_CVX);

    if (low_cvx_comp == 0) {
      low_cvx_zero ++;
      ilow ++;
    }
    else if(low_cvx_comp == 1)
      low_cvx_pos ++;
    else if(low_cvx_comp == -1)
      low_cvx_neg ++;
  }
     
     /*----------------
       upp convex hull
       ---------------*/

  for(i = 0; i < upp_cvx_hull; i++){
    upp_cvx_comp = lf2comp(d_cvxu[i], 0., EDITOR_ZERO_CVX);
    
    if(upp_cvx_comp == 0){
      upp_cvx_zero ++;
      iupp ++;
    }
    else if(upp_cvx_comp == 1)
      upp_cvx_pos ++;
    else if(upp_cvx_comp == -1)
      upp_cvx_neg ++;
  }
     
     /*--------------------------------------------------------------------
       convex hull axis intersection
       -------------------------------------------------------------------*/
     
     /* check if convex are flat and zero */
  
  if(low_cvx_zero == low_cvx_hull && upp_cvx_zero == low_cvx_hull){
    *Tmin = 0.0;
    *Tmax = 0.0;
    return_value = EDITOR_NO_ROOT;
  }           /* fg is positive */
  else if(dcvxl_s_comp == 1 && low_cvx_neg == 0 && dcvxl_e_comp != 0 &&
	  dcvxu_s_comp != 0 && dcvxu_e_comp != 0){
    *Tmin = 0.0;
    *Tmax = 0.0;
    return_value = EDITOR_NO_ROOT;
  }           /* fg is negative */
  else if(dcvxu_s_comp == -1 && upp_cvx_pos == 0 && dcvxl_s_comp != 0 &&
	  dcvxl_e_comp != 0 && dcvxu_e_comp != 0){
    *Tmin = 0.0;
    *Tmax = 0.0;
    return_value = EDITOR_NO_ROOT;
  }    /*check if convex hull touches at left */
  else if (dcvxl_s_comp == 0 && dcvxu_s_comp == 0 &&
	   ((low_cvx_pos == low_cvx_hull-1 &&
	     upp_cvx_pos ==
	     upp_cvx_hull-1)  ||
	    (low_cvx_neg == low_cvx_hull-1 &&
	     upp_cvx_neg == upp_cvx_hull-1))){
    *Tmin = t_cvxl[0];
    *Tmax = t_cvxl[0];
    return_value = EDITOR_CONVEX_HULL_INTERSECT;
  }     /*check if convex hull touches at right */
  else if (dcvxl_e_comp == 0 && dcvxl_e_comp == 0 &&
	   ((low_cvx_pos == low_cvx_hull-1 &&
	     upp_cvx_pos ==
	     upp_cvx_hull-1)  ||
	    (low_cvx_neg == low_cvx_hull-1 &&
	     upp_cvx_neg == upp_cvx_hull-1))){
    *Tmin = t_cvxl[low_cvx_hull-1];
    *Tmax = t_cvxl[low_cvx_hull-1];
    return_value = EDITOR_CONVEX_HULL_INTERSECT;
  } /*Calculate intersections between convex hull and axis*/
  else if((low_cvx_pos > 0 && low_cvx_neg > 0) ||
	  (upp_cvx_pos > 0 && upp_cvx_neg > 0)){
    
    low_intsct = 0;
    intsct = calc_intsct(low_cvx_hull, t_cvxl, d_cvxl, upp_cvx_hull,
			 t_cvxu, d_cvxu, t_param, &low_intsct);
	  
    if(intsct == 1){
      if(low_intsct == 1){
	if(dcvxl_s_comp <= 0){ 
	  *Tmin = t_cvxl[0];
	  *Tmax = t_param[0];
	}
	else{
	  *Tmin = t_param[0];
	  *Tmax = t_cvxl[low_cvx_hull-1];
	}
      }
      else{
	if(dcvxu_s_comp >= 0){
	  *Tmin = t_cvxl[0];
	  *Tmax = t_param[0];
	}
	else{
	  *Tmin = t_param[0];
	  *Tmax = t_cvxl[low_cvx_hull-1];
	}
      }
    }
    else {
      *Tmin = 1.0e+10;
      for(i=0;
	  i<intsct; i++){
	if (t_param[i] < (*Tmin))
	  *Tmin = t_param[i];
      }
      
      *Tmax = -1.0e+10;
      for (i=0;
	   i<intsct; i++) {
	if (t_param[i] > (*Tmax))
	  *Tmax = t_param[i];
      }
    }
    return_value = EDITOR_CONVEX_HULL_INTERSECT;
  }
  else {
    *Tmin = knot_s;
    *Tmax = knot_e;

    return_value = EDITOR_CONVEX_HULL_INTERSECT;
  }

  return(return_value);
}

/*---------------------------------------------------------------------------
  Calculate intersections between convex hull and axis
  -------------------------------------------------------------------------*/

int calc_intsct(int low_cvx_hull, double *t_cvxl, double *d_cvxl,
		int upp_cvx_hull, double *t_cvxu, double *d_cvxu,
		double *t_param, int *low_intsct)
{
  int i, intsct, low_cvx_prev, low_cvx_comp, upp_cvx_prev, upp_cvx_comp;
  double t1, d1, t2, d2, slope;

  intsct= 0;

     /* lower convex hull */

  for (i=1; i<low_cvx_hull; i++) {
    t1 = t_cvxl[i-1];
    d1 = d_cvxl[i-1];
    t2 = t_cvxl[i];
    d2 = d_cvxl[i];
    
    low_cvx_comp = lf2comp(d2, 0., EDITOR_ZERO_CVX);
    low_cvx_prev = lf2comp(d1, 0., EDITOR_ZERO_CVX);
	  
    if (low_cvx_comp * low_cvx_prev <0) {
      slope = (d2-d1)/(t2-t1);
      t_param[intsct] = -(d1 - slope*t1)/slope;
      (*low_intsct)++;
      intsct++;
    }
    else if (low_cvx_comp == 0 && i != low_cvx_hull-1) {
      t_param[intsct] = t2;
      (*low_intsct)++;
      intsct++;
    }
  }
     
     /* upper convex hull */
     
  for (i=1; i<upp_cvx_hull; i++) {
    t1 = t_cvxu[i-1];
    d1 = d_cvxu[i-1];
    t2 = t_cvxu[i];
    d2 = d_cvxu[i];

    upp_cvx_comp = lf2comp(d2, 0., EDITOR_ZERO_CVX);
    upp_cvx_prev = lf2comp(d1, 0., EDITOR_ZERO_CVX);
    
    if (upp_cvx_comp * upp_cvx_prev <0) {
      slope = (d2-d1)/(t2-t1);
      t_param[intsct] = -(d1 - slope*t1)/slope;
      intsct++;
    }
    else if (upp_cvx_comp == 0 && i != upp_cvx_hull-1) {
      t_param[intsct] = t2;
      intsct++;
    }
  }
     
  return(intsct);
}
  
/*--------------------------------------------------------------------------
   Calculates the lower/upper (id=0,1) limit of the convex hull,
   corresponding to a set of n cooplanar points with x,y coordinates.
   NOTICE that the x-coordinates must be given in monotonic order.
   The routine returns the number of convex-hull points;
   the points are returned in xc,yc.
  -------------------------------------------------------------------------*/

int calc_convex_hull(int id, int n, double *x, double *y, double *xc,
		     double *yc)
{
  int i, j, k, iold;
  double2 slope, slope_min, slope_max;

  iold=0;   
  k=0;  
  yc[k]=y[k];  
  xc[k]=x[k];

  switch(id) {
  case 0:	/* Lower-limit of convex hull */
    while (iold < n-1) {
      slope_min = 1.0e+30;
      for (j=iold+1; j<n; j++) {
	slope = (y[j]-yc[k]) / (x[j]-xc[k]);
	if (slope < slope_min)  {
	  i = j;  
	  slope_min = slope;
	}
      }
      k++;  
      yc[k] = y[i];  
      xc[k] = x[i];  
      iold = i;
    }
    break;
  case 1:	/* Upper-limit of convex hull */
    while (iold < n-1) {
      slope_max = -1.0e+30;
      for (j=iold+1; j<n; j++) {
	slope = (y[j]-yc[k]) / (x[j]-xc[k]);
	if (slope > slope_max) {
	  i=j;
	  slope_max = slope;
	}
      }
      k++;  
      yc[k] = y[i]; 
      xc[k] = x[i]; 
      iold = i;
    }
    break;
  default:	/* Error */
    printf("Argument error in calc_convex_hull()\n");
    exit(1);
  }

  return(k+1);
}

/*--------------------------------------------------------------------------
   Given two intervals [a1,b1] and [a2,b2], it calculates the common part,
   returned in [a,b].
   If there is no common part, the routine returns FALSE (0), otherwise
   it returns TRUE (1).
   ------------------------------------------------------------------------*/

int common_interval(double a1, double b1, double a2, double b2,
		    double *a, double *b)
{
  double tmin, tmax, eps = EDITOR_ZEROD;
  int compare;

/* Check arguments */
  if (a1>b1+eps  ||  a2>b2+eps) {
    printf("Error in common_interval() [%lf %lf] [%lf %lf]\n",a1,b1,a2,b2);
    exit(1); 
  }
  
  tmin = MAX(a1,a2);	  tmax = MIN(b1,b2);

  compare = lf2comp(tmin, tmax, eps);

  if (compare == 1)
    return(0);    /* No common part */
  else {
    *a = tmin;   *b = tmax;
    return(1);
  }
}

/* Copyright (C) 1993 Massachusetts Institute of Technology, Cambridge, MA
   All rights reserved
*/

/* ray.c */

/* This code finds intersection points of B-spline and ray */
/* Based on Chun-Y Hu: ~chun/master/hidline/pclass/rayintersect.c */

/* arr_vec()
   BsplineToBezier()
   close_enough()
   deriv_sign()
   FilterNegative()
   implicitize_ray()
   int_number()
   loop_rayintersect()
   next_u()
   point_bspl()
   signchange_ex()
   solve_int()
   solve_u
   split_bspl()
   split_it()
   vec_arr()
   vec_neq()
*/

#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include "gen.h"
#include "bspl.h"
#include "editor.h"

int loop_rayintersect(ParCurv **egeom, int nCurves, double *ray, double *point)
{
  ParCurv **bezier;
  double rmin;
  short i, j, ni, nb;

  rmin = 2 * (egeom[0]->order - 1) * PRECISION;
  
  for (j=ni=0; j<nCurves; j++) {
    nb = BsplineToBezier(egeom[j], &bezier);
    for (i=0; i<nb; i++)
      ni += int_number(bezier[i], ray , point, rmin);

    for (i=0; i<nb; i++)
      free_egeom(bezier[i]);
    free_garray1((char *)bezier);
  }
  
  return ni;
}

/*-------------------------------------------------------------------*/
/* This function is to find the number of the intersection point of  */
/* a ray and a piece of a b-spline.                                   */
/*-------------------------------------------------------------------*/
int int_number(ParCurv *egeom, double *ray, double *point, double rmin)
{
  ParCurv *kid1, *kid2; /* pointers for the two split pieces */
  double implicit[3];   /* coefficients x,y,c of implicit function of ray */
  double *coef, u;
  short i, n, nsign;

  coef = dbl_array1(egeom->ncontpts);
  
  implicitize_ray(point, ray, implicit);
     
  for (i=0; i<egeom->ncontpts; i++)
    coef[i] = implicit[0]*egeom->contpts[i]->x +
              implicit[1]*egeom->contpts[i]->y +
	      implicit[2];
 
  nsign = signchange_ex(coef, egeom->ncontpts, rmin);
     
  if (nsign == 0) {
    free_darray1(coef);
    return 0;
  }
  else if (nsign == 1 && deriv_sign(coef, egeom->ncontpts, rmin)) {
    u = solve_int(egeom, coef); /* record the u */
    free_darray1(coef);

    if (FilterNegative(u, egeom, point) == 1)
      return 1;
    else
      return 0;
  }

  free_darray1(coef);
  kid1 = egeomalloc1(egeom->order, egeom->ncontpts);
  kid2 = egeomalloc1(egeom->order, egeom->ncontpts);
	  
    /* kids have the same order and ncontpts with fathers' */
	  
  split_bspl(egeom, kid1, kid2);

  n = int_number(kid1, ray, point, rmin) +
      int_number(kid2, ray, point, rmin);
    
  free_egeom(kid1);
  free_egeom(kid2);

  return n;
}

/*-----------------------------------------------------*/
/* This function is to implicitize the function of ray */
/*-----------------------------------------------------*/
void implicitize_ray(double *point, double *ray, double *implicit)
{
  implicit[0]=  ray[1];
  implicit[1]= -ray[0];
  implicit[2]=  ray[0]*point[1] - ray[1]*point[0];
}
						  
/*----------------------------------------------------*/
/* This function is to find the number of signchanges */
/*----------------------------------------------------*/
int signchange_ex(double *coeff, int ncontpts, double rmin)
{
  short i, nsign = 0, sign, sign2;

  if (coeff[0] < -rmin)       /* establish sign of coeff[0] */
    sign = -1;
  else if (coeff[0] > rmin)
    sign = 1;
  else
    sign = 0;

  sign2 = sign;

  for (i = 1; i < ncontpts; i++) {
    if (coeff[i] < -rmin)      
      sign = -1;
    else if (coeff[i] > rmin)
      sign = 1;
    else
      sign = 0;
      
    if (sign != sign2)    /* keep track of changes */
      nsign++;
      sign = sign2;
    }

  return nsign;
}

/*----------------------------------------------------------------------*/
/* this function is to see if there is sign change of the derivative of */
/* a bspl coefficient                                                   */
/*----------------------------------------------------------------------*/
int deriv_sign(double *coef, int ncontpts, double rmin)
{
  double *deriv;
  short i, nsign;

  deriv = dbl_array1(ncontpts - 1 );

  for (i=0; i<ncontpts-1; i++)
    deriv[i] = coef[i+1] - coef[i];

  nsign = signchange_ex(deriv, ncontpts-1, rmin);

  free_darray1(deriv);

  if (nsign == 0)
    return 1;

  return 0;
}

/*------------------------------------------------------------------------*/
/* This function is to generate an imaginary b-spline curve and solve the */
/* unique solution of the polynominal                                     */
/*------------------------------------------------------------------------*/
double solve_int(ParCurv *egeom, double *coeff)
{
  ParCurv *sgeom;
  double u;
  short i;

  sgeom = egeomalloc1(egeom->order, egeom->ncontpts);
  
  for (i=0; i<egeom->kmem; i++ )
    sgeom->knots[i] = egeom->knots[i];

  for (i=0; i<egeom->ncontpts; i++) {
    sgeom->contpts[i]->x = coeff[i];
    sgeom->contpts[i]->y = 0.0;
    sgeom->contpts[i]->z = 0.0;
    sgeom->contpts[i]->w = 1.0;
  }

  u = solve_u(sgeom, egeom->knots[0], egeom->knots[egeom->kmem-1]);
  free_egeom(sgeom);

  return u;                  
}

/*-------------------------------------------------------------------*/
/* This function is to solve the u value with the coefficient of the */
/* b-spline basis function                                           */
/*-------------------------------------------------------------------*/
double solve_u(ParCurv *egeom, double start_pr, double end_pr)
{
  double next;

  if (start_pr > end_pr)
    printf("It dosen't converge \n");

  if (close_enough(point_bspl(egeom, start_pr)))
    return start_pr;

  else if (close_enough(point_bspl(egeom, end_pr)))
    return end_pr;

  else if (close_enough(point_bspl(egeom, next = next_u(egeom, start_pr))))
    return next;

  return solve_u(egeom, next, end_pr);
}

/*------------------------------------------------*/
/* This function is to check if u is close enough */
/*------------------------------------------------*/
int close_enough(double f)
{
  if (fabs(f) <= 1.0e-8)
    return 1;
  else
    return 0;
}

/*---------------------------------------------*/
/* This function is to find next approximate u */
/*---------------------------------------------*/
double next_u(ParCurv *egeom, double u)
{
  vector *r;
  double next, start_knot, end_knot, x0, x1;

  r = evalderivbsp(egeom, u, 0);
  x0 = r->x;
  vectfree(r);

  r = evalderivbsp(egeom, u, 1);
  x1 = r->x;
  vectfree(r);

  next = u - x0/x1;

  if (next < (start_knot = egeom->knots[0]))
    return start_knot;

  else if (next > (end_knot = egeom->knots[egeom->kmem-1]))
    return end_knot;

  return(next);
}

/*-----------------------------------------------------------*/
/* This function is to compute the point of an integral bspl */
/*-----------------------------------------------------------*/
double point_bspl(ParCurv *egeom, double u)
{
  vector *r;
  double pr;

  r = evalderivbsp(egeom, u, 0);
  pr = r->x;
  vectfree(r);
  
  return (pr);
}

/*--------------------------------------*/
/* Filter out values less than point[0] */
/*--------------------------------------*/
int FilterNegative(double u, ParCurv *egeom, double *point)
{
  vector *r;
  double x;

  r = rbspeval(egeom, u, 0);
  x = r->x;
  vectfree(r);

  if (x - point[0] > 0.0)
    return 1;

  return 0;
}

/*---------------------------------------------------*/
/* This function is to split bspline into two pieces */
/*---------------------------------------------------*/
void split_bspl(ParCurv *egeom, ParCurv *kid1, ParCurv *kid2)
{
  ParCurv *sgeom ;
  double ip; /* integer part */
  double inserted_u;
  double **dadarray, **kidarray;
  int i, nknots;

  nknots = 2*egeom->order + egeom->ncontpts;
                                               
  modf((egeom->order + egeom->ncontpts)/2.0, &ip); /* ip is the middle point */
       						   /* of the new knot index  */
  inserted_u = (egeom->knots[0] + egeom->knots[egeom->kmem-1]) / 2.0; 

  sgeom = egeomalloc1(egeom->order, nknots - egeom->order);
  /* child has the same order, but egeom->order more ncontpts */

  for (i=0; i<(int)ip; i++)
    sgeom->knots[i] = egeom->knots[i];

  for (i=(int)ip; i<(int)ip + egeom->order; i++)
    sgeom->knots[i] = (egeom->knots[0] + egeom->knots[egeom->kmem-1]) / 2.0;
  /* that will add new knot in the middle point */

  for (i=(short)ip + egeom->order; i<nknots; i++)
    sgeom->knots[i] = egeom->knots[i - egeom->order];

  dadarray = dbl_array2(egeom->ncontpts, 4);
  kidarray = dbl_array2(egeom->ncontpts + egeom->order, 4);

  vec_arr(egeom->ncontpts, egeom->contpts, dadarray);

  curve_oslo1(egeom->order, egeom->ncontpts, nknots, egeom->knots,
	      sgeom->knots, dadarray, kidarray);

  arr_vec(sgeom->ncontpts, kidarray, sgeom->contpts);

  free_darray2(dadarray);
  free_darray2(kidarray);

  split_it(inserted_u, sgeom, kid1, kid2);

  free_egeom(sgeom);
}

/*------------------------------------------------------*/
/* This function is to convert VECTOR to a double array */
/*------------------------------------------------------*/
void vec_arr(int ncontpts, vector **contpts, double **array)
{
  short i;

  for (i=0; i<ncontpts; i++) {
    array[i][0] = contpts[i]->x;
    array[i][1] = contpts[i]->y;
    array[i][2] = contpts[i]->z;
    array[i][3] = contpts[i]->w;
  }
}

/*------------------------------------------------------*/
/* This function is to convert double array to a VECTOR */
/*------------------------------------------------------*/
void arr_vec(int ncontpts, double **array, vector **contpts)
{
  short i;

  for (i=0; i<ncontpts; i++) {
    contpts[i]->x = array[i][0];
    contpts[i]->y = array[i][1];
    contpts[i]->z = array[i][2];
    contpts[i]->w = array[i][3];
  }
}

/*------------------------------------------------------------------------*/
/* This function is to split a b-spline which has been subdivided already */
/* into two pieces                                                        */
/*------------------------------------------------------------------------*/
void split_it(double u, ParCurv *father, ParCurv *kid1, ParCurv *kid2)
{
  short i, j;
  int dadnk, kidnk;   /* the number of kids' knot vector */

  kidnk = kid1->order + kid1->ncontpts;
  dadnk = father->order + father->ncontpts;

  /* split the knot vector in the u */
  for (i=j=0; father->knots[i]<=u; i++,j++)   /* copy the knot vector */
    kid1->knots[i] = father->knots[i];

  for (i=0; i<kidnk; i++)
    kid2->knots[i] = father->knots[i + j - father->order];      

  /* split the contpts where the two consecutive contpts are the same */
  for (i=j=0; vec_neq(father->contpts[i], father->contpts[i+1]); i++,j++)
    copyvector( father->contpts[i] , kid1->contpts[i] );
  
  copyvector(father->contpts[j], kid1->contpts[j]);

  for (i=j+1; i<father->ncontpts; i++)
    copyvector(father->contpts[i], kid2->contpts[i-j-1]);
}

/*---------------------------------------------------------------*/
/* This function is to compare two vectors to see if they are != */
/*---------------------------------------------------------------*/
int vec_neq(vector *vector1, vector *vector2)
{
  if (vector1->x == vector2->x && vector1->y == vector2->y && 
      vector1->z == vector2->z && vector1->w == vector2->w)
    return 0;

  return 1;
}

/*---------------------------------------------*/
/* Decompose B-spline curve into bezier curves */
/*---------------------------------------------*/
int BsplineToBezier(ParCurv *egeom, ParCurv ***bezier)
{
  ParCurv *egm;
  vector *r;
  double *knots, **inpoly, **outpoly, tang;
  int nKnots, nPts;
  short i, j, nb, nc, order;

  order = egeom->order;
  nc    = egeom->ncontpts;

  nKnots = checknot(egeom);

  knots   = dbl_array1(nKnots);
  prodknot(egeom, knots);

  inpoly  = dbl_array2(nKnots, 4);
  outpoly = dbl_array2(nKnots, 4);

  for(i=0; i<nc; i++) {
    inpoly[i][0] = egeom->contpts[i]->x;
    inpoly[i][1] = egeom->contpts[i]->y;
    inpoly[i][2] = egeom->contpts[i]->z;
    inpoly[i][3] = egeom->contpts[i]->w;
  }
  curve_oslo1(egeom->order, egeom->ncontpts, nKnots, egeom->knots, knots,
	      inpoly, outpoly);
  free_darray2(inpoly);

  nb = (nKnots-order)/order;
  (*bezier) = (ParCurv **)gen_array1(nb, sizeof(ParCurv *));

  for (j=0; j*order < nKnots-order; j++) {
    (*bezier)[j] = egeomalloc1(order, order);
    (*bezier)[j]->order    = order;
    (*bezier)[j]->ncontpts = order;

    for(i=0; i<order; i++)
      (*bezier)[j]->knots[i] = 0.0;

    for(i=order; i<2*order; i++)
      (*bezier)[j]->knots[i] = 1.0;

    for(i=0; i<(*bezier)[j]->order; i++) {
      (*bezier)[j]->contpts[i]->x = outpoly[i + j*order][0];
      (*bezier)[j]->contpts[i]->y = outpoly[i + j*order][1];
      (*bezier)[j]->contpts[i]->z = outpoly[i + j*order][2];
      (*bezier)[j]->contpts[i]->w = outpoly[i + j*order][3];
    }
  }
  free_darray2(outpoly);

  return nb;
}

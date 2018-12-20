/* Copyright (C) 1992 Massachusetts Institute of Technology, Cambridge, MA
   All rights reserved
*/

/* decasteljau_curve.c */

/* build_w_table()
   decasteljau_curve()
   decasteljau_curve_contpts()
   linear_interpolate_vector()
   linear_interpolate_w()
*/

#include <stdio.h>
#include "gen.h"
#include "editor.h"

/*----------------------------------------------------------------*/
/*  This program subdivides a rational bezier curve               */ 
/*  into 2 rational bezier curves using decasteljau algorithm.    */
/*----------------------------------------------------------------*/

void decasteljau_curve(ParCurv *egeome, double t, ParCurv **egeome_1,
		       ParCurv **egeome_2)
{
  double a, b;

  a =  egeome->knots[0];
  b = egeome->knots[egeome->kmem-1];

  *egeome_1 = egeomalloc1(egeome->order, egeome->ncontpts);
  *egeome_2 = egeomalloc1(egeome->order, egeome->ncontpts);

  bezier_knots(egeome->order, a, t, (*egeome_1)->knots);
  bezier_knots(egeome->order, t, b, (*egeome_2)->knots);

  decasteljau_curve_contpts(egeome->contpts, t, egeome->order, 
			    (*egeome_1)->contpts, (*egeome_2)->contpts, a, b);
}

/*--------------------------------------------------------*/
/* This program subdivides rational bezier curve  control */ 
/* points into 2 ratinal bezier curve control points.     */
/*--------------------------------------------------------*/

void decasteljau_curve_contpts(vector **contpts, double t, int n,
			       vector **contpts_1, vector **contpts_2,
			       double a, double b)
{
  vector ***nxn_table;
  int level; /* level */
  int i,j;
  double **w_table;
  vector zero_vec;

  zero_vec.x= zero_vec.y= zero_vec.z= 0.0;
  zero_vec.w=1.0;

  nxn_table = vec_array2(n,n);
  w_table = build_w_table(contpts, t, n, a, b);

  for(i=0; i<n; i++)
    copyvector(contpts[i], nxn_table[i][0]); 

  for(i=0 ; i<n ; i++) {
    nxn_table[i][0]->x = nxn_table[i][0]->x / nxn_table[i][0]->w;
    nxn_table[i][0]->y = nxn_table[i][0]->y / nxn_table[i][0]->w;
    nxn_table[i][0]->z = nxn_table[i][0]->z / nxn_table[i][0]->w;
    nxn_table[i][0]->w = 1.0;
  }

  level = 1;

  for(i=1; i<n; i++) {
    level = i;
    for (j=i; j < n; j++)
      linear_interpolate_vector(t, nxn_table[j][level], w_table[j][level],
				nxn_table[j-1][level-1],
				w_table[j-1][level-1], nxn_table[j][level-1],
				w_table[j][level-1], a, b);
  }
     
  for(i=0; i<n; i++) 	 
    for (j=i; j < n; j++) {
      nxn_table[j][i]->w = (w_table[j][i]) * (nxn_table[j][i]->w);
      nxn_table[j][i]->x = (w_table[j][i]) * (nxn_table[j][i]->x);
      nxn_table[j][i]->y = (w_table[j][i]) * (nxn_table[j][i]->y);
      nxn_table[j][i]->z = (w_table[j][i]) * (nxn_table[j][i]->z);
    }
  for(i=0; i<n; i++) 
    copyvector(nxn_table[i][i], contpts_1[i]);

  for(i=0; i<n; i++) 
    copyvector( nxn_table[n-1][n-i-1], contpts_2[i]);

  free_varray2(nxn_table, n , n);
  free_darray2(w_table);
}

void linear_interpolate_vector(double t, vector *result, double w1,
			       vector *vec1, double w2, vector *vec2,
			       double w3, double a, double b)
{
  vector temp1, temp2;

  scale_vect1((b - t)/(b - a)*w2/w1, vec1, &temp1);
  scale_vect1((t - a)/(b - a)*w3/w1, vec2, &temp2);
  add_vect1(&temp1, &temp2, result);
}

double **build_w_table(vector **contpts, double t, int n, double a, double b)
{
  int level; /* level */
  int i,j;
  double **w_table;
  
  w_table = dbl_array2(n, n);

  for(i=0; i<n; i++)
    w_table[i][0] = contpts[i]->w; 
  level = 1;

  for(i=1; i<n; i++) {
    level = i;
    for (j=i; j < n; j++)
      linear_interpolate_w(t, &w_table[j][level], w_table[j-1][level-1],
			   w_table[j][level-1], a, b);
  }
     
  return w_table;
}

void linear_interpolate_w(double t, double *result, double vec1,
			  double vec2, double a, double b)
{
  *result = (b - t)/(b - a)*vec1 + (t - a)/(b - a)*vec2;
}

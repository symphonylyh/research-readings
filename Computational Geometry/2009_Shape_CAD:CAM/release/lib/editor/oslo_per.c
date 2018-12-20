/************************************************************************
 *                                                                      *
                        Copyright (C) 1992 by
        Massachusetts Institute of Technology, Cambridge, MA
                         All rights reserved
 *                                                                      *
 ************************************************************************/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "gen.h"
#include "bspl.h"
#include "editor.h"

#define MACHPREC 1.1e-13

void curve_oslo_per(short order, int ncontpts, short q, double *tau,
		    double *t, double **inper, double **outper)

     /* q      number of knots in the new knot vector          */
     /* inper  array of original control points, [ncontpts][4] */
     /* tau    vector of original knots, [ncontpts+1]          */
     /* t      vector of new augmented knots, [q]              */
     /* outper array of new control points,  [q-1][4]          */
{
  short i, j; 
  double *t1, **inpoly, **outpoly;

  t1 = dbl_array1(q);

  inpoly = dbl_array2(ncontpts + 2*order, 4);
  outpoly = dbl_array2(q + 2*order, 4);

  for(i=0; i<ncontpts; i++) 
    for (j=0; j<4; j++) 
      outper[i][j] = inper[i][j];
   
  i = j = 0; 
  
  while (fabs(t[i]-tau[i]) < MACHPREC) 
    i++;
  if (i<order) {
    do t1[j] = t[j]; 
    while (t[j++] <= tau[order-1] - MACHPREC);
    for (i=order; i<=ncontpts; i++)
      t1[j++] = tau[i];
    curve_osloper_lt(order, &ncontpts, j,tau, t1, inpoly, outpoly, outper);
    if (j<q) 
      curve_osloper_rt(order, &ncontpts, q, tau, t, inpoly, outpoly, outper);
  }
  else
    curve_osloper_rt(order, &ncontpts, q, tau, t, inpoly, outpoly, outper);

  free_darray1(t1);
  free_darray2(inpoly);
  free_darray2(outpoly);
}

void curve_osloper_rt(short order, int *ncontpts, short q, double *knots,
		      double *t, double **inpoly, double **outpoly,
		      double **outper)

     /* q        number of knots in the new knot vector          */
     /* **inpoly array of original control points, [ncontpts][3] */
     /* knots    vector of original knots, [ncontpts+order]      */
     /* t        vector of new augmented knots, [q]              */
     /* outpoly  array of new control points,  [q-1][3]          */
{
  short i, j;
  double *t1, *tau;

  t1 = dbl_array1(q+order);
  tau = dbl_array1(q+order); 

  for (i=0; i<=(*ncontpts); i++)
    tau[i] = knots[i] ;
  for (i=1; i<order; i++) 
    tau[(*ncontpts)+i] = tau[(*ncontpts)+i-1] + knots[i]-knots[i-1];
  for (i=0; i<q; i++)
    t1[i] = t[i];
  for (i=1; i<order; i++)
    t1[q+i-1] = t1[q+i-2] + t[i]-t[i-1];

  for (i=0; i<(*ncontpts); i++)
    for (j=0; j<4; j++)
      inpoly[i][j] = outper[i][j];

  i = q+order-1;
  curve_oslo2(order, *ncontpts, i, tau, t1, inpoly, outpoly);

  for (i=0; i<q; i++)
    knots[i] = t[i];
  
  *ncontpts = q-1; 
  for (i=0; i< *ncontpts; i++) 
    for (j=0; j<4; j++)
      outper[i][j] = outpoly[i][j];

  free_darray1(t1);
  free_darray1(tau);
}

void curve_osloper_lt(short order, int *ncontpts, short q, double *knots,
		      double *t, double **inpoly, double **outpoly,
		      double **outper)
{
  short i, j, k, m;
  double *t1, *tau;

  t1 = dbl_array1(q+order);
  tau = dbl_array1(q+order); 

  for (i=0; i<= *ncontpts; i++)
    tau[i+order-1] = knots[i];
  for (i=order-1; i>0; i--) {
    j = *ncontpts + i + 1 - order; 
    tau[i-1] = tau[i]-knots[j]+knots[j-1];
  }
  for (i=0; i<q; i++)
    t1[i+order-1] = t[i] ;
  for (i=order-1; i>0; i--) {
    j = q + i - order; 
    t1[i-1] = t1[i]-t[j] + t[j-1];
  }

  for (i=0; i< *ncontpts; i++) {
    j = (i+order-1) % (*ncontpts);
    for (m=0; m<4; m++) 
      inpoly[j][m] = outper[i][m];
  }

  i = q+order-1; 
  curve_oslo2(order, *ncontpts, i, tau, t1, inpoly, outpoly);
  
  for (i=0; i<q; i++)
    knots[i] = t[i];

  *ncontpts = q-1;

  for (i=0; i< *ncontpts; i++) {
    j = (i+order-1) % (*ncontpts);
    for (m=0; m<4; m++) 
      outper[i][m] = outpoly[j][m];
  }
  free_darray1(t1);
  free_darray1(tau);
}

void curve_oslo2(short order, int ncontpts, short q, double *tau,
		 double *t, double **inpoly, double **outpoly)
{
  double2 **au;

  au = dbl_array2(q-1, q);

  trans_per(order, ncontpts, q, tau, t, au);

  matrixmu(q-1, ncontpts, 4, q, 4, 4, &au[0][0],&inpoly[0][0],&outpoly[0][0]);
  free_darray2(au);
}

void trans_per(short k, short n, short q, double *tau, double *t, 
	       double **a)
{
  short i, j, l, mu;
  double **alpha;
  alpha = dbl_array2(q, k+1);

  for(i=0; i<q-1; i++)
    for(j=0; j<q; j++)
      a[i][j] = 0.0;
  
  for(j=0; j<q-k; j++) {
    mu = find(k+n, tau, t[j]); 

    for (i=0; i<q; i++)   /* set all alpha to zero */
      for (l=0; l<k+1; l++)
	alpha[i][l] = 0.0;

    disc_bsp(mu, j, k, t, tau, alpha);
    
    for(i=MAX(mu-k+1, 0); i<=mu; i++) {
      a[j][i] = alpha[i][k];
    }
  }
  free_darray2(alpha);
}

void disc_bsp(short mu, short j, short k, double *t, double *tau, 
	      double **alpha)
{
  short r, i, mu2;
  double beta_1, beta, d1, d2, tj;

  alpha[mu][1] = 1.0;

  mu2 = mu;

  for(r=1; r<k; r++) {
    beta_1 = 0.0;
    tj = t[j+r];

    for(i=MAX(mu2,1); i<=mu; i++) {
      d1 = tj - tau[i];
      d2 = tau[i+r] - tj;
      beta = alpha[i][r] / (d1 + d2);
      alpha[i-1][r+1] = d2*beta + beta_1;
      beta_1 = d1*beta;
    }
    alpha[mu][r+1] = beta_1;
    mu2 = mu2 - 1;
  }
  if(mu == 0)
    alpha[mu][k] = 1.0;
}	

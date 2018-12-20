/************************************************************************
			    Copyright (C) 1995
     Massachusetts Institute of Technology, Cambridge, Massachusetts.
			   All rights reserved.

     This portion of source code was prepared in the Ocean Engineering
     Design Laboratory by Jingfang Zhou.
     Modified by S. L. Abrams, 13 June 1995
 ************************************************************************/

#include <iostream>
#include <cmath>
#include "conv.h"
#ifdef USE_INTERVAL
#include "interval.h"
typedef Interval real;
#else
#ifdef USE_RAT
#include "rat.h"
typedef Rat real;
#else
typedef double real;
#endif
#endif

using namespace std;

void convPolyBox(real **bp, short **ordlist, real **intv, int numvar,
		 int numeq)
{
  for (int i=0; i<numeq; i++)
    convOneEq(bp[i], ordlist[i], intv, numvar);
}

#ifdef USE_RAT
static Int binomial (int n, int k)   /* Binomial */
{
  Int result = 1;
  if (k > n) return (0);
  for (int i=1; i <= k; i++)
    result *= (Int)(n-i+1);
  for (i=1; i <= k; i++)
    result /= (Int)i;
  return (result);
}
#else
static real binomial (int n, int k)   /* Binomial */
{
  int i;
  real result = 1;
  if (k > n) return (0);
  for (i=1; i <= k; i++)
    result *= (real)(n-i+1);
  for (i=1; i <= k; i++)
    result /= (real)i;
  return (result);
}
#endif

void convOneEq(real *&bp, short *ordlist, real **intv0, int numvar)
{
  int i,j;
  int cosize=1;
  real **intv = new real*[numvar];
  for (i=0; i<numvar; i++){
    intv[numvar-i-1] = new real[2];
    for (j=0; j<2; j++)
      intv[numvar-i-1][j] = intv0[i][j];
  }
    
  for (i=0; i<numvar; i++)
    cosize *= (ordlist[i]);
  
  real *new_bp = new real[cosize], *tmp_bp;
  
  for (i=0; i < numvar; i++) {
    int j,k,l,m;		/* Loop counters */
    int lastl=1, skip=1;	/* Loop variables */
    int lastj=ordlist[i];
    for (j=0; j < i; j++)	/* Set skip */
      lastl *= ordlist[j];
    for (j=i+1; j < numvar; j++)
      skip *= ordlist[j];	/* Set lastl */
    for (l=0; l<lastl; l++) {
      for (k=0; k < skip; k++) {
	real *vect = new real[lastj];	/* Allocate multiplication vector */
	real *p = &bp[k + cosize * l / lastl];	     /* Set p */
	real *q = &new_bp[k + cosize * l / lastl];   /* Set q */
	for (j=0; j < lastj; j++)
	  vect[j] = p[j*skip];	/* Assign vect[] */
	for (j=0; j < lastj; j++) {
	  real sum = 0;         /* The summation */
	  for (m=0; m < lastj-j; m++) {
#ifdef USE_RAT
	    real t1 = (intv[i][0]) ^ (long)m;
	    real t2(binomial(m+j, m), Int(1));
#else
	    real t1 = pow(intv[i][0], m);
	    real t2 = binomial(m+j, m);
#endif

	    t1 *= t2;
	    sum += t1 * vect[m+j];
	  }
	  real t1 = intv[i][1]-intv[i][0];
#ifdef USE_RAT
	  t1 ^= (long)j;
#else
	  t1 = pow(t1, j);
#endif
	  q[j*skip] = t1 * sum;	/* Reassign into array */
	}
	delete [] vect;         /* redundant delete[lastj] vect */
      }
    }
    tmp_bp = bp;
    bp = new_bp;                /* Keep work from last step */
    new_bp = tmp_bp;
  }

  delete [] new_bp;
  for (i=0; i<numvar; i++)
    delete [] intv[i];
  delete [] intv;
}

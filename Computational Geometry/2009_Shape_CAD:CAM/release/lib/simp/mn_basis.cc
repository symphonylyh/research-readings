// Copyright (C) Massachusetts Institute of Technology, 1995
// All rights reserved

#include <cmath>
#include "multinom.h"

using namespace std;

#ifdef USE_RAT
static Int binomial (int n, int k)

/* Binomial */

{
  int i;
  Int result = 1;
  if (k > n) return (0);
  for (i=1; i <= k; i++)
    result *= (Int)(n-i+1);
  for (i=1; i <= k; i++)
    result /= (Int)i;
  return (result);
}

static Rat coeff (int j, int k, int n)

/* The conversion coefficient lambda (j,k) of order n */

{
  Rat R(binomial(k,j),
	binomial(n,j));
  return (R);
}

static Rat Coeff (int j, int k, int n)

/* The conversion coefficient lambda (j,k) of order n */

{
  Int pow;   Int one = 1;  Int mone = -1;

  if( (k-j)%2 ) pow = mone;
  else pow = one; 
  Int iR = binomial(n,k)*binomial(k,j);
  iR = iR * pow;
  Rat R( iR, one );
  return (R);
}

#else
static real binomial (int n, int k)

/* Binomial */

{
  int i;
  real result = 1;
  if (k > n) return (0);
  for (i=1; i <= k; i++)
    result *= (n-i+1);
  for (i=1; i <= k; i++)
    result /= i;
  return (result);
}

static real coeff (int j, int k, int n)

/* The conversion coefficient lambda (j,k) of order n */

{
  real R = binomial(k,j) / binomial(n,j);
  return (R);
}


static real Coeff (int j, int k, int n)

{
  real pow;
  if( (k-j)%2 ) pow = -1.;
  else pow = 1.; 
  real R = pow*binomial(n,k)*binomial(k,j);
  return (R);
}
#endif



multinomial *multinomial::monotobern (multinomial *mn_out=0) const

/* Convert a multinomial from monomial to Bernstein basis. */

{
  multinomial *mn_new;		/* New multinomial */
  const multinomial *in = this;
  if (!mn_out)
    mn_new = new multinomial
      (ordlist);
  else
    mn_new = mn_out;
  int i;
  for (i=0; i < in->ndim(); i++) {
    int j,k,l,m;		/* Loop counters */
    int lastl=1, skip=1;	/* Loop variables */
    int lastj=in->ordlist[i]+1;
    for (j=0; j < i; j++)	/* Set skip */
      lastl *= in->ordlist[j] + 1;
    for (j=i+1; j < in->ndim(); j++)
      skip *=			/* Set lastl */
	in->ordlist[j] + 1;
    for (l=0; l < lastl; l++) {
      for (k=0; k < skip; k++) {
	real *vect =		/* Allocate multiplication vector */
	  new real[lastj];
	real *p =		/* Set p */
	  &bp[k + in->cosize *
		  l / lastl];
	real *q =		/* Set q */
	  &mn_new->bp[k + in->cosize *
		   l / lastl];

/*
	cout << "Initial p is at arraypos " << k + in->cosize*l/lastl <<
	  " skipsize is " << skip << endl;
*/

	for (j=0; j < lastj; j++)
	  vect[j] = p[j*skip];	/* Assign vect[] */
	for (j=0; j < lastj; j++) {
	  real sum = 0.;	 /* The summation */
	  for (m=0; m < lastj; m++) {
	    sum +=
	      coeff(m,j,in->ordlist[i]) * vect[m];
	  }
	  q[j*skip] = sum;	/* Reassign into array */
	}
	delete[] vect;  /* redundant delete[lastj] vect */
      }
    }
    in = mn_new;		/* Keep work from last step */
  }
  return (mn_new);
}



multinomial *multinomial::berntomono (multinomial *mn_out=0) const

/* Convert a multinomial from Bernstein to monomial basis. */

{
  multinomial *mn_new;		/* New multinomial */
  const multinomial *in = this;
  if (!mn_out)
    mn_new = new multinomial
      (ordlist);
  else
    mn_new = mn_out;
  int i;
  for (i=0; i < in->ndim(); i++) {
    int j,k,l,m;		/* Loop counters */
    int lastl=1, skip=1;	/* Loop variables */
    int lastj=in->ordlist[i]+1;
    for (j=0; j < i; j++)	/* Set skip */
      lastl *= in->ordlist[j] + 1;
    for (j=i+1; j < in->ndim(); j++)
      skip *=			/* Set lastl */
	in->ordlist[j] + 1;
    for (l=0; l < lastl; l++) {
      for (k=0; k < skip; k++) {
	real *vect =		/* Allocate multiplication vector */
	  new real[lastj];
	real *p =		/* Set p */
	  &bp[k + in->cosize *
		  l / lastl];
	real *q =		/* Set q */
	  &mn_new->bp[k + in->cosize *
		   l / lastl];

/*
	cout << "Initial p is at arraypos " << k + in->cosize*l/lastl <<
	  " skipsize is " << skip << endl;
*/

	for (j=0; j < lastj; j++)
	  vect[j] = p[j*skip];	/* Assign vect[] */
	for (j=0; j < lastj; j++) {
	  real sum = 0.;	 /* The summation */
	  for (m=0; m < lastj; m++) {
	    sum +=
	      Coeff(m,j,in->ordlist[i]) * vect[m];
	  }
	  q[j*skip] = sum;	/* Reassign into array */
	}
	delete[] vect;  /* redundant delete[lastj] vect */
      }
    }
    in = mn_new;		/* Keep work from last step */
  }
  return (mn_new);
}


























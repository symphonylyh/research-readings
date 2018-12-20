/****************************************************************************
  mn_sub.cc            last edit 10.2.92
 ***************************************************************************/

#include "lib1.h"
#include "multinom.h"
using namespace std;

/* Subdivision of an explicit Bernstein multinomial */

multinomial *multinomial::sub
  (short index,                 /* Coordinate direction in which to sub. */
   real a, real b,              /* The interval in [0,1] in which to sub. */
   multinomial *mn_out=0) const /* Optional output multinomial */

/* Subdivides an explicit Bernstein polynomial in the direction
   'index.' This operation may be VERY expensive, as the looping
   below should indicate. The routine allocates its own new multinomial,
   to avoid overwriting the old. Returns the new multinomial on
   exit, or 0 if allocation fails. Note: 'index' runs from 0 to ndim-1 on
   input. The optional argument mn_out is provided to allow us to overwrite
   an already allocated multinomial. */

{
  multinomial *mn_new;          /* New multinomial */
  int i,j,k,l;                  /* Loop counters */
  int skip=1, step, lastl=1;    /* Loop terminators and incrementers */
  short biksize =               /* Set size of array */
    ordlist [index] + 1;
  real_array2 bik               /* Coefficient array */
    (biksize, biksize);
  real *p;                      /* Pointer to new coefficient block */
  real *q;                      /* Pointer to old coefficient block */

//  cout << " in sub " << endl;
//  cout << " a= " << a << " b= " << b << endl;



#ifndef USE_RAT
  if (b != 0.)
#else
  if (b != Rat(0,1))
#endif
    a  = a / b;                     /* Effective value of a in [0,b] */
  
//  cout << " a= " << a << " b= " << b << endl; 

  real abar ;
#ifndef USE_RAT
  abar =  1.-a;	        /* Just 1-a */
#else
  abar = Rat(1,1)-a;
#endif
  real bbar;
#ifndef USE_RAT
  bbar =  1.-b;	        /* 1-b */
#else
  bbar =  Rat(1,1)-b;
#endif
//  cout << " index = " << index << endl;
//  cout << " a= " << a << " b= " << b << endl; 
//  cout << " abar = " << abar << " bbar = " <<bbar << endl; 

  if (!mn_out)                  /* mn_out does not exist */
    mn_new = new multinomial    /* Allocate the new multinomial */
      (ordlist);
  else
    mn_new = mn_out;
  p = mn_new->bp;               /* Initialize p and q */
  q = bp;
  for (i=0; i < index; i++)
    lastl *= ordlist[i] + 1;    /* Set lastl */
  for (i=index+1; i < ndim(); i++)
    skip *= ordlist[i] + 1;     /* Set skip */

/* Begin subdivision. The variable p is used to traverse the coefficient
   block. */

  for (l=0; l < lastl; l++) {   /* Outer loop */
    for (j=0; j < skip; j++, p++, q++) {   /* Inner loop */

/* For this iteration, get all ordlist[index]+1 coefficients taken by
   varying the index'th coefficient and holding the rest fixed */

      for (step=0, i=0; i <= ordlist[index]; i++, step += skip)
	bik[i][0] =             /* Put the coefficients in */
	  q[step];
      for (k=1; k <= ordlist[index]; k++){   /* Actual subdivison */
	for (i=k; i <= ordlist[index]; i++){
	  bik [i][k] =
	    bbar * bik [i-1][k-1] +
	      b * bik [i][k-1]; /* Now from 0 to b */
//	    cout << "bik[" << i << "][" << k << "]=" << bik[i][k] <<endl;
     }
   }
      for (i=0; i <= ordlist[index]; i++)
	bik [i][0] = bik [i][i];   /* Get set for second subdivision */
      for (k=1; k <= ordlist[index]; k++)  { /* Second subdivision */
	for (i=k; i <= ordlist[index]; i++){
	  bik [i][k] =
	    abar * bik [i-1][k-1] +
	      a * bik [i][k-1];   /* Now from a to b */
//	    cout << "bik[" << i << "][" << k << "]=" << bik[i][k] <<endl;
      }
   }
      for (step=0, i=0; i <= ordlist[index]; i++, step += skip)
	 p[step] =              /* Put in the subdivided values */
	   bik [ordlist[index]][ordlist[index]-i];
    }
    p += skip * ordlist[index];   /* Increment p and q */
    q += skip * ordlist[index];
  }
//  cout << " mn_new = " << *mn_new << endl;

  return (mn_new);              /* And return */
}

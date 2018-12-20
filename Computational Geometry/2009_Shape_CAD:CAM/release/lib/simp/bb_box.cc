/********************************************************************
  bb_box.cc          last edit 10.8.92
 *******************************************************************/
#include "bbox.h"
#include "interval.h"

#ifdef USE_RAT
static real crit(8,10);
#else
static real crit=0.8;
#endif

using namespace std;

/* Functions for supporting the Bbox structure */

Bbox::Bbox(int nd) : a(nd), b(nd)

/* Constructor function for the box. Make it equal, by default, to
   [0,1]^n */

{
  int i;                        /* Loop counter */
  for (i=0; i < nd; i++) {
#ifndef USE_RAT
    a[i] = 0.0;
    b[i] = 1.0;
#else
    a[i] = Rat(0,1);
    b[i] = Rat(1,1);
#endif
  }
  n = nd;
  stat = 0;
}

Bbox::~Bbox() {}                    /* Destructor */

Bbox::Bbox(Bbox &B) : a(B.n), b(B.n)

/* Construct a box from another box. */

{
  int i;
  for (i=0; i < B.n; i++) {
    a[i] = B.a[i];
    b[i] = B.b[i];
  }
  n = B.n;
  stat = B.stat;
}

void Bbox::copy (Bbox &B)

/* Copy box B into current box. Assume that number of dimensions is
   the same. */

{
  int i;
  for (i=0; i < B.n; i++) {
    a[i] = B.a[i];
    b[i] = B.b[i];
  }
  stat = B.stat;
}

void Bbox::scale (Bbox &B2)

/* Scale the current box with respect to a second box. In other words,
   if [0,1]^n is mapped to B2, what is our box mapped to? */

{
  int i;                        /* Loop counter */

#ifdef USE_INTERVAL

  for (i=0; i < n; i++) {
       real eps = B2.b[i] - B2.a[i];
       a [i] *= eps;
       a [i] += B2.a[i];
       a [i] = a[i].lower();
       b [i] *= eps;
       b [i] += B2.a[i];
       b [i] = b[i].upper();
  }
#else

  for (i=0; i < n; i++) {
       real eps = B2.b[i] - B2.a[i];
       a [i] *= eps;
       a [i] += B2.a[i];
       
       b [i] *= eps;
       b [i] += B2.a[i];
       
  }
#endif
}

int Bbox::small (real eps)

/* Return 1 if every side of the box is <= eps, otherwise return 0. */

{
  int i;                        /* Loop counter */
  for (i=0; i < n; i++)
    if ((b[i] - a[i]) > eps)
      return 0;
  return 1;
}



void Bbox::subcheck (Bbox &B, int_array &check, int_array &ret)

/* Return an integer array of n elements. If the ith element of the array is
   1, binary subdivision in that direction is required. Otherwise, the
   element is 0. The status flag is set to 0 if no binary subdivision is
   required, otherwise to 1. Whether or not binary subdivision is required
   is determined by the constant crit, which is defined (hardwired) above. */

{
  int i=0;                      /* Loop counter */
  int last=-1;
  stat = 0;                     /* Initial value of status */
  for (i=0; i < n; i++) {
    if (!check[i]) {
      ret[i] = 0;
      continue;
    }
    if ((b[i]-a[i]) * crit <=
	(B.b[i] - B.a[i])) {
      stat = 1;
      if (last != -1) {
	if (B.b[i] - B.a[i] >
	    B.b[last] - B.a[last]) {
	  ret[i] = 1;
     //	  cout << ' ' << i << ' ';
	  last = i;
	} else ret[i] = 0;
      } else {
	last = i;
	ret[i] = 1;
     //	cout << ' ' << i << ' ';
      }
    } else
      ret[i] = 0;
  }
}

ostream &operator << (ostream &stream, Bbox &B)

/* A stream inserter for boxes */

{
     int i;
     real low;
     real upper;
     for (i=0; i < B.n; i++) {

	  low = B.a[i];
	  upper = B.b[i];

	  stream << "[" << low << "," <<
	       upper << "]";
	  if (i != B.n-1)
	       stream << "x";
     }
     return (stream);
}








#include "multinom.h"

#include <iostream>
using namespace std;

#ifndef NO
#define NO 0
#endif 

#ifndef YES
#define YES 1
#endif
/* Special allocation for the real_array given an ordlist */

real_array::real_array (sht_array olist)

{
  int i,n=1;                
  for (i=0; i < olist.size(); i++)
    n *= (olist[i] + 1);
  init (n);
}



mn_array::mn_array (mn_array &mna)

/* Copy constructor for mn_array */

{ 
  init (mna.size());
  for (int i=0; i < n; i++)
    ((multinomial **) p) [i] =
      new multinomial (mna[i]);
}

mn_array::~mn_array ()

/* Destructor */

{
  for (int i=0; i < n; i++)
    delete ((multinomial **) p) [i];
  dealloc ();
}

int mn_array::every_points_not_far_away_from_axis_by_eps(real re)
{
     for (int i=0; i < size(); i++)
	  if ((*this)[i].every_points_not_far_away_from_axis_by_eps(re) == NO)
	      return NO;
     return YES;
}



int mn_array::convex_hull_cross_every_axis()
{
     for (int i=0; i < size(); i++)
	  if ((*this)[i].convex_hull_cross_every_axis() == NO)
	      return NO;
     return YES;
}

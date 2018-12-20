// Copyright (C) 1995 Massachusetts Institute of Technology, Cambridge, MA
// All rights reserved

#ifndef MULTINOM_H		/* Header file for simultaneous poly. */
#define MULTINOM_H

//#include <stdlib.h>
#include "lib1.h"
#include "compiler.h"

/* The following structure is a generalized structure for multinomials.
   The internal representation consists of
     ordlist = An array indicating the degree of the multinomial
	       in each direction
     bp = "Base pointer" to the coefficient block 
	  (i.e. start of coefficients)
     cosize = the length of bp (# of coefficients)
   These fields are all "const" -- once they're set, we don't want anything
   short of nuclear war to change them.

   Notice that there is nothing that requires the multinomial to be expressed 
   in a specific basis. Bernstein or monomial are both acceptable. */

class multinomial {
public:

  sht_array ordlist;      /* List of orders */

  real_array bp;          /* Base pointer to coefficients */

  int cosize;                   /* Size of coefficient block */

  multinomial (const sht_array &); 
  multinomial (multinomial &mn);

  ~multinomial();
  multinomial *sub              /* Subdivide the multinomial */
    (short, real, real,
     multinomial *) const;
  multinomial *monotobern	/* Do a basis conversion */
    (multinomial *) const;
  multinomial *berntomono	/* Do a basis conversion */
    (multinomial *) const;
  int ndim () const             /* Data abstraction */
    { return ordlist.size(); }
  
  int every_points_not_far_away_from_axis_by_eps(real);
  int convex_hull_cross_every_axis();

  friend std::istream &operator      /* Input */
    >> (std::istream&, multinomial*&);
  friend std::ostream &operator      /* Output */
    << (std::ostream&, multinomial&);
};

/* Our next class object is an array of multinomials, which we will use
   frequently. It inherits characteristics from the generic_array structure.
   It is structured as an array of multinomial pointers so that
   initialization is done individually. However, the copy constructor
   copies everything, not just pointers. Furthermore, the destructor frees
   everything. */

class mn_array : public generic_array<multinomial*> { /* Array of multinomials */
protected:
  mn_array () {}                /* Default constructor -- for inheritance */
  void init (int siz)           /* Initializer */
    { p = alloc (siz); n = siz; }
public:
  mn_array (int siz)            /* Create the mn_array, pointers only */
    { init (siz); }
  mn_array (mn_array &mna);	/* Copy constructor */
  ~mn_array ();			/* Destructor -- kill everything */
   multinomial &operator[] (int i) /* Indexer */
     {  return *p[i];}
   //   { return *(((multinomial **) p)[i]); }
  int every_points_not_far_away_from_axis_by_eps(real);
  int convex_hull_cross_every_axis();

  operator multinomial** ()     /* Conversion to array of pointers */
  { return (multinomial **) p; }
};

#endif

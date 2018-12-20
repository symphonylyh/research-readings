// Copyright (C) 1995 Massachusetts Institute of Technology, Cambridge, MA
// All rights reserved

#ifndef COMPILER_H		/* Compilation-specific stuff */        
#define COMPILER_H

#include <iostream.h>
#include <iomanip.h>
#include "lib1.h"

#define TRUE 1
#define FALSE 0
#define MAYBE -1

#ifdef USE_INTERVAL
#include "interval.h"
typedef Interval real;

class real_array : public Interval_array {
 public:
     real_array();
     real_array (int);
     real_array (sht_array);
     real_array (real_array& );
     ~real_array();
     real_array operator=(real_array);
     
  operator real* () const     /* Conversion operator */
    { return ((real *) p); }

  friend int if_overlap(real_array&, real_array&);
  friend ostream& operator<<( ostream& , real_array&);

  friend rlist_element;
};
typedef Interval_array2 real_array2;
#else
#ifdef USE_RAT
#include "rat.h"
typedef Rat real;

class real_array : public Rat_array {
 public:
     real_array();
     real_array (int);
     real_array (sht_array);
     real_array (real_array& );
     ~real_array();
     real_array operator=(real_array);
     
  operator real* () const     /* Conversion operator */
    { return ((real *) p); }

  friend int if_overlap(real_array&, real_array&);
  friend ostream& operator<<( ostream& , real_array&);

  friend rlist_element;
};
typedef Rat_array2 real_array2;
#else
typedef double real;

class real_array : public dbl_array {
public:  
  real_array();
  real_array (int siz) : dbl_array(siz) {}
  real_array (sht_array);
  real_array (real_array& );
  ~real_array();
  real_array operator=(real_array &ra)  
  {
       if (n == ra.n) memcpy (p, (real *) ra.p, n*sizeof(real));
       else {
	    resize(ra.n);
	    memcpy (p, (real *) ra.p, n*sizeof(real));
       }
       return *this;
  }
  
  operator real* () const     /* Conversion operator */
    { return ((real *) p); }

  friend int if_overlap(real_array&, real_array&);
  friend ostream& operator<<( ostream& , real_array&);
  friend rlist_element;
};
typedef double_array2 real_array2;
#endif
#endif

typedef real* real_pointer;
#endif

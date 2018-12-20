// Copyright (C) 1995 Massachusetts Institute of Technology, Cambridge, MA
// All rights reserved

#ifndef _RATPOLY_H_		/* Header file for simultaneous rat. poly. */
#define _RATPOLY_H_

//#include <stdlib.h>
#include "rat.h"
#include "lib1.h"
#include "compiler.h"

#define SMALL 1e-05
#define RP_TOL_NUM 1
#define RP_TOL_DEN 100

/* The following structure is a generalized structure for multinomials. At
   the topmost level, ndim is set to the dimension of the polynomial, order
   is set to the order (or degree) of the Bernstein polynomials in the
   first coordinate direction. The next field in the union is set to point
   to 'order' different multinomials, each of degree ndim-1. This process
   is continued recursively until we have 'order1' * 'order2' * 'order3' *
   ... * 'order(ndim-1)' different polynomials of dimension 1 (where
   order1, order2, etc. denote the orders of the multinomial in each of the
   coordinate directions). At this bottom level, the union field used is
   the coeffs field, which points to an array of coefficients. NOTE: The
   size of each array is order+1, because it is indexed from 0 .. order. 

   The three fields bp, cosize, and ordlist are only set in the "root"
   multinomial. Because of this, we waste quite a bit of space, but we
   save evaluation time and code. */

class multinomial {
public:

  short ndim;			/* Number of dimensions */
  Rat *bp;			/* Base pointer to coefficients */
  int cosize;			/* Size of coefficient block */
  short *ordlist;		/* List of orders */

  multinomial () {};		/* Constructor with no args */
  multinomial (short, short *);	/* Constructor, given ndim, ordlist */
  ~multinomial ();		/* Destructor */
  Rat *find (short *);		/* Address of a coefficient */
  Rat *check (short,Rat&,Rat&);	/* Check for intersection interval */
  multinomial *copy ();		/* Return a copy of the multinomial */
  Rat eval (Rat *);		/* Evaluate a multinomial */
  multinomial *partial (short);	/* Partial derivative */
  multinomial *sub		/* Subdivde the multinomial */
    (short, Rat &, Rat &,
     multinomial *);
  multinomial *monotobern
    (multinomial *) const;
  friend std::istream &operator	/* Input */
    >> (std::istream&, multinomial*&);
  friend std::ostream &operator	/* Output */
    << (std::ostream&, multinomial&);
};

struct stackel {
  Rat *a;
  Rat *b;
  multinomial **mn_list;
  int i0;
};

class Rstack {
  stackel *base;		/* Base of the stack */
  int maxsiz;			/* Maximum size it can attain */
  int n;			/* Number of elements */

public:
  Rstack();			/* Constructor */
  ~Rstack();			/* Destructor */
  void push(Rat *, Rat *, multinomial **, int);
  int pop(Rat *&, Rat *&, multinomial **&, int &);
};

struct Point {
  Rat x;
  Rat y;
};

struct sort_key {
  Rat dx;
  Rat dy;
  char flag;
  int index;
};				/* A sorting structure for the hull alg. */

typedef struct LIST_EL {	/* A linked-list element */
                        real *u;   /* N-dimensional coordinate */
			struct LIST_EL *next;
		       } list_element;

class list {
public:
  list_element *head;
  list_element *tail;
  friend list *concat (list*, list*);
};

Rat *Rat_array (int);
Rat **Rat_array_2 (int, int);
void free_Rat_array (Rat *, int);
void free_Rat_array_2 (Rat **, int, int);
void free_Point_array (Point *, int);
void free_sort (sort_key *, int);
int rp_install (multinomial **mn_list);
list *rp_minimize (Rat *a, Rat *b);
list *rp_inter (multinomial **mn_list, Rat &, int);

void rp_ginit (int);
void rp_draw_hull (Point *, int);
void rp_show_graph (Rat &, Rat &, Rat &, Rat &, int);
void rp_show_graph_2 (Rat &, Rat &);
void rp_reset_graph (int);

extern int rp_graph;

inline Rat *Rat_array (int len)

/* Allocate a group of Rats using calloc. Saves calling constructor
   functions. */

{
  return ((Rat *) calloc (len, sizeof (Rat)));
}

inline Point *Point_array (int len)

/* Allocate a group of Points using calloc to avoid calling Rat
   constructor. */

{
  return ((Point *) calloc (len, sizeof (Point)));
}

inline Point *Prealloc (Point *r, int len)

/* Call realloc on a point array; assumes that we are reducing the
   size and the last few elements don't need to be destroyed. */

{
  return ((Point *) realloc (r, len * sizeof (Point)));
}

inline sort_key *sort_array (int len)

/* Allocate a group of sort_keys using calloc to avoid calling Rat
   constructor. */

{
  return ((sort_key *) calloc (len, sizeof (sort_key)));
}

inline Rstack::Rstack ()

/* Constructor for stacks */

{
  base = (stackel *)		/* Allocate the base */
    malloc (sizeof (stackel));
  maxsiz = 1;			/* Initial max size */
  n = 0;			/* Empty */
}

inline Rstack::~Rstack ()

/* Destructor */

{
  free ((void *) base);
}

inline void Rstack::push (Rat *a, Rat *b, multinomial **mn_list, int i0)

/* Push elements onto the stack */

{
  if (n >= maxsiz) {
    maxsiz *= 2;		/* Resize it */
    base = (stackel *)
      realloc (base,
	       maxsiz * sizeof (stackel));
  }
  base[n].a = a;			/* Store the arguments */
  base[n].b = b;
  base[n].mn_list = mn_list;
  base[n].i0 = i0;
  n++;
}

inline int Rstack::pop (Rat *&a, Rat *&b, multinomial **&mn_list, int &i0)

/* Pop elements off of the stack */

{
  if (n==0)
    return (0);
  else {
    n--;
    a = base[n].a;
    b = base[n].b;
    mn_list = base[n].mn_list;
    i0 = base[n].i0;
    return (1);
  }
}

#endif

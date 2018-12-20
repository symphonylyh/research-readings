// Copyright (C) 1995 Massachusetts Institute of Technology, Cambridge, MA
// All rights reserved

#ifndef BBOX_H			/* Header file for boxes */
#define BBOX_H

#include "lib1.h"
#include "compiler.h"

#ifndef TRUE
#define TRUE 1
#define FALSE 0
#endif

class Bbox {

public:

  real_array a;                 /* Lower bounds */
  real_array b;                 /* Upper bounds */
  int n;                        /* Dimension of box */
  int stat;                     /* Status flag */
  Bbox(int);                    /* Constructor -- makes [0,1]^n */
  Bbox(Bbox&);                  /* Essentially a copying function */
  ~Bbox();                    /* Destructor */
  int flag();                   /* Status flag checker */
  void copy(Bbox &);            /* Copier */
  void scale(Bbox &);           /* Scale box w.r.t. another box */
  void swap(Bbox &);            /* Swap info w/ another box */
  int small(real);              /* Box small enough? */
  void subcheck                 /* Subdivision check */
    (Bbox &, int_array &, int_array &);
  friend std::ostream &operator      /* Stream inserter */
    << (std::ostream&, Bbox&);
};

inline int Bbox::flag()

/* Return the status flag */

{
  return stat;
}

#endif

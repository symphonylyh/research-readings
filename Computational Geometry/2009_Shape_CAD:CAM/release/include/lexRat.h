// Copyright (C) 1995 Massachusetts Institute of Technology, Cambridge, MA
// All rights reserved

#ifndef LEXPP_H
#define LEXPP_H

#ifdef USE_INTERVAL

#include "interval.h"
typedef Interval real;
#include "pp_solver.h"
#include "simpoly.h"

extern real atorv();

extern "C" { 
  int lexCode        (mn_array *&mn_list);
}
extern int lexConv   (char *, real **&, short **&, int &, int &);

#else
#ifdef USE_RAT

#include "ratpoly.h"
#include "pp_solver.h"
extern Rat &atofv    ();
extern Rat &atorv    ();

extern "C" { 
  int lexCode           (multinomial **&mn_list);
 }
extern void rconv    (multinomial *, multinomial *);
extern int lexConv   (char *, Rat **&, short **&, int &, int &);

#else

typedef double real;
#include "pp_solver.h"
#include "simpoly.h"

extern double atorv  ();

extern "C" { 
  int lexCode        (mn_array *&mn_list);
}
extern int lexConv   (char *, real **&, short **&, int &, int &);

#endif
#endif

extern int atop      (int , int*);
extern int get_total (int, int*);
extern int rsort     (int, int*, int*);
extern int sort      (int , int*, int*);

#endif	

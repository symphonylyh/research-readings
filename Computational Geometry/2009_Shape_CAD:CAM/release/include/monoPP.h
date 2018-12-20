// Copyright (C) 1995 Massachusetts Institute of Technology, Cambridge, MA
// All rights reserved

#ifndef MONOPP_H
#define MONOPP_H

#ifndef USE_RAT

#include "simpoly.h"
extern void monoToBern (real **, short **, int, int);
extern void bernToMono (real **, short **, int, int);

#else

#include "ratpoly.h"
extern void rconv      (multinomial *, multinomial *);
extern void monoToBern (Rat **, short **, int, int);

#endif

#endif

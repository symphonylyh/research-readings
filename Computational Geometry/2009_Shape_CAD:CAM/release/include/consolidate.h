// Copyright (C) Massachusetts Institute of Technology, 1995
// All rights reserved

#ifndef CONSOLIDATE_H
#define CONSOLIDATE_H

#ifdef USE_INTERVAL
#include "interval.h"
typedef Interval real;
#else
typedef double real;
#endif

void consolidate (real epsCon, real eps, real **bp, short **ordlists,
		  int numeq, int numvar, real **&roots, int &numroots);

#endif

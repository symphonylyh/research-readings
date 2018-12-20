// Copyright (C) Massachusetts Institute of Technology, 1995
// All rights reserved

#ifndef CONV_H
#define CONV_H

#ifdef USE_INTERVAL
#include "interval.h"
#define CONV_H_REAL Interval
#else
#ifdef USE_RAT
#include "rat.h"
#define CONV_H_REAL Rat
#else
#define CONV_H_REAL double
#endif
#endif

void convOneEq   (CONV_H_REAL *&, short *,  CONV_H_REAL **, int);
void convPolyBox (CONV_H_REAL **, short **, CONV_H_REAL **, int, int);

#endif

/* Copyright (C) Massachusetts Institute of Technology, 1995
   All rights reserved
*/

/* mdist.h
*/

#ifndef MDIST_H
#define MDIST_H

#include <stdio.h>
#include "gen.h"

#ifdef __cplusplus
extern "C" {
#endif

double curvToCurv       (ParCurv *, ParCurv *, double, double *, double *);
double curvToSurf       (ParCurv *, ParSurf *, double, double *, double *,
			 double *);
double distance2        (double, double, double, double, double, double);
int    getNextKnot      (double *, int);
void   mDistCurv        (ParCurv *, ListCurv *, double, double **, int, int,
			 int, double);
void   mDistSurf        (ParSurf *, ListCurv *, double, double **, int, int);
void   minMaxBox        (double *, double *, double *, double *);
double pointToCurv      (double, double, double, ParCurv *, double, double *);
double pointToSurf      (double, double, double, ParSurf *, double, double *,
			 double *);
double surfToSurf       (ParSurf *, ParSurf *, double, double *, double *,
			 double *, double *);
void   updateCurvBox    (vector **, int, double *);
void   updateSurfBox    (vector ***, int, int, double *);
void   writeMinDistFile (FILE *, ListCurv *, double **, int *, int, char **,
			 int, char *, int);

#ifdef __cplusplus
}
#endif

#endif


// Copyright (C) Massachusetts Institute of Technology, 1995
// All rights reserved

// mdist_cc.h

#ifndef MDIST_CC_H
#define MDIST_CC_H

#include "geom.h"
#include "gen.h"
#include "multinom.h"
#include "simpoly.h"

const int MDIST_PC = 1;
const int MDIST_PS = 2;
const int MDIST_CC = 2;
const int MDIST_CS = 3;
const int MDIST_SS = 4;

mn_array *convertToMonomial (geom **);
geom     *dist_bez_cc       (ParCurv *, ParCurv *);
geom     *dist_bez_cs       (ParCurv *, ParSurf *);
geom     *dist_bez_pc       (double, double, double, ParCurv *);
geom     *dist_bez_ps       (double, double, double, ParSurf *);
geom     *dist_bez_ss       (ParSurf *, ParSurf *);
geom    **distObj           (geom *);
void      obj               (geom **);
double    robustMinDistance (double, double, double, ParCurv *, double,
			     double *);
double    robustMinDistance (double, double, double, ParSurf *, double,
			     double *, double *);
double    robustMinDistance (ParCurv *, ParCurv *, double, double *, double *);
double    robustMinDistance (ParCurv *, ParSurf *, double, double *, double *,
			     double *);
double    robustMinDistance (ParSurf *, ParSurf *, double, double *, double *,
			     double *, double *);
rootlist *solveMinDistance  (geom *, double);

extern "C" {
  double  pointToCurv       (double, double, double, ParCurv *, double,
			     double *); 
  void    pointToCurv_n     (double, double, double, ParCurv *, double,
			     double*, double *); 
  double  pointToSurf       (double, double, double, ParSurf *, double,
			     double *, double *);
  void    pointToSurf_n       (double, double, double, ParSurf *, double,
			     double *, double *, double *);
  double curvToCurv         (ParCurv *, ParCurv *, double, double *,
			     double *);
  void   curvToCurv_n         (ParCurv *, ParCurv *, double, double *,
			     double *, double *);
  double curvToSurf         (ParCurv *, ParSurf *, double, double *,
			     double *, double *);
  void   curvToSurf_n         (ParCurv *, ParSurf *, double, double *,
			     double *, double *, double *);
  double surfToSurf         (ParSurf *, ParSurf *, double, double *,
			     double *, double *, double *);
  void   surfToSurf_n         (ParSurf *, ParSurf *, double, double *,
			     double *, double *, double *, double *);
}

#endif

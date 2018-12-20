/*---------------------------------------------------------*/
/* Copyright(c) 1994 Massachusetts Institute of Technology */
/* All rights reserved                                     */
/* by Chun-Yi Hu        July 26 1994                       */
/*---------------------------------------------------------*/
/* This routine define a new class of interval  which has  */
/* more features and its rounding process is more loyal to */
/* the original data than the old one.                     */

#ifndef REFINED_INTERVAL_H
#define REFINED_INTERVAL_H

#include "interval.h"

class refined_interval: public Interval {
 public:
     refined_interval();
     refined_interval (double);
     refined_interval (double, double);
     refined_interval (Interval);
     ~refined_interval ();
     int in_interval(double);
     int contain(double);
};

#endif

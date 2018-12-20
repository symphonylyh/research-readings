
/*---------------------------------------------------------*/
/* Copyright(c) 1994 Massachusetts Institute of Technology */
/* by Chun-Yi Hu        July 26 1994                       */
/*---------------------------------------------------------*/

#include "refined_interval.h"


refined_interval::refined_interval()
: Interval()
{}


refined_interval::refined_interval (double d)
: Interval(d)
{}

refined_interval::refined_interval (double lo, double up)
: Interval(lo,up)
{}

refined_interval::refined_interval (Interval inter)
: Interval(inter)
{}

refined_interval::~refined_interval ()
{}

/* this routine is to decide whether the num is in the interval. */
int refined_interval::in_interval(double num)
{

     if ( get_low() <= num && num <= get_upp())
	  return TRUE;
     else 
	  return FALSE;
}

int refined_interval::contain(double num)
{
     if ( get_low() <= num && num <= get_upp())
	  return TRUE;
     else 
	  return FALSE;

}

/************************************************************************
 *									*
			    Copyright (C) 1992
     Massachusetts Institute of Technology, Cambridge, Massachusetts.
			   All rights reserved.

     This portion of source code was prepared in the Ocean Engineering
     Design Laboratory by Seamus T. Tuohy on 1/25/90.
 *									*
 ************************************************************************/
/*                                                                      *
          routine: normalize()

	  This function normalizes a knot vector such that it
          ranges from 0 to 1. The knots must be positive and in 
	  ascending order. If the last knot in the array is found
          to be zero, or very small ( < 10E-8 ) the routine returns
          a value of negative one (-1), else a value of one (1) is
          returned.

 *                                                                      *
 ************************************************************************
 *                                                                      *
                      INPUT
      NAME       TYPE       DESCRIPTION
      knots     double      array of knots ( i.e. the knot vector )
      num        int        number of knots ( length of array knots )

 ************************************************************************
                      OUTPUT
      NAME       TYPE       DESCRIPTION
      knots     double      upon return the array contains the normalized
                            knot vector.
 *                                                                      *
 ************************************************************************
 *                                                                      *
                         HISTORY
    1/25/90  Hatched
    1/29/90  corrected zeroimg out to start at largest and go to smallest
    1/30/90  transferred from personal library to praxlib
    7/15/92  added to new editor library
 *                                                                      *
 ************************************************************************/

#include <math.h>
#include "editor.h"

short knot_normalize(double *knots, short num)
{
  short i;

  for(i=num-1; i>=0; i--)
    knots[i] = knots[i] - knots[0];    /* zero out the knots */

  if (knots[num-1] < ZERO*0.01)        /* check to see if last knot is */
    return(-1);                        /* small or zero or negative */

  for(i=0; i<num; i++)
    knots[i] = knots[i]/knots[num-1];  /* normalize the knots wrt last knot */

  return(1);
}

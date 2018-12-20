/************************************************************************
 *									*
		Copyright (C) 1989 by Bradley A. Moran
			 All rights reserved.

    The  author hereby grants to  MIT  permission to reproduce and  to
    distribute  copies of  this source  code  document in whole or  in
    part.
 *									*
 ************************************************************************/

# include <malloc.h>
# include <gen.h>

/* Function: avec02b()
 * Purpose: Allocate a vector structure
 * Method: The vector structure is defined in the header file "gen.h".
 *         It holds a single homogeneous coordinate with fields x,y,z,w.
 *         We use the function gen_array1 to allocate dynamic memory
 *         from the system heap.
 *         The vector is normalized with w = 1, in other words, x,y,z
 *         are a cartesian point.
 * Return: The address of the allocated vector structure
 * Note: This function is generally not used. Use vectalloc() instead.
 */
/* Functions referenced by avec02b() are:
 *  gen_array1()
 */

vector *avec02b (void)		/* allocate & initialize vector */

{
     vector *vec;

     vec = (vector *) gen_array1(1, sizeof(vector));

     if (!vec)			/* print error message */
	  perror ("avec02b");
     else
	  vec->w = 1.0;		/* normalize homogeneous coordinate */
     return (vec);		/* could be NULL, check returned value */
}

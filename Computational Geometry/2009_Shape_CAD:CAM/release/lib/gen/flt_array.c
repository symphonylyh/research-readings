/************************************************************************
 *									*
		Copyright (C) 1989 by Bradley A. Moran
			 All rights reserved.

    The  author hereby grants to  MIT  permission to reproduce and  to
    distribute  copies of  this source  code  document in whole or  in
    part.
 *									*
 ************************************************************************/

#include <stdio.h>
# include <malloc.h>
# include <gen.h>

/* Function: flt_array1()
 * Purpose: Allocate a one-dimensional float array
 * Method: Use the general purpose allocator gen_array1 to allocate the
 *         array.
 * Arguments:
 *  nel - the size of the array
 * Return: The address of the allocated array.
 */
/* Functions referenced by flt_array1() are:
 *  errormsg()
 *  gen_array1()
 */
/* Functions that reference flt_array1() are:
 *  DrawTrimCurv()
 *  StoreTrimValues()
 */

float *flt_array1 (unsigned nel)
{
  float *v;
  char line[256];

  v = (float *) gen_array1 (nel, (unsigned) sizeof(float));

  if (!v) {
    sprintf(line, "allocation failure in flt_array1(), requesting %d\n%s",
	    nel, MemoryStatus());
    errormsg(8, line);
  }

  return(v);
}

/* Function: flt_array2()
 * Purpose: Allocate a two-dimensional float array
 * Method: Use the general purpose allocator gen_array2 to allocate the
 *         array.
 * Arguments:
 *  nrow - the number of rows
 *  ncol - the number of columns
 * Return: The address of the allocated array.
 */
/* Functions referenced by flt_array2() are:
 *  errormsg()
 *  gen_array2()
 */
/* Functions that reference flt_array2() are:
 *  DrawSurfCurv()
 *  DrawSurfHeight()
 *  StoreCosValues()
 *  StoreCurvValues()
 *  StoreTrimValues()
 */

float **flt_array2 (unsigned nrow, unsigned ncol)

{
  float **m;
  char line[256];

  m = (float **) gen_array2 (nrow, ncol, (unsigned) sizeof(float));

  if (!m) {
    sprintf(line, "allocation failure in flt_array2(), requesting %dx%d\n%s",
	      nrow, ncol, MemoryStatus());
    errormsg (9, line);
  }

  return(m);
}

/* Function: free_farray1()
 * Purpose: Deallocate one-dimensional float array
 * Method: Use the general purpose memory deallocator free_garray1 to
 *         deallocate the memory.
 * Arguments:
 *  v - the address of the float array
 */
/* Functions referenced by free_farray1() are:
 *  free_garray1()
 */
/* Functions that reference free_farray1() are:
 *  DrawTrimCurv()
 *  StoreTrimValues()
 */

void free_farray1(float *v)
{
  free_garray1((char *)v);
}

/* Function: free_farray1()
 * Purpose: Deallocate two-dimensional float array
 * Method: Use the general purpose memory deallocator free_garray1 to
 *         deallocate the memory.
 * Arguments:
 *  m - the address of the array
 */
/* Functions referenced by free_farray2() are:
 *  free_garray1()
 */
/* Functions that reference free_farray2() are:
 *  DeleteCos()
 *  DeleteCosPts()
 *  DeleteCurv()
 *  DeleteCurvPts()
 *  DeleteTrim()
 *  DeleteTrimPts()
 *  DrawSurfCurv()
 *  DrawSurfHeight()
 *  TrimSurfaceCB()
 */

void free_farray2(float **m)
{
  free_garray1((char *)&m[0][0]);
  free_garray1((char *)m);
}

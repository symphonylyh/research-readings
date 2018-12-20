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
#include <malloc.h>
#include <gen.h>

/* Function: ptr_array1()
 * Purpose: Allocate a one-dimensional array of char pointers
 * Method: Use the general purpose allocator gen_array1() to allocate the
 *         array.
 * Arguments:
 *  nel - the size of the array
 * Return: The address of the allocated array.
 */
/* Functions referenced by ptr_array1() are:
 *  errormsg()
 *  gen_array1()
 */
/* Functions that reference ptr_array1() are:
 *  copyegeom()
 *  reflectegeom()
 *  transformegeom()
 */

char **ptr_array1 (unsigned nel)
{
  char **ptr, line[256];

  ptr = (char **) gen_array1 (nel, (unsigned) sizeof(char *));

  if (!ptr) {
    sprintf(line, "allocation failure in ptr_array1(), requesting %d\n%s",
	    nel, MemoryStatus());
    errormsg (25, line);
  }

  return (ptr);
}

/* Function: ptr_array2()
 * Purpose: Allocate a two-dimensional array of char pointers
 * Method: Use the general purpose allocator gen_array1() to allocate the
 *         array.
 * Arguments:
 *  nrow - the number of rows
 *  ncol - the number of columns
 * Return: The address of the allocated array.
 */
/* Functions referenced by ptr_array2() are:
 *  errormsg()
 *  gen_array2()
 */
/* Functions that reference ptr_array2() are:
 *  copyfgeom()
 *  reflectfgeom()
 *  transformfgeom()
 */

char ***ptr_array2 (unsigned nrow, unsigned ncol)
{
  char ***ptr, line[256];

  ptr = (char ***) gen_array2 (nrow, ncol, (unsigned) sizeof(char *));

  if (!ptr) {
    sprintf(line, "allocation failure 1 in ptr_array2(), requesting %dx%d\n%s",
	      nrow, ncol, MemoryStatus());
    errormsg (22, line);
  }

  return (ptr);
}

/* Function: free_parray1()
 * Purpose: Deallocate one-dimensional array of char pointers
 * Method: Use the general purpose memory deallocator free_garray1() to
 *         deallocate the memory.
 * Arguments:
 *  ptr - the address of the array
 */
/* Functions referenced by free_parray1() are:
 *  free_garray1()
 */

void free_parray1 (char **ptr)	 /* free a 1D array of pointers */
{
  free_garray1((char *)ptr);
}

/* Function: free_parray2()
 * Purpose: Deallocate two-dimensional array of char pointers
 * Method: Use the general purpose memory deallocator free_garray1() to
 *         deallocate the memory.
 * Arguments:
 *  ptr - the address of the array
 */
/* Functions referenced by free_parray2() are:
 *  free_garray1()
 */

void free_parray2 (char ***ptr)	 /* free a 2D array of pointers */
{
  free_garray1((char *)ptr[0]);
  free_garray1((char *)ptr);
}

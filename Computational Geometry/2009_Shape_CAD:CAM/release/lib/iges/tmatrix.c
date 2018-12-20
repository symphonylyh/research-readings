/* Copyright (C) 1993 Massachusetts Institute of Technology, Cambridge, MA
   All rights reserved
*/

/* tmatrix.c */

/* alloc_tmatrix()
 * free_tmatrix()
*/

#include <stdio.h>
#include <malloc.h>
#include "iges.h"

/********* alloc_tmatrix() *********
* 1     Purpose
* 
*       Allocate an IGES Type124 (Transformation Matrix) structure.
* 
* 2     Specification
* 
*       #include "iges.h"
*       Type124 *alloc_tmatrix(void);
* 
* 3     Description
* 
*       This function allocates a data structure to contain the information of
*       an IGES Type124 (Transformation Matrix) entity.
* 
* 4     References
* 
*       [1]   Digital Representation for Communication of Product Definition
*             Data, US PRO/IPO-100, Initial Graphics Exchange Specification
* 	    (IGES) 5.2, IGES/PDES Organization, U.S. Product Data Association,
* 	    Fairfax, VA, November 1993.
* 
* 6     Return Values, Error Indicators and Warnings
* 
*       If insufficient dynamic memory is available, the function returns the
*       value NULL.
* 
* 8     Further Comments
* 
*       Free the structure with free_tmatrix().
* 
* Functions that reference alloc_tmatrix() are:
*  ReadIges124()
*/

Type124 *alloc_tmatrix(void)
{
  Type124 *mat;

  mat = (Type124 *)calloc(1, sizeof(Type124));
  return mat;
}

/********* free_tmatrix() *********
* 1     Purpose
* 
*       Deallocate an IGES Type124 (Transformation Matrix) structure.
* 
* 2     Specification
* 
*       #include "iges.h"
*       void free_tmatrix(Type124 *mat);
* 
* 3     Description
* 
*       This function deallocates an IGES Type124 (Transformation Matrix)
*       structure.
* 
* 4     References
* 
*       [1]   Digital Representation for Communication of Product Definition
*             Data, US PRO/IPO-100, Initial Graphics Exchange Specification
* 	    (IGES) 5.2, IGES/PDES Organization, U.S. Product Data Association,
* 	    Fairfax, VA, November 1993.
* 
* 5     Parameters
* 
*           1.  Type124 * mat
*               On entry:  the address of the structure.
* 
* 8     Further Comments
* 
*       The structure should be allocated with alloc_tmatrix().
* 
* Functions that reference free_tmatrix() are:
*  OpenIgesFileCB()
*/

void free_tmatrix(Type124 *mat)
{
  free(mat);
}

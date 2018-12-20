/* Copyright (C) 1992 Massachusetts Institute of Technology, Cambridge, MA
   All rights reserved
*/

/* copious.c */

/* alloc_copious()
 * free_copious()
*/

#include <malloc.h>
#include "iges.h"

/********* alloc_copious() *********
* 1     Purpose
* 
*       Allocate an IGES Type106 (Copious Data) structure.
* 
* 2     Specification
* 
*       #include "iges.h"
*       Type106 *allocate_copious(int ip, int n);
* 
* 3     Description
* 
*       This function allocates a data structure to contain the information of
*       an IGES Type106 (Copious Data) entity.
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
*           1.  int ip
*               On entry:  flag specifying the interpretation of the coordinate
* 	      data.
* 
*                ip    Interpretation
*                ---   ---------------------------------------
*                 1    Coordinate pairs (x,y)
*                 2    Coordinate triples (x,y,z)
*                 3    Coordinate sextuples (x,y,z,i,j,k)
* 
*           2.  int n
*               On entry:  the number of data coordinates contained in the
* 	      structure.
* 
* 6     Return Values, Error Indicators and Warnings
* 
*       If insufficient dynamic memory is available, the function returns the
*       value NULL.
* 
* 8     Further Comments
* 
*       Free the data stucture with free_copious().
* 
* Functions referenced by alloc_copious() are:
*  dbl_array1()
*  gen_array1()
*
* Functions that reference alloc_copious() are:
*  ReadIges106()
*  SaveIgesCos()
*  SaveIgesCurv()
*  SaveIgesGrid()
*  SaveIgesList()
*  SaveIgesMinDist()
*  SaveIgesSurf()
*  SaveIgesUv()
*/

Type106 *alloc_copious(int ip, int n)
{
  Type106 *copious;

  if (copious = (Type106 *)gen_array1(1, sizeof(Type106))) {
    copious->n = n;
    switch (ip) {
    case 1:    /* form 1 and 11: (x,y) pairs */
    case 11:
      copious->ip = 1;
      copious->pts.pair.x = dbl_array1(copious->n);
      copious->pts.pair.y = dbl_array1(copious->n);
      break;
    case 2:    /* form 2 and 12: (x,y,z) triples */
    case 12:
      copious->ip = 2;
      copious->pts.triple.x = dbl_array1(copious->n);
      copious->pts.triple.y = dbl_array1(copious->n);
      copious->pts.triple.z = dbl_array1(copious->n);
      break;
    case 3:    /* form 3 and 13: (x,y,z,i,j,k) sextuples */
    case 13:
      copious->ip = 3;
      copious->pts.sextuple.x = dbl_array1(copious->n);
      copious->pts.sextuple.y = dbl_array1(copious->n);
      copious->pts.sextuple.z = dbl_array1(copious->n);
      copious->pts.sextuple.i = dbl_array1(copious->n);
      copious->pts.sextuple.j = dbl_array1(copious->n);
      copious->pts.sextuple.k = dbl_array1(copious->n);
      break;
    }
  }
  return (copious);
}

/********* free_copious() *********
* 1     Purpose
* 
*       Deallocate IGES Type106 (Copious Data) structure.
* 
* 2     Specification
* 
*       #include "iges.h"
*       void free_copious(Type106 *copious);
* 
* 3     Description
* 
*       This function deallocates an IGES Type106 (Copious Data) structure.
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
*           1.  Type106 * copious
*               On entry:  the address of the structure.
* 
* 8     Further Comments
* 
*       The structure should be allocated with alloc_copious().
* 
* Functions referenced by free_copious() are:
*  free_darray1()
*  free_garray1()
*
* Functions that reference free_copious() are:
*  FreeIgesDirList()
*  ReadIgesList()
*  ReadIgesUv()
*  SaveIgesCos()
*  SaveIgesCurv()
*  SaveIgesGrid()
*  SaveIgesList()
*  SaveIgesMinDist()
*  SaveIgesSurf()
*  SaveIgesUv()
*/

void free_copious(Type106 *copious)
{
  switch (copious->ip) {
  case 1:    /* form 1 and 11: (x,y) pairs */
  case 11:
    free_darray1(copious->pts.pair.x);
    free_darray1(copious->pts.pair.y);
    break;
  case 2:    /* form 2 and 12: (x,y,z) triples */
  case 12:
    free_darray1(copious->pts.triple.x);
    free_darray1(copious->pts.triple.y);
    free_darray1(copious->pts.triple.y);
    break;
  case 3:    /* form 3 and 13: (x,y,z,i,j,k) sextuples */
  case 13:
    free_darray1(copious->pts.sextuple.x);
    free_darray1(copious->pts.sextuple.y);
    free_darray1(copious->pts.sextuple.y);
    free_darray1(copious->pts.sextuple.i);
    free_darray1(copious->pts.sextuple.j);
    free_darray1(copious->pts.sextuple.k);
    break;
  }
  free_garray1((char *)copious);
}

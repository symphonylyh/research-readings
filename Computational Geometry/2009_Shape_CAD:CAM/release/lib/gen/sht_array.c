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

/* Function: sht_array1()
 * Purpose: Allocate a one-dimensional short array
 * Method: Use the general purpose allocator gen_array1() to allocate the
 *         array.
 * Arguments:
 *  nel - the size of the array
 * Return: The address of the allocated array.
 */
/* Functions referenced by sht_array1() are:
 *  errormsg()
 *  gen_array1()
 */
/* Functions that reference sht_array1() are:
 *  AlphaFacetCB()
 *  camber_surface()
 *  ConstrainedCB()
 *  ContourFacets()
 *  ConvertToCos()
 *  CopyFacet()
 *  CopyMinDist()
 *  Delaunay()
 *  EvalFacet()
 *  InitFullFacet()
 *  IsAbove()
 *  main()
 *  MakeTriangle()
 *  MeasuredLineCB()
 *  MeasuredPoints()
 *  MinDistance()
 *  MinDistCurvCB()
 *  NewGroupCB()
 *  OrthoDistance()
 *  ProbeRadiusCB()
 *  PruneFacets()
 *  ReadDeslabFacet()
 *  ReadDeslabHull()
 *  ReadDeslabLocal()
 *  ReadDeslabMinDist()
 *  ReflectFacet()
 *  RobustCurvCB()
 *  RobustDistCB()
 *  SampleFacetGridCB()
 *  SubDistance()
 *  toFacet()
 *  TraceContour()
 *  TransformFacet()
 *  TriangulateStructured()
 *  TriangulateUnstructured()
 *  TrimToFacetCB()
 *  Unconloc()
 *  Unconstrained2dCB()
 */

short *sht_array1 (unsigned nel)
{
  short *v;
  char line[256];

  v = (short *) gen_array1 (nel, (unsigned) sizeof(short));

  if (!v) {
    sprintf(line, "allocation failure in sht_array1(), requesting %d\n%s",
	    nel, MemoryStatus());
    errormsg (5, line);
  }

  return (v);
}

/* Function: sht_array2()
 * Purpose: Allocate a two-dimensional short array
 * Method: Use the general purpose allocator gen_array1() to allocate the
 *         array.
 * Arguments:
 *  nrow - the number of rows
 *  ncol - the number of columns
 * Return: The address of the allocated array.
 */
/* Functions referenced by sht_array2() are:
 *  errormsg()
 *  gen_array2()
 */
/* Functions that reference sht_array2() are:
 *  AlphaFacetCB()
 *  AlphaShape()
 *  ContourFacets()
 *  ReadDeslabHull()
 */

short **sht_array2 (unsigned nrow, unsigned ncol)
{
  short **m;
  char line[256];

  m = (short **) gen_array2 (nrow, ncol, (unsigned) sizeof(short));

  if (!m) {
    sprintf(line, "allocation failure in sht_array2(), requesting %dx%d\n%s",
	      nrow, ncol, MemoryStatus());
    errormsg (6, line);
  }

  return (m);
}

/* Function: free_sarray1()
 * Purpose: Deallocate one-dimensional short array
 * Method: Use the general purpose memory deallocator free_garray1() to
 *         deallocate the memory.
 * Arguments:
 *  v - the address of the array
 */
/* Functions referenced by free_sarray1() are:
 *  free_garray1()
 */
/* Functions that reference free_sarray1() are:
 *  AlphaFacetCB()
 *  camber_surface()
 *  ContourFacets()
 *  ContourFacetsCB()
 *  ConvertToCos()
 *  DeleteFacet()
 *  DeleteHull()
 *  DeleteList()
 *  DeleteTrim()
 *  DeleteTrimPts()
 *  EvalFacetOkCB()
 *  InitFullFacet()
 *  IsAbove()
 *  main()
 *  MakeTriangle()
 *  MeasuredLineCB()
 *  MeasuredPoints()
 *  NewGroupCB()
 *  OrthoDistance()
 *  PruneFacets()
 *  SampleFacetGridCB()
 *  TraceContour()
 *  TriangulateStructured()
 *  TrimSurfaceCB()
 */

void free_sarray1 (short *v)
{
  free_garray1((char *)v);
}

/* Function: free_sarray2()
 * Purpose: Deallocate two-dimensional short array
 * Method: Use the general purpose memory deallocator free_garray1() to
 *         deallocate the memory.
 * Arguments:
 *  m - the address of the array
 */
/* Functions referenced by free_sarray2() are:
 *  free_garray1()
 */
/* Functions that reference free_sarray2() are:
 *  AlphaFacetCB()
 *  AlphaShape()
 *  ContourFacetsCB()
 *  DeleteHull()
 */

void free_sarray2 (short **m)
{
  free_garray1((char *)&m[0][0]);
  free_garray1((char *)m);
}

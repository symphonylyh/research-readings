/************************************************************************
 *									*
		Copyright (C) 1989 by Bradley A. Moran
			 All rights reserved.

    The  author hereby grants to  MIT  permission to reproduce and  to
    distribute  copies of  this source  code  document in whole or  in
    part.
 *									*
 ************************************************************************/

#include <malloc.h>
#include "gen.h"

/* Function: gen_array1()
 * Purpose: Allocate a one-dimensional array with arbitrary sized elements
 * Method: Use the standard library function calloc() to allocate the
 *         array from the system heap.
 * Arguments:
 *  nel - the size of the array
 *  elsize - the size of each element of the array in bytes
 * Return: The address of the allocated array.
 */
/* Functions that reference gen_array1() are:
 *  AllocIgesDirectory()
 *  AllocIgesGlobal()
 *  alloc_copious()
 *  avec02b()
 *  bit_array1()
 *  BlendSurfCB()
 *  BsplineToBezier()
 *  camber_surface()
 *  CheckCosContinuity()
 *  CheckCurvContinuity()
 *  CheckSurfContinuity()
 *  ContourFacets()
 *  ContourFacetsCB()
 *  copyegeom()
 *  copyegeom_leo()
 *  copyfgeom()
 *  CopyGridSurf()
 *  CopyListCurv()
 *  CopyParUv()
 *  CopyTrimSurf()
 *  dbl_array1()
 *  dbl_array3()
 *  DecodeIgesDirectory()
 *  DecodeIgesGlobal()
 *  DecodeIgesTerminate()
 *  egeomalloc1()
 *  egeomalloc2()
 *  fgeomalloc1()
 *  fgeomalloc2()
 *  flt_array1()
 *  GenCylCB()
 *  GetIgesDirList()
 *  InitializeStore()
 *  InsertTriangle()
 *  int_array1()
 *  InvertTrimSurf()
 *  JoinCoss()
 *  JoinCurvs()
 *  JoinSurfs()
 *  LoftSurfCB()
 *  main()
 *  merge_can_off()
 *  OpenIgesFileCB()
 *  ParseIgesRecord()
 *  ParSurf_approx()
 *  ParSurf_to_SurfRev()
 *  ParSurf_to_TabCyl()
 *  pgeomalloc()
 *  PowSurf_to_ParSurf_loft()
 *  ptr_array1()
 *  ReadGridSurf()
 *  ReadIges110()
 *  ReadIges118()
 *  ReadIges120()
 *  ReadIges122()
 *  ReadIges142()
 *  ReadIges144()
 *  ReadListCurv()
 *  ReadParCurv_Per()
 *  ReadParSurf_Peru()
 *  ReadParUv()
 *  ReadTrimSurf()
 *  reflectegeom()
 *  reflectfgeom()
 *  ReflectGridSurf()
 *  ReflectListCurv()
 *  ReflectTrimSurf()
 *  RobustCurvCB()
 *  RobustDistCB()
 *  sgeomalloc()
 *  sht_array1()
 *  SurfSplitAddUCB()
 *  SurfSplitAddVCB()
 *  SurfSplitRemoveUCB()
 *  SurfSplitRemoveVCB()
 *  transformegeom()
 *  transformfgeom()
 *  TransformGridSurf()
 *  TransformListCurv()
 *  TransformTrimSurf()
 *  vectalloc()
 *  vec_array1()
 *  vec_array2()
 */

char *gen_array1 (unsigned nel, unsigned elsize)
{
  char *v;
  struct mallinfo mi;

  v = calloc(nel, elsize);

  return v;
}

/* Function: gen_array2()
 * Purpose: Allocate a two-dimensional array with arbitrary sized elements
 * Method: Use the standard library function calloc() to allocate the
 *         array from the system heap.
 * Arguments:
 *  nrow - the number of rows
 *  ncol - the number of columns
 *  elsize - the size of each element of the array in bytes
 * Return: The address of the allocated array.
 */
/* Functions that reference gen_array2() are:
 *  BsplineToBezierSurf()
 *  CheckCosContinuity()
 *  CheckCurvContinuity()
 *  CheckSurfContinuity()
 *  compare_patches()
 *  dbl_array2()
 *  DrawSurfOffset2()
 *  DrawSurfTwoSided()
 *  fgeom_array2()
 *  flt_array2()
 *  int_array2()
 *  JoinCoss()
 *  JoinCurvs()
 *  JoinSurfs()
 *  ptr_array2()
 *  sht_array2()
 *  SplitSurface()
 *  SubDistance()
 *  tolerance_region()
 *  vec_array2()
 */

char **gen_array2 (unsigned nrow, unsigned ncol, unsigned elsize)
{
  register char **mptr;		/* matrix pointer */
  struct mallinfo mi;

				/* allocate the row-pointers first */
  mptr = (char **) malloc (nrow * sizeof(char *));

  if (mptr) {			/* error check on initial allocation */

    				/* try to allocate memory block */
    *mptr = (char *) calloc (nrow * ncol, elsize);

    if (*mptr) {		/* error check on large memory block */

      register int i;
				/* perform initialization */
      for (i = 1; i < nrow; ++i)
	mptr[i] = (*mptr + i * ncol * elsize);

    } else {			/* large block failed, free row */
				/* pointers to prevent net growth */
      free ((char *) mptr);
      mptr = (char **) 0;	/* reset return value to nil for error */
				/* signal to calling routine */
    }
  }

  return (mptr);
}

/* Function: free_garray1()
 * Purpose: Deallocate one-dimensional array
 * Method: Use the standard library function free() to deallocate the memory.
 * Arguments:
 *  v - the address of the array
 */
/* Functions that reference free_garray1() are:
 *  BsplineToBezierCurvCB()
 *  CheckCosContinuity()
 *  CheckCurvContinuity()
 *  CheckSurfContinuity()
 *  ContourFacetsCB()
 *  DeleteSurf()
 *  FreeIgesDirectory()
 *  FreeIgesGlobal()
 *  free_barray1()
 *  free_copious()
 *  free_darray1()
 *  free_darray2()
 *  free_darray3()
 *  free_egeom()
 *  free_egeom_array1()
 *  free_farray1()
 *  free_farray2()
 *  free_fgeom()
 *  free_iarray1()
 *  free_iarray2()
 *  free_parray1()
 *  free_parray2()
 *  free_pgeom()
 *  free_sarray1()
 *  free_sarray2()
 *  free_sgeom()
 *  free_varray1()
 *  free_varray2()
 *  GetNextInteger()
 *  GetNextReal()
 *  JoinCoss()
 *  JoinCurvs()
 *  JoinSurfs()
 *  loop_rayintersect()
 *  main()
 *  mDistCurv()
 *  OpenIgesFileCB()
 *  RobustCurvCB()
 *  SplitSurfCancelCB()
 *  SplitSurfOkCB()
 *  SurfSplitRemoveUCB()
 *  SurfSplitRemoveVCB()
 *  vectfree()
 */

void free_garray1(char *v)
{
  free(v);
}

/* Function: free_garray2()
 * Purpose: Deallocate two-dimensional array
 * Method: Use the standard library function free() to deallocate the memory.
 * Arguments:
 *  m - the address of the array
 */
/* Functions that reference free_garray2() are:
 *  BsplineToBezierSurfCB()
 *  CheckCosContinuity()
 *  CheckCurvContinuity()
 *  CheckSurfContinuity()
 *  DrawSurfOffset2()
 *  DrawSurfTwoSided()
 *  free_fgeom_array2()
 *  JoinCoss()
 *  JoinCurvs()
 *  JoinSurfs()
 *  mDistSurf()
 *  SplitSurface()
 *  SubDistance()
 */

void free_garray2(char **a)
{
  free(&a[0][0]);
  free(a);
}

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

/* Function: vec_array1()
 * Purpose: Allocate a one-dimensional vector array
 * Method: Use the general purpose allocator gen_array1() to allocate the
 *         array.
 * Arguments:
 *  nel - the size of the array
 * Return: The address of the allocated array.
 */
/* Functions referenced by vec_array1() are:
 *  errormsg()
 *  gen_array1()
 */
/* Functions that reference vec_array1() are:
 *  boehm_curve()
 *  boehm_surface()
 *  CalcIsophotes()
 *  CalcReflectionLines()
 *  CamberBisect()
 *  CamberBisect3d()
 *  ConstrainedCB()
 *  ConvertToCos()
 *  cyl_frenet_tr()
 *  DrawIsophotes()
 *  DrawReflectionLines()
 *  egeomalloc1()
 *  fcn()
 *  frenet_tr()
 *  raise_byone()
 *  ReadIgesLine()
 *  ReadParCurv()
 *  ReadParCurv_Per()
 *  reflectegeom()
 *  solve_bc()
 *  tol_sample_offsets()
 *  Unconloc()
 */

vector **vec_array1 (unsigned nel)
{
  register int i;
  register vector **v;
  char line[256];

  v = (vector **) gen_array1 (nel, (unsigned) sizeof(vector *));

  if (!v) {
    sprintf(line, "allocation failure 1 in vec_array1(), requesting %d\n%s",
	    nel, MemoryStatus());
    errormsg (14, line);
  }

  for (i = 0; i < nel; i++) {
    v[i] = (vector *) gen_array1 (1, sizeof(vector));

    if (!v[i]) {
      sprintf(line, "allocation failure 2 in vec_array1(), requesting %d\n%s",
	      nel, MemoryStatus());
      errormsg (14, line);
    }

    v[i]->w = 1.0;
  }
  return (v);
}

/* Function: vec_array2()
 * Purpose: Allocate a two-dimensional vector array
 * Method: Use the general purpose allocator gen_array1() to allocate the
 *         array.
 * Arguments:
 *  nrow - the number of rows
 *  ncol - the number of columns
 * Return: The address of the allocated array.
 */
/* Functions referenced by vec_array2() are:
 *  errormsg()
 *  gen_array1()
 *  gen_array2()
 */
/* Functions that reference vec_array2() are:
 *  BezierContpts()
 *  copyfgeom()
 *  decasteljau_curve_contpts()
 *  fgeomalloc1()
 *  MinDistance()
 *  ParCurv_to_PowCurv()
 *  pgeomalloc()
 *  ReadParSurf()
 *  ReadParSurf_Peru()
 *  reflectfgeom()
 *  sgeomalloc()
 *  SubDistance()
 *  transformfgeom()
 */

vector ***vec_array2 (unsigned nrow, unsigned ncol)
{
  register vector ***m;
  register int    i, j;
  char line[256];

  m = (vector ***) gen_array2 (nrow, ncol,
			       (unsigned) sizeof(vector **));

  if (!m) {
    sprintf(line, "allocation failure 1 in vec_array2(), requesting %dx%d\n%s",
	      nrow, ncol, MemoryStatus());
    errormsg (15, line);
  }

  for (i = 0; i < nrow; i++)
    for (j = 0; j < ncol; j++) {

      m[i][j] = (vector *) gen_array1 (1, sizeof(vector));

      if (!m[i][j]) {
	sprintf(line, "allocation failure 3 in vec_array2(), requesting %dx%d\n%s",
		  nrow, ncol, MemoryStatus());
	errormsg (16, line);
      }

      m[i][j]->w = 1.0;
    }
  return (m);
}

/* Function: free_varray1()
 * Purpose: Deallocate one-dimensional vector array
 * Method: Deallocate each vector using vectfree() and then use the
 *         general purpose memory deallocator free_garray1() to
 *         deallocate the memory.
 * Arguments:
 *  m - the address of the array
 *  nel - the size of the array
 */
/* Functions referenced by free_varray1() are:
 *  free_garray1()
 *  vectfree()
 */
/* Functions that reference free_varray1() are:
 *  approx_fnbc()
 *  approx_fn_per()
 *  boehm_curve()
 *  boehm_surface()
 *  CalcIsophotes()
 *  CalcReflectionLines()
 *  CamberBisect()
 *  CamberBisect3d()
 *  ConstrainedCB()
 *  copyegeom()
 *  cross_link()
 *  cyl_generatrix()
 *  DrawIsophotes()
 *  DrawReflectionLines()
 *  fcn()
 *  fit_curve()
 *  free_egeom()
 *  generatrix()
 *  raise_byone()
 *  ReadIgesLine()
 *  ReadParCurv()
 *  ReadParCurv_Per()
 *  reflectegeom()
 *  solve_bc()
 *  tol_sample_offsets()
 *  transformegeom()
 *  Unconloc()
 */

void free_varray1 (vector **v, unsigned nel)		
{
  register int i;
  for (i = nel - 1; i >= 0; i--)
    vectfree(v[i]);
  free_garray1((char *)v);
}

/* Function: free_varray2()
 * Purpose: Deallocate two-dimensional vector array
 * Method: Deallocate each vector using vectfree() and then use the
 *         general purpose memory deallocator free_garray1() to
 *         deallocate the memory.
 * Arguments:
 *  m - the address of the array
 *  nrow - the number of rows
 *  ncol - the number of columns
 */
/* Functions referenced by free_varray2() are:
 *  free_garray1()
 *  vectfree()
 */
/* Functions that reference free_varray2() are:
 *  copyfgeom()
 *  decasteljau_curve_contpts()
 *  FreeBezierContpts()
 *  free_fgeom()
 *  free_pgeom()
 *  free_sgeom()
 *  MinDistance()
 *  ParCurv_to_PowCurv()
 *  ReadParSurf()
 *  ReadParSurf_Peru()
 *  reflectfgeom()
 *  SubDistance()
 *  transformfgeom()
 */

void free_varray2 (vector ***m, unsigned nrow, unsigned ncol)
{
  register int i, j;
  for (i = nrow - 1; i >= 0; i--)
  for (j = ncol - 1; j >= 0; j--)
    vectfree(m[i][j]);
  free_garray1((char *)m[0]);
  free_garray1((char *)m);
}

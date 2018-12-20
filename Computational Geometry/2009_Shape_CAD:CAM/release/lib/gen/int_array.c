/************************************************************************
 *									*
		Copyright (C) 1989 by Bradley A. Moran
			 All rights reserved.

    The  author hereby grants to  MIT  permission to reproduce and  to
    distribute  copies of  this source  code  document in whole or  in
    part.
 *									*
 ************************************************************************/

#include <gen.h>
#include <malloc.h>
#include <stdio.h>

/* Function: int_array1()
 * Purpose: Allocate a one-dimensional int array
 * Method: Use the general purpose allocator gen_array1() to allocate the
 *         array.
 * Arguments:
 *  nel - the size of the array
 * Return: The address of the allocated array.
 */
/* Functions referenced by int_array1() are:
 *  errormsg()
 *  gen_array1()
 */
/* Functions that reference int_array1() are:
 *  addsections_canal()
 *  approx_cubic_surf()
 *  build_offset()
 *  build_offsets_int()
 *  calc_gram_percurv()
 *  camber_surface()
 *  CheckBox()
 *  CheckCosContinuity()
 *  CheckCurvContinuity()
 *  CheckSurfContinuity()
 *  ConstrainedCB()
 *  ConvexHullSurf()
 *  convexhull_test()
 *  con_localize()
 *  CopyMinDist()
 *  curve_fit()
 *  cyl_sample_gencyl()
 *  Delaunay()
 *  evalbsp_per()
 *  evalderivbsp_per()
 *  evalderivsurf_per()
 *  evalrsurf_per()
 *  evalsurf_per()
 *  FindUniqueKnots()
 *  get_error()
 *  get_error_cv()
 *  integral_build_offset()
 *  int_sample_gencyl()
 *  JoinCoss()
 *  JoinCurvs()
 *  JoinSurfs()
 *  knot_rem_and_pert()
 *  knot_rem_cv()
 *  LoftSurfCB()
 *  loft_integral()
 *  loft_rational()
 *  main()
 *  MeasuredPoints()
 *  merge_can_off()
 *  merge_tol_edges()
 *  MinDistance()
 *  MinDistCurvCB()
 *  NewGroupCB()
 *  n_matrix_per()
 *  OpenIgesFileCB()
 *  OrthoDistance()
 *  ParSurf_approx()
 *  perturb_knots_cv()
 *  ProbeRadiusCB()
 *  ReadDeslabLocal()
 *  ReadDeslabMinDist()
 *  RobustCurvCB()
 *  RobustDistCB()
 *  sample_csurf()
 *  sample_gencyl()
 *  SaveIgesTrim()
 *  scattered_fit()
 *  solve_gram()
 *  solve_gram_vector()
 *  soslo()
 *  StoreTrimValues()
 *  SubDistance()
 *  toFacet()
 *  tolerance_region()
 *  tol_loft_rational()
 *  Unconloc()
 *  Unconstrained2dCB()
 */

int *int_array1 (unsigned nel)
{
  int *v;
  char line[256];

  v = (int *) gen_array1 (nel, (unsigned) sizeof(int));

  if (!v) {
    sprintf(line, "allocation failure in int_array1(), requesting %d\n%s",
	    nel, MemoryStatus());
    errormsg (5, line);
  }

  return (v);
}

/* Function: int_array2()
 * Purpose: Allocate a two-dimensional int array
 * Method: Use the general purpose allocator gen_array1() to allocate the
 *         array.
 * Arguments:
 *  nrow - the number of rows
 *  ncol - the number of columns
 * Return: The address of the allocated array.
 */
/* Functions referenced by int_array2() are:
 *  errormsg()
 *  gen_array2()
 */
/* Functions that reference int_array2() are:
 *  addsections_canal()
 *  integral_sample_offset()
 *  int_sample_gencyl()
 *  monotobern()
 *  sample_offset()
 *  tol_sample_offsets()
 */

int **int_array2 (unsigned nrow, unsigned ncol)
{
  int **m;
  char line[256];

  m = (int **) gen_array2 (nrow, ncol, (unsigned) sizeof(int));

  if (!m) {
    sprintf(line, "allocation failure in int_array2(), requesting %dx%d\n%s",
	    nrow, ncol, MemoryStatus());
    errormsg (6, line);
  }

  return(m);
}

/* Function: free_iarray1()
 * Purpose: Deallocate one-dimensional int array
 * Method: Use the general purpose memory deallocator free_garray1() to
 *         deallocate the memory.
 * Arguments:
 *  v - the address of the array
 */
/* Functions referenced by free_iarray1() are:
 *  free_garray1()
 */
/* Functions that reference free_iarray1() are:
 *  addsections_canal()
 *  approx_cubic_surf()
 *  build_offset()
 *  build_offsets_int()
 *  calc_gram_percurv()
 *  camber_surface()
 *  CheckBox()
 *  CheckCosContinuity()
 *  CheckCurvContinuity()
 *  CheckSurfContinuity()
 *  ConstrainedCB()
 *  ConvexHullSurf()
 *  convexhull_test()
 *  con_localize()
 *  curve_fit()
 *  cyl_sample_gencyl()
 *  Delaunay()
 *  DeleteList()
 *  evalbsp_per()
 *  evalderivbsp_per()
 *  evalderivsurf_per()
 *  evalrsurf_per()
 *  evalsurf_per()
 *  FindUniqueKnots()
 *  get_error()
 *  get_error_cv()
 *  integral_build_offset()
 *  int_sample_gencyl()
 *  JoinCoss()
 *  JoinCurvs()
 *  JoinSurfs()
 *  knot_rem_and_pert()
 *  knot_rem_cv()
 *  LoftSurfCB()
 *  loft_integral()
 *  loft_rational()
 *  main()
 *  MeasuredPoints()
 *  merge_can_off()
 *  merge_tol_edges()
 *  MinDistance()
 *  MinDistCurvCB()
 *  NewGroupCB()
 *  n_matrix_per()
 *  OpenIgesFileCB()
 *  OrthoDistance()
 *  ParSurf_approx()
 *  perturb_knots_cv()
 *  RobustCurvCB()
 *  RobustDistCB()
 *  sample_csurf()
 *  sample_gencyl()
 *  SaveIgesTrim()
 *  scattered_fit()
 *  solve_gram()
 *  solve_gram_vector()
 *  soslo()
 *  StoreTrimValues()
 *  SubDistance()
 *  toFacet()
 *  tolerance_region()
 *  tol_loft_rational()
 *  Unconloc()
 *  Unconstrained2dCB()
 */

void free_iarray1 (int *v)
{
  free_garray1((char *)v);
}

/* Function: free_iarray2()
 * Purpose: Deallocate two-dimensional int array
 * Method: Use the general purpose memory deallocator free_garray1() to
 *         deallocate the memory.
 * Arguments:
 *  m - the address of the array
 */
/* Functions referenced by free_iarray2() are:
 *  free_garray1()
 */
/* Functions that reference free_iarray2() are:
 *  addsections_canal()
 *  integral_sample_offset()
 *  int_sample_gencyl()
 *  monotobern()
 *  sample_offset()
 *  tol_sample_offsets()
 */

void free_iarray2(int **m)
{
  free_garray1((char *)&m[0][0]);
  free_garray1((char *)m);
}

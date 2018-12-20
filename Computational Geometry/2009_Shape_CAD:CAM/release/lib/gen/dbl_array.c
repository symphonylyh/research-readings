/************************************************************************
 *									*
		Copyright (C) 1989 by Bradley A. Moran
			 All rights reserved.

    The  author hereby grants to  MIT  permission to reproduce and  to
    distribute  copies of  this source  code  document in whole or  in
    part.
 *									*
 ************************************************************************/

#include <string.h>
# include <malloc.h>
# include <gen.h>

/* Function: dbl_array1()
 * Purpose: Allocate a one-dimensional double array
 * Method: Use the general purpose allocator gen_array1 to allocate the
 *         array.
 * Arguments:
 *  nel - the size of the array
 * Return: The address of the allocated array.
 */
/* Functions referenced by dbl_array1() are:
 *  errormsg()
 *  gen_array1()
 */
/* Functions that reference dbl_array1() are:
 *  alloc_copious()
 *  ApproxCurvCB()
 *  approx_cubic_curv()
 *  approx_cubic_curv1()
 *  approx_cubic_surf()
 *  approx_fnbc()
 *  approx_fnbc_knot()
 *  approx_fn_per()
 *  approx_user_knots()
 *  avg_arc_length()
 *  bernrs1()
 *  BlendSurfCB()
 *  boehm_curve()
 *  boehm_surface()
 *  bone_intersect()
 *  bone_intersect3()
 *  BsplineToBezier()
 *  build_offset()
 *  build_offsets_int()
 *  calc_disc()
 *  calc_disc_curv()
 *  calc_hartley()
 *  calc_hartley_per()
 *  calc_interval_egeom()
 *  calc_knots_foley()
 *  calc_knots_hartley()
 *  calc_knots_per()
 *  calc_par_foley()
 *  calc_par_hartley()
 *  CamberBisect()
 *  CamberBisect3d()
 *  camber_line_2d()
 *  camber_line_3d()
 *  camber_surface()
 *  CheckCosContinuity()
 *  CheckCurvContinuity()
 *  CheckSurfContinuity()
 *  check_one_dir()
 *  ConstrainedCB()
 *  ConvertToCos()
 *  convexbox()
 *  ConvexHullSurf()
 *  ConvexHullTest()
 *  convexhull_test()
 *  con_localize()
 *  copyegeom()
 *  copyfgeom()
 *  CopyMinDist()
 *  CopyUniqueKnots()
 *  critical_bone()
 *  critical_bone3()
 *  CurveOrientation()
 *  curve_fit()
 *  curve_osloper_lt()
 *  curve_osloper_rt()
 *  curve_oslo_per()
 *  Delaunay()
 *  del_dbl_knots()
 *  deriv_sign()
 *  EditKnotsCos()
 *  EditKnotsCurv()
 *  egeomalloc1()
 *  egeomalloc2()
 *  fgeomalloc1()
 *  fgeomalloc2()
 *  findroot()
 *  FindUniqueKnots()
 *  FitCosCB()
 *  FitCurvCB()
 *  FitSurfCB()
 *  fit_curve()
 *  fn_camber()
 *  fn_camber3()
 *  fn_thckns()
 *  fn_thckns3()
 *  fundamentals()
 *  GenCylCB()
 *  global_pos_error()
 *  integral_build_offset()
 *  integral_sample_offset()
 *  interp_ortho()
 *  inter_curv()
 *  inter_percurv()
 *  inter_surf()
 *  int_curve()
 *  int_number()
 *  JoinCoss()
 *  JoinCurvs()
 *  JoinSurfs()
 *  knot_factor0()
 *  knot_factor1()
 *  LoftSurfCB()
 *  loft_integral()
 *  main()
 *  MeasuredPoints()
 *  MergeKnots()
 *  merge_can_off()
 *  MinDistance()
 *  MinDistCurvCB()
 *  monotobern()
 *  nbasisd()
 *  nbasisd_per()
 *  offset_curve()
 *  OrthoDistance()
 *  OrthoProjCosCB()
 *  OrthoProjCurvCB()
 *  ParCurv_to_PowCurv()
 *  ParSurf_approx()
 *  partan_dir()
 *  pgeomalloc()
 *  PowCurv_to_ParCurv()
 *  PowSurf_to_ParSurf_int()
 *  PowSurf_to_ParSurf_loft()
 *  ProbeRadiusCB()
 *  project_curve_onsurf()
 *  project_curve_on_surf()
 *  raise_surf()
 *  ReadDeslabLocal()
 *  ReadDeslabMinDist()
 *  ReadParCurv()
 *  ReadParCurv_Per()
 *  ReadParSurf()
 *  ReadParSurf_Peru()
 *  reflectegeom()
 *  reflectfgeom()
 *  renew_knots()
 *  rk4()
 *  RobustCurvCB()
 *  RobustDistCB()
 *  rootisolate()
 *  sample()
 *  sample_blend()
 *  sample_csurf()
 *  sample_mod_surf()
 *  sample_offset()
 *  sample_per()
 *  sample_proj()
 *  sample_surf_test()
 *  ScatteredSurfCB()
 *  scattered_fit()
 *  sgeomalloc()
 *  solve_gram_percurv()
 *  solve_gram_pervector()
 *  soslo()
 *  splitpoly2()
 *  StoreCosValues()
 *  StoreCurvValues()
 *  StoreTrimValues()
 *  SubDistance()
 *  SubdivideCos()
 *  SubdivideCurv()
 *  SurfSplitAddUCB()
 *  SurfSplitAddVCB()
 *  SurfSplitRemoveUCB()
 *  SurfSplitRemoveVCB()
 *  tol_sample_offsets()
 *  transformegeom()
 *  transformfgeom()
 *  trim_start_par()
 *  trim_start_par3()
 *  Unconloc()
 *  Unconstrained2dCB()
 *  unify_knots()
 *  unify_knots_cv()
 */

double *dbl_array1 (unsigned nel)

/* allocate a nel element array of doubles */

{
  double *dbl;
  char line[256];

  dbl = (double *) gen_array1 (nel, (unsigned) sizeof(double));

  if (!dbl) {
    sprintf(line, "allocation failure in dbl_array1(), requesting %d\n%s",
	    nel, MemoryStatus());
    errormsg (11, line);
  }

  return (dbl);
}

/* Function: dbl_array2()
 * Purpose: Allocate a two-dimensional double array
 * Method: Use the general purpose allocator gen_array2 to allocate the
 *         array.
 * Arguments:
 *  nrow - the number of rows
 *  ncol - the number of columns
 * Return: The address of the allocated array.
 */
/* Functions referenced by dbl_array2() are:
 *  errormsg()
 *  gen_array2()
 */
/* Functions that reference dbl_array2() are:
 *  addpoints()
 *  AlphaFacetCB()
 *  ApproxCurvCB()
 *  approx_cubic_curv()
 *  approx_cubic_curv1()
 *  approx_cubic_surf()
 *  approx_fnbc()
 *  approx_fnbc_knot()
 *  approx_fn_per()
 *  approx_user_knots()
 *  bernrs1()
 *  berntomono()
 *  BezierKnots()
 *  BsplineToBezier()
 *  build_offset()
 *  build_offsets_int()
 *  build_w_table()
 *  calc_gram_curv()
 *  calc_gram_percurv()
 *  calc_gram_surf()
 *  camber_surface()
 *  compute_lim()
 *  compute_transfm()
 *  ConstrainedCB()
 *  ContourFacets()
 *  ConvexHullTest()
 *  convexhull_test()
 *  con_localize()
 *  CopyFacet()
 *  CopyListCurv()
 *  CopyMinDist()
 *  CopyParUv()
 *  CopyVect()
 *  CurveOrientation()
 *  curve_oslo1()
 *  curve_oslo2()
 *  curve_oslo_per()
 *  cyl_frenet_tr()
 *  cyl_generatrix()
 *  Delaunay()
 *  DrawCurvFunction()
 *  EditKnotsCos()
 *  EditKnotsCurv()
 *  EditUvCB()
 *  evalbsp()
 *  evalbsp_per()
 *  evalderivbsp()
 *  evalderivbsp_per()
 *  evalderivsurf()
 *  evalderivsurf_per()
 *  evalrsurf()
 *  evalrsurf_per()
 *  evalsurf()
 *  evalsurf_per()
 *  FacetHullCB()
 *  FindBot()
 *  FindLft()
 *  FindRht()
 *  FindTop()
 *  FitCosCB()
 *  generatrix()
 *  GridToListCB()
 *  integral_build_offset()
 *  inter_curv()
 *  inter_curv_leo()
 *  inter_percurv()
 *  inter_surf()
 *  int_curve()
 *  IsoLoopIntersect()
 *  ListToUvPtsCB()
 *  LoftSurfCB()
 *  loft_integral()
 *  loft_rational()
 *  main()
 *  mDistCurv()
 *  MeasuredPoints()
 *  MergeKnots()
 *  MinDistance()
 *  MinDistCurvCB()
 *  MinDistToList()
 *  monotobern()
 *  n_matrix_per()
 *  OffsetSurfCB()
 *  offset_curve()
 *  offset_surf()
 *  OrthoDistance()
 *  ParCurv_iso()
 *  partan_dir()
 *  partial()
 *  PostScriptCurvFunction()
 *  PowCurv_to_ParCurv()
 *  PowSurf_to_ParSurf_int()
 *  ProbeRadiusCB()
 *  project_curve_on_surf()
 *  raise_surf()
 *  ReadDeslabFacet()
 *  ReadDeslabHull()
 *  ReadDeslabLocal()
 *  ReadDeslabMinDist()
 *  ReadDeslabVect()
 *  ReadIgesList()
 *  ReadIgesUv()
 *  ReadListCurv()
 *  ReadParUv()
 *  recover_cv()
 *  ReflectPts()
 *  RobustCurvCB()
 *  RobustDistCB()
 *  rotate_section()
 *  SampleCurvCB()
 *  SampleSurfListCB()
 *  solve_gram()
 *  solve_gram_percurv()
 *  solve_gram_pervector()
 *  solve_gram_vector()
 *  soslo()
 *  splitpoly()
 *  split_bspl()
 *  StoreSurfValues()
 *  SubDistance()
 *  SubdivideCos()
 *  SubdivideCurv()
 *  subdiv_cv()
 *  tolerance_region()
 *  tol_loft_rational()
 *  transform()
 *  TransformIgesArray()
 *  TransformIgesVect()
 *  TransformListCurv()
 *  TransformPts()
 *  trans_per()
 *  TriangulateStructured()
 *  TriangulateUnstructured()
 *  TrimToFacetCB()
 *  Unconloc()
 *  Unconstrained2dCB()
 *  UvOkCB()
 *  UvPointsCB()
 *  UvPtsToListCB()
 *  zevalderivsurf()
 */

double **dbl_array2 (unsigned nrow, unsigned ncol)

/* nrow by ncol array of doubles */

{
  double **m;
  char line[256];

  m = (double **) gen_array2 (nrow, ncol, (unsigned)
			       sizeof(double));

  if (!m) {
    sprintf(line, "allocation failure in dbl_array2(), requesting %dx%d\n%s",
	      nrow, ncol, MemoryStatus());
    errormsg (12, line);
  }

  return (m);
}

/* Function: dbl_array3()
 * Purpose: Allocate a three-dimensional double array
 * Method: Use the general purpose allocator gen_array1 to allocate the
 *         array.
 * Arguments:
 *  nrow - the number of rows
 *  ncol - the number of columns
 *  depth - the number of layers
 * Return: The address of the allocated array.
 */
/* Functions referenced by dbl_array3() are:
 *  errormsg()
 *  gen_array1()
 */
/* Functions that reference dbl_array3() are:
 *  CopyGridSurf()
 *  evalderivbsp()
 *  evalderivbsp_per()
 *  evalderivsurf()
 *  evalderivsurf_per()
 *  evaldersurf()
 *  evaldersurf_per()
 *  main()
 *  mDistSurf()
 *  ReadGridSurf()
 *  SampleFacetGridCB()
 *  SampleSurfGridCB()
 *  TransformGridSurf()
 *  zevalderivsurf()
 *  zevaldsurf()
 */

double ***dbl_array3 (unsigned nrow, unsigned ncol, unsigned depth)

/* nrow by ncol by depth array of doubles */

{
  double ***t;
  int i, j;
  char *c1, *c2;
  char line[256];

  t = (double ***) gen_array1 (1, nrow*sizeof(double **));
  if (!t) {
    sprintf(line, "allocation failure 1 in dbl_array3(), requesting %dx%dx%d\n%s",
	      nrow, ncol, depth, MemoryStatus());
    errormsg (14, line);
  }
  c1 = gen_array1(1, nrow*ncol*sizeof(double *));
  c2 = (char *) gen_array1 (nrow*ncol*depth, sizeof(double));
  if (!c1) {
    sprintf(line, "allocation failure 2 in dbl_array3(), requesting %dx%dx%d\n%s",
	      nrow, ncol, depth, MemoryStatus());
    errormsg (15, line);
  }
  if (!c2) {
    sprintf(line, "allocation failure 3 in dbl_array3(), requesting %dx%dx%d/n%s",
	      nrow, ncol, depth, MemoryStatus());
    errormsg (16, line);
  }
  for (i = 0; i < nrow; i++) {	/* initialize the row pointers */
    t[i] = (double **) (c1 + i*ncol*sizeof(double *));
    for (j = 0; j < ncol; j++)	/* initialize the row and col pointers */
      t[i][j] = (double *) (c2 + (i*ncol + j)*depth*sizeof(double));
  }
  return (t);
}

/* Function: free_darray1()
 * Purpose: Deallocate a one-dimensional double array
 * Method: Use the general purpose deallocator free_garray1() to 
 *         deallocate the array.
 * Arguments:
 *  v - the address of the array
 */
/* Functions referenced by free_darray1() are:
 *  free_garray1()
 */
/* Functions that reference free_darray1() are:
 *  ApproxCurvCB()
 *  approx_cubic_curv()
 *  approx_cubic_curv1()
 *  approx_cubic_surf()
 *  approx_fnbc()
 *  approx_fnbc_knot()
 *  approx_fn_per()
 *  approx_user_knots()
 *  bernrs1()
 *  BlendSurfCB()
 *  boehm_curve()
 *  boehm_surface()
 *  bone_intersect()
 *  bone_intersect3()
 *  BsplineToBezierKnots()
 *  BsplineToBezierSurf()
 *  build_offset()
 *  build_offsets_int()
 *  calc_disc()
 *  calc_disc_curv()
 *  calc_hartley()
 *  calc_hartley_per()
 *  calc_interval_egeom()
 *  calc_knots_foley()
 *  calc_knots_hartley()
 *  calc_knots_per()
 *  calc_par_foley()
 *  calc_par_hartley()
 *  CamberBisect()
 *  CamberBisect3d()
 *  camber_line_2d()
 *  camber_line_3d()
 *  camber_surface()
 *  CheckCosContinuity()
 *  CheckCurvContinuity()
 *  CheckSurfContinuity()
 *  check_one_dir()
 *  ConstrainedCB()
 *  ConvertToCos()
 *  convexbox()
 *  ConvexHullSurf()
 *  ConvexHullTest()
 *  convexhull_test()
 *  con_localize()
 *  copyegeom()
 *  copyfgeom()
 *  critical_bone()
 *  critical_bone3()
 *  CurveOrientation()
 *  curve_fit()
 *  curve_osloper_lt()
 *  curve_osloper_rt()
 *  curve_oslo_per()
 *  curv_cont()
 *  Delaunay()
 *  DeleteCos()
 *  DeleteCosPts()
 *  DeleteCurv()
 *  DeleteCurvPts()
 *  DeleteList()
 *  DeleteSurf()
 *  DeleteTrim()
 *  DeleteTrimPts()
 *  del_dbl_knots()
 *  deriv_sign()
 *  EditKnotsCos()
 *  EditKnotsCurv()
 *  findroot()
 *  FitCosCB()
 *  FitCurvCB()
 *  FitSurfCB()
 *  fit_curve()
 *  fn_camber()
 *  fn_camber3()
 *  fn_thckns()
 *  fn_thckns3()
 *  free_copious()
 *  free_egeom()
 *  free_fgeom()
 *  free_pgeom()
 *  free_sgeom()
 *  GenCylCB()
 *  global_pos_error()
 *  integral_build_offset()
 *  integral_sample_offset()
 *  inter_curv()
 *  inter_percurv()
 *  inter_surf()
 *  int_curve()
 *  int_number()
 *  JoinCoss()
 *  JoinCurvs()
 *  JoinSurfs()
 *  LoftSurfCB()
 *  loft_integral()
 *  main()
 *  MeasuredPoints()
 *  MergeKnots()
 *  merge_can_off()
 *  MinDistance()
 *  MinDistCurvCB()
 *  monotobern()
 *  nbasisd()
 *  nbasisd_per()
 *  offset_curve()
 *  OrthoDistance()
 *  OrthoProjCosCB()
 *  OrthoProjCurvCB()
 *  ParCurv_to_PowCurv()
 *  ParSurf_approx()
 *  partan_dir()
 *  PowCurv_to_ParCurv()
 *  PowSurf_to_ParSurf_int()
 *  PowSurf_to_ParSurf_loft()
 *  project_curve_onsurf()
 *  project_curve_on_surf()
 *  raise_surf()
 *  ReadParCurv()
 *  ReadParCurv_Per()
 *  ReadParSurf()
 *  ReadParSurf_Peru()
 *  reflectegeom()
 *  reflectfgeom()
 *  renew_knots()
 *  rk4()
 *  RobustCurvCB()
 *  RobustDistCB()
 *  rootisolate()
 *  sample()
 *  sample_blend()
 *  sample_csurf()
 *  sample_mod_surf()
 *  sample_offset()
 *  sample_per()
 *  sample_proj()
 *  sample_surf_test()
 *  ScatteredSurfCB()
 *  scattered_fit()
 *  solve_gram_percurv()
 *  solve_gram_pervector()
 *  soslo()
 *  splitpoly2()
 *  SplitSurfCancelCB()
 *  SplitSurfOkCB()
 *  SubDistance()
 *  SubdivideCos()
 *  SubdivideCurv()
 *  SurfSplitAddUCB()
 *  SurfSplitAddVCB()
 *  SurfSplitRemoveUCB()
 *  SurfSplitRemoveVCB()
 *  tol_sample_offsets()
 *  transformegeom()
 *  transformfgeom()
 *  TrimSurfaceCB()
 *  Unconloc()
 *  Unconstrained2dCB()
 *  unify_knots()
 *  unify_knots_cv()
 */

void free_darray1 (double *v)

{
  free_garray1((char *)v);
}

/* Function: free_darray2()
 * Purpose: Deallocate a two-dimensional double array
 * Method: Use the general purpose deallocator free_garray1() to 
 *         deallocate the array.
 * Arguments:
 *  m - the address of the array
 */
/* Functions referenced by free_darray2() are:
 *  free_garray1()
 */
/* Functions that reference free_darray2() are:
 *  addpoints()
 *  AlphaFacetCB()
 *  ApproxCurvCB()
 *  approx_cubic_curv()
 *  approx_cubic_curv1()
 *  approx_cubic_surf()
 *  approx_fnbc()
 *  approx_fnbc_knot()
 *  approx_user_knots()
 *  bernrs1()
 *  berntomono()
 *  BsplineToBezier()
 *  BsplineToBezierSurf()
 *  build_offset()
 *  build_offsets_int()
 *  calc_gram_curv()
 *  calc_gram_percurv()
 *  calc_gram_surf()
 *  camber_surface()
 *  compute_lim()
 *  ConstrainedCB()
 *  ContourFacets()
 *  ConvexHullTest()
 *  convexhull_test()
 *  con_confun()
 *  con_localize()
 *  CurveOrientation()
 *  curve_oslo1()
 *  curve_oslo2()
 *  curve_oslo_per()
 *  cyl_generatrix()
 *  decasteljau_curve_contpts()
 *  Delaunay()
 *  DeleteFacet()
 *  DeleteHull()
 *  DeleteList()
 *  DeleteSurf()
 *  DeleteSurfPts()
 *  DeleteUv()
 *  DeleteVect()
 *  DrawCurvFunction()
 *  EditKnotsCos()
 *  EditKnotsCurv()
 *  EditUvCB()
 *  evalbsp()
 *  evalbsp_per()
 *  evalderivbsp()
 *  evalderivbsp_per()
 *  evalderivsurf()
 *  evalderivsurf_per()
 *  evalrsurf()
 *  evalrsurf_per()
 *  evalsurf()
 *  evalsurf_per()
 *  FacetHullCB()
 *  FitCosCB()
 *  generatrix()
 *  integral_build_offset()
 *  inter_curv()
 *  inter_curv_leo()
 *  inter_percurv()
 *  inter_surf()
 *  int_curve()
 *  IsoLoopIntersect()
 *  localize_diagnostic()
 *  localize_sumsq()
 *  localize_sumsq_opt()
 *  LoftSurfCB()
 *  loft_integral()
 *  loft_rational()
 *  main()
 *  mDistCurv()
 *  MeasuredPoints()
 *  MergeKnots()
 *  MinDistance()
 *  MinDistCurvCB()
 *  monotobern()
 *  n_matrix_per()
 *  OffsetSurfCB()
 *  offset_curve()
 *  OrthoDistance()
 *  ParCurv_iso()
 *  partan_dir()
 *  PostScriptCurvFunction()
 *  PowCurv_to_ParCurv()
 *  PowSurf_to_ParSurf_int()
 *  raise_surf()
 *  recover_cv()
 *  RobustCurvCB()
 *  RobustDistCB()
 *  rotate_section()
 *  solve_gram()
 *  solve_gram_percurv()
 *  solve_gram_pervector()
 *  solve_gram_vector()
 *  soslo()
 *  splitpoly()
 *  split_bspl()
 *  SubDistance()
 *  SubdivideCos()
 *  SubdivideCurv()
 *  subdiv_cv()
 *  tolerance_region()
 *  tol_loft_rational()
 *  transform()
 *  transformegeom()
 *  transformfgeom()
 *  TransformGridSurf()
 *  TransformIgesArray()
 *  TransformIgesVect()
 *  TransformListCurv()
 *  TransformPts()
 *  trans_per()
 *  Unconloc()
 *  Unconstrained2dCB()
 *  UvOkCB()
 *  zevalderivsurf()
 */

void free_darray2 (double **m)

{
  free_garray1((char *)&m[0][0]);
  free_garray1((char *)m);
}

/* Function: free_darray3()
 * Purpose: Deallocate a three-dimensional double array
 * Method: Use the general purpose deallocator free_garray1() to 
 *         deallocate the array.
 * Arguments:
 *  t - the address of the array
 */
/* Functions referenced by free_darray3() are:
 *  free_garray1()
 */
/* Functions that reference free_darray3() are:
 *  DeleteGrid()
 *  evalderivbsp()
 *  evalderivbsp_per()
 *  evalderivsurf()
 *  evaldersurf()
 *  evaldersurf_per()
 *  main()
 *  mDistSurf()
 *  zevalderivsurf()
 *  zevaldsurf()
 */

void free_darray3 (double ***t)

{
  free_garray1((char *)&t[0][0][0]);
  free_garray1((char *)&t[0][0]);
  free_garray1((char *)t);
}

/************************************************************************
 *									*
			 Copyright (C) 1989, 1990
	   Massachusetts Institute of Technology, Cambridge, MA
			    All rights reserved
 *									*
 ************************************************************************/
/************************************************************************
 *									*
			   Modification History

     16 Jul 90 - allocate new knot array to be same as old->?kmem
     23 Apr 89 - changed # of knots copied, fixed fgeom
     01 Feb 89 - written by Bradley A. Moran
 *									*
 ************************************************************************/

#include <stdlib.h>
#include <memory.h>
#include "gen.h"

/* Function: copyegeom()
 * Purpose: Copy a NURBS curve structure
 * Method: The NURBS curve structure ParCurv is defined in the header
 *         file "gen.h".
 * Arguments:
 *  old - address of the existing structure
 *  new - address of an existing structure into which old is copied.
 *        If this is equal to NULL a new structure is allocated
 * Return: Address of the structure containing the copy
 */
/* Functions referenced by copyegeom() are:
 *  copyvector()
 *  dbl_array1()
 *  free_darray1()
 *  free_varray1()
 *  gen_array1()
 *  ptr_array1()
 */
/* Functions that reference copyegeom() are:
 *  addsections_canal()
 *  ApproxCurvCB()
 *  BsplineToBezierCurvCB()
 *  CamberLineCB()
 *  camber_surface()
 *  CheckCosContinuity()
 *  CheckCurvContinuity()
 *  clip_bezier_curve()
 *  ContourFacetsCB()
 *  ConvexHullTest()
 *  convexhull_test()
 *  CopyCos()
 *  CopyCurv()
 *  CopyTrimSurf()
 *  CurveOrientation()
 *  cyl_sample_gencyl()
 *  EditKnotsCos()
 *  EditKnotsCurv()
 *  FairCos()
 *  FairCurv()
 *  int_sample_gencyl()
 *  InvertCurve()
 *  InvertTrimSurf()
 *  JoinCoss()
 *  JoinCurvs()
 *  knot_rem_cv()
 *  LoftSurfCB()
 *  MeasuredLineCB()
 *  MergeKnots()
 *  offset_curve()
 *  ParSurf_approx()
 *  raise_byone()
 *  ReadIgesRuledSurf()
 *  ReadIgesSurfRev()
 *  ReadIgesTabCyl()
 *  RecursiveInt()
 *  ReflectTrimSurf()
 *  RevCurvParam()
 *  RuledSurfCB()
 *  sample_csurf()
 *  sample_gencyl()
 *  SubdivideCos()
 *  SubdivideCurv()
 *  SurfRev_to_ParSurf()
 *  TransformTrimSurf()
 *  Unconstrained2dCB()
 */

ParCurv *copyegeom (register ParCurv *old, register ParCurv *new)

/* Copy an egeom; allocate new if nil */

{
     register int i, nbytes;

/*
 * Check to see if the target pointer is nil.  If it is, it needs to be
 * allocated.
 */
     if (!new) {
	  new = (ParCurv *) gen_array1(1, sizeof(ParCurv));

/*
 * Upon successful allocation, copy all contents of the original and then
 * allocate the knot vector and the control points.  Copy the knots and
 * control points.
 */
	  if (new) {
	       (void) memcpy ((char *) new, (char *) old, sizeof(ParCurv));
	       new->knots = dbl_array1 ((unsigned) old->kmem);
	       new->contpts = (vector **) ptr_array1 ((unsigned)
						      old->pmem);
	  }
	  else {
	       perror ("copyegeom");
	       exit (1);
	  }
     }
/*
 * If the target structure was previously allocated, copy information item
 * by item.  Check the size of the knot vector and control point array.
 * Increase the size of either if it is too small.  Copy the knots and
 * control points.
 */
     else {
	  if (new->kmem < old->kmem) {
	       free_darray1 (new->knots);
	       new->knots = dbl_array1 ((unsigned) old->kmem);
	       new->kmem = old->kmem;
	  }
	  if (new->pmem < old->pmem) {
	       free_varray1 (new->contpts, (unsigned) new->ncontpts);
	       new->contpts = (vector **) ptr_array1 ((unsigned) old->pmem);
	       new->pmem = old->pmem;
	  }
	  new->type = old->type;
	  new->ncontpts = old->ncontpts;
	  new->order = old->order;
     }
     nbytes = sizeof(double) * old->kmem;
     (void) memcpy ((char *) new->knots, (char *) old->knots, nbytes);

     for (i = 0; i < old->ncontpts; i++)
	  new->contpts[i] = copyvector (old->contpts[i], new->contpts[i]);

     return new;
}

/* Function: copyfgeom()
 * Purpose: Copy a NURBS surface structure
 * Method: The NURBS surface structure ParSurf is defined in the header
 *         file "gen.h".
 * Arguments:
 *  old - address of the existing structure
 *  new - address of an existing structure into which old is copied.
 *        If this is equal to NULL a new structure is allocated
 * Return: Address of the structure containing the copy
 */
/* Functions referenced by copyfgeom() are:
 *  copyvector()
 *  dbl_array1()
 *  free_darray1()
 *  free_varray2()
 *  gen_array1()
 *  ptr_array2()
 *  vec_array2()
 */
/* Functions that reference copyfgeom() are:
 *  addpoints_surf()
 *  ApproxSurfCB()
 *  BsplineToBezierSurf()
 *  BsplineToBezierSurfCB()
 *  build_offset()
 *  camber_surface()
 *  CheckSurfContinuity()
 *  check_surf_param()
 *  ConstrainedCB()
 *  ContourFacetsCB()
 *  ConvexHullSurf()
 *  convexhull_surf()
 *  CopySurf()
 *  CopyTrimSurf()
 *  EditKnotsSurf()
 *  FairSurf()
 *  InvertSurface()
 *  InvertTrimSurf()
 *  JoinSurfs()
 *  knot_rem_and_pert()
 *  merge_can_off()
 *  ParSurf_approx()
 *  PowSurf_to_ParSurf_loft()
 *  renew_knots()
 *  RevSurfUParam()
 *  RevSurfVParam()
 *  SplitSurf()
 *  SplitSurface()
 *  SubdivSurfNxM()
 *  SwitchSurfParam()
 *  TransposeSurface()
 *  TrimSurfaceCB()
 *  Unconloc()
 */

ParSurf *copyfgeom (ParSurf *old, ParSurf *new)

{
     register int i, j;

     if (!new) {
	  new = (ParSurf *) gen_array1(1, sizeof(ParSurf));

	  if (new) {
	       (void) memcpy ((char *) new, (char *) old, sizeof(ParSurf));
	       new->uknots = dbl_array1 ((unsigned) old->ukmem);
	       new->vknots = dbl_array1 ((unsigned) old->vkmem);
	       new->contpts = vec_array2 ((unsigned) old->upmem,
					  (unsigned) old->vpmem);
	  }
	  else {
	       perror ("copyfgeom");
	       exit (1);
	  }
     }
     else {
	  if (new->ukmem < old->ukmem) {
	       free_darray1 (new->uknots);
	       new->uknots = dbl_array1 ((unsigned) old->ukmem);
	       new->ukmem = old->ukmem;
	  }
	  if (new->vkmem < old->vkmem) {
	       free_darray1 (new->vknots);
	       new->vknots = dbl_array1 ((unsigned) old->vkmem);
	       new->vkmem = old->vkmem;
	  }
	  if (new->vpmem < old->vpmem || new->upmem < old->upmem) {
	       free_varray2 (new->contpts, (unsigned) new->ucontpts,
			     (unsigned) new->vcontpts);
	       new->contpts = (vector ***)
		    ptr_array2 ((unsigned) old->upmem, (unsigned)
				old->vpmem);
	       new->upmem = old->upmem, new->vpmem = old->vpmem;
	  }
	  new->type = old->type;
	  new->uorder = old->uorder, new->ucontpts = old->ucontpts;
	  new->vorder = old->vorder, new->vcontpts = old->vcontpts;
     }
     (void) memcpy ((char *) new->uknots, (char *) old->uknots,
		    (int) (old->ukmem * sizeof(double)));
     (void) memcpy ((char *) new->vknots, (char *) old->vknots,
		    (int) (old->vkmem * sizeof(double)));

     for (i = 0; i < old->ucontpts; i++)
	  for (j = 0; j < old->vcontpts; j++)
	       new->contpts[i][j] = copyvector (old->contpts[i][j],
						new->contpts[i][j]);
     return new;
}

/* Function: copyvector()
 * Purpose: Copy a vector structure
 * Method: The vector structure is defined in the header file "gen.h".
 * Arguments:
 *  old - address of the existing structure
 *  new - address of an existing structure into which old is copied.
 *        If this is equal to NULL a new structure is allocated.
 * Return: Address of the structure containing the copy
 */
/* Functions referenced by copyvector() are:
 *  vectalloc()
 */
/* Functions that reference copyvector() are:
 *  add_curv()
 *  add_surf()
 *  approx_cubic_curv()
 *  approx_cubic_surf()
 *  approx_fnbc()
 *  approx_fnbc_knot()
 *  approx_fn_per()
 *  approx_user_knots()
 *  BezierContpts()
 *  boehm_curve()
 *  boehm_surface()
 *  calc_hartley_per()
 *  calc_knots_per()
 *  calc_par_hartley()
 *  calc_points_curv()
 *  Camber2dTo3d()
 *  Camber2d_3d()
 *  CamberBisect()
 *  CamberBisect3d()
 *  camber_3d()
 *  camber_surf_3d()
 *  check_surf_param()
 *  circle_gen()
 *  cone_gen()
 *  construct_bspl()
 *  construct_jbspl()
 *  con_confun()
 *  copyegeom()
 *  copyegeom_leo()
 *  copyfgeom()
 *  copy_contpts()
 *  corner_offset()
 *  cross_link()
 *  cylinder_gen()
 *  cyl_frenet_tr()
 *  decasteljau_curve_contpts()
 *  eval_fnbc()
 *  extract_cv()
 *  extract_edge()
 *  extract_patch()
 *  fair_knot()
 *  fair_per_knot()
 *  fn()
 *  fn_camber()
 *  get_cv_removal()
 *  get_iso_removal()
 *  inter_curv()
 *  inter_curv_leo()
 *  inter_percurv()
 *  inter_surf()
 *  line_gen()
 *  localize_diagnostic()
 *  localize_sumsq()
 *  localize_sumsq_opt()
 *  MakeCamberLine3d()
 *  matr_m_vect()
 *  MergeKnots()
 *  merge_tol_edges()
 *  mult4x4()
 *  opt_param()
 *  opt_param_per()
 *  ParCurv_hodograph()
 *  ParCurv_to_PowCurv()
 *  ParSurf_to_PowSurf()
 *  plane_gen()
 *  PowCurv_to_ParCurv()
 *  PowSurf_to_ParSurf_int()
 *  raise_byone()
 *  ReadParCurv()
 *  ReadParCurv_Per()
 *  ReadParSurf()
 *  ReadParSurf_Peru()
 *  RevCurvParam()
 *  RevSurfUParam()
 *  RevSurfVParam()
 *  RulSurf_to_ParSurf()
 *  solve_bc()
 *  sphere_gen()
 *  SplitSurf()
 *  split_it()
 *  store_edge()
 *  store_vertices()
 *  subbezier()
 *  SubDistance()
 *  subdivbs()
 *  SubdivideCos()
 *  SubdivideCurv()
 *  SurfRev_to_ParSurf()
 *  swapsurf()
 *  SwitchSurfParam()
 *  torus_gen()
 *  unitvector()
 *  unitvector1()
 *  vec_rot()
 *  vec_rot_fix()
 */

vector *copyvector (vector *orig, vector *copy)

{
     if (!copy)			/* allocate target if nil */
	  copy = vectalloc ();
     (void) memcpy ((char *) copy, (char *) orig, (int) sizeof(vector));
     return copy;
}

/************************************************************************
 *									*
			Copyright (C) 1989 by
	Massachusetts Institute of Technology, Cambridge, MA
			 All rights reserved
 *									*
 ************************************************************************/
/************************************************************************
 *									*
			 Modification History

     06 Feb 90 -- made free_egeom dispatch on kmem and pmem
     30 Mar 89 -- changed free_egeom to use free_?array1
 *									*
 ************************************************************************/

#include <malloc.h>
#include <gen.h>

/* Function: egeomalloc1()
 * Purpose: Allocate a NURBS curve structure
 * Method:  The NURBS curve structure is defined in the header file
 *          "gen.h".
 *          The allocated structure is large enough to hold a curve
 *          of the specified order and control points.
 * Arguments:
 *  order - the order of the curve
 *  ncontpts - the number of control points of the curve
 * Return: The address of the allocated structure
 */
/* Functions referenced by egeomalloc1() are:
 *  dbl_array1()
 *  errormsg()
 *  gen_array1()
 *  vec_array1()
 */
/* Functions that reference egeomalloc1() are:
 *  ApproxCurvCB()
 *  bernrs1()
 *  BsplineToBezier()
 *  calc_disc_surf()
 *  camber_surface()
 *  circle_gen()
 *  clip_bezier_curve()
 *  compare_curves()
 *  construct_bspl()
 *  ConvertToCos()
 *  ConvexHullTest()
 *  convexhull_test()
 *  corner_offset()
 *  CurveOrientation()
 *  curve_fit()
 *  cyl_generatrix_lk()
 *  decasteljau_curve()
 *  EditKnotsCos()
 *  EditKnotsCurv()
 *  extract_cv()
 *  fair_cubic_surf()
 *  FitCosCB()
 *  FitCurvCB()
 *  generatrix_lk()
 *  global_pos_error()
 *  int_curve()
 *  int_number()
 *  knot_rem_cv()
 *  line_gen()
 *  link_alloc()
 *  MakeCamberLine()
 *  MakeCamberLine3d()
 *  MakeThickness()
 *  MakeThickness3d()
 *  MergeKnots()
 *  OrthoDistance()
 *  ParCurv_hodograph()
 *  ParCurv_iso()
 *  ParSurf_approx()
 *  ParSurf_to_SurfRev()
 *  ParSurf_to_TabCyl()
 *  PowCurv_to_ParCurv()
 *  raise_byone()
 *  ReadIges126()
 *  ReadParCurv()
 *  recover_cv()
 *  repeat_offset()
 *  sample_blend()
 *  sample_surf_test()
 *  solve_int()
 *  split_bspl()
 *  SubdivideCos()
 *  SubdivideCurv()
 */

				/* compact egeom */
ParCurv *egeomalloc1 (int order, int ncontpts)

{
     ParCurv *egm;
     char line[256];

     egm = (ParCurv *) gen_array1 (1, sizeof(ParCurv));

     if (!egm) {
       sprintf(line,
	       "allocation failure in egeomalloc1(), requesting %d %d\n%s",
	       order, ncontpts, MemoryStatus);
	  errormsg (0, line);
     }

     egm->knots = dbl_array1(order + ncontpts);   /* knot vector */
     egm->contpts = vec_array1(ncontpts);         /* control points array */

     egm->order = order, egm->ncontpts = ncontpts;
     egm->kmem = order + ncontpts;   /* size of knot vector */
     egm->pmem = ncontpts;           /* size of control points array */

     return (egm);
}

/* Function: egeomalloc2()
 * Purpose: Allocate a NURBS curve structure
 * Method:  The NURBS curve structure is defined in the header file
 *          "gen.h".
 *          The allocated structure is large enough to hold a curve
 *          of up to the specified maximum order and control points.
 * Arguments:
 *  order - the order of the curve
 *  ncontpts - the number of control points of the curve
 * Return: The address of the allocated structure
 */
/* Functions referenced by egeomalloc2() are:
 *  dbl_array1()
 *  errormsg()
 *  gen_array1()
 */
/* Functions that reference egeomalloc2() are:
 *  approx_fnbc()
 *  approx_fnbc_knot()
 *  approx_fn_per()
 *  camber_surface()
 *  ConvexHullSurf()
 *  convexhull_surf()
 *  EditKnotsCos()
 *  EditKnotsCurv()
 *  fit_curve()
 *  LoftSurfCB()
 *  merge_can_off()
 *  offset_curve()
 *  ParSurf_approx()
 *  raise_surf()
 *  sample_csurf()
 *  SubdivideCos()
 *  SubdivideCurv()
 */

				/* non-compact egeom */
ParCurv *egeomalloc2 (int maxorder, int maxcontpts)

{
     ParCurv *egm;
     char *c, line[256];

     egm = (ParCurv *) gen_array1 (1, sizeof(ParCurv));

     if (!egm) {
       sprintf(line,
	       "allocation failure 1 in egeomalloc2(), requesting %d %d\n%s",
	       maxorder, maxcontpts, MemoryStatus());
       errormsg (0, line);
     }

     egm->knots = dbl_array1 (maxorder + maxcontpts);   /* knot vector */
     egm->contpts = (vector **) gen_array1 (maxcontpts,
					    sizeof(vector *));
     if (!egm->contpts) {
       sprintf(line,
	       "allocation failure 2 in egeomalloc2(), requesting %d %d\n%s",
	       maxorder, maxcontpts, MemoryStatus());
       errormsg (20, line);
     }

     egm->kmem = maxorder + maxcontpts;  /* maximum size of knot vector */
     egm->pmem = maxcontpts;             /* maximum size of control points */

     return (egm);
}

/* Function: free_egeom()
 * Purpose: Deallocate a NURBS curve structure
 * Method: First deallocate the knot vector and control points arrays,
 *         then use the general purpose deallocator free_garray1 to
 *         deallocate the curve structure
 * Arguments:
 *  egm - address of the curve structure
 */
/* Functions referenced by free_egeom() are:
 *  free_darray1()
 *  free_garray1()
 *  free_varray1()
 */
/* Functions that reference free_egeom() are:
 *  addsections_canal()
 *  approx_fnbc()
 *  approx_fnbc_knot()
 *  approx_fn_per()
 *  avg_arc_length()
 *  bernrs1()
 *  BlendSurfCB()
 *  BsplineToBezierCurvCB()
 *  calc_disc_surf()
 *  CamberLineCB()
 *  camber_surface()
 *  clip_bezier_curve()
 *  compare_curves()
 *  ContourFacetsCB()
 *  ConvexHullSurf()
 *  ConvexHullTest()
 *  convexhull_surf()
 *  convexhull_test()
 *  corner_offset()
 *  CurveOrientation()
 *  cyl_sample_gencyl()
 *  DeleteCos()
 *  DeleteCurv()
 *  DeleteTrim()
 *  EditKnotsCos()
 *  EditKnotsCurv()
 *  FairCos()
 *  FairCurv()
 *  fair_cubic_surf()
 *  fit_curve()
 *  FreeIgesDirList()
 *  free_earray1()
 *  free_earray2()
 *  free_egeom_array1()
 *  get_error_cv()
 *  global_pos_error()
 *  intersect_curve_to_axis()
 *  int_number()
 *  int_sample_gencyl()
 *  knot_rem_cv()
 *  LoftSurfCB()
 *  loop_rayintersect()
 *  main()
 *  mDistCurv()
 *  MeasuredLineCB()
 *  offset_curve()
 *  OrthoDistance()
 *  ParCurv_approx()
 *  ParSurf_approx()
 *  PowSurf_to_ParSurf_loft()
 *  raise_byone()
 *  raise_surf()
 *  ReadDeslabLocal()
 *  ReadDeslabMinDist()
 *  ReadIgesRuledSurf()
 *  ReadIgesSurfRev()
 *  ReadIgesTabCyl()
 *  RecursiveInt()
 *  RevCurvParam()
 *  RobustDistCB()
 *  RuledSurfCB()
 *  RulSurf_to_ParSurf()
 *  sample_blend()
 *  sample_csurf()
 *  sample_gencyl()
 *  sample_mod_surf()
 *  sample_surf_test()
 *  SaveIgesCos()
 *  SaveIgesSurf()
 *  SaveIgesTrim()
 *  solve_int()
 *  SpanwiseIsU()
 *  split_bspl()
 *  SubdivideCos()
 *  SubdivideCurv()
 *  SurfRev_to_ParSurf()
 *  TabCyl_to_ParSurf()
 *  Unconstrained2dCB()
 */

void free_egeom (ParCurv *egm)
{
  if (egm->pmem)                 /* delete control points array */
    free_varray1 (egm->contpts, (unsigned) egm->ncontpts);
  
  if (egm->kmem)                 /* delete knot vector */
    free_darray1 (egm->knots);

  free_garray1((char *)egm);     /* delete structure */
}

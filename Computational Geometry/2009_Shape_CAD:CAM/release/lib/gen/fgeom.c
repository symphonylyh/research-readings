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
 *									*
 ************************************************************************/

#include <malloc.h>
#include "gen.h"

/* Function: fgeomalloc1()
 * Purpose: Allocate a NURBS surface structure
 * Method:  The NURBS surface structure is defined in the header file
 *          "gen.h".
 *          The allocated structure is large enough to hold a surface
 *          of the specified order and control points.
 * Arguments:
 *  uord - the u order of the surface
 *  vord - the v order of the surface
 *  upts - the number of u control points of the surface
 *  vpts - the number of v control points of the surface
 * Return: The address of the allocated structure
 */
/* Functions referenced by fgeomalloc1() are:
 *  dbl_array1()
 *  errormsg()
 *  gen_array1()
 *  vec_array2()
 */
/* Functions that reference fgeomalloc1() are:
 *  addpoints_surf()
 *  ApproxSurfCB()
 *  BlendSurfCB()
 *  BsplineToBezierSurf()
 *  check_one_dir()
 *  check_surf_param()
 *  cone_gen()
 *  construct_jbspl()
 *  corner_offset()
 *  cylinder_gen()
 *  EditKnotsSurf()
 *  extract_patch()
 *  FitSurfCB()
 *  knot_rem_and_pert()
 *  loft_integral()
 *  loft_rational()
 *  merge_can_off()
 *  ParSurf_approx()
 *  patch_with_com_knots()
 *  plane_gen()
 *  PowSurf_to_ParSurf()
 *  ReadIges128()
 *  ReadParSurf()
 *  recover_patch()
 *  renew_knots()
 *  RulSurf_to_ParSurf()
 *  scattered_fit()
 *  sphere_gen()
 *  SplitSurf()
 *  SurfRev_to_ParSurf()
 *  swapsurf()
 *  SwitchSurfParam()
 *  TabCyl_to_ParSurf()
 *  tol_loft_rational()
 *  torus_gen()
 */

ParSurf *fgeomalloc1 (int uord, int vord, int upts, int vpts)
{
  ParSurf *fgm;
  char line[256];

  fgm = (ParSurf *) gen_array1 (1, sizeof(ParSurf));
  if (!fgm) {
    sprintf(line,
	    "allocation failure in fgeomalloc1(), requesting %d %d %d %d\n%s",
	    uord, vord, upts, vpts, MemoryStatus()); 
    errormsg (0, line);
  }

  fgm->uknots = dbl_array1 (uord + upts);   /* u knot vector */
  fgm->vknots = dbl_array1 (vord + vpts);   /* v knot vector */
  fgm->contpts = vec_array2 (upts, vpts);   /* control points array */

  fgm->uorder = uord;         /* u order */
  fgm->ucontpts = upts;       /* number of u control points */
  fgm->vorder = vord;         /* v order */
  fgm->vcontpts = vpts;       /* number of v control points */

  fgm->ukmem = uord + upts;   /* size of u knot vector */
  fgm->upmem = upts;          /* size of u dimension of control points array */
  fgm->vkmem = vord + vpts;   /* size of v knot vector */
  fgm->vpmem = vpts;          /* size of v dimension of control points array */

  return (fgm);
}

/* Function: fgeomalloc2()
 * Purpose: Allocate a NURBS surface structure
 * Method:  The NURBS surface structure is defined in the header file
 *          "gen.h".
 *          The allocated structure is large enough to hold a surface
 *          of up to the specified maximum order and control points.
 * Arguments:
 *  maxuorder - the maximum u order of the surface
 *  maxvorder - the maximum v order of the surface
 *  maxupts - the maximum number of u control points of the surface
 *  maxvpts - the maximum number of v control points of the surface
 * Return: The address of the allocated structure
 */
/* Functions referenced by fgeomalloc2() are:
 *  dbl_array1()
 *  errormsg()
 *  gen_array1()
 */
/* Functions that reference fgeomalloc2() are:
 *  ConvexHullSurf()
 *  convexhull_surf()
 *  integral_offset()
 *  integral_offset_surf()
 *  merge_can_off()
 *  rational_offset()
 *  subbezier()
 */

ParSurf *fgeomalloc2 (int maxuorder, int maxvorder, int maxupts,
		      int maxvpts)
{
  ParSurf *fgm;
  int i;
  char *c, line[256];

  fgm = (ParSurf *) gen_array1 (1, sizeof(ParSurf));
  if (!fgm) {
    sprintf(line,
	    "allocation failure 1 in fgeomalloc2(), requesting %d %d %d %d\n%s",
	    maxuorder, maxvorder, maxupts, maxvpts, MemoryStatus());
    errormsg (0, line);
  }

  fgm->uknots = dbl_array1 (maxuorder + maxupts);   /* u knot vector */
  fgm->vknots = dbl_array1 (maxvorder + maxvpts);   /* v knot vector */
  fgm->contpts = (vector ***)gen_array1(1, maxupts*sizeof(vector **));

  if (!fgm->contpts) {
    sprintf(line,
	     "allocation failure 2 in fgeomalloc2(), requesting %d %d %d %d\n%s",
	    maxuorder, maxvorder, maxupts, maxvpts, MemoryStatus());
    errormsg (15, line);
  }

  c = (char *) gen_array1 (maxupts*maxvpts, sizeof(vector *));
  if (!c) {
    sprintf(line,
	    "allocation failure 3 in fgeomalloc2(), requesting %d %d %d %d\n%s",
	    maxuorder, maxvorder, maxupts, maxvpts, MemoryStatus());
    errormsg (16, line);
  }

  fgm->ukmem = maxuorder + maxupts;   /* maximum size of u knot vector */
  fgm->upmem = maxupts;       /* max size of control points array in u */
  fgm->vkmem = maxvorder + maxvpts;   /* maximum size of v knot vector */
  fgm->vpmem = maxvpts;       /* max size of control points array in v */

  for (i = 0; i < maxupts; i++)
    fgm->contpts[i] = (vector **) (c + i*maxvpts*sizeof(vector *));

  return (fgm);
}

/* Function: free_fgeom()
 * Purpose: Deallocate a NURBS surface structure
 * Method: First deallocate the knot vector and control points arrays,
 *         then use the general purpose deallocator free_garray1() to
 *         deallocate the surface structure
 * Arguments:
 *  fgm - address of the surface structure
 */
/* Functions referenced by free_fgeom() are:
 *  free_darray1()
 *  free_garray1()
 *  free_varray2()
 */
/* Functions that reference free_fgeom() are:
 *  addpoints_surf()
 *  addsections_canal()
 *  BlendSurfCB()
 *  BsplineToBezierSurf()
 *  BsplineToBezierSurfCB()
 *  build_offset()
 *  CamberLineCB()
 *  camber_surface()
 *  canal_to_int()
 *  canal_to_rat()
 *  check_one_dir()
 *  ConstrainedCB()
 *  ConvexHullSurf()
 *  convexhull_surf()
 *  cyl_canal_to_rat()
 *  DeleteSurf()
 *  DeleteTrim()
 *  EditKnotsSurf()
 *  FairSurf()
 *  FreeIgesDirList()
 *  free_fgeom_array2()
 *  get_error()
 *  main()
 *  mDistSurf()
 *  MeasuredLineCB()
 *  merge_can_off()
 *  ParSurf_approx()
 *  PowSurf_to_ParSurf_loft()
 *  ReadDeslabLocal()
 *  ReadDeslabMinDist()
 *  renew_knots()
 *  RevSurfUParam()
 *  RevSurfVParam()
 *  sample_csurf()
 *  SplitSurf()
 *  SplitSurface()
 *  subbezier()
 *  SubDistance()
 *  SubdivSurfNxM()
 *  SurfRev_to_ParSurf()
 *  SwitchSurfParam()
 *  tolerance_region()
 *  TrimSurfaceCB()
 *  Unconloc()
 *  unify_knots()
 */

void free_fgeom (ParSurf *fgm)	 /* free the memory of an fgeom */
{
  if (fgm->upmem && fgm->vpmem)
    free_varray2 (fgm->contpts, (unsigned) fgm->ucontpts, (unsigned)
		  fgm->vcontpts);
  if (fgm->ukmem)
    free_darray1 (fgm->uknots);
  if (fgm->vkmem)
    free_darray1 (fgm->vknots);
  free_garray1((char *)fgm);
}

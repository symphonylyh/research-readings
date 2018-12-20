/************************************************************************
 *									*
		Copyright (C) 1989 by Bradley A. Moran
			 All rights reserved.

    The  author hereby grants to  MIT  permission to reproduce and  to
    distribute  copies of  this source  code  document in whole or  in
    part.
 *									*
 ************************************************************************/

# include <math.h>
# include <memory.h>
# include <stdio.h>
# include <gen.h>

#define MACHPREC 1.1102230246251567E-13	  /* changed on 25 JUN 90 */

/* Function: add_vect()
 * Purpose: Add two vectors
 * Method: A new vector structure is allocated and holds the sum
 *         of the two vectors.
 * Arguments:
 *  v1 - address of the first vector
 *  v2 - address of the second vector
 * Return: The address of the newly allocated vector structure containing
 *         the sum.
 */
/* Functions referenced by add_vect() are:
 *  vectalloc()
 */
/* Functions that reference add_vect() are:
 *  build_offsets_int()
 *  CamberBisect()
 *  CamberBisect3d()
 *  CamberLineCB()
 *  camber_line_2d()
 *  camber_line_3d()
 *  cross_link()
 *  CylinderTo3D()
 *  DrawSurfPivot()
 *  eval_fnbc()
 *  eval_vert()
 *  funccurve()
 *  integral_build_offset()
 *  integral_sample_offset()
 *  normaldu()
 *  sample_offset()
 *  tol_sample_offsets()
 */

vector *add_vect (vector *v1, vector *v2)
{
  vector *v;
  double2 w1, w2;

  v = vectalloc ();

  w1 = v1->w;
  w2 = v2->w;
  v->x = v1->x/w1 + v2->x/w2;   /* convert from homogeneous */
  v->y = v1->y/w1 + v2->y/w2;   /* to cartesian coordinates */
  v->z = v1->z/w1 + v2->z/w2;   /* by dividing by y */
  v->w = 1.0;

  return (v);
}

/* Function: add_vec1()
 * Purpose: Add two vectors
 * Method: Place the sum of two vectors in a third vector.
 * Arguments:
 *  v1 - address of the first vector
 *  v2 - address of the second vector
 *  v - address of the vector to hold the differece
 */
/* Functions that reference add_vect1() are:
 *  axis_and_terminate()
 *  build_offset()
 *  camber_line_2d()
 *  camber_line_3d()
 *  check_cv()
 *  check_patch_u()
 *  cone_gen()
 *  construct_bspl()
 *  construct_jbspl()
 *  corner_offset()
 *  cylinder_gen()
 *  DrawSurfPivot()
 *  eval_fn()
 *  eval_fnbc()
 *  eval_vert()
 *  fair_knot()
 *  fair_per_knot()
 *  find_dir()
 *  fn_camber()
 *  fn_camber3()
 *  fun_camber()
 *  fun_camber3()
 *  geod()
 *  geod_par()
 *  g_camber()
 *  g_camber3()
 *  lead_chord()
 *  linear_interpolate_vector()
 *  matr_m_vect()
 *  merge_tol_edges()
 *  offnorm()
 *  offnorm2d()
 *  plane_gen()
 *  PowCurv_error()
 *  PowCurv_eval()
 *  PowSurf_error()
 *  PowSurf_eval()
 *  PowSurf_iso()
 *  raise_byone()
 *  sample_mod_surf()
 *  solve_bc()
 *  sphere_gen()
 *  torus_gen()
 */

void add_vect1 (vector *v1, vector *v2, vector *v)
{
  double2 w1, w2;

  w1 = v1->w;
  w2 = v2->w;
  v->x = v1->x/w1 + v2->x/w2;   /* convert from homogeneous to  */
  v->y = v1->y/w1 + v2->y/w2;   /* cartesian coordinates */
  v->z = v1->z/w1 + v2->z/w2;   /* by dividing by y */
  v->w = 1.0;
}

/* Function: cross()
 * Purpose: Cross product of two vectors
 * Method: Allocate a new vector structure to hold the cross product
 *         of two vectors
 * Arguments:
 *  v1 - address of the first vector
 *  v2 - address of the second vector
 * Return: The address of a newly allocated vector structure containing
 *         the cross product of the two input vectors.
 */
/* Functions referenced by cross() are:
 *  vectalloc()
 */
/* Functions that reference cross() are:
 *  cross_link()
 *  cyl_frenet_tr()
 *  DrawCurvEvaluate()
 *  DrawCurvTorsion()
 *  DrawSurfEvaluate()
 *  DrawSurfPivot()
 *  EvaluateCurv()
 *  EvaluateSurf()
 *  fcn()
 *  fgeodpar()
 *  fn()
 *  frenet_tr()
 *  funccurve()
 *  GetCurvValues()
 *  GetSurfValues()
 *  iso_curv0()
 *  iso_curv1()
 *  max_curve_curvature()
 *  normalcurv()
 *  normaldu()
 *  normalsurf()
 *  normalsurf1()
 *  offnorm2d()
 *  parallel()
 *  PostScriptCurvEvaluate()
 *  PostScriptCurvTorsion()
 *  PostScriptSurfEvaluate()
 *  pretransmat()
 *  rbspcur()
 *  SurfPivotUpdate()
 *  surf_normal()
 *  transmat()
 *  vectang()
 */

vector *cross (vector *v1, vector *v2)
{
  vector *v;
  double2 w1, w2, w12;

  v = vectalloc ();

  w1 = v1->w;
  w2 = v2->w;
  w12 = w1*w2;
  v->x = (v1->y*v2->z - v2->y*v1->z)/w12;   /* convert from homogeneous */
  v->y = (v1->z*v2->x - v2->z*v1->x)/w12;   /* to cartesian coordinates */
  v->z = (v1->x*v2->y - v2->x*v1->y)/w12;   /* by dividing by y */

  return (v);
}

/* Function: cross1()
 * Purpose: Cross product of two vectors
 * Method: Place the cross product of two vectors in a third vector
 * Arguments:
 *  v1 - address of the first vector
 *  v2 - address of the second vector
 *  v - address of the vector to contain the cross product
 */
/* Functions that reference cross1() are:
 *  axis_and_terminate()
 *  build_offset()
 *  CheckCurvPlanar()
 *  circle_gen()
 *  cone_gen()
 *  cross_link()
 *  curv_dev()
 *  cylinder_gen()
 *  fcn()
 *  fgeodpar()
 *  offnorm()
 *  PowSurf_error()
 */

void cross1 (vector *v1, vector *v2, vector *v)
{
  double2 xyzw[4], w1, w2, w12;

  w1 = v1->w;
  w2 = v2->w;
  w12 = w1*w2;
  xyzw[0] = (v1->y*v2->z - v2->y*v1->z)/w12;   /* convert from homogeneous */
  xyzw[1] = (v1->z*v2->x - v2->z*v1->x)/w12;   /* to cartesian coordinates */
  xyzw[2] = (v1->x*v2->y - v2->x*v1->y)/w12;   /* by dividing by y */
  xyzw[3] = 1.0;
  memcpy ((char *) v, (char *) xyzw, sizeof(vector));
}

/* Function: distance()
 * Purpose: Euclidean distance between two vectors
 * Method: Return the Euclidean distance of two vectors.
 * Arguments:
 *  v1 - address of the first vector
 *  v2 - address of the second vector
 * Return: The Euclidean (square root of the sum of the squared differences)
 *         of the two input vectors
 */
/* Functions that reference distance() are:
 *  calc_hartley()
 *  calc_hartley_per()
 *  calc_knots_chord()
 *  calc_knots_hartley()
 *  calc_knots_per()
 *  CamberBisect()
 *  CamberBisect3d()
 *  camber_line_2d()
 *  camber_line_3d()
 *  camber_surface()
 *  checkoff()
 *  check_surf_param()
 *  compute_diffeq()
 *  compute_ortho_proj()
 *  ConvexHullSurf()
 *  ConvexHullTest()
 *  convexhull_surf()
 *  convexhull_test()
 *  cross_link()
 *  find_deviation()
 *  find_dev_surf()
 *  global_pos_error()
 *  integral_sample_offset()
 *  main()
 *  MakeCamberLine()
 *  MakeCamberLine3d()
 *  MinDistance()
 *  MinDistCurvCB()
 *  position_err()
 *  PowCurv_error()
 *  PowParCurv_compare()
 *  PowParSurf_compare()
 *  PowSurf_error()
 *  sample()
 *  sample_mod_surf()
 *  sample_offset()
 *  sample_per()
 *  sample_proj()
 *  SubDistance()
 *  tol_sample_offsets()
 */

double2	distance (vector *v1, vector *v2)
{
  double2 rt, s, x, y, z, w1, w2;

  w1 = v1->w;
  w2 = v2->w;
  x = v1->x/w1 - v2->x/w2;   /* convert from homogeneous */
  y = v1->y/w1 - v2->y/w2;   /* to cartesian coordinates */
  z = v1->z/w1 - v2->z/w2;   /* by dividing by y */
  s = x*x + y*y + z*z;       /* sum of the squares */
  if (s != 0.0)
    rt = SQRT (s);           /* square root of the sum */
  else
    rt = s;   /****** sqrt(0.0) gives NaN ******/

  return (rt);
}

/* Function: distanc1e()
 * Purpose: Euclidean distance between two vectors
 * Method: Return the Euclidean distance of two vectors.  The second
 *         vector is specified by three scalar values, not by a vector
 *         structure.
 * Arguments:
 *  v1 - address of the first vector
 *  x - x component of the second vector
 *  y - y component of the second vector
 *  z - z component of the second vector
 * Return: The Euclidean (square root of the sum of the squared differences)
 *         of the two input vectors
 */

double distance1 (vector *v1, double x, double y, double z)
{
  double dx, dy, dz, rt, s;

  dx = v1->x/v1->w - x;        /* convert from homogeneous to */
  dy = v1->y/v1->w - y;        /* cartesian coordinates by */
  dz = v1->z/v1->w - z;        /* dividing by w */
  s = dx*dx + dy*dy + dz*dz;   /* sum of the squares */
  if (s != 0.0)
    rt = SQRT (s);             /* square root of the sum */
  else
    rt = s;   /****** sqrt(0.0) gives NaN ******/
  return (rt);
}

/* Function: dot()
 * Purpose: Dot product of two vectors
 * Arguments:
 *  v1 - address of the first vector
 *  v2 - address of the second vector
 * Return: Dot product (inner product) of the two vectors
 */
/* Functions that reference dot() are:
 *  axis_and_terminate()
 *  CalcIsophotes()
 *  CalcReflectionLines()
 *  calc_knots_foley()
 *  calc_par_foley()
 *  camber_line_2d()
 *  camber_line_3d()
 *  CheckCurvPlanar()
 *  compute_diffeq()
 *  compute_ortho_proj()
 *  con_confun()
 *  corner_offset()
 *  curvature()
 *  curvature1()
 *  cyl_frenet_tr()
 *  deriv_dev()
 *  deriv_dev_per()
 *  DrawCurvTorsion()
 *  EvaluateCurv()
 *  EvaluateSurf()
 *  fcn()
 *  fgeodpar()
 *  find_init()
 *  fn()
 *  footpoint_fcn()
 *  fundamentals()
 *  geod()
 *  geod_par()
 *  GetCurvValues()
 *  GetSurfValues()
 *  init_ortho_proj()
 *  init_proj()
 *  localize_diagnostic()
 *  localize_sumsq()
 *  localize_sumsq_2d()
 *  localize_sumsq_opt()
 *  main()
 *  MinDistance()
 *  MinDistCurvCB()
 *  normalcurv()
 *  normaldu()
 *  offnorm2d()
 *  PostScriptCurvTorsion()
 *  PowCurv_error()
 *  PowSurf_error()
 *  pretransmat()
 *  signed_distance()
 *  SubDistance()
 *  SurfPivotUpdate()
 *  transmat()
 *  vectang()
 *  vect_colinear()
 *  vec_angle()
 */

double2 dot (vector *v1, vector *v2)
{
     double2 dp, w12;

     w12 = v1->w*v2->w;
     dp = (v1->x*v2->x + v1->y*v2->y + v1->z*v2->z)/w12;

     return (dp);
}

/* Function: glue_vector()
 * Purpose: Initialize a vector
 * Method: Allocate a new vector structure and initialize its three
 *         scalar fields, x, y, and z.  Recall that vectalloc()
 *         initializes the w field to 1.0.
 * Arguments:
 *  x - the x component of the vector
 *  y - the y component of the vector
 *  z - the z component of the vector
 * Return: The address of the newly allocated vector structure that
 *         is initialized to the specified values.
 */
/* Functions referenced by glue_vector() are:
 *  vectalloc()
 */
/* Functions that reference glue_vector() are:
 *  funccurve()
 *  pretransmat()
 *  transmat()
 */

vector *glue_vector (double2 x, double2 y, double2 z)
{
  vector *v;

  v = vectalloc (); 

  v->x = x;
  v->y = y;
  v->z = z;

  return (v);
}

/* Function: glue_vector1()
 * Purpose: Initialize a vector
 * Method: Initialize its three scalar fields, x, y, and z, of a
 *         vector structure.  Recall that vectalloc() initializes
 *         the w field to 1.0.
 * Arguments:
 *  x - the x component of the vector
 *  y - the y component of the vector
 *  z - the z component of the vector
 *  v - address of the vector structure that is initialized to the
 *      specified values.
 */
/* Functions that reference glue_vector1() are:
 *  matr_m_vect()
 */

void glue_vector1 (double2 x, double2 y, double2 z, vector *v)
{
  v->x = x, v->y = y, v->z = z, v->w = 1.0;
}

/* Function: mag()
 * Purpose: Magnitude of vector
 * Arguments:
 *  v - address of the vector structure
 * Return: The magnitude (square root of the sum of the squares) of
 *         the vector
 */
/* Functions that reference mag() are:
 *  arc_length()
 *  calc_knots_foley()
 *  calc_par_foley()
 *  CheckCurvClosed()
 *  CheckCurvPlanar()
 *  CheckSurfClosedU()
 *  CheckSurfClosedV()
 *  check_cv()
 *  check_patch_u()
 *  circle_gen()
 *  cone_gen()
 *  curvature()
 *  curvature1()
 *  curv_cont()
 *  curv_dev()
 *  cylinder_gen()
 *  cyl_frenet_tr()
 *  cyl_sample_gencyl()
 *  deriv_dev()
 *  deriv_dev_per()
 *  DrawCurvTorsion()
 *  EvaluateCurv()
 *  EvaluateSurf()
 *  fcn()
 *  fgeodpar()
 *  find_dir()
 *  fn()
 *  GetCurvValues()
 *  GetSurfValues()
 *  int_sample_gencyl()
 *  iso_curv0()
 *  iso_curv1()
 *  knot_holzle()
 *  knot_holzle_fn()
 *  max_chord()
 *  max_curve_curvature()
 *  normalcurv()
 *  normalsurf()
 *  normalsurf1()
 *  offnorm2d()
 *  PostScriptCurvTorsion()
 *  PowCurv_error()
 *  pretransmat()
 *  rbspcur()
 *  sample_gencyl()
 *  SurfPivotUpdate()
 *  surf_normal()
 *  transmat()
 *  vectang()
 *  vector_dircos()
 *  vect_colinear()
 *  vec_angle()
 */

double2	mag (vector *v)
{
  double2 len, x, y, z;

  x = v->x/v->w;   /* convert from homogeneous to cartesian coordinates */
  y = v->y/v->w;   /* by dividing by w */
  z = v->z/v->w;
  len = SQRT (x*x + y*y + z*z);   /* square root of sum of squares */

  return (len);
}

/* Function: scale4()
 * Purpose: Scale a vector
 * Method: Scale all fields, x, y, and z, of the vector
 * Arguments:
 *  a - scale factor
 *  v - address of the vector structure
 */
/* Functions referenced by scale4() are:
 */
/* Functions that reference scale4() are:
 *  axis_and_terminate()
 */

void scale4 (double2 a, vector *v)
{
  v->x = a*v->x;
  v->y = a*v->y;
  v->z = a*v->z;
}

/* Function: scale_vect()
 * Purpose: Scale a vector
 * Method: Allocate a vector structure and fill it with a scaled
 *         copy of the fields, x, y, and z, of the input vector.
 *         The input vector is unchanged.
 * Arguments:
 *  a - scale factor
 *  v1 - address of the input vector
 * Return: The address of the newly allocated vector structure
 *         containing a scaled copy of the input vector
 */
/* Functions referenced by scale_vect() are:
 *  vectalloc()
 */
/* Functions that reference scale_vect() are:
 *  bias()
 *  bone_intersect()
 *  bone_intersect3()
 *  DrawSurfPivot()
 *  eval_fn()
 *  eval_fnbc()
 *  fn_camber()
 *  fn_camber3()
 *  funccurve()
 *  fun_camber()
 *  fun_camber3()
 *  g_camber()
 *  g_camber3()
 *  iso_curv0()
 *  iso_curv1()
 *  lead_chord()
 *  PowCurv_eval()
 *  PowSurf_eval()
 */

vector *scale_vect (double2 a, vector *v1)
{
  vector *v;

  v = vectalloc ();

  v->x = a*v1->x;
  v->y = a*v1->y;
  v->z = a*v1->z;
  v->w = v1->w;
  
  return (v);
}

/* Function: scale_vec1t()
 * Purpose: Scale a vector
 * Method: Scale the fields, x, y, and z, of the input vector.
 *         The input vector is unchanged.
 * Arguments:
 *  a - scale factor
 *  v1 - address of the input vector
 *  v - address of the vector  containing a scaled copy of the input vector
 */
/* Functions that reference scale_vect1() are:
 *  bias()
 *  build_offset()
 *  build_offsets_int()
 *  CamberBisect()
 *  CamberBisect3d()
 *  CamberLineCB()
 *  camber_line_2d()
 *  camber_line_3d()
 *  circle_gen()
 *  cone_gen()
 *  construct_bspl()
 *  construct_jbspl()
 *  convexhull_test()
 *  corner_offset()
 *  CylinderTo3D()
 *  cylinder_gen()
 *  cyl_frenet_tr()
 *  DrawSurfPivot()
 *  EvaluateCurv()
 *  eval_curve_bounded()
 *  eval_fn()
 *  eval_fnbc()
 *  eval_surface_bounded()
 *  eval_vert()
 *  fair_knot()
 *  fair_per_knot()
 *  fcn()
 *  fdist()
 *  fdist_curv()
 *  fgeodpar()
 *  find_dir()
 *  fn_camber()
 *  fn_camber3()
 *  geod()
 *  geod_par()
 *  GetCurvValues()
 *  hodograph_surf()
 *  integral_build_offset()
 *  integral_sample_offset()
 *  iso_curv0()
 *  iso_curv1()
 *  linear_interpolate_vector()
 *  matr_m_vect()
 *  merge_tol_edges()
 *  normalcurv()
 *  normaldu()
 *  offnorm()
 *  offnorm2d()
 *  ParCurv_hodograph()
 *  ParCurv_to_PowCurv()
 *  ParSurf_to_PowSurf()
 *  PowCurv_error()
 *  PowCurv_eval()
 *  PowSurf_error()
 *  PowSurf_eval()
 *  PowSurf_iso()
 *  raise_byone()
 *  sample_mod_surf()
 *  sample_offset()
 *  sphere_gen()
 *  tol_sample_offsets()
 *  torus_gen()
 */

void scale_vect1 (double2 a, vector *v1, vector *v)
{
  v->x = a*v1->x;
  v->y = a*v1->y;
  v->z = a*v1->z;
  v->w = v1->w;
}

/* Function: sethomogeq()
 * Purpose: Scale a vector homo
 * Method: Scale all four the fields, x, y, z, and w, of a vector.
 * Arguments:
 *  v - address of the vector
 *  homo - scale factor
 */
/* Functions that reference sethomogeq() are:
 *  compute_diffeq()
 *  compute_ortho_proj()
 *  cross_link()
 *  repeat_offset()
 *  SurfRev_to_ParSurf()
 */

void sethomogeq (vector *v, double2 homo)
{
  v->x *= homo;
  v->y *= homo;
  v->z *= homo;
  v->w *= homo;
}

/* Function: sub_vect()
 * Purpose: Subtract two vectors
 * Method: A new vector structure is allocated and holds the differance
 *         of the two vectors.
 * Arguments:
 *  v1 - address of the first vector
 *  v2 - address of the second vector
 * Return: The address of the newly allocated vector structure containing
 *         the difference.
 */
/* Functions referenced by sub_vect() are:
 *  vectalloc()
 */
/* Functions that reference sub_vect() are:
 *  bias()
 *  Bisector()
 *  build_offsets_int()
 *  CamberBisect()
 *  CamberBisect3d()
 *  camber_line_2d()
 *  camber_line_3d()
 *  check_cv()
 *  check_patch_u()
 *  circle_gen()
 *  cone_gen()
 *  cross_link()
 *  cylinder_gen()
 *  EvaluateCurv()
 *  eval_vert()
 *  fn_camber()
 *  fn_camber3()
 *  GetCurvValues()
 *  init_proj()
 *  iso_curv0()
 *  iso_curv1()
 *  main()
 *  max_chord()
 *  MinDistance()
 *  MinDistCurvCB()
 *  normalcurv()
 *  normaldu()
 *  offnorm2d()
 *  repeat_offset()
 *  SubDistance()
 *  tol_sample_offsets()
 *  torus_gen()
 */

vector *sub_vect (vector *v1, vector *v2)
{
  vector *v;
  double2 w1, w2;
  v = vectalloc ();
  w1 = v1->w, w2 = v2->w;
  v->x = v1->x/w1 - v2->x/w2, v->y = v1->y/w1 - v2->y/w2;
  v->z = v1->z/w1 - v2->z/w2;
  return (v);
}

/* Function: sub_vect1()
 * Purpose: Subtract two vectors
 * Method: Place the difference of two vectors in a third vector.
 * Arguments:
 *  v1 - address of the first vector
 *  v2 - address of the second vector
 *  v - address of the vector to hold the differece
 */
/* Functions that reference sub_vect1() are:
 *  arc_length()
 *  axis_and_terminate()
 *  bone_intersect()
 *  bone_intersect3()
 *  calc_knots_foley()
 *  calc_par_foley()
 *  CheckCurvClosed()
 *  CheckCurvPlanar()
 *  CheckSurfClosedU()
 *  CheckSurfClosedV()
 *  check_cv()
 *  check_patch_u()
 *  convexhull_test()
 *  corner_offset()
 *  cyl_frenet_tr()
 *  cyl_sample_gencyl()
 *  eval_curve_bounded()
 *  eval_surface_bounded()
 *  fair_knot()
 *  fair_per_knot()
 *  fdist()
 *  fdist_curv()
 *  fn_camber()
 *  fn_camber3()
 *  hodograph_surf()
 *  int_sample_gencyl()
 *  iso_curv0()
 *  iso_curv1()
 *  LinearCosNormal()
 *  localize_diagnostic()
 *  localize_sumsq()
 *  localize_sumsq_2d()
 *  localize_sumsq_opt()
 *  ParCurv_hodograph()
 *  plane_gen()
 *  sample_gencyl()
 *  solve_bc()
 *  trim_start_par()
 *  trim_start_par3()
 */

void sub_vect1 (vector *v1, vector *v2, vector *v)
{
  double2 w1, w2;
  w1 = v1->w, w2 = v2->w;
  v->x = v1->x/w1 - v2->x/w2, v->y = v1->y/w1 - v2->y/w2;
  v->z = v1->z/w1 - v2->z/w2, v->w = 1.0;
}

/* Function: unitvector()
 * Purpose: Normalize vector
 * Method: Allocate a new vector structure and fill it with a copy of the
 *         input vector normalized to unit length.
 * Arguments:
 *  v - address of the input vector
 * Return: The address of the newly allocated vector structure 
 *         containing a copy of the input vector normalized to unit
 *         length.
 */
/* Functions referenced by unitvector() are:
 *  copyvector()
 *  vectalloc()
 */
/* Functions that reference unitvector() are:
 *  cyl_frenet_tr()
 *  frenet_tr()
 *  funccurve()
 *  pretransmat()
 *  repeat_offset()
 *  transmat()
 */

vector *unitvector (vector *v)
{
  double2 length;
  vector *unit;

  unit = vectalloc ();

  /* calculate the length of the vector */
  length = SQRT (v->x*v->x + v->y*v->y + v->z*v->z) / v->w;
  if (FABS(length) < MACHPREC) {  /* if zero length, just copy */
    fputs ("input vector is zero - result not unit\n\r", stderr);  
    copyvector (v, unit);
  }
  else {                          /* normalize to unit length */
    unit->x = v->x/length;
    unit->y = v->y/length;
    unit->z = v->z/length;
  }

  return (unit);
}

/* Function: unitvector1()
 * Purpose: Normalize vector
 * Method: Normalize a copy of the input vector to unit length.
 * Arguments:
 *  v - address of the input vector
 *  unit - address of the vector containing a copy of the input vector
 *         normalized to unit length.
 */
/* Functions referenced by unitvector1() are:
 *  copyvector()
 */
/* Functions that reference unitvector1() are:
 *  axis_and_terminate()
 *  Bisector()
 *  build_offset()
 *  CalcReflectionLines()
 *  calc_disc_curv()
 *  camber_line_2d()
 *  camber_line_3d()
 *  circle_gen()
 *  cone_gen()
 *  corner_offset()
 *  cross_link()
 *  cylinder_gen()
 *  DrawFacetEvaluateMap()
 *  DrawSurfEvaluate()
 *  DrawSurfPivot()
 *  EvaluateCurv()
 *  EvaluateSurf()
 *  fcn()
 *  fn()
 *  GetCurvValues()
 *  GetSurfValues()
 *  init_proj()
 *  LinearCosNormal()
 *  normalcurv()
 *  normalsurf()
 *  normalsurf1()
 *  offnorm()
 *  offnorm2d()
 *  PostScriptFacetEvaluateMap()
 *  PostScriptSurfEvaluate()
 *  PowSurf_error()
 *  SurfPivotUpdate()
 *  surf_normal()
 *  torus_gen()
 *  trim_start_par()
 *  trim_start_par3()
 */

void unitvector1 (vector *v, vector *unit)
{
  double2 length;

  /* calculate the length of the vector */
  length = SQRT (v->x*v->x + v->y*v->y + v->z*v->z) / v->w;

  if (FABS(length) < MACHPREC) {   /* if zero length, just copy */
    fputs ("input vector is zero - result not unit\n\r", stderr);
    copyvector (v, unit);
  }
  else {                           /* normalize to unit length */
    unit->x = v->x/length;
    unit->y = v->y/length;
    unit->z = v->z/length;
    unit->w = 1.0;
  }
}

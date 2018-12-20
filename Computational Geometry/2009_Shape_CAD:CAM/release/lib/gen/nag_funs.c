/************************************************************************
 *									*
			    Copyright (C) 1990
     Massachusetts Institute of Technology, Cambridge, Massachusetts.
			   All rights reserved.

     This portion of source code was prepared in the Ocean Engineering
     Design Laboratory by Bradley A. Moran.
 *									*
 ************************************************************************/

/* These are actual pointers to user-supplied routines */

/* Many NAG routines require functions with the constant names
 * funct1_(), funct2_(), and hess2_()
 * Because the same function name is used by different NAG routines
 * it is necessary to have the actual function func1_() use a
 * function pointer to call the appropriate function
 */

/* pfunct1 - a function pointer to a function that will be called
 *           by funct1_()
 * pfunct2 - a function pointer to a function that will be called
 *           by funct2_()
 * phess2 - a function pointer to a function that will be called
 *           by hess2_()
 * pfsave1 - saved copy of the current value of pfunct1
 * pfsave2 - saved copy of the current value of pfunct2
 */

static void (*pfunct1) (int *, double *, double *);
static void (*pfunct2) (int *, double *, double *, double *);
static void (*phess2) (int *, double *, double *, int *, double *);
static void (*pfsave1) (int *, double *, double *);
static void (*pfsave2) (int *, double *, double *, double *);

/* Function: set_fun1_ptr()
 * Purpose: Set the function pointer used by func1_()
 * Arguments:
 *  userfun1 - function pointer
 */
/* Functions that reference set_fun1_ptr() are:
 *  ExtractLeadingEdge2D()
 *  ExtractLeadingEdge3D()
 *  find_dist_curv()
 *  MaxCamber()
 *  MaxCamber3d()
 *  maxcur()
 *  maxcur_csurf()
 *  max_camber()
 *  max_camber3()
 */

void set_fun1_ptr (void (*userfun1)(int *, double *, double *))
{
     pfunct1 = userfun1;
}

/* Function: set_fun2_ptr()
 * Purpose: Set the function pointer used by func2_()
 * Arguments:
 *  userfun2 - function pointer
 */
/* Functions that reference set_fun2_ptr() are:
 *  find_surf_dist()
 *  find_surf_dist1()
 *  find_surf_dist2()
 *  lead_chord()
 *  localize_sparse_2d()
 *  localize_sumsq_opt()
 */

void set_fun2_ptr (void (*userfun2)(int *, double *, double *, double *))

/* initializes the external pointer */

{
     pfunct2 = userfun2;
}

/* Function: save_fun1_ptr()
 * Purpose: Save the current value of the function pointer pfunct1
 */
/* Functions that reference save_fun1_ptr() are:
 *  find_main_trim_curve()
 */

void save_fun1_ptr(void)
{
  pfsave1 = pfunct1;
}

/* Function: save_fun2_ptr()
 * Purpose: Save the current value of the function pointer pfunct2
 */
/* Functions that reference save_fun2_ptr() are:
 *  footpoint()
 *  footpoint_curve()
 */

void save_fun2_ptr(void)
{
  pfsave2 = pfunct2;
}

/* Function: restore_fun1_ptr()
 * Purpose: Restore the saved value of the function pointer pfunct2
 */
/* Functions that reference restore_fun1_ptr() are:
 *  find_min_trim_curv()
 */

void restore_fun1_ptr(void)
{
  pfunct1 = pfsave1;
}

/* Function: restore_fun2_ptr()
 * Purpose: Restore the saved value of the function pointer pfunct2
 */
/* Functions that reference restore_fun2_ptr() are:
 *  footpoint()
 *  footpoint_curve()
 */

void restore_fun2_ptr(void)
{
  pfunct2 = pfsave2;
}

/* Function: set_hes2_ptr()
 * Purpose: Set the function pointer used by hess2_()
 * Arguments:
 *  userhes2 - function pointer
 */

void set_hes2_ptr (void (*userhes2)(int *, double *, double *, int *,
				    double *))
{
     phess2 = userhes2;
}

/* Function: funct1_()
 * Purpose: Function called by NAG routines
 * Method: Use the function pointer to invoke the appropriate function.
 * Arguments: The arguments are defined by the invoking NAG routine
 */

void funct1_ (int *niv, double *xiv, double *fx)
{
     (*pfunct1) (niv, xiv, fx);
}

/* Function: funct2_()
 * Purpose: Function called by NAG routines
 * Method: Use the function pointer to invoke the appropriate function.
 * Arguments: The arguments are defined by the invoking NAG routine
 */

void funct2_ (int *niv, double *xiv, double *fx, double *gx)
{
     (*pfunct2) (niv, xiv, fx, gx);
}

/* Function: hess2_()
 * Purpose: Function called by NAG routines
 * Method: Use the function pointer to invoke the appropriate function.
 * Arguments: The arguments are defined by the invoking NAG routine
 */

void hess2_ (int *niv, double *hxiv, double *hlc, int *lh, double *hdc)
{
     (*phess2) (niv, hxiv, hlc, lh, hdc);
}

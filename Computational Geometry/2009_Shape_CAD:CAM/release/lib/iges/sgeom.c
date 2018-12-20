/* Copyright (C) 1992 Massachusetts Institute of Technology, Cambridge, MA
   All rights reserved
*/

/* sgeom.c */
/* based on: /usr/deslab/epraxiteles/lib/gen/sgeom.c */

/* free_sgeom()
 * PowSurf_eval()
 * sgeomalloc()
*/

#include <malloc.h>
#include "iges.h"
#include "bspl.h"

/********* free_sgeom() *********
* 1     Purpose
* 
*       Deallocate parametric spline surface structure.
* 
* 2     Specification
* 
*       #include "iges.h"
*       void free_sgeom(PowSurf *sgeom);
* 
* 3     Description
* 
*       This function deallocates an parametric spline surface structure.
* 
* 4     References
* 
*       [1]   Digital Representation for Communication of Product Definition
*             Data, US PRO/IPO-100, Initial Graphics Exchange Specification
* 	    (IGES) 5.2, IGES/PDES Organization, U.S. Product Data Association,
* 	    Fairfax, VA, November 1993.
* 
* 5     Parameters
* 
*           1.  PowSurf * sgeom
*               On entry:  the address of the structure.
* 
* 8     Further Comments
* 
*       The structure should be allocated with sgeomalloc().
* 
* Functions referenced by free_sgeom() are:
*  free_darray1()
*  free_garray1()
*  free_varray2()
*
* Functions that reference free_sgeom() are:
*  ReadIgesPowerSurf()
*  SaveIgesSurf()
*/

void free_sgeom(PowSurf *sgeom)
{
  int i;

/************************************************************************
        free_sgeom() frees the memory associated with a 
	PowSurf structure
 ************************************************************************/
  
  if (sgeom->upmem && sgeom->vpmem)
    for (i=0; i<16; i++)
      free_varray2(sgeom->contpts[i], sgeom->usegmts, sgeom->vsegmts);
  
  if (sgeom->ukmem)
    free_darray1(sgeom->uknots);
  if (sgeom->vkmem)
    free_darray1(sgeom->vknots);
  
  free_garray1((char *)sgeom);
}

/********* PowSurf_eval() *********
* 1     Purpose
* 
*       Evaluate parametric spline surface.
* 
* 2     Specification
* 
*       #include "iges.h"
*       vector *PowSurf_eval(PowCuv *sgeom, double u, double v, int uderiv,
*                            int vderiv);
* 
* 3     Description
* 
*       This function evaluates the uderivth ; vderivth  derivative of a
*       parametric spline surface.
* 
* 4     References
* 
*       [1]   Digital Representation for Communication of Product Definition
*             Data, US PRO/IPO-100, Initial Graphics Exchange Specification
* 	    (IGES) 5.2, IGES/PDES Organization, U.S. Product Data Association,
* 	    Fairfax, VA, November 1993.
* 
* 5     Parameters
* 
*           1.  PowSurf * sgeom
*               On entry:  the address of a structure containing the parametric
* 	      spline surface.
* 
*           2.  double u
*               On entry:  the u parameter value at which the surface is
* 	      evaluated.
* 
*           3.  double v
*               On entry:  the v parameter value at which the surface is
* 	      evaluated.
* 
*           4.  int uderiv
*               On entry:  flag specifying which partial derivative in the u
* 	      direction to evaluate:
* 	      if 0, then return the zero-th derivative (or position);
* 	      if 1, return the first partial derivative; and so on.
* 
*           5.  int vderiv
*               On entry:  flag specifying which partial derivative in the v
* 	      direction to evaluate:
* 	      if 0, then return the zero-th derivative (or position);
* 	      if 1, return the first partial derivative; and so on.
* 
* 6     Return Values, Error Indicators and Warnings
* 
*       The  address  of  a  vector  structure  containing  the  evaluated
*       derivative  of  the  parametric spline surface.
* 
* 8     Further Comments
* 
*       The returned vector structure is newly allocated by this function.  It
*       should be explicitly deallocated using vectfree().
* 
* Functions referenced by PowSurf_eval() are:
*  add_vect1()
*  find()
*  scale_vect()
*  scale_vect1()
*
* Functions that reference PowSurf_eval() are:
*  calc_points_surf()
*  PowParSurf_compare()
*/

vector *PowSurf_eval(PowSurf *sgeom, double u, double v, int uderiv,
		     int vderiv)
{
  vector vt0, vt1, vt2, vt3, *val;
  double s, t;
  int j, k;

/************************************************************************
         PowSurf_eval() returns a pointer to a vector containing
	 the ideriv(th) derivative of the power basis Surface, sgeom,
	 evaluated at the parametric values t and w.
 ************************************************************************/

  j = find(sgeom->usegmts + 1, sgeom->uknots, u);
  if (j >= sgeom->usegmts)
    j = sgeom->usegmts - 1;

  k = find (sgeom->vsegmts + 1, sgeom->vknots, v);
  if (k >= sgeom->vsegmts)
    k = sgeom->vsegmts - 1;

  s = u - sgeom->uknots[j];
  t = v - sgeom->vknots[k];

  scale_vect1(s, sgeom->contpts[15][j][k], &vt3);
  add_vect1(&vt3, sgeom->contpts[14][j][k], &vt3);
  scale_vect1(s, &vt3, &vt3);
  add_vect1(&vt3, sgeom->contpts[13][j][k], &vt3);
  scale_vect1(s, &vt3, &vt3);
  add_vect1(&vt3, sgeom->contpts[12][j][k], &vt3);

  scale_vect1(s, sgeom->contpts[11][j][k], &vt2);
  add_vect1(&vt2, sgeom->contpts[10][j][k], &vt2);
  scale_vect1(s, &vt2, &vt2);
  add_vect1(&vt2, sgeom->contpts[9][j][k], &vt2);
  scale_vect1(s, &vt2, &vt2);
  add_vect1(&vt2, sgeom->contpts[8][j][k], &vt2);

  scale_vect1(s, sgeom->contpts[7][j][k], &vt1);
  add_vect1(&vt1, sgeom->contpts[6][j][k], &vt1);
  scale_vect1(s, &vt1, &vt1);
  add_vect1(&vt1, sgeom->contpts[5][j][k], &vt1);
  scale_vect1(s, &vt1, &vt1);
  add_vect1(&vt1, sgeom->contpts[4][j][k], &vt1);

  scale_vect1(s, sgeom->contpts[3][j][k], &vt0);
  add_vect1(&vt0, sgeom->contpts[2][j][k], &vt0);
  scale_vect1(s, &vt0, &vt0);
  add_vect1(&vt0, sgeom->contpts[1][j][k], &vt0);
  scale_vect1(s, &vt0, &vt0);
  add_vect1(&vt0, sgeom->contpts[0][j][k], &vt0);

  val = scale_vect (t, &vt3);
  add_vect1(&vt2, val, val);
  scale_vect1(t, val, val);
  add_vect1(&vt1, val, val);
  scale_vect1(t, val, val);
  add_vect1(&vt0, val, val);

  return (val);
}

/********* sgeomalloc() *********
* 1     Purpose
* 
*       Allocate a parametric spline surface structure.
* 
* 2     Specification
* 
*       #include "iges.h"
*       PowSurf * sgeomalloc(int usegmts, int vsegmts);
* 
* 3     Description
* 
*       This function allocates a data structure that will contain a parametric
*       spline surface.
* 
* 4     References
* 
*       [1]   Digital Representation for Communication of Product Definition
*             Data, US PRO/IPO-100, Initial Graphics Exchange Specification
* 	    (IGES) 5.2, IGES/PDES Organization, U.S. Product Data Association,
* 	    Fairfax, VA, November 1993.
* 
* 5     Parameters
* 
*           1.  int usegmts
*               On entry:  the number of polynomial segments in the u direction
* 	      of the spline.
* 
*           2.  int vsegmts
*               On entry:  the number of polynomial segments in the v direction
* 	      of the spline.
* 
* 6     Return Values, Error Indicators and Warnings
* 
*       The address of the structure is returned.  If insufficient dynamic
*       memory is available, the value NULL is returned.
* 
* 8     Further Comments
* 
*       The structure should be deallocated with free_sgeom().
* 
* Functions referenced by sgeomalloc() are:
*  dbl_array1()
*  errormsg()
*  gen_array1()
*  vec_array2()
*
* Functions that reference sgeomalloc() are:
*  ParSurf_to_PowSurf()
*  ReadIges114()
*/

PowSurf *sgeomalloc(int usegmts, int vsegmts)
{
  PowSurf *sgeom;
  int i;

/************************************************************************
        sgeomalloc() returns a pointer to a PowSurf structure 
	dimensioned by usegmts and vsegmts
 ************************************************************************/

  if (!(sgeom = (PowSurf *)gen_array1(1, sizeof(PowSurf))))
    errormsg (0, "allocation failure in sgeomalloc()");

  sgeom->type = PSurfacePow;
  sgeom->uknots = dbl_array1(usegmts+1);
  sgeom->vknots = dbl_array1(vsegmts+1);

  for (i=0; i<16; i++)
    sgeom->contpts[i] = vec_array2(usegmts, vsegmts);
  
  sgeom->usegmts = usegmts;
  sgeom->vsegmts = vsegmts;
  sgeom->ukmem = usegmts + 1;
  sgeom->upmem = usegmts;
  sgeom->vkmem = vsegmts + 1;
  sgeom->vpmem = vsegmts;
  
  return (sgeom);
}

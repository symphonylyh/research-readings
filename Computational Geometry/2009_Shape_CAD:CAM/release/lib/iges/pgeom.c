/* Copyright (C) 1992 Massachusetts Institute of Technology, Cambridge, MA
   All rights reserved
*/

/* pgeom.c */
/* based on: /usr/deslab/epraxiteles/lib/gen/pgeom.c */

/* free_pgeom()
 * pgeomalloc()
 * PowCurv_eval()
*/

#include <stdio.h>
#include <malloc.h>
#include "iges.h"
#include "bspl.h"

/********* free_pgeom() *********
* 1     Purpose
* 
*       Deallocate parametric spline curve structure.
* 
* 2     Specification
* 
*       #include "iges.h"
*       void free_pgeom(PowCurv *pgeom);
* 
* 3     Description
* 
*       This function deallocates an parametric spline curve structure.
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
*           1.  PowCurv * pgeom
*               On entry:  the address of the structure.
* 
* 8     Further Comments
* 
*       The structure should be allocated with pgeomalloc().
* 
* Functions referenced by free_pgeom() are:
*  free_darray1()
*  free_garray1()
*  free_varray2()
*
* Functions that reference free_pgeom() are:
*  PowSurf_to_ParSurf_loft()
*  ReadIgesPowerCurv()
*  SaveIgesCos()
*  SaveIgesCurv()
*  SaveIgesSurf()
*/

void free_pgeom(PowCurv *pgeom)
{
/************************************************************************
        free_pgeom() frees the memory associated with a 
	PowCurv structure
 ************************************************************************/

  if (pgeom->pmem)
    free_varray2(pgeom->contpts, pgeom->nsegmts, 4);
  
  if (pgeom->kmem)
    free_darray1(pgeom->knots);
  
  free_garray1((char *)pgeom);
}

/********* pgeomalloc() *********
* 1     Purpose
* 
*       Allocate a parametric spline curve structure.
* 
* 2     Specification
* 
*       #include "iges.h"
*       PowCurv * pgeomallc(int nsegmts);
* 
* 3     Description
* 
*       This function allocates a data structure that will contain a parametric
*       spline curve.
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
*           1.  int nsegmts
*               On entry:  the number of polynomial segments in the spline.
* 
* 6     Return Values, Error Indicators and Warnings
* 
*       The address of the structure is returned.  If insufficient dynamic
*       memory is available, the value NULL is returned.
* 
* 8     Further Comments
* 
*       The structure should be deallocated with free_pgeom().
* 
* Functions referenced by pgeomalloc() are:
*  dbl_array1()
*  errormsg()
*  gen_array1()
*  vec_array2()
*
* Functions that reference pgeomalloc() are:
*  ParCurv_to_PowCurv()
*  PowSurf_iso()
*  PowSurf_to_ParSurf_loft()
*  ReadIges112()
*/

PowCurv *pgeomalloc(int nsegmts)
{
  PowCurv *pgeom;

/************************************************************************
        pgeomalloc() returns a pointer to a PowCurv structure 
	dimensioned by nsegmts
 ************************************************************************/
  
  if (!(pgeom = (PowCurv *)gen_array1(1, sizeof(PowCurv))))
    errormsg (0, "allocation failure in pgeomalloc()");

  pgeom->type = PCurvePow;
  pgeom->knots = dbl_array1(1 + nsegmts);
  pgeom->contpts = vec_array2(nsegmts, 4);
  
  pgeom->nsegmts = nsegmts;
  pgeom->kmem = 1 + nsegmts;
  pgeom->pmem = nsegmts;
  
  return (pgeom);
}

/********* PowCurv_eval() *********
* 1     Purpose
* 
*       Evaluate parametric spline curve.
* 
* 2     Specification
* 
*       #include "iges.h"
*       vector *PowCurv_eval(PowCuv *pgeom, double t, int ideriv);
* 
* 3     Description
* 
*       This function evaluates the iderivth  derivative of a parametric
*       spline curve.
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
*           1.  PowCurv * pgeom
*               On entry:  the address of a structure containing the parametric
* 	      spline curve.
* 
*           2.  double t
*               On entry:  the parameter value at which the curve is evaluated.
* 
*           3.  int ideriv
*               On entry:  flag specifying which derivative to evaluate:
* 	      if 0, then return the zero-th derivative (or position);
* 	      if 1, return the first derivative; and so on.
* 
* 6     Return Values, Error Indicators and Warnings
* 
*       The  address  of  a  vector  structure  containing  the  evaluated
*       derivative  of  the  parametric spline curve.
* 
* 8     Further Comments
* 
*       The returned vector structure is newly allocated by this function.
*       It should be explicitly deallocated using vectfree().
* 
* Functions referenced by PowCurv_eval() are:
*  add_vect1()
*  find()
*  scale_vect()
*  scale_vect1()
*
* Functions that reference PowCurv_eval() are:
*  calc_points_curv()
*  PowParCurv_compare()
*  WriteIges112()
*/

vector *PowCurv_eval(PowCurv *pgeom, double t, int ideriv)
{
  vector *v;
  double s;
  int j;

/************************************************************************
         PowCurv_eval() returns a pointer to a vector containing
	 the ideriv(th) derivativeof the power basis curve, pgeom,
	 evaluated ate the parametric value t.
 ************************************************************************/

  j = find(pgeom->nsegmts + 1, pgeom->knots, t);
  if (j >= pgeom->nsegmts)
    j = pgeom->nsegmts - 1;

  s = t - pgeom->knots[j];
  v = scale_vect(s, pgeom->contpts[j][3]);
  add_vect1(v, pgeom->contpts[j][2], v);
  scale_vect1(s, v, v);
  add_vect1(v, pgeom->contpts[j][1], v);
  scale_vect1(s, v, v);
  add_vect1(v, pgeom->contpts[j][0], v);

  return(v);
}

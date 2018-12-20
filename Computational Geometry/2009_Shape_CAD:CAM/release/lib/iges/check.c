/* Copyright (C) 1992 Massachusetts Institute of Technology, Cambridge, MA
   All rights reserved
*/

/* check.c */

/* CheckCurv2d()
 * CheckCurvClosed()
 * CheckCurvIntegral()
 * CheckCurvPlanar()
 * CheckIgesFile()
 * CheckSurfClosedU()
 * CheckSurfCloseV()
 * CheckSurfIntegral()
*/

#include <stdio.h>
#include <math.h>
#include "iges.h"

/********* CheckCurv2d() ********* 1 Purpose
* 
*       Check if NURBS curve is 2d.
* 
* 2     Specification
* 
*       #include "iges.h"
*       int CheckCurv2d(ParCurv *egeom);
* 
* 3     Description
* 
*       This function determines if the input NURBS curve is 2d.
* 
* 5     Parameters
* 
*           1.  ParCurv * egeom
*               On entry:  the pointer to a structure containing NURBS curve.
* 
* 6     Return Values, Error Indicators and Warnings
* 
*       If the curve is 2d, the value 1 is returned; otherwise, 0.
* 
* 8     Further Comments
* 
*       The curve is considered 2d if all of the z components of the
*       control points are in the range 0.0+/-10.0E-10.
* 
* Functions that reference CheckCurv2d() are:
* ReadDeslabCurv()
*/

int CheckCurv2d(ParCurv *egeom)
{
  int i;

  for (i=0; i<egeom->ncontpts; i++)                    /* if all z's are 0 */
    if (fabs(egeom->contpts[i]->z) > PRAX_IGES_ZERO)   /* curve is 2d */
      return (0);
  return (1);
}

/******** CheckCurvClosed() *********
* 1     Purpose
* 
*       Check if NURBS curve is closed.
* 
* 2     Specification
* 
*       #include "iges.h"
*       int CheckCurvClosed(ParCurv *egeom);
* 
* 3     Description
* 
*       This function determines if the input NURBS curve is closed (periodic).
* 
* 5     Parameters
* 
*           1.  ParCurv * egeom
*               On entry:  the pointer to a structure containing NURBS curve.
* 
* 6     Return Values, Error Indicators and Warnings
* 
*       If the curve is closed, the value 1 is returned; otherwise, 0.
* 
* 8     Further Comments
* 
*       The curve is considered closed if distance between the first and last
*       control points is less than 10.0-10.
* 
* Functions referenced by CheckCurvClosed() are:
*  mag()
*  sub_vect1()
*
* Functions that reference CheckCurvClosed() are:
*  WriteIges126()
*/

int CheckCurvClosed(ParCurv *egeom)
{
  int i;
  vector du;

  /* if starting and ending positions are equal, curve is closed */
  sub_vect1(egeom->contpts[0], egeom->contpts[egeom->ncontpts-1], &du);
  if (mag(&du) >= PRAX_IGES_ZERO)
    return (0);

  return (1);
}

/********* CheckCurvIntegral() *********
* 1     Purpose
* 
*       Check if NURBS curve is integral.
* 
* 2     Specification
* 
*       #include "iges.h"
*       int CheckCurvIntegral(ParCurv *egeom);
* 
* 3     Description
* 
*       This function determines if the input NURBS curve is integral.
* 
* 5     Parameters
* 
*           1.  ParCurv * egeom
*               On entry:  the pointer to a structure containing NURBS curve.
* 
* 6     Return Values, Error Indicators and Warnings
* 
*       If the curve is integral, the value 1 is returned; otherwise, 0.
* 
* 8     Further Comments
* 
*       The curve is considered integral if all of the w components of the
*       control points are in the range 1.0+/-10.0E-10.
* 
* Functions that reference CheckCurvIntegral() are:
*  WriteIges126()
*/

int CheckCurvIntegral(ParCurv *egeom)
{
  int i;

  for (i=0; i<egeom->ncontpts; i++)          /* if all w's are 1, then curve */
    if (fabs(egeom->contpts[i]->w - 1.0) > PRAX_IGES_ZERO)   /*  is integral */
      return (0);
  return (1);
}

/********* CheckCurvPlanar() *********
* 1     Purpose
* 
*       Check if NURBS curve is planar.
* 
* 2     Specification
* 
*       #include "iges.h"
*       int CheckCurvPlanar(ParCurv *egeom, double *x, double *y, double *z);
* 
* 3     Description
* 
*       This function determines if the input NURBS curve is planar.
* 
* 5     Parameters
* 
*           1.  ParCurv * egeom
*               On entry:  the pointer to a structure containing NURBS curve.
* 
*           2.  double * x
*               On exit:  if the curve is planar, the address of the x
* 	      component of the unit normal vector of the plane in which
* 	      the curve is defined; otherwise, undefined.
* 
*           3.  double * y
*               On exit:  if the curve is planar, the address of the y
* 	      component of the unit normal vector of the plane in which
* 	      the curve is defined; otherwise, undefined.
* 
*           4.  double * z
*               On exit:  if the curve is planar, the address of the z
* 	      component of the unit normal vector of the plane in which
* 	      the curve is defined; otherwise, undefined.
* 
* 6     Return Values, Error Indicators and Warnings
* 
*       If the curve is planar, the value 1 is returned; otherwise, 0.
* 
* 8     Further Comments
* 
*       The  curve is considered planar if cosine of the angle between a
*       segment of the control polygon and the normal of a surface
*       containing all other segments of the control polyon is
*       in the range of 1.0+/-10.0E-8 radians.
* 
* Functions referenced by CheckCurvPlanar() are:
*  cross1()
*  dot()
*  mag()
*  sub_vect1()
*
* Functions that reference CheckCurvPlanar() are:
*  WriteIges126()
*/

int CheckCurvPlanar(ParCurv *egeom, double *x, double *y, double *z)
{
  vector v1, v2, norm;
  double m;
  int i, r_val = 1;

  for (i=1; i<egeom->ncontpts-2; i++) {
    sub_vect1(egeom->contpts[i], egeom->contpts[0], &v1);
    sub_vect1(egeom->contpts[i+1], egeom->contpts[0], &v2);
    cross1 (&v1, &v2, &norm);
    if ((m = mag (&norm)) > 0.001) {   /* normal is the cross product */
      norm.w = m;                      /* of adjacent segments of the */
      *x = norm.x / norm.w;            /* control polygon */
      *y = norm.y / norm.w;
      *z = norm.z / norm.w;	       
      break;
    }
  }
  for (i=1; i<egeom->ncontpts; i++) {
    sub_vect1 (egeom->contpts[i], egeom->contpts[0], &v1);
    if (fabs(dot(&v1, &norm)) > 1.0e-08) {  /* if all segments of control */
      r_val = 0;                            /* polygon are perpindicular to */
      break;                                /* normal, curve is planar */
    }
  }
  return (r_val);
}

/********* CheckIgesFile() *********
*1     Purpose
*
*      Check if file is in IGES format.
*
*2     Specification
*
*      #include "iges.h"
*      int CheckIgesFile(FILE *fp, char *record);
*
*3     Description
*
*      This function examines the Start section to determine if the file
*      is an IGES file, and if so, if it is in the normal (non-compressed,
*      non-binary) format.
*
*4     References
*
*      [1]   Digital Representation for Communication of Product
*            Definition Data, US PRO/IPO-100, Initial Graphics Exchange
*	    Specification (IGES) 5.2, IGES/PDES Organization, U.S.
*	    Product Data Association, Fairfax, VA, November 1993.
*
*5     Parameters
*
*          1.  FILE * fp
*              On entry:  file pointer of the open file.
*
*          2.  char * record
*              On entry:  the address of a character string containing
*	      the first 80 character line of the file.
*
*6     Return Values, Error Indicators and Warnings
*
*      This function returns one of the following values:
*
*   Value Mmemonic                   Interpretation
*   ----- ------------------------   -----------------------------------------
*     0   _                          IGES in non-compressed, non-binary format
*    -1   PRAX_ERR_IGES_BINARY       IGES binrary format
*    -2   PRAX_ERR_IGES_COMPRESSED   IGES compressed format
*    -3   PRAX_ERR_NOT_IGES          Not an IGES format
*
* Functions referenced by CheckIgesFile() are:
*  ReadNextIgesRecord()
*
* Functions that reference CheckIgesFile() are:
*  OpenIgesFileCB()
*/

int CheckIgesFile(FILE *fp, char *record)
{
  char code;

  if ((code = ReadNextIgesRecord(fp, record, 0)) != EOF) {
    if (code == PRAX_IGES_BINARY)
      return (PRAX_ERR_IGES_BINARY);
    else if (code == PRAX_IGES_COMPRESSED)
      return (PRAX_ERR_IGES_COMPRESSED);
    else if (code != PRAX_IGES_START)
      return (PRAX_ERR_NOT_IGES);
    else
      return (0);
  }
  else
    return (PRAX_ERR_EOF);
}

/********* CheckSurfClosedU() *********
* 1     Purpose
* 
*       Check if NURBS surface is closed in the u direction.
* 
* 2     Specification
* 
*       #include "iges.h"
*       int CheckSurfClosedU(ParSurf *fgeom);
* 
* 3     Description
* 
*       This function determines if the input NURBS surface is closed
*       (periodic) in the u parametric direction.
* 
* 5     Parameters
* 
*           1.  ParSurf * fgeom
*               On entry:  the pointer to a structure containing NURBS surface.
* 
* 6     Return Values, Error Indicators and Warnings
* 
*       If the surface is closed, the value 1 is returned; otherwise, 0.
* 
* 8     Further Comments
* 
*       The surface is considered closed if the distances between analogous
*       members of the first and last row of control points is less than
*       10.0E-10.
* 
* Functions referenced by CheckSurfClosedU() are:
*  mag()
*  sub_vect1()
*
* Functions that reference CheckSurfClosedU() are:
*  WriteIges128()
*/

int CheckSurfClosedU(ParSurf *fgeom)
{
  int i;
  vector du;

  for(i=0; i<fgeom->vcontpts; i++) {
    sub_vect1(fgeom->contpts[0][i], fgeom->contpts[fgeom->ucontpts-1][i], &du);
    if (mag(&du) >= PRAX_IGES_ZERO) /* if u=0 & u=1 are position continuous */
      return (0);                   /* surface is closed in u */
  }
  return (1);
}

/********* CheckSurfClosedV() *********
* 1     Purpose
* 
*       Check if NURBS surface is closed in the v direction.
* 
* 2     Specification
* 
*       #include "iges.h"
*       int CheckSurfClosedV(ParSurf *fgeom);
* 
* 3     Description
* 
*       This function determines if the input NURBS surface is closed
*       (periodic) in the v parametric direction.
* 
* 5     Parameters
* 
*           1.  ParSurf * fgeom
*               On entry:  the pointer to a structure containing NURBS surface.
* 
* 6     Return Values, Error Indicators and Warnings
* 
*       If the surface is closed, the value 1 is returned; otherwise, 0.
* 
* 8     Further Comments
* 
*       The surface is considered closed if the distances between analogous
*       members of the first and last column of control points is less than
*       10.0E-10.
* 
* Functions referenced by CheckSurfClosedV() are:
*  mag()
*  sub_vect1()
*
* Functions that reference CheckSurfClosedV() are:
*  WriteIges128()
*/

int CheckSurfClosedV(ParSurf *fgeom)
{
  int i;
  vector dv;

  for(i=0; i<fgeom->ucontpts; i++) {
    sub_vect1(fgeom->contpts[0][i], fgeom->contpts[i][fgeom->vcontpts-1], &dv);
    if (mag(&dv) >= PRAX_IGES_ZERO) /* if v=0 & v=1 are position continuous */
      return (0);                   /* surface is closed in v */
  }
  return(1);
}

/********* CheckSurfIntegral() *********
* 1     Purpose
* 
*       Check if NURBS surface is integral.
* 
* 2     Specification
* 
*       #include "iges.h"
*       int CheckSurfIntegral(ParSurf *fgeom);
* 
* 3     Description
* 
*       This function determines if the input NURBS surface is integral.
* 
* 5     Parameters
* 
*           1.  ParSurf * fgeom
*               On entry:  the pointer to a structure containing NURBS surface.
* 
* 6     Return Values, Error Indicators and Warnings
* 
*       If the surface is integral, the value 1 is returned; otherwise, 0.
* 
* 8     Further Comments
* 
*       The surface is considered integral if all of the z components of the
*       control points are in the range 1.0+/-10.0E-10.
* 
* Functions that reference CheckSurfIntegral() are:
*  WriteIges128()
*/

int CheckSurfIntegral(ParSurf *fgeom)
{
  int i, j;

  for (i=0; i<fgeom->ucontpts; i++)     /* if all w's are 0, the surface */
    for (j=0; j<fgeom->vcontpts; j++)   /* if integral */
      if (fabs(fgeom->contpts[i][j]->w - 1.0) > PRAX_IGES_ZERO)
	return (0);
  return (1);
}

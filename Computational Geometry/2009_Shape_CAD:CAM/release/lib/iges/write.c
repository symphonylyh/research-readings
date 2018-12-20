/* Copyright (C) 1992 Massachusetts Institute of Technology, Cambridge, MA
   All rights reserved
*/

/* write.c */

/* WriteIges106()
 * WriteIges110()
 * WriteIges112()
 * WriteIges114()
 * WriteIges118()
 * WriteIges120()
 * WriteIges122()
 * WriteIges126()
 * WriteIges128()
 * WriteIges142()
 * WriteIges144()
*/

#include <stdio.h>
#include <malloc.h>
#include "iges.h"

/********* WriteIges106() *********
* 1     Purpose
* 
*       Write an IGES Type106 (Copious Data) entity.
* 
* 2     Specification
* 
*       #include "iges.h"
*       int  WriteIges106(FILE *fp, Type106 *copious, char *eop, char *eor,
*                         int iDirectory, int *iParameter);
* 
* 3     Description
* 
*       This function writes a Type106 (Copious Data) entity to an IGES file.
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
*           1.  FILE * fp
*               On entry:  the file structure of the IGES file.
* 
* 
*           2.  Type106 * copious
*               On entry:  the address of a structure containing the Type106
* 	      (Copious Data) entity.
* 
*           3.  char eop
*               On entry:  the parameter delimiter character.
* 
*           4.  char eor
*               On entry:  the record delimiter character.
* 
*           5.  int iDirectory
*               On entry:  the line number of the first line of this entity's
* 	      Directory entry.
* 
*           6.  int * iParameter
*               On  entry:  the address of a variable containing the current
* 	      line number,  before the write.
* 
*               On exit:  the address of a variable containing the current line
* 	      number, after the write.
* 
* 6     Return Values, Error Indicators and Warnings
* 
*       The number of lines used to write the entity to the file is returned.
* 
* Functions referenced by WriteIges106() are:
*  AddNextInteger()
*  AddNextReal()
*
* Functions that reference WriteIges106() are:
*  SaveIgesCos()
*  SaveIgesCurv()
*  SaveIgesGrid()
*  SaveIgesList()
*  SaveIgesMinDist()
*  SaveIgesSurf()
*  SaveIgesUv()
*/

int WriteIges106(FILE *fp, Type106 *copious, char eop, char eor,
		   int iDirectory, int *iParameter)
{
  char record[73];
  int i, iCount, iLast;

  iCount = *iParameter;
  record[0] = '\0';
  AddNextInteger(fp, 106, eop, record, iDirectory, PRAX_IGES_PARAMETER,
		 iParameter, 0);
  AddNextInteger(fp, copious->ip, eop, record, iDirectory, PRAX_IGES_PARAMETER,
		 iParameter, 0);
  AddNextInteger(fp, copious->n, eop, record, iDirectory, PRAX_IGES_PARAMETER,
		 iParameter, 0);
  switch (copious->ip) {
  case 1:
  case 11:
    AddNextReal(fp, copious->pts.pair.zt, eop, record, iDirectory,
		PRAX_IGES_PARAMETER, iParameter, 0);
    for (i=0; i<copious->n; i++) {
      iLast = (i < copious->n - 1 ? 0 : 1);
      AddNextReal(fp, copious->pts.pair.x[i], eop, record, iDirectory,
		  PRAX_IGES_PARAMETER, iParameter, 0);
      AddNextReal(fp, copious->pts.pair.y[i], (iLast ? eor : eop), record,
		  iDirectory, PRAX_IGES_PARAMETER, iParameter, (iLast ? 1 :0));
    }
    break;
  case 2:
  case 12:
    for (i=0; i<copious->n; i++) {
      iLast = (i < copious->n - 1 ? 0 : 1);
      AddNextReal(fp, copious->pts.triple.x[i], eop, record, iDirectory,
		  PRAX_IGES_PARAMETER, iParameter, 0);
      AddNextReal(fp, copious->pts.triple.y[i], eop, record, iDirectory,
		  PRAX_IGES_PARAMETER, iParameter, 0);
      AddNextReal(fp, copious->pts.triple.z[i], (iLast ? eor : eop), record,
		  iDirectory, PRAX_IGES_PARAMETER, iParameter, (iLast ? 1 :0));
    }
    break;
  case 3:
  case 13:
    for (i=0; i<copious->n; i++) {
      iLast = (i < copious->n - 1 ? 0 : 1);
      AddNextReal(fp, copious->pts.sextuple.x[i], eop, record, iDirectory,
		  PRAX_IGES_PARAMETER, iParameter, 0);
      AddNextReal(fp, copious->pts.sextuple.y[i], eop, record, iDirectory,
		  PRAX_IGES_PARAMETER, iParameter, 0);
      AddNextReal(fp, copious->pts.sextuple.z[i], eop, record, iDirectory,
		  PRAX_IGES_PARAMETER, iParameter, 0);
      AddNextReal(fp, copious->pts.sextuple.i[i], eop, record, iDirectory,
		  PRAX_IGES_PARAMETER, iParameter, 0);
      AddNextReal(fp, copious->pts.sextuple.j[i], eop, record, iDirectory,
		  PRAX_IGES_PARAMETER, iParameter, 0);
      AddNextReal(fp, copious->pts.sextuple.k[i], (iLast ? eor : eop), record,
		  iDirectory, PRAX_IGES_PARAMETER, iParameter, (iLast ? 1 :0));
    }
    break;
  }
  return (*iParameter - iCount);
}

/********* WriteIges110() *********
* 1     Purpose
* 
*       Write an IGES Type110 (Line) entity.
* 
* 2     Specification
* 
*       #include "iges.h"
*       int WriteIges110(FILE *fp, Type110 *line, char *eop, char *eor,
*                        int iDirectory, int *iParameter);
* 
* 3     Description
* 
*       This function writes a Type110 (Line) entity to an IGES file.
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
*           1.  FILE * fp
*               On entry:  the file structure of the IGES file.
* 
* 
*           2.  Type110 * line
*               On entry:  the address of a structure containing the Type110
* 	      (Line) entity.
* 
*           3.  char eop
*               On entry:  the parameter delimiter character.
* 
*           4.  char eor
*               On entry:  the record delimiter character.
* 
*           5.  int iDirectory
*               On entry:  the line number of the first line of this entity's
* 	      Directory entry.
* 
*           6.  int * iParameter
*               On  entry:  the address of a variable containing the current
* 	      line number,  before the write.
* 
*               On exit:  the address of a variable containing the current line
* 	      number, after the write.
* 
* 6     Return Values, Error Indicators and Warnings
* 
*       The number of lines used to write the entity to the file is returned.
* 
* Functions referenced by WriteIges110() are:
*  AddNextInteger()
*  AddNextReal()
*
* Functions that reference WriteIges110() are:
*  SaveIgesSurf()
*/

int WriteIges110(FILE *fp, Type110 *line, char eop, char eor,
		   int iDirectory, int *iParameter)
{
  char record[73];
  int iCount;

  iCount = *iParameter;
  record[0] = '\0';
  AddNextInteger(fp, 110, eop, record, iDirectory, PRAX_IGES_PARAMETER,
		 iParameter, 0);
  AddNextReal(fp, line->x1, eop, record, iDirectory, PRAX_IGES_PARAMETER,
	       iParameter, 0);
  AddNextReal(fp, line->y1, eop, record, iDirectory, PRAX_IGES_PARAMETER,
	       iParameter, 0);
  AddNextReal(fp, line->z1, eop, record, iDirectory, PRAX_IGES_PARAMETER,
	       iParameter, 0);
  AddNextReal(fp, line->x2, eop, record, iDirectory, PRAX_IGES_PARAMETER,
	       iParameter, 0);
  AddNextReal(fp, line->y2, eop, record, iDirectory, PRAX_IGES_PARAMETER,
	       iParameter, 0);
  AddNextReal(fp, line->z2, eor, record, iDirectory, PRAX_IGES_PARAMETER,
	       iParameter, 1);

  return (*iParameter - iCount);
}

/********* WriteIges112() *********
* 1     Purpose
* 
*       Write an IGES Type112 (Parametric Spline Curve) entity.
* 
* 2     Specification
* 
*       #include "iges.h"
*       int  WriteIges112(FILE *fp, PowCurv *pgeom, char *eop, char *eor,
*                         int iDirectory, int *iParameter);
* 
* 3     Description
* 
*       This function writes a Type112 (Parametric Spline Curve) entity to an
*       IGES file.
* 
* 4     References
* 
* 
*       [1]   Digital Representation for Communication of Product Definition
*             Data, US PRO/IPO-100, Initial Graphics Exchange Specification
* 	    (IGES) 5.2, IGES/PDES Organization, U.S. Product Data Association,
* 	    Fairfax, VA, November 1993.
* 
* 5     Parameters
* 
* 
*           1.  FILE * fp
*               On entry:  the file structure of the IGES file.
* 
*           2.  PowCurv * pgeom
*               On  entry:   the  address  of  a  structure  containing  the
* 	      Type112  (Parametric  Spline Curve) entity.
* 
*           3.  char eop
*               On entry:  the parameter delimiter character.
* 
*           4.  char eor
*               On entry:  the record delimiter character.
* 
*           5.  int iDirectory
*               On entry:  the line number of the first line of this entity's
* 	      Directory entry.
* 
*           6.  int * iParameter
*               On  entry:  the address of a variable containing the current
* 	      line number,  before the write.
* 
*               On exit:  the address of a variable containing the current line
* 	      number, after the write.
* 
* 6     Return Values, Error Indicators and Warnings
* 
*       The number of lines used to write the entity to the file is returned.
* 
* Functions referenced by WriteIges112() are:
*  AddNextInteger()
*  AddNextReal()
*  PowCurv_eval()
*  vectfree()
*
* Functions that reference WriteIges112() are:
*  SaveIgesCos()
*  SaveIgesCurv()
*  SaveIgesSurf()
*/

int WriteIges112(FILE *fp, PowCurv *pgeom, char eop, char eor,
		   int iDirectory, int *iParameter)
{
  static double fact[] = {1.0, 1.0, 2.0, 6.0};
  vector *tp[4];
  char record[73];
  int i, iCount, j;
  
  iCount = *iParameter;
  record[0] = '\0';
  AddNextInteger(fp, 112, eop, record, iDirectory, PRAX_IGES_PARAMETER,
		 iParameter, 0);
  AddNextInteger(fp, pgeom->order-1, eop, record, iDirectory,
		 PRAX_IGES_PARAMETER, iParameter, 0);
  AddNextInteger(fp, pgeom->order-2, eop, record, iDirectory,
		 PRAX_IGES_PARAMETER, iParameter, 0);
  AddNextInteger(fp, 3, eop, record, iDirectory,
		 PRAX_IGES_PARAMETER, iParameter, 0);
  AddNextInteger(fp, pgeom->nsegmts, eop, record, iDirectory,
		 PRAX_IGES_PARAMETER, iParameter, 0);

  for (i=0; i<=pgeom->nsegmts; ++i)    /* breakpoints */
    AddNextReal(fp, pgeom->knots[i], eop, record, iDirectory,
		PRAX_IGES_PARAMETER, iParameter, 0);
     
  for (i=0; i<pgeom->nsegmts; ++i) {   /* coefficients */
    for (j=0; j<4; ++j)
      AddNextReal(fp, pgeom->contpts[i][j]->x, eop, record, iDirectory,
		  PRAX_IGES_PARAMETER, iParameter, 0);
    for (j=0; j<4; ++j)
      AddNextReal(fp, pgeom->contpts[i][j]->y, eop, record, iDirectory,
		  PRAX_IGES_PARAMETER, iParameter, 0);
    for (j=0; j<4; ++j)
      AddNextReal(fp, pgeom->contpts[i][j]->z, eop, record, iDirectory,
		  PRAX_IGES_PARAMETER, iParameter, 0);
  }
  for (i=0; i<4; ++i) {   /* Taylor series for x */
    tp[i] = PowCurv_eval (pgeom, pgeom->knots[pgeom->nsegmts], i);
    AddNextReal(fp, tp[i]->x/fact[i], eop, record, iDirectory,
		PRAX_IGES_PARAMETER, iParameter, 0);
  }
  for (i=0; i<4; i++)
    AddNextReal(fp, tp[i]->y/fact[i], eop, record, iDirectory,
		PRAX_IGES_PARAMETER, iParameter, 0);
  for (i=0; i<4; i++) {
    AddNextReal(fp, tp[i]->z/fact[i], (i < 3 ? eop : eor), record, iDirectory,
		PRAX_IGES_PARAMETER, iParameter, (i < 3 ? 0 : 1));
    vectfree(tp[i]);
  }
     
  return (*iParameter - iCount);
}

/********* WriteIges114() *********
* 1     Purpose
* 
*       Write an IGES Type114 (Parametric Spline Surface) entity.
* 
* 2     Specification
* 
*       #include "iges.h"
*       int  WriteIges114(FILE *fp, PowSurf *sgeom, char *eop, char *eor,
*                         int  iDirectory,  int *iParameter);
* 
* 3     Description
* 
*       This function writes a Type114 (Parametric Spline Surface) entity to
*       an IGES file.
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
*           1.  FILE * fp
*               On entry:  the file structure of the IGES file.
* 
*           2.  PowSurf * sgeom
*               On  entry:   the  address  of  a  structure  containing  the
* 	      Type114  (Parametric  Spline Surface) entity.
* 
*           3.  char eop
*               On entry:  the parameter delimiter character.
* 
*           4.  char eor
*               On entry:  the record delimiter character.
* 
*           5.  int iDirectory
*               On entry:  the line number of the first line of this entity's
* 	      Directory entry.
* 
*           6.  int * iParameter
*               On  entry:  the address of a variable containing the current
* 	      line number,  before the write.
* 
*               On exit:  the address of a variable containing the current
* 	      line number, after the write.
* 
* 6     Return Values, Error Indicators and Warnings
* 
*       The number of lines used to write the entity to the file is returned.
* 
* Functions referenced by WriteIges114() are:
*  AddNextInteger()
*  AddNextReal()
*
* Functions that reference WriteIges114() are:
*  SaveIgesSurf()
*/

int WriteIges114(FILE *fp, PowSurf *sgeom, char eop, char eor,
		   int iDirectory, int *iParameter)
{
  char record[73];
  int i, j, k, iCount, iLast = 0;

  iCount = *iParameter;
  record[0] = '\0';
  AddNextInteger(fp, 114, eop, record, iDirectory, PRAX_IGES_PARAMETER,
		 iParameter, 0);
  AddNextInteger(fp, 3, eop, record, iDirectory, PRAX_IGES_PARAMETER,
		 iParameter, 0);
  AddNextInteger(fp, 0, eop, record, iDirectory, PRAX_IGES_PARAMETER,
		 iParameter, 0);
  AddNextInteger(fp, sgeom->usegmts, eop, record, iDirectory,
		 PRAX_IGES_PARAMETER, iParameter, 0);
  AddNextInteger(fp, sgeom->vsegmts, eop, record, iDirectory,
		 PRAX_IGES_PARAMETER, iParameter, 0);

  for (i=0; i<=sgeom->usegmts; i++)
    AddNextReal(fp, sgeom->uknots[i], eop, record, iDirectory,
		PRAX_IGES_PARAMETER, iParameter, 0);
  for (i=0; i<=sgeom->vsegmts; i++)
    AddNextReal(fp, sgeom->vknots[i], eop, record, iDirectory,
		PRAX_IGES_PARAMETER, iParameter, 0);

  for (i=0; i<sgeom->usegmts; i++) {
    for (j=0; j<sgeom->vsegmts; j++) {
      for (k=0; k<16; k++)		/* X coefficients */
	AddNextReal(fp, sgeom->contpts[k][i][j]->x, eop, record, iDirectory,
		    PRAX_IGES_PARAMETER, iParameter, 0);
      for (k=0; k<16; k++)		/* Y coefficients */
	AddNextReal(fp, sgeom->contpts[k][i][j]->y, eop, record, iDirectory,
		    PRAX_IGES_PARAMETER, iParameter, 0);
      for (k=0; k<16; k++)		/* Z coefficients */
	AddNextReal(fp, sgeom->contpts[k][i][j]->z, eop, record, iDirectory,
		    PRAX_IGES_PARAMETER, iParameter, 0);
    }
    for (k=0; k<48; k++) {		/* 48 placeholders */
      iLast = (i < sgeom->usegmts-1 || k < 47 ? 0 : 1);
      AddNextReal(fp, 0.0, (iLast ? eor : eop), record, iDirectory,
		  PRAX_IGES_PARAMETER, iParameter, (iLast ? 1 : 0));
    }
  }  
  return (*iParameter - iCount);
}

/********* WriteIges118() *********
* 1     Purpose
* 
*       Write an IGES Type118 (Ruled Surface) entity.
* 
* 2     Specification
* 
*       #include "iges.h"
*       int WriteIges118(FILE *fp, Type118 *ruled, char *eop, char *eor,
*                        int iDirectory, int *iParameter);
* 
* 3     Description
* 
*       This function writes a Type118 (Ruled Surface) entity to an IGES file.
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
*           1.  FILE * fp
*               On entry:  the file structure of the IGES file.
* 
*           2.  Type118 * ruled
*               On entry:  the address of a structure containing the Type118
* 	      (Ruled Surface) entity.
* 
*           3.  char eop
*               On entry:  the parameter delimiter character.
* 
*           4.  char eor
*               On entry:  the record delimiter character.
* 
*           5.  int iDirectory
*               On entry:  the line number of the first line of this entity's
* 	      Directory entry.
* 
*           6.  int * iParameter
*               On  entry:  the address of a variable containing the current
* 	      line number,  before the write.
* 
*               On exit:  the address of a variable containing the current line
* 	      number, after the write.
* 
* 6     Return Values, Error Indicators and Warnings
* 
*       The number of lines used to write the entity to the file is returned.
* 
* Functions referenced by WriteIges118() are:
*  AddNextInteger()
*
* Functions that reference WriteIges118() are:
*  SaveIgesSurf()
*/

int WriteIges118(FILE *fp, Type118 *ruled, char eop, char eor,
		   int iDirectory, int *iParameter)
{
  char record[73];
  int iCount;

  iCount = *iParameter;
  record[0] = '\0';
  AddNextInteger(fp, 118, eop, record, iDirectory, PRAX_IGES_PARAMETER,
		 iParameter, 0);
  AddNextInteger(fp, ruled->de1, eop, record, iDirectory, PRAX_IGES_PARAMETER,
		 iParameter, 0);
  AddNextInteger(fp, ruled->de2, eop, record, iDirectory, PRAX_IGES_PARAMETER,
		 iParameter, 0);
  AddNextInteger(fp, ruled->dirflg, eop, record, iDirectory,
		 PRAX_IGES_PARAMETER, iParameter, 0);
  AddNextInteger(fp, ruled->devflg, eor, record, iDirectory,
		 PRAX_IGES_PARAMETER, iParameter, 1);

  return (*iParameter - iCount);
}

/********* WriteIges120() *********
* 1     Purpose
* 
*       Write an IGES Type120 (Surface Of Revolution) entity.
* 
* 2     Specification
* 
*       #include "iges.h"
*       int  WriteIges120(FILE  *fp,  Type120  *surfrev,  char  *eop,
*                         char  *eor,  int  iDirectory,  int *iParameter);
* 
* 3     Description
* 
*       This function writes a Type120 (Surface Of Revolution) entity to an
*       IGES file.
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
*           1.  FILE * fp
*               On entry:  the file structure of the IGES file.
* 
*           2.  Type120 * surfrev
*               On entry:  the address of a structure containing the Type120
* 	      (Surface Of Revolution) entity.
* 
*           3.  char eop
*               On entry:  the parameter delimiter character.
* 
*           4.  char eor
*               On entry:  the record delimiter character.
* 
*           5.  int iDirectory
*               On entry:  the line number of the first line of this entity's
* 	      Directory entry.
* 
*           6.  int * iParameter
*               On  entry:  the address of a variable containing the current
* 	      line number,  before the write.
* 
*               On exit:  the address of a variable containing the current
* 	      line number, after the write.
* 
* 6     Return Values, Error Indicators and Warnings
* 
*       The number of lines used to write the entity to the file is returned.
* 
* Functions referenced by WriteIges120() are:
*  AddNextInteger()
*  AddNextReal()
*
* Functions that reference WriteIges120() are:
*  SaveIgesSurf()
*/

int WriteIges120(FILE *fp, Type120 *surf, char eop, char eor,
		   int iDirectory, int *iParameter)
{
  char record[73];
  int iCount;

  iCount = *iParameter;
  record[0] = '\0';
  AddNextInteger(fp, 120, eop, record, iDirectory, PRAX_IGES_PARAMETER,
		 iParameter, 0);
  AddNextInteger(fp, surf->de1, eop, record, iDirectory, PRAX_IGES_PARAMETER,
		 iParameter, 0);
  AddNextInteger(fp, surf->de2, eop, record, iDirectory, PRAX_IGES_PARAMETER,
		 iParameter, 0);
  AddNextReal(fp, surf->start, eop, record, iDirectory,
	       PRAX_IGES_PARAMETER, iParameter, 0);
  AddNextReal(fp, surf->term, eor, record, iDirectory,
	       PRAX_IGES_PARAMETER, iParameter, 1);

  return (*iParameter - iCount);
}

/********* WriteIges122() *********
* 1     Purpose
* 
*       Write an IGES Type122 (Tabulated Cylinder) entity.
* 
* 2     Specification
* 
*       #include "iges.h"
*       int WriteIges122(FILE *fp, Type122 *tabcyl, char *eop, char *eor,
*                        int iDirectory, int *iParameter);
* 
* 3     Description
* 
*       This function writes a Type122 (Tabulated Cylinder) entity to an IGES
*       file.
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
*           1.  FILE * fp
*               On entry:  the file structure of the IGES file.
* 
*           2.  Type122 * tabcyl
*               On  entry:  the address of a structure containing the Type122
* 	      (Tabulated Cylinder) entity.
* 
*           3.  char eop
*               On entry:  the parameter delimiter character.
* 
*           4.  char eor
*               On entry:  the record delimiter character.
* 
*           5.  int iDirectory
*               On entry:  the line number of the first line of this entity's
* 	      Directory entry.
* 
*           6.  int * iParameter
*               On  entry:  the address of a variable containing the current
* 	      line number,  before the write.
* 
*               On exit:  the address of a variable containing the current line
* 	      number, after the write.
* 
* 6     Return Values, Error Indicators and Warnings
* 
*       The number of lines used to write the entity to the file is returned.
* 
* Functions referenced by WriteIges122() are:
*  AddNextInteger()
*  AddNextReal()
*
* Functions that reference WriteIges122() are:
*  SaveIgesSurf()
*/

int WriteIges122(FILE *fp, Type122 *surf, char eop, char eor,
		   int iDirectory, int *iParameter)
{
  char record[73];
  int iCount;

  iCount = *iParameter;
  record[0] = '\0';
  AddNextInteger(fp, 122, eop, record, iDirectory, PRAX_IGES_PARAMETER,
		 iParameter, 0);
  AddNextInteger(fp, surf->de, eop, record, iDirectory, PRAX_IGES_PARAMETER,
		 iParameter, 0);
  AddNextReal(fp, surf->lx, eop, record, iDirectory, PRAX_IGES_PARAMETER,
	       iParameter, 0);
  AddNextReal(fp, surf->ly, eop, record, iDirectory, PRAX_IGES_PARAMETER,
	       iParameter, 0);
  AddNextReal(fp, surf->lz, eor, record, iDirectory, PRAX_IGES_PARAMETER,
	       iParameter, 1);

  return (*iParameter - iCount);
}

/********* WriteIges126() *********
* 1     Purpose
* 
*       Write an IGES Type126 (NURBS Curve) entity.
* 
* 2     Specification
* 
*       #include "iges.h"
*       int WriteIges126(FILE *fp, ParSurf *egeom, char *eop, char *eor,
*                        int iDirectory, int *iParameter);
* 
* 3     Description
* 
*       This function writes a Type126 (NURBS Curve) entity to an IGES file.
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
*           1.  FILE * fp
*               On entry:  the file structure of the IGES file.
* 
*           2.  ParCurv * egeom
*               On entry:  the address of a structure containing the Type126
* 	      (NURBS Curve) entity.
* 
*           3.  char eop
*               On entry:  the parameter delimiter character.
* 
*           4.  char eor
*               On entry:  the record delimiter character.
* 
*           5.  int iDirectory
*               On entry:  the line number of the first line of this entity's
* 	      Directory entry.
* 
*           6.  int * iParameter
*               On  entry:  the address of a variable containing the current
* 	      line number,  before the write.
* 
*               On exit:  the address of a variable containing the current
* 	      line number, after the write.
* 
* 6     Return Values, Error Indicators and Warnings
* 
*       The number of lines used to write the entity to the file is returned.
* 
* Functions referenced by WriteIges126() are:
*  AddNextInteger()
*  AddNextReal()
*  CheckCurvClosed()
*  CheckCurvIntegral()
*  CheckCurvPlanar()
*
* Functions that reference WriteIges126() are:
*  SaveIgesCos()
*  SaveIgesCurv()
*  SaveIgesSurf()
*  SaveIgesTrim()
*/

int WriteIges126(FILE *fp, ParCurv *egeom, char eop, char eor,
		   int iDirectory, int *iParameter)
{
  char record[73];
  double nx, ny, nz;
  int i, iCount, iPlanar = 0;

  iCount = *iParameter;
  record[0] = '\0';
  AddNextInteger(fp, 126, eop, record, iDirectory, PRAX_IGES_PARAMETER,
		 iParameter, 0);
  AddNextInteger(fp, egeom->ncontpts-1, eop, record, iDirectory,
		 PRAX_IGES_PARAMETER, iParameter, 0);
  AddNextInteger(fp, egeom->order-1, eop, record, iDirectory,
		 PRAX_IGES_PARAMETER, iParameter, 0);
  if (!(iPlanar = CheckCurvPlanar(egeom, &nx, &ny, &nz))) {
    nx = 0.0;
    ny = 0.0;
    nz = 0.0;
  }
  AddNextInteger(fp, iPlanar, eop, record, iDirectory, PRAX_IGES_PARAMETER,
		 iParameter, 0);
  AddNextInteger(fp, CheckCurvClosed(egeom), eop, record, iDirectory,
		 PRAX_IGES_PARAMETER, iParameter, 0);
  AddNextInteger(fp, CheckCurvIntegral(egeom), eop, record, iDirectory,
		 PRAX_IGES_PARAMETER, iParameter, 0);
  AddNextInteger(fp, 0, eop, record, iDirectory, PRAX_IGES_PARAMETER,
		 iParameter, 0);

  for (i=0; i<egeom->order + egeom->ncontpts; i++)
    AddNextReal(fp, egeom->knots[i], eop, record, iDirectory,
		PRAX_IGES_PARAMETER, iParameter, 0);

  for (i=0; i<egeom->ncontpts; i++)
    AddNextReal(fp, egeom->contpts[i]->w, eop, record, iDirectory,
		PRAX_IGES_PARAMETER, iParameter, 0);

  /* Praxiteles stores (x,y,z,w) */
  /* IGES want to see (x/w, y/w, z/w, w) */

  for (i=0; i<egeom->ncontpts; i++) {
    AddNextReal(fp, egeom->contpts[i]->x/egeom->contpts[i]->w,
		eop, record, iDirectory, PRAX_IGES_PARAMETER, iParameter, 0);
    AddNextReal(fp, egeom->contpts[i]->y/egeom->contpts[i]->w,
		eop, record, iDirectory, PRAX_IGES_PARAMETER, iParameter, 0);
    AddNextReal(fp, egeom->contpts[i]->z/egeom->contpts[i]->w,
		eop, record, iDirectory, PRAX_IGES_PARAMETER, iParameter, 0);
    }
  AddNextReal(fp, egeom->knots[0], eop, record, iDirectory,
	      PRAX_IGES_PARAMETER, iParameter, 0);
  AddNextReal(fp, egeom->knots[egeom->ncontpts], eop, record, iDirectory,
	      PRAX_IGES_PARAMETER, iParameter, 0);
  AddNextReal(fp, nx, eop, record, iDirectory, PRAX_IGES_PARAMETER,
	      iParameter, 0);
  AddNextReal(fp, ny, eop, record, iDirectory, PRAX_IGES_PARAMETER,
	      iParameter, 0);
  AddNextReal(fp, nz, eor, record, iDirectory, PRAX_IGES_PARAMETER,
	      iParameter, 1);

  return (*iParameter - iCount);
}

/********* WriteIges128() *********
* 1     Purpose
* 
*       Write an IGES Type128 (NURBS Surface) entity.
* 
* 2     Specification
* 
*       #include "iges.h"
*       int WriteIges128(FILE *fp, ParSurf *fgeom, char *eop, char *eor,
*       int iDirectory, int *iParameter);
* 
* 3     Description
* 
*       This function writes a Type128 (NURBS Surface) entity to an IGES file.
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
*           1.  FILE * fp
*               On entry:  the file structure of the IGES file.
* 
*           2.  ParSurf * fgeom
*               On entry: the address of a structure containing the Type128
* 	      (NURBS Surface) entity.
* 
*           3.  char eop
*               On entry:  the parameter delimiter character.
* 
*           4.  char eor
*               On entry:  the record delimiter character.
* 
*           5.  int iDirectory
*               On entry:  the line number of the first line of this entity's
* 	      Directory entry.
* 
*           6.  int * iParameter
*               On  entry:  the address of a variable containing the current
* 	      line number,  before the write.
* 
*               On exit:  the address of a variable containing the current line
* 	      number, after the write.
* 
* 6     Return Values, Error Indicators and Warnings
* 
*       The number of lines used to write the entity to the file is returned.
* 
* Functions referenced by WriteIges128() are:
*  AddNextInteger()
*  AddNextReal()
*  CheckSurfClosedU()
*  CheckSurfClosedV()
*  CheckSurfIntegral()
*
* Functions that reference WriteIges128() are:
*  SaveIgesCos()
*  SaveIgesSurf()
*  SaveIgesTrim()
*/

int WriteIges128(FILE *fp, ParSurf *fgeom, char eop, char eor,
		  int iDirectory, int *iParameter)
{
  char record[73];
  int i, iCount, j;

  iCount = *iParameter;
  record[0] = '\0';
  AddNextInteger(fp, 128, eop, record, iDirectory, PRAX_IGES_PARAMETER,
		 iParameter, 0);
  AddNextInteger(fp, fgeom->ucontpts-1, eop, record, iDirectory,
		 PRAX_IGES_PARAMETER, iParameter, 0);
  AddNextInteger(fp, fgeom->vcontpts-1, eop, record, iDirectory,
		 PRAX_IGES_PARAMETER, iParameter, 0);
  AddNextInteger(fp, fgeom->uorder-1, eop, record, iDirectory,
		 PRAX_IGES_PARAMETER, iParameter, 0);
  AddNextInteger(fp, fgeom->vorder-1, eop, record, iDirectory,
		 PRAX_IGES_PARAMETER, iParameter, 0);
  AddNextInteger(fp, CheckSurfClosedU(fgeom), eop, record, iDirectory,
		 PRAX_IGES_PARAMETER, iParameter, 0);
  AddNextInteger(fp, CheckSurfClosedV(fgeom), eop, record, iDirectory,
		 PRAX_IGES_PARAMETER, iParameter, 0);
  AddNextInteger(fp, CheckSurfIntegral(fgeom), eop, record, iDirectory,
		 PRAX_IGES_PARAMETER, iParameter, 0);
  AddNextInteger(fp, 0, eop, record, iDirectory, PRAX_IGES_PARAMETER,
		 iParameter, 0);
  AddNextInteger(fp, 0, eop, record, iDirectory, PRAX_IGES_PARAMETER,
		 iParameter, 0);

  for (i=0; i<fgeom->uorder + fgeom->ucontpts; i++)
    AddNextReal(fp, fgeom->uknots[i], eop, record, iDirectory,
		PRAX_IGES_PARAMETER, iParameter, 0);

  for (i=0; i<fgeom->vorder + fgeom->vcontpts; i++)
    AddNextReal(fp, fgeom->vknots[i], eop, record, iDirectory,
		PRAX_IGES_PARAMETER, iParameter, 0);

  for (j=0; j<fgeom->vcontpts; j++)
    for (i=0; i<fgeom->ucontpts; i++)
      AddNextReal(fp, fgeom->contpts[i][j]->w, eop, record, iDirectory,
		  PRAX_IGES_PARAMETER, iParameter, 0);

  /* Praxiteles stores (x,y,z,w) */
  /* IGES want to see (x/w, y/w, z/w, w) */

  for (j=0; j<fgeom->vcontpts; j++)
    for (i=0; i<fgeom->ucontpts; i++) {
      AddNextReal(fp, fgeom->contpts[i][j]->x/fgeom->contpts[i][j]->w,
		  eop, record, iDirectory, PRAX_IGES_PARAMETER, iParameter, 0);
      AddNextReal(fp, fgeom->contpts[i][j]->y/fgeom->contpts[i][j]->w,
		  eop, record, iDirectory, PRAX_IGES_PARAMETER, iParameter, 0);
      AddNextReal(fp, fgeom->contpts[i][j]->z/fgeom->contpts[i][j]->w,
		  eop, record, iDirectory, PRAX_IGES_PARAMETER, iParameter, 0);
    }
  AddNextReal(fp, fgeom->uknots[fgeom->uorder-1], eop, record, iDirectory,
	      PRAX_IGES_PARAMETER, iParameter, 0);
  AddNextReal(fp, fgeom->uknots[fgeom->ucontpts], eop, record, iDirectory,
	      PRAX_IGES_PARAMETER, iParameter, 0);
  AddNextReal(fp, fgeom->vknots[fgeom->vorder-1], eop, record, iDirectory,
	      PRAX_IGES_PARAMETER, iParameter, 0);
  AddNextReal(fp, fgeom->vknots[fgeom->vcontpts], eor, record, iDirectory,
	      PRAX_IGES_PARAMETER, iParameter, 1);

  return (*iParameter - iCount);
}

/********* WriteIges142() *********
* 1     Purpose
* 
*       Write an IGES Type142 (Curve on Parametric Surface) entity.
* 
* 2     Specification
* 
*       #include "iges.h"
*       int WriteIges142(FILE *fp, Type142 *cos, char *eop, char *eor,
*                        int iDirectory, int *iParameter);
* 
* 3     Description
* 
*       This function writes a Type142 (Curve on Parametric Surface) entity
*       to an IGES file.
* 
* 4     References
* 
* 
*       [1]   Digital Representation for Communication of Product Definition
*             Data, US PRO/IPO-100, Initial Graphics Exchange Specification
* 	    (IGES) 5.2, IGES/PDES Organization, U.S. Product Data Association,
* 	    Fairfax, VA, November 1993.
* 
* 5     Parameters
* 
*           1.  FILE * fp
*               On entry:  the file structure of the IGES file.
* 
*           2.  Type142 * cos
*               On entry:  the address of a structure containing the Type142
* 	      (Curve on Parametric Surface) entity.
* 
*           3.  char eop
*               On entry:  the parameter delimiter character.
* 
*           4.  char eor
*               On entry:  the record delimiter character.
* 
*           5.  int iDirectory
*               On entry:  the line number of the first line of this entity's
* 	      Directory entry.
* 
*           6.  int * iParameter
*               On  entry:  the address of a variable containing the current
* 	      line number,  before the write.
* 
*               On exit:  the address of a variable containing the current line
* 	      number, after the write.
* 
* 6     Return Values, Error Indicators and Warnings
* 
*       The number of lines used to write the entity to the file is returned.
* 
* Functions referenced by WriteIges142() are:
*  AddNextInteger()
*
* Functions that reference WriteIges142() are:
*  SaveIgesCos()
*  SaveIgesTrim()
*/

int WriteIges142(FILE *fp, Type142 *cos, char eop, char eor,
		   int iDirectory, int *iParameter)
{
  char record[73];
  int iCount;

  iCount = *iParameter;
  record[0] = '\0';
  AddNextInteger(fp, 142, eop, record, iDirectory, PRAX_IGES_PARAMETER,
		 iParameter, 0);
  AddNextInteger(fp, cos->crtn, eop, record, iDirectory, PRAX_IGES_PARAMETER,
		 iParameter, 0);
  AddNextInteger(fp, cos->sptr, eop, record, iDirectory, PRAX_IGES_PARAMETER,
		 iParameter, 0);
  AddNextInteger(fp, cos->bptr, eop, record, iDirectory, PRAX_IGES_PARAMETER,
		 iParameter, 0);
  AddNextInteger(fp, cos->cptr, eop, record, iDirectory, PRAX_IGES_PARAMETER,
		 iParameter, 0);
  AddNextInteger(fp, cos->pref, eor, record, iDirectory, PRAX_IGES_PARAMETER,
		 iParameter, 1);

  return (*iParameter - iCount);
}

/********* WriteIges144() *********
* 1     Purpose
* 
*       Write an IGES Type144 (Trimmed (Parametric) Surface) entity.
* 
* 2     Specification
* 
*       #include "iges.h"
*       int WriteIges144(FILE *fp, Type144 *trim, char *eop, char *eor,
*                        int iDirectory, int *iParameter);
* 
* 3     Description
* 
*       This function writes a Type144 (Trimmed (Parametric) Surface) entity
*       to an IGES file.
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
*           1.  FILE * fp
*               On entry:  the file structure of the IGES file.
* 
*           2.  Type144 * trim
*               On entry:  the address of a structure containing the Type144
* 	      (Trimmed (Parametric) Surface) entity.
* 
*           3.  char eop
*               On entry:  the parameter delimiter character.
* 
*           4.  char eor
*               On entry:  the record delimiter character.
* 
*           5.  int iDirectory
*               On entry:  the line number of the first line of this entity's
* 	      Directory entry.
* 
*           6.  int * iParameter
*               On  entry:  the address of a variable containing the current
* 	      line number,  before the write.
* 
*               On exit:  the address of a variable containing the current line
* 	      number, after the write.
* 
* 6     Return Values, Error Indicators and Warnings
* 
*       The number of lines used to write the entity to the file is returned.
* 
* Functions referenced by WriteIges144() are:
*  AddNextInteger()
*
* Functions that reference WriteIges144() are:
*  SaveIgesTrim()
*/

int WriteIges144(FILE *fp, Type144 *trim, char eop, char eor,
		   int iDirectory, int *iParameter)
{
  char record[73];
  int i, iCount;

  iCount = *iParameter;
  record[0] = '\0';
  AddNextInteger(fp, 144, eop, record, iDirectory, PRAX_IGES_PARAMETER,
		 iParameter, 0);
  AddNextInteger(fp, trim->pts, eop, record, iDirectory, PRAX_IGES_PARAMETER,
		 iParameter, 0);
  AddNextInteger(fp, trim->n1, eop, record, iDirectory, PRAX_IGES_PARAMETER,
		 iParameter, 0);
  AddNextInteger(fp, trim->n2, eop, record, iDirectory, PRAX_IGES_PARAMETER,
		 iParameter, 0);
  AddNextInteger(fp, trim->pt0, eop, record, iDirectory, PRAX_IGES_PARAMETER,
		 iParameter, trim->n2 > 0 ? 0 : 1);
  for (i=0; i<trim->n2-1; i++)
    AddNextInteger(fp, trim->pti[i], eop, record, iDirectory,
		   PRAX_IGES_PARAMETER, iParameter, 0);
  if (trim->n2 > 0)
    AddNextInteger(fp, trim->pti[trim->n2-1], eor, record, iDirectory,
		   PRAX_IGES_PARAMETER, iParameter, 1);

  return (*iParameter - iCount);
}

/* Copyright (C) 1997 Massachusetts Institute of Technology, Cambridge, MA
   All rights reserved
*/

/* read.c */

/* ReadIges100()
 * ReadIges102()
 * ReadIges102Curv()
 * ReadIges104()
 * ReadIges106()
 * ReadIges108()
 * ReadIges110()
 * ReadIges112()
 * ReadIges114()
 * ReadIges116()
 * ReadIges118()
 * ReadIges120()
 * ReadIges122()
 * ReadIges124()
 * ReadIges126()
 * ReadIges128()
 * ReadIges132()
 * ReadIges142()
 * ReadIges144()
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h>
#include <math.h>
#include "iges.h"
#include "bspl.h"
#include "editor.h"
char *strdup(const char *);

Type100 *ReadIges100(FILE *fp, char *record, int *ip, char eop, char eor)
{
  Type100 *circle;
  double val;
  int ival;

  GetNextInteger(fp, record, ip, eop, eor, 64, &ival);
  if (ival != 100)
    return (NULL);

  if (circle = (Type100 *)gen_array1(1, sizeof(Type100))) {
    GetNextReal(fp, record, ip, eop, eor, 64, &val);
    circle->zt = val;
    GetNextReal(fp, record, ip, eop, eor, 64, &val);
    circle->x1 = val;
    GetNextReal(fp, record, ip, eop, eor, 64, &val);
    circle->y1 = val;
    GetNextReal(fp, record, ip, eop, eor, 64, &val);
    circle->x2 = val;
    GetNextReal(fp, record, ip, eop, eor, 64, &val);
    circle->y2 = val;
    GetNextReal(fp, record, ip, eop, eor, 64, &val);
    circle->x3 = val;
    GetNextReal(fp, record, ip, eop, eor, 64, &val);
    circle->y3 = val;
  }
  return (circle);
}

Type102 *ReadIges102(FILE *fp, char *record, int *ip, char eop, char eor)
{
  Type102 *composite;
  int i, ival;

  GetNextInteger(fp, record, ip, eop, eor, 64, &ival);
  if (ival != 102)
    return (NULL);

  if (composite = (Type102 *)gen_array1(1, sizeof(Type102))) {
    GetNextInteger(fp, record, ip, eop, eor, 64, &ival);
    composite->n = ival;
    composite->de = int_array1(composite->n);
    for (i=0; i<composite->n; i++) {
      GetNextInteger(fp, record, ip, eop, eor, 64, &ival);
      composite->de[i] = ival;
    }
  }
  return (composite);
}

Type104 *ReadIges104(FILE *fp, char *record, int *ip, char eop, char eor)
{
  Type104 *conic;
  double val;
  int ival;

  GetNextInteger(fp, record, ip, eop, eor, 64, &ival);
  if (ival != 104)
    return (NULL);

  if (conic = (Type104 *)gen_array1(1, sizeof(Type104))) {
    GetNextReal(fp, record, ip, eop, eor, 64, &val);
    conic->a = val;
    GetNextReal(fp, record, ip, eop, eor, 64, &val);
    conic->b = val;
    GetNextReal(fp, record, ip, eop, eor, 64, &val);
    conic->c = val;
    GetNextReal(fp, record, ip, eop, eor, 64, &val);
    conic->d = val;
    GetNextReal(fp, record, ip, eop, eor, 64, &val);
    conic->e = val;
    GetNextReal(fp, record, ip, eop, eor, 64, &val);
    conic->f = val;
    GetNextReal(fp, record, ip, eop, eor, 64, &val);
    conic->zt = val;
    GetNextReal(fp, record, ip, eop, eor, 64, &val);
    conic->x1 = val;
    GetNextReal(fp, record, ip, eop, eor, 64, &val);
    conic->y1 = val;
    GetNextReal(fp, record, ip, eop, eor, 64, &val);
    conic->x2 = val;
    GetNextReal(fp, record, ip, eop, eor, 64, &val);
    conic->y2 = val;
  }
  return (conic);
}

Type108 *ReadIges108(FILE *fp, char *record, int *ip, char eop, char eor)
{
  Type108 *plane;
  double val;
  int ival;

  GetNextInteger(fp, record, ip, eop, eor, 64, &ival);
  if (ival != 108)
    return (NULL);

  if (plane = (Type108 *)gen_array1(1, sizeof(Type108))) {
    GetNextReal(fp, record, ip, eop, eor, 64, &val);
    plane->a = val;
    GetNextReal(fp, record, ip, eop, eor, 64, &val);
    plane->b = val;
    GetNextReal(fp, record, ip, eop, eor, 64, &val);
    plane->c = val;
    GetNextReal(fp, record, ip, eop, eor, 64, &val);
    plane->d = val;
    GetNextInteger(fp, record, ip, eop, eor, 64, &ival);
    plane->ptr = ival;
    GetNextReal(fp, record, ip, eop, eor, 64, &val);
    plane->x = val;
    GetNextReal(fp, record, ip, eop, eor, 64, &val);
    plane->y = val;
    GetNextReal(fp, record, ip, eop, eor, 64, &val);
    plane->z = val;
    GetNextReal(fp, record, ip, eop, eor, 64, &val);
    plane->size = val;
  }
  return (plane);
}

/********* ReadIges106() *********
* 1     Purpose
* 
*       Read an IGES Type106 (Copious Data) entity.
* 
* 2     Specification
* 
*       #include "iges.h"
*       Type106 *ReadIges106(FILE *fp, char *record, int *ip, char *eop,
*                            char *eor);
* 
* 3     Description
* 
*       This function reads a Type106 (Copious Data) entity from an IGES file.
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
*           2.  char * record
*               On  entry:  the  address  of  a  character  string  containing
* 	      the  current  line,  before  the read.
* 
*               On exit:  the address of a character string containing the
* 	      current line, after the read.
* 
*           3.  int * ip
*               On entry:  the address of a variable containing the current
* 	      character position, before the read.
* 
*               On exit:  the address of a variable containing the current
* 	      character position, after the read.
* 
*           4.  char eop
*               On entry:  the parameter delimiter character.
* 
*           5.  char eor
*               On entry:  the record delimiter character.
* 
* 6     Return Values, Error Indicators and Warnings
* 
*       The address of the structure containing the Type106 (Copious Data) 
*       entity is returned.
* 
* 8     Further Comments
* 
*       The returned structure should be deallocated with free_copious().
* 
* Functions referenced by ReadIges106() are:
*  alloc_copious()
*  GetNextInteger()
*  GetNextReal()
*
* Functions that reference ReadIges106() are:
*  ReadIgesList()
*  ReadIgesUv()
*/

Type106 *ReadIges106(FILE *fp, char *record, int *ip, char eop, char eor)
{
  Type106 *copious;
  int ipp, n;
  int i;

  GetNextInteger(fp, record, ip, eop, eor, 64, &n);
  if (n != 106)
    return (NULL);

  GetNextInteger(fp, record, ip, eop, eor, 64, &ipp);
  GetNextInteger(fp, record, ip, eop, eor, 64, &n);
  if (copious = alloc_copious(ipp, n)) {
    switch (copious->ip) {
    case 1:    /* form 1 or 11: (x,y) pairs */
    case 11:
      GetNextReal(fp, record, ip, eop, eor, 64, &copious->pts.pair.zt);
      for (i=0; i<copious->n; i++) {
	GetNextReal(fp, record, ip, eop, eor, 64, &copious->pts.pair.x[i]);
	GetNextReal(fp, record, ip, eop, eor, 64, &copious->pts.pair.y[i]);
      }
      break;
    case 2:    /* form 2 or 12: (x,y,z) triples */
    case 12:
      for (i=0; i<copious->n; i++) {
	GetNextReal(fp, record, ip, eop, eor, 64, &copious->pts.triple.x[i]);
	GetNextReal(fp, record, ip, eop, eor, 64, &copious->pts.triple.y[i]);
	GetNextReal(fp, record, ip, eop, eor, 64, &copious->pts.triple.z[i]);
      }
      break;
    case 3:    /* form 3 or 13: (x,y,z,i,j,k) sextuples */
    case 13:
      for (i=0; i<copious->n; i++) {
	GetNextReal(fp, record, ip, eop, eor, 64, &copious->pts.sextuple.x[i]);
	GetNextReal(fp, record, ip, eop, eor, 64, &copious->pts.sextuple.y[i]);
	GetNextReal(fp, record, ip, eop, eor, 64, &copious->pts.sextuple.z[i]);
	GetNextReal(fp, record, ip, eop, eor, 64, &copious->pts.sextuple.i[i]);
	GetNextReal(fp, record, ip, eop, eor, 64, &copious->pts.sextuple.j[i]);
	GetNextReal(fp, record, ip, eop, eor, 64, &copious->pts.sextuple.k[i]);
      }
      break;
    }
  }
  return (copious);
}

/********* ReadIges110() *********
* 1     Purpose
* 
*       Read an IGES Type110 (Line) entity.
* 
* 2     Specification
* 
*       #include "iges.h"
*       Type110 *ReadIges110(FILE *fp, char *record, int *ip, char *eop,
*                            char *eor);
* 
* 3     Description
* 
*       This function reads aIGES Type110 (Line) entity from an IGES file.
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
*           2.  char * record
*               On  entry:  the  address  of  a  character  string  containing
* 	      the  current  line,  before  the read.
* 
*               On exit:  the address of a character string containing the
* 	      current line, after the read.
* 
*           3.  int * ip
*               On entry:  the address of a variable containing the current
* 	      character position, before the read.
* 
*               On exit:  the address of a variable containing the current
* 	      character position, after the read.
* 
*           4.  char eop
*               On entry:  the parameter delimiter character.
* 
*           5.  char eor
*               On entry:  the record delimiter character.
* 
* 6     Return Values, Error Indicators and Warnings
* 
*       The address of the structure containing the Type110 (Line) entity is
*       returned.
* 
* Functions referenced by ReadIges110() are:
*  gen_array1()
*  GetNextInteger()
*  GetNextReal()
*
* Functions that reference ReadIges110() are:
*  ReadIgesLine()
*/

Type110 *ReadIges110(FILE *fp, char *record, int *ip, char eop, char eor)
{
  Type110 *line;
  double val;
  int ival;

  GetNextInteger(fp, record, ip, eop, eor, 64, &ival);
  if (ival != 110)
    return (NULL);

  if (line = (Type110 *)gen_array1(1, sizeof(Type110))) {
    GetNextReal(fp, record, ip, eop, eor, 64, &val);
    line->x1 = val;
    GetNextReal(fp, record, ip, eop, eor, 64, &val);
    line->y1 = val;
    GetNextReal(fp, record, ip, eop, eor, 64, &val);
    line->z1 = val;
    GetNextReal(fp, record, ip, eop, eor, 64, &val);
    line->x2 = val;
    GetNextReal(fp, record, ip, eop, eor, 64, &val);
    line->y2 = val;
    GetNextReal(fp, record, ip, eop, eor, 64, &val);
    line->z2 = val;
  }
  return (line);
}

/********* ReadIges112() *********
* 1     Purpose
* 
*       Read an IGES Type112 (Parametric Spline Curve) entity.
* 
* 2     Specification
* 
*       #include "iges.h"
*       PowCurv *ReadIges112(FILE *fp, char *record, int *ip, char *eop,
*                            char *eor);
* 
* 3     Description
* 
*       This function reads a Type112 (Parametric Spline Curve) entity from an
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
*           2.  char * record
*               On  entry:  the  address  of  a  character  string  containing
* 	      the  current  line,  before  the read.
* 
*               On exit:  the address of a character string containing the
* 	      current line, after the read.
* 
*           3.  int * ip
*               On entry:  the address of a variable containing the current
* 	      character position, before the read.
* 
*               On exit:  the address of a variable containing the current
* 	      character position, after the read.
* 
*           4.  char eop
*               On entry:  the parameter delimiter character.
* 
*           5.  char eor
*               On entry:  the record delimiter character.
* 
* 6     Return Values, Error Indicators and Warnings
* 
*       The address of the structure containing the Type112 (Parametric Spline
*       Curve) entity is returned.
* 
* 8     Further Comments
* 
*       The returned structure should be deallocated with free_pgeom().
* 
* Functions referenced by ReadIges112() are:
*  GetNextInteger()
*  GetNextReal()
*  pgeomalloc()
*
* Functions that reference ReadIges112() are:
*  ReadIgesPowerCurv()
*/

PowCurv *ReadIges112(FILE *fp, char *record, int *ip, char eop, char eor)
{
  PowCurv *pgeom;
  int ctype, h, n, ndim;
  int i;
  
  GetNextInteger(fp, record, ip, eop, eor, 64, &ctype);
  if (ctype != 112)
    return (NULL);

  GetNextInteger(fp, record, ip, eop, eor, 64, &ctype);
  GetNextInteger(fp, record, ip, eop, eor, 64, &h);
  GetNextInteger(fp, record, ip, eop, eor, 64, &ndim);
  GetNextInteger(fp, record, ip, eop, eor, 64, &n);
  if (pgeom = pgeomalloc(n)) {
    pgeom->nsegmts = n;
    pgeom->order = 4;

    for (i=0; i<=n; i++)
      GetNextReal(fp, record, ip, eop, eor, 64, &pgeom->knots[i]);

    for (i=0; i<n; i++) {
      GetNextReal(fp, record, ip, eop, eor, 64, &pgeom->contpts[i][0]->x);
      GetNextReal(fp, record, ip, eop, eor, 64, &pgeom->contpts[i][1]->x);
      GetNextReal(fp, record, ip, eop, eor, 64, &pgeom->contpts[i][2]->x);
      GetNextReal(fp, record, ip, eop, eor, 64, &pgeom->contpts[i][3]->x);

      GetNextReal(fp, record, ip, eop, eor, 64, &pgeom->contpts[i][0]->y);
      GetNextReal(fp, record, ip, eop, eor, 64, &pgeom->contpts[i][1]->y);
      GetNextReal(fp, record, ip, eop, eor, 64, &pgeom->contpts[i][2]->y);
      GetNextReal(fp, record, ip, eop, eor, 64, &pgeom->contpts[i][3]->y);

      GetNextReal(fp, record, ip, eop, eor, 64, &pgeom->contpts[i][0]->z);
      GetNextReal(fp, record, ip, eop, eor, 64, &pgeom->contpts[i][1]->z);
      GetNextReal(fp, record, ip, eop, eor, 64, &pgeom->contpts[i][2]->z);
      GetNextReal(fp, record, ip, eop, eor, 64, &pgeom->contpts[i][3]->z);
    }
  }
  return (pgeom);
}

/********* ReadIges114() *********
* 1     Purpose
* 
*       Read an IGES Type114 (Parametric Spline Curve) entity.
* 
* 2     Specification
* 
*       #include "iges.h"
*       PowSurf *ReadIges114(FILE *fp, char *record, int *ip, char *eop,
*                            char *eor);
* 
* 3     Description
* 
*       This function reads a Type114 (Parametric Spline Surface) entity from
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
*           2.  char * record
*               On  entry:  the  address  of  a  character  string  containing
* 	      the  current  line,  before  the read.
* 
*               On exit:  the address of a character string containing the
* 	      current line, after the read.
* 
*           3.  int * ip
*               On entry:  the address of a variable containing the current
* 	      character position, before the read.
* 
*               On exit:  the address of a variable containing the current
* 	      character position, after the read.
* 
*           4.  char eop
*               On entry:  the parameter delimiter character.
* 
*           5.  char eor
*               On entry:  the record delimiter character.
* 
* 6     Return Values, Error Indicators and Warnings
* 
*       The address of the structure containing the Type114 (Parametric Spline
*       Surface) entity is returned.
* 
* 8     Further Comments
* 
*       The returned structure should be deallocated with free_geom().
* 
* Functions referenced by ReadIges114() are:
*  GetNextInteger()
*  GetNextReal()
*  sgeomalloc()
*
* Functions that reference ReadIges114() are:
*  ReadIgesPowerSurf()
*/

PowSurf *ReadIges114(FILE *fp, char *record, int *ip, char eop, char eor)
{
  PowSurf *sgeom;
  double dummy;
  int ctype, m, n, ptype;
  int i, j, k;

  GetNextInteger(fp, record, ip, eop, eor, 64, &ctype);
  if (ctype != 114)
    return (NULL);

  GetNextInteger(fp, record, ip, eop, eor, 64, &ctype);
  GetNextInteger(fp, record, ip, eop, eor, 64, &ptype);
  GetNextInteger(fp, record, ip, eop, eor, 64, &m);
  GetNextInteger(fp, record, ip, eop, eor, 64, &n);

  if (sgeom = sgeomalloc(m, n)) {
    sgeom->uorder = 4;
    sgeom->vorder = 4;
    sgeom->usegmts = m;
    sgeom->vsegmts = n;

    for (i=0; i<=m; i++)
      GetNextReal(fp, record, ip, eop, eor, 64, &sgeom->uknots[i]);
    for (i=0; i<=n; i++)
      GetNextReal(fp, record, ip, eop, eor, 64, &sgeom->vknots[i]);

    for (i=0; i<m; i++) {
      for (j=0; j<n; j++) {
	for (k=0; k<16; k++)
	  GetNextReal(fp, record, ip, eop, eor, 64,
		      &sgeom->contpts[k][i][j]->x);
	for (k=0; k<16; k++)
	  GetNextReal(fp, record, ip, eop, eor, 64,
		      &sgeom->contpts[k][i][j]->y);
	for (k=0; k<16; k++)
	  GetNextReal(fp, record, ip, eop, eor, 64,
		      &sgeom->contpts[k][i][j]->z);
      }
      for (j=0; j<48; j++)
	GetNextReal(fp, record, ip, eop, eor, 64, &dummy);
    }
  }
  return (sgeom);
}

Type116 *ReadIges116(FILE *fp, char *record, int *ip, char eop, char eor)
{
  Type116 *point;
  double val;
  int ival;

  GetNextInteger(fp, record, ip, eop, eor, 64, &ival);
  if (ival != 116)
    return (NULL);

  if (point = (Type116 *)gen_array1(1, sizeof(Type116))) {
    GetNextReal(fp, record, ip, eop, eor, 64, &val);
    point->x = val;
    GetNextReal(fp, record, ip, eop, eor, 64, &val);
    point->y = val;
    GetNextReal(fp, record, ip, eop, eor, 64, &val);
    point->z = val;
    GetNextInteger(fp, record, ip, eop, eor, 64, &ival);
    point->ptr = ival;
  }
  return (point);
}

/********* ReadIges118() *********
* 1     Purpose
* 
*       Read an IGES Type118 (Ruled Surface) entity.
* 
* 2     Specification
* 
*       #include "iges.h"
*       Type118 *ReadIges118(FILE *fp, char *record, int *ip, char *eop,
*       char *eor);
* 
* 3     Description
* 
*       This function reads a Type118 (Ruled Surface) entity from an IGES file.
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
*           2.  char * record
*               On  entry:  the  address  of  a  character  string  containing
* 	      the  current  line,  before  the read.
* 
*               On exit:  the address of a character string containing the
* 	      current line, after the read.
* 
*           3.  int * ip
*               On entry:  the address of a variable containing the current
* 	      character position, before the read.
* 
*               On exit:  the address of a variable containing the current
* 	      character position, after the read.
* 
*           4.  char eop
*               On entry:  the parameter delimiter character.
* 
*           5.  char eor
*               On entry:  the record delimiter character.
* 
* 6     Return Values, Error Indicators and Warnings
* 
*       The address of the structure containing the Type118 (Ruled Surface)
*       entity is returned.
* 
* Functions referenced by ReadIges118() are:
*  gen_array1()
*  GetNextInteger()
*
* Functions that reference ReadIges118() are:
*  ReadIgesRuledSurf()
*/

Type118 *ReadIges118(FILE *fp, char *record, int *ip, char eop, char eor)
{
  Type118 *ruled;
  int ival;

  GetNextInteger(fp, record, ip, eop, eor, 64, &ival);
  if (ival != 118)
    return (NULL);

  if (ruled = (Type118 *)gen_array1(1, sizeof(Type118))) {
    GetNextInteger(fp, record, ip, eop, eor, 64, &ival);
    ruled->de1 = ival;
    GetNextInteger(fp, record, ip, eop, eor, 64, &ival);
    ruled->de2 = ival;
    GetNextInteger(fp, record, ip, eop, eor, 64, &ival);
    ruled->dirflg = ival;
    GetNextInteger(fp, record, ip, eop, eor, 64, &ival);
    ruled->devflg = ival;
  }
  return (ruled);
}

/********* ReadIges120() *********
* 1     Purpose
* 
*       Read an IGES Type120 (Surface of Revolution) entity.
* 
* 2     Specification
* 
*       #include "iges.h"
*       Type120 *ReadIges120(FILE *fp, char *record, int *ip, char *eop,
*                            char *eor);
* 
* 3     Description
* 
*       This function reads a Type120 (Surface of Revolution) entity from an
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
*           2.  char * record
*               On  entry:  the  address  of  a  character  string  containing
* 	      the  current  line,  before  the read.
* 
*               On exit:  the address of a character string containing the
* 	      current line, after the read.
* 
*           3.  int * ip
*               On entry:  the address of a variable containing the current
* 	      character position, before the read.
* 
*               On exit:  the address of a variable containing the current
* 	      character position, after the read.
* 
*           4.  char eop
*               On entry:  the parameter delimiter character.
* 
*           5.  char eor
*               On entry:  the record delimiter character.
* 
* 6     Return Values, Error Indicators and Warnings
* 
*       The  address  of  the  structure  containing  the  Type120  (Surface
*       of  Revolution)  entity  is returned.
* 
* Functions referenced by ReadIges120() are:
*  gen_array1()
*  GetNextInteger()
*  GetNextReal()
*
* Functions that reference ReadIges120() are:
*  ReadIgesSurfRev()
*/

Type120 *ReadIges120(FILE *fp, char *record, int *ip, char eop, char eor)
{
  Type120 *surf;
  double val;
  int ival;

  GetNextInteger(fp, record, ip, eop, eor, 64, &ival);
  if (ival != 120)
    return (NULL);

  if (surf = (Type120 *)gen_array1(1, sizeof(Type120))) {
    GetNextInteger(fp, record, ip, eop, eor, 64, &ival);
    surf->de1 = ival;
    GetNextInteger(fp, record, ip, eop, eor, 64, &ival);
    surf->de2 = ival;
    GetNextReal(fp, record, ip, eop, eor, 64, &val);
    surf->start = val;
    GetNextReal(fp, record, ip, eop, eor, 64, &val);
    surf->term = val;
  }
  return (surf);
}

/********* ReadIges122() *********
* 1     Purpose
* 
*       Read an IGES Type122 (Tabulated Cylinder) entity.
* 
* 2     Specification
* 
*       #include "iges.h"
*       Type122 *ReadIges122(FILE *fp, char *record, int *ip, char *eop,
*       char *eor);
* 
* 3     Description
* 
*       This function reads a Type122 (Tabulated Cylinder) entity from an IGES
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
*           2.  char * record
*               On  entry:  the  address  of  a  character  string  containing
* 	      the  current  line,  before  the read.
* 
*               On exit:  the address of a character string containing the
* 	      current line, after the read.
* 
*           3.  int * ip
*               On entry:  the address of a variable containing the current
* 	      character position, before the read.
* 
*               On exit:  the address of a variable containing the current
* 	      character position, after the read.
* 
*           4.  char eop
*               On entry:  the parameter delimiter character.
* 
*           5.  char eor
*               On entry:  the record delimiter character.
* 
* 6     Return Values, Error Indicators and Warnings
* 
*       The  address  of  the  structure  containing  the  Type122  (Tabulated
*       Cylinder)  entity  is  returned.
* 
* Functions referenced by ReadIges122() are:
*  gen_array1()
*  GetNextInteger()
*  GetNextReal()
*
* Functions that reference ReadIges122() are:
*  ReadIgesTabCyl()
*/

Type122 *ReadIges122(FILE *fp, char *record, int *ip, char eop, char eor)
{
  Type122 *surf;
  double val;
  int ival;

  GetNextInteger(fp, record, ip, eop, eor, 64, &ival);
  if (ival != 122)
    return (NULL);

  if (surf = (Type122 *)gen_array1(1, sizeof(Type122))) {
    GetNextInteger(fp, record, ip, eop, eor, 64, &ival);
    surf->de = ival;
    GetNextReal(fp, record, ip, eop, eor, 64, &val);
    surf->lx = val;
    GetNextReal(fp, record, ip, eop, eor, 64, &val);
    surf->ly = val;
    GetNextReal(fp, record, ip, eop, eor, 64, &val);
    surf->lz = val;
  }
  return (surf);
}

/********* ReadIges124() *********
* 1     Purpose
* 
*       Read an IGES Type124 (Transformation Matrix) entity.
* 
* 2     Specification
* 
*       #include "iges.h"
*       Type124 *ReadIges124(FILE *fp, char *record, int *ip, char *eop,
*                            char *eor);
* 
* 3     Description
* 
*       This function reads a Type124 (Transformation Matrix) entity from an
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
*           2.  char * record
*               On  entry:  the  address  of  a  character  string  containing
* 	      the  current  line,  before  the read.
* 
*               On exit:  the address of a character string containing the
* 	      current line, after the read.
* 
*           3.  int * ip
*               On entry:  the address of a variable containing the current
* 	      character position, before the read.
* 
*               On exit:  the address of a variable containing the current
* 	      character position, after the read.
* 
*           4.  char eop
*               On entry:  the parameter delimiter character.
* 
*           5.  char eor
*               On entry:  the record delimiter character.
* 
* 6     Return Values, Error Indicators and Warnings
* 
*       The  address  of  the  structure  containing  the  Type124
*       (Transformation  Matrix)  entity  isreturned.
* 
* 8     Further Comments
* 
*       The returned structure should be deallocated with free_tmatrix().
* 
* Functions referenced by ReadIges124() are:
*  alloc_tmatrix()
*  GetNextInteger()
*  GetNextReal()
*
* Functions that reference ReadIges124() are:
*  OpenIgesFileCB()
*/

Type124 *ReadIges124(FILE *fp, char *record, int *ip, char eop, char eor)
{
  Type124 *mat;
  int ctype;
  short i, j;

  GetNextInteger(fp, record, ip, eop, eor, 64, &ctype);
  if (ctype != 124)
    return (NULL);

  if (mat = alloc_tmatrix())
    for (i=0; i<3; i++) {
      for (j=0; j<3; j++)
	GetNextReal(fp, record, ip, eop, eor, 64, &mat->r[i][j]);
      GetNextReal(fp, record, ip, eop, eor, 64, &mat->t[i]);
    }
  return mat;
}

/********* ReadIges126() *********
* 1     Purpose
* 
*       Read an IGES Type126 (NURBS Curve) entity.
* 
* 2     Specification
* 
*       #include "iges.h"
*       ParCurv *ReadIges126(FILE *fp, char *record, int *ip, char *eop,
*                            char *eor);
* 
* 3     Description
* 
*       This function reads a Type126 (NURBS Curve) entity from an IGES file.
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
*           2.  char * record
*               On  entry:  the  address  of  a  character  string  containing
* 	      the  current  line,  before  the read.
* 
*               On exit:  the address of a character string containing the
* 	      current line, after the read.
* 
*           3.  int * ip
*               On entry:  the address of a variable containing the current
* 	      character position, before the read.
* 
*               On exit:  the address of a variable containing the current
* 	      character position, after the read.
* 
*           4.  char eop
*               On entry:  the parameter delimiter character.
* 
*           5.  char eor
*               On entry:  the record delimiter character.
* 
* 6     Return Values, Error Indicators and Warnings
* 
*       The address of the structure containing the Type126 (NURBS Curve) entity
*       is returned.
* 
* 8     Further Comments
* 
*       The returned structure should be deallocated with free_egeom().
* 
* Functions referenced by ReadIges126() are:
*  egeomalloc1()
*  GetNextInteger()
*  GetNextReal()
*
* Functions that reference ReadIges126() are:
*  ReadIgesCos()
*  ReadIgesCurv()
*  ReadIgesTrim()
*/

ParCurv *ReadIges126(FILE *fp, char *record, int *ip, char eop, char eor)
{
  ParCurv *egeom;
  double dummy;
  int degree, index;
  int i;
  
  GetNextInteger(fp, record, ip, eop, eor, 64, &index);
  if (index != 126)
    return (NULL);

  GetNextInteger(fp, record, ip, eop, eor, 64, &index);
  GetNextInteger(fp, record, ip, eop, eor, 64, &degree);
  if (egeom = egeomalloc1(degree+1, index+1)) {
    egeom->type = PCurveOpen;

    GetNextInteger(fp, record, ip, eop, eor, 64, &index); /* ignore next */
    GetNextInteger(fp, record, ip, eop, eor, 64, &index); /* 3 integers */
    GetNextInteger(fp, record, ip, eop, eor, 64, &index);

    GetNextInteger(fp, record, ip, eop, eor, 64, &index);
    if (index)
      egeom->type = PCurvePer;

    for (i=0; i<egeom->order + egeom->ncontpts; i++)
      GetNextReal(fp, record, ip, eop, eor, 64, &egeom->knots[i]);

    for (i=0; i<egeom->ncontpts; i++)
      GetNextReal(fp, record, ip, eop, eor, 64, &egeom->contpts[i]->w);

    for (i=0; i<egeom->ncontpts; i++) {
      GetNextReal(fp, record, ip, eop, eor, 64, &egeom->contpts[i]->x);
      GetNextReal(fp, record, ip, eop, eor, 64, &egeom->contpts[i]->y);
      GetNextReal(fp, record, ip, eop, eor, 64, &egeom->contpts[i]->z);

/* IGES stores control points as (xw,yw,zw,w) */
/* Praxiteles wants to see (x,y,z,w) */

      egeom->contpts[i]->x *= egeom->contpts[i]->w;
      egeom->contpts[i]->y *= egeom->contpts[i]->w;
      egeom->contpts[i]->z *= egeom->contpts[i]->w;
    }

    for (i=0; i<5; i++)                                  /* ignore last */
      GetNextReal(fp, record, ip, eop, eor, 64, &dummy); /* five reals */
  }
  return (egeom);
}

/********* ReadIges128() *********
* 1     Purpose
* 
*       Read an IGES Type128 (NURBS Surface) entity.
* 
* 2     Specification
* 
*       #include "iges.h"
*       ParSurf *ReadIges128(FILE *fp, char *record, int *ip, char *eop,
*                            char *eor);
* 
* 3     Description
* 
*       This function reads a Type128 (NURBS Surface) entity from an IGES file.
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
*           2.  char * record
*               On  entry:  the  address  of  a  character  string  containing
* 	      the  current  line,  before  the read.
* 
*               On exit:  the address of a character string containing the
* 	      current line, after the read.
* 
*           3.  int * ip
*               On entry:  the address of a variable containing the current
* 	      character position, before the read.
* 
*               On exit:  the address of a variable containing the current
* 	      character position, after the read.
* 
*           4.  char eop
*               On entry:  the parameter delimiter character.
* 
*           5.  char eor
*               On entry:  the record delimiter character.
* 
* 6     Return Values, Error Indicators and Warnings
* 
*       The address of the structure containing the Type128 (NURBS Surface)
*       entity is returned.
* 
* 8     Further Comments
* 
*       The returned structure should be deallocated with free_fgeom().
* 
* Functions referenced by ReadIges128() are:
*  fgeomalloc1()
*  GetNextInteger()
*  GetNextReal()
*
* Functions that reference ReadIges128() are:
*  ReadIgesSurf()
*  ReadIgesTrim()
*/

ParSurf *ReadIges128(FILE *fp, char *record, int *ip, char eop, char eor)
{
  ParSurf *fgeom;
  double dummy;
  int udegree, vdegree, uindex, vindex;
  int i, j;

  GetNextInteger(fp, record, ip, eop, eor, 64, &vindex);
  if (vindex != 128)
    return (NULL);

  GetNextInteger(fp, record, ip, eop, eor, 64, &uindex);
  GetNextInteger(fp, record, ip, eop, eor, 64, &vindex);
  GetNextInteger(fp, record, ip, eop, eor, 64, &udegree);
  GetNextInteger(fp, record, ip, eop, eor, 64, &vdegree);

  if (fgeom = fgeomalloc1(udegree+1, vdegree+1, uindex+1, vindex+1)) {
    fgeom->type = PSurfaceOpen;

    GetNextInteger(fp, record, ip, eop, eor, 64, &vindex); /* ignore next */
    GetNextInteger(fp, record, ip, eop, eor, 64, &vindex); /* five integers */
    GetNextInteger(fp, record, ip, eop, eor, 64, &vindex);
    GetNextInteger(fp, record, ip, eop, eor, 64, &vindex);
    GetNextInteger(fp, record, ip, eop, eor, 64, &vindex);

    for (i=0; i<fgeom->uorder + fgeom->ucontpts; i++)
      GetNextReal(fp, record, ip, eop, eor, 64, &fgeom->uknots[i]);

    for (i=0; i<fgeom->vorder + fgeom->vcontpts; i++)
      GetNextReal(fp, record, ip, eop, eor, 64, &fgeom->vknots[i]);

    for (j=0; j<fgeom->vcontpts; j++)
      for (i=0; i<fgeom->ucontpts; i++)
	GetNextReal(fp, record, ip, eop, eor, 64, &fgeom->contpts[i][j]->w);

    for (j=0; j<fgeom->vcontpts; j++)
      for (i=0; i<fgeom->ucontpts; i++) {
	GetNextReal(fp, record, ip, eop, eor, 64, &fgeom->contpts[i][j]->x);
	GetNextReal(fp, record, ip, eop, eor, 64, &fgeom->contpts[i][j]->y);
	GetNextReal(fp, record, ip, eop, eor, 64, &fgeom->contpts[i][j]->z);

/* IGES stores control points as (xw,yw,zw,w) */
/* Praxiteles wants to see (x,y,z,w) */

	fgeom->contpts[i][j]->x *= fgeom->contpts[i][j]->w;
	fgeom->contpts[i][j]->y *= fgeom->contpts[i][j]->w;
	fgeom->contpts[i][j]->z *= fgeom->contpts[i][j]->w;
      }

    for (i=0; i<4; i++)                                  /* ignore last */
      GetNextReal(fp, record, ip, eop, eor, 64, &dummy); /* four reals */
  }
  return (fgeom);
}

Type132 *ReadIges132(FILE *fp, char *record, int *ip, char eop, char eor)
{
  Type132 *connect;
  int ival;

  GetNextInteger(fp, record, ip, eop, eor, 64, &ival);
  if (ival != 132)
    return (NULL);

  connect = (Type132 *)calloc(1, sizeof(Type132));

  GetNextReal(fp, record, ip, eop, eor, 64, &connect->x);
  GetNextReal(fp, record, ip, eop, eor, 64, &connect->y);
  GetNextReal(fp, record, ip, eop, eor, 64, &connect->z);

  GetNextInteger(fp, record, ip, eop, eor, 64, &connect->ptr);

  GetNextInteger(fp, record, ip, eop, eor, 64, &connect->tf);
  GetNextInteger(fp, record, ip, eop, eor, 64, &connect->ff);

  GetNextString(fp, record, ip, eop, eor, 64, &connect->cid);
  GetNextInteger(fp, record, ip, eop, eor, 64, &connect->pttcid);

  GetNextString(fp, record, ip, eop, eor, 64, &connect->cfn);
  GetNextInteger(fp, record, ip, eop, eor, 64, &connect->pttcfn);

  GetNextInteger(fp, record, ip, eop, eor, 64, &connect->cpid);
  GetNextInteger(fp, record, ip, eop, eor, 64, &connect->fc);
  GetNextInteger(fp, record, ip, eop, eor, 64, &connect->sf);

  GetNextInteger(fp, record, ip, eop, eor, 64, &connect->psfi);

  return (connect);
}

/********* ReadIges142() *********
* 1     Purpose
* 
*       Read an IGES Type142 (Curve on Parametric Surface) entity.
* 
* 2     Specification
* 
*       #include "iges.h"
*       Type142 *ReadIges142(FILE *fp, char *record, int *ip, char *eop,
*                            char *eor);
* 
* 3     Description
* 
*       This function reads a Type142 (Curve on Parametric Surface) entity
*       from an IGES file.
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
*           2.  char * record
*               On  entry:  the  address  of  a  character  string  containing
* 	      the  current  line,  before  the read.
* 
*               On exit:  the address of a character string containing the
* 	      current line, after the read.
* 
*           3.  int * ip
*               On entry:  the address of a variable containing the current
* 	      character position, before the read.
* 
*               On exit:  the address of a variable containing the current
* 	      character position, after the read.
* 
*           4.  char eop
*               On entry:  the parameter delimiter character.
* 
*           5.  char eor
*               On entry:  the record delimiter character.
* 
* 6     Return Values, Error Indicators and Warnings
* 
*       The address of the structure containing the Type142 (Curve on
*       Parametric Surface) entity is returned.
* 
* Functions referenced by ReadIges142() are:
*  gen_array1()
*  GetNextInteger()
*
* Functions that reference ReadIges142() are:
*  ReadIgesCos()
*  ReadIgesTrim()
*/

Type142 *ReadIges142(FILE *fp, char *record, int *ip, char eop, char eor)
{
  Type142 *cs;
  int ival;

  GetNextInteger(fp, record, ip, eop, eor, 64, &ival);
  if (ival != 142)
    return (NULL);

  if (cs = (Type142 *)gen_array1(1, sizeof(Type142))) {
    GetNextInteger(fp, record, ip, eop, eor, 64, &ival);
    cs->crtn = ival;
    GetNextInteger(fp, record, ip, eop, eor, 64, &ival);
    cs->sptr = ival;
    GetNextInteger(fp, record, ip, eop, eor, 64, &ival);
    cs->bptr = ival;
    GetNextInteger(fp, record, ip, eop, eor, 64, &ival);
    cs->cptr = ival;
    GetNextInteger(fp, record, ip, eop, eor, 64, &ival);
    cs->pref = ival;
  }
  return (cs);
}

/********* ReadIges144() *********
* 1     Purpose
* 
*       Read an IGES Type144 (Trimmed (Parametric) Surface) entity.
* 
* 2     Specification
* 
*       #include "iges.h"
*       Type144 *ReadIges144(FILE *fp, char *record, int *ip, char *eop,
*                            char *eor);
* 
* 3     Description
* 
*       This function reads a Type144 (Trimmed (Parametric) Surface) entity
*       from an IGES file.
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
*           2.  char * record
*               On  entry:  the  address  of  a  character  string  containing
* 	      the  current  line,  before  the read.
* 
*               On exit:  the address of a character string containing the
* 	      current line, after the read.
* 
*           3.  int * ip
*               On entry:  the address of a variable containing the current
* 	      character position, before the read.
* 
*               On exit:  the address of a variable containing the current
* 	      character position, after the read.
* 
*           4.  char eop
*               On entry:  the parameter delimiter character.
* 
*           5.  char eor
*               On entry:  the record delimiter character.
* 
* 6     Return Values, Error Indicators and Warnings
* 
*       The address of the structure containing the Type144 (Curve on
*       Parametric Surface) entity is returned.
* 
* Functions referenced by ReadIges144() are:
*  gen_array1()
*  GetNextInteger()
*
* Functions that reference ReadIges144() are:
*  ReadIgesTrim()
*/

Type144 *ReadIges144(FILE *fp, char *record, int *ip, char eop, char eor)
{
  Type144 *cs;
  int i, ival;

  GetNextInteger(fp, record, ip, eop, eor, 64, &ival);
  if (ival != 144)
    return (NULL);

  if (cs = (Type144 *)gen_array1(1, sizeof(Type144))) {
    GetNextInteger(fp, record, ip, eop, eor, 64, &ival);
    cs->pts = ival;
    GetNextInteger(fp, record, ip, eop, eor, 64, &ival);
    cs->n1 = ival;
    GetNextInteger(fp, record, ip, eop, eor, 64, &ival);
    cs->n2 = ival;
    GetNextInteger(fp, record, ip, eop, eor, 64, &ival);
    cs->pt0 = ival;
    cs->pti = (int *)gen_array1(cs->n2, sizeof(int));
    for (i=0; i<cs->n2; i++) {
      GetNextInteger(fp, record, ip, eop, eor, 64, &ival);
      cs->pti[i] = ival;
    }
  }
  return (cs);
}

ParCurv *ReadIges102Curv(FILE *fp, char *record, int *ip, char eop, char eor,
			 Type102 *composite)
{
  ParCurv *egeom, **egeoms;
  PowCurv *pgeom;
  Type100 *arc;
  Type104 *conic;
  Type110 *lin;
  Type116 *point;
  Type132 *connect;
  Type124 *xm = NULL;
  vector **v;
  double errors[2], err, factor = 0.1, maxError, PERT, r;
  char line[80], token[9];
  int ent, i, j, k, join_key = 2, n, nc = 1, parm, POS_FLAG = 0, xform;

  egeoms = (ParCurv **)gen_array1(composite->n, sizeof(ParCurv *));

  token[8] = '\0';
  for (i=j=0; i<composite->n; i++) {
    SkipToThis(fp, record, ip, composite->de[i], 'D');
    sscanf(record, "%8d%8d%8c%8c%8c%8c%8c", &ent, &parm, line, line, line,
	   line, token);
    xform = atoi(token);

    if (xform) {
      SkipToThis(fp, record, ip, xform, 'D'); /* directory entry of matrix */
      sscanf(&record[8], "%8c", token);     /* pointer to parameter record */
      k = atoi(token);
      SkipToThis(fp, record, ip, k, 'P');   /* go to parameter record */
      if (xm)
	free_tmatrix(xm);
      xm = ReadIges124(fp, record, ip, eop, eor);   /* read the matrix */
    }

    switch (ent) {
    case 100:        /* circular arc */
      SkipToThis(fp, record, ip, parm, 'P');

      arc = ReadIges100(fp, record, ip, eop, eor);

      v = vec_array1(3);
      v[0]->x = arc->x1;
      v[0]->y = arc->y1;
      v[0]->z = arc->zt;
      v[1]->x = arc->x2;
      v[1]->y = arc->y2;
      v[1]->z = arc->zt;
      v[2]->x = arc->x3;
      v[2]->y = arc->y3;
      v[2]->z = arc->zt;
      r = sqrt((v[1]->x-v[0]->x)*(v[1]->x-v[0]->x) +
	       (v[1]->y-v[0]->y)*(v[1]->y-v[0]->y));
      egeoms[j] = arc_gen(v, r);
      free_varray1(v, 3);
      free(arc);

      if (xform)
	TransformIgesVect(xm, egeoms[j]->ncontpts, egeoms[j]->contpts);
      j++;
      break;

    case 104:        /* conic arc */
printf("\n>>>   conic arc\n\n");fflush(stdout);
      SkipToThis(fp, record, ip, parm, 'P');

      conic = ReadIges104(fp, record, ip, eop, eor);

      /* for the time being, (crudely) approximate the arc with a
       * linear curve between the starting and ending points */

      egeoms[j] = egeomalloc1(2, 2);
      egeoms[j]->knots[0] = egeoms[j]->knots[1] = 0.0;
      egeoms[j]->knots[2] = egeoms[j]->knots[3] = 1.0;
      egeoms[j]->contpts[0]->x = conic->x1;
      egeoms[j]->contpts[0]->y = conic->y1;
      egeoms[j]->contpts[0]->z = conic->zt;
      egeoms[j]->contpts[1]->x = conic->x2;
      egeoms[j]->contpts[1]->y = conic->y2;
      egeoms[j]->contpts[1]->z = conic->zt;

      free(conic);

      if (xform)
	TransformIgesVect(xm, egeoms[j]->ncontpts, egeoms[j]->contpts);
      j++;
      break;

    case 110:        /* line */
      SkipToThis(fp, record, ip, parm, 'P');

      lin = ReadIges110(fp, record, ip, eop, eor);
      v = vec_array1(2);
      v[0]->x = lin->x1;
      v[0]->y = lin->y1;
      v[0]->z = lin->z1;
      v[0]->w = 1.0;
      v[1]->x = lin->x2;
      v[1]->y = lin->y2;
      v[1]->z = lin->z2;
      v[1]->w = 1.0;
      egeoms[j] = line_gen(v);
      free_varray1(v, 2);
      free(lin);

      if (xform)
	TransformIgesVect(xm, egeoms[j]->ncontpts, egeoms[j]->contpts);
      j++;
      break;

    case 112:        /* parametric spline */
      SkipToThis(fp, record, ip, parm, 'P');

      pgeom = ReadIges112(fp, record, ip, eop, eor);
      PowCurv_error(pgeom, errors);
      egeoms[j] = PowCurv_to_ParCurv(pgeom, NULL, &maxError, nc);
      free_pgeom(pgeom);

      if (xform)
	TransformIgesVect(xm, egeoms[j]->ncontpts, egeoms[j]->contpts);
      j++;
      break;

    case 116:        /* point */
      SkipToThis(fp, record, ip, parm, 'P');

      point = ReadIges116(fp, record, ip, eop, eor);
      free(point);
      break;

    case 126:        /* NURBS */
      SkipToThis(fp, record, ip, parm, 'P');

      egeoms[j] = ReadIges126(fp, record, ip, eop, eor);

      if (xform)
	TransformIgesVect(xm, egeoms[j]->ncontpts, egeoms[j]->contpts);
      j++;
      break;

    case 132:        /* connect point */
      SkipToThis(fp, record, ip, parm, 'P');

      connect = ReadIges132(fp, record, ip, eop, eor);
      free(connect->cid);
      free(connect->cfn);
      free(connect);
      break;

    default:
printf(">>>   unknown Entity Type: %d\n", ent);fflush(stdout);	  
    }
  }
  free(composite->de);
  free(composite);
  if (xm)
    free_tmatrix(xm);

  if (j > 1) {
    egeom = construct_bspl(egeoms, j, POS_FLAG, join_key);
    PERT = max_pert_cv(egeoms, j, egeom, err, factor);
    knot_rem_cv(egeoms, j, egeom);
    perturb_knots_cv(egeoms, j, egeom, PERT);
  }
  else
    egeom = copyegeom(egeoms[0], NULL);

  for (i=0; i<j; i++)
    free_egeom(egeoms[i]);
  free_garray1((char *)egeoms);

  return egeom;
}

/* Copyright (C) 1992 Massachusetts Institute of Technology, Cambridge, MA
   All rights reserved
*/

/* print.c */

/* PrintIgesDirectory()
 * PrintIgesGlobal()
 * PrintIgesParamter()
 * PrintIgesTerminate()
*/

#include <stdio.h>
#include "iges.h"

/********* PrintIgesDirectory() *********
* 1     Purpose
* 
*       Print IGES file Directory structure.
* 
* 2     Specification
* 
*       #include "iges.h"
*       void PrintIgesDirectory(DirStruct *dir, int level);
* 
* 3     Description
* 
*       This function prints to the standard output unit the contents of a
*       structure containing the information from an IGES file Directory
*        section.
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
*           1.  DirStruct * dir
*               On entry:  the address of the Directory structure.
*
*           2.  int level
*               On entry:  if 0, print short directory entry; if 1, print
*                          the complete entry.
*/
 
 void PrintIgesDirectory(DirStruct *dir, int level)
 {
   static char *blankStatus[] = {
    "Visible", "Blanked"
  };
  static char *subordinateSwitch[] = {
    "Independent", "Physically dependent", "Logically dependent",
    "Physically and logically dependent"
  };
  static char *entityUseFlag[] = {
    "Geometry", "Annotation", "Definition", "Other", "Logical/Positional",
    "2D parametric"
  };
  static char *hierarchy[] = {
    "Global top down", "Global defer", "Use hierarchy property"
  };

  printf(" DE %d, Entity type: %d", dir->DE, dir->type);
  switch (dir->type) {
  case 100:
    printf(" (Circular arc)\n");
    break;
  case 102:
    printf(" (Composite curve)\n");
    break;
  case 106:
    printf(" (Copious data)\n");
    break;
  case 108:
    printf(" (Plane)\n");
    break;
  case 110:
    printf(" (Line)\n");
    break;
  case 112:
    printf(" (Parametric spline curve)\n");
    break;
  case 114:
    printf(" (Parametric spline surface)\n");
    break;
  case 116:
    printf(" (Point)\n");
    break;
  case 118:
    printf(" (Ruled surface)\n");
    break;
  case 120:
    printf(" (Surface of revolution)\n");
    break;
  case 122:
    printf(" (Tabulated cylinder)\n");
    break;
  case 124:
    printf(" (Tranformation matrix)\n");
    break;
  case 126:
    printf(" (NURBS curve)\n");
    break;
  case 128:
    printf(" (NURBS surface)\n");
    break;
  case 142:
    printf(" (Curve on parametric surface)\n");
    break;
  case 144:
    printf(" (Trimmed (parametric) surface)\n");
    break;
  default:
    printf(" (Unknown)\n");
  }
  printf("  Parameter data:              %d\n", dir->parameter);
  if (level) {
    printf("  Structure:                   %d\n", dir->structure);
    printf("  Line font pattern:           %d\n", dir->pattern);
    printf("  Level:                       %d\n", dir->level);
    printf("  View:                        %d\n", dir->view);
    printf("  Transformation matrix:       %d\n", dir->matrix);
    printf("  Label display associativity: %d\n", dir->associativity);
    printf("  Blank status:                \"%s\"\n", blankStatus[dir->blank]);
    printf("  Subordinate entity switch:   \"%s\"\n",
	   subordinateSwitch[dir->subordinate]);
    printf("  Entity use flag:             \"%s\"\n", entityUseFlag[dir->use]);
    printf("  Hierarchy:                   \"%s\"\n",
	   hierarchy[dir->hierarchy]);
    printf("  Line weight:                 %d\n", dir->weight);
    printf("  Color:                       %d\n", dir->color);
  }
  printf("  Parameter line count:        %d\n", dir->count);
  printf("  Form number:                 %d\n", dir->form);
  printf("  Entity label:                \"%s\"\n", dir->label);
  printf("  Subscript number:            %d\n", dir->subscript);
}

/********* PrintIgesGlobal() *********
* 1     Purpose
* 
*       Print IGES file Global section.
* 
* 2     Specification
* 
*       #include "iges.h"
*       void PrintIgesGlobal(GlobalStruct *global);
* 
* 3     Description
* 
*       This function writes the contents of an IGES file Global structure to
*       the standard output unit.
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
*           1.  GlobalStruct * global
*               On entry:  the address of the structure containing the Global
* 	      section.
*/

void PrintIgesGlobal(GlobalStruct *global)
{
  static char *unitFlag[] = {
    "Inches", "Millimeters", "<See units>", "Feet", "Miles", "Meters",
    "Kilometers", "Mils", "Microns", "centimeters", "Microinches"
  };
  static char *version[] = {
    "1.0", "ANSI Y14.26M-1981", "2.0", "3.0", "ASME/ANSI Y14.26M-1987",
    "4.0", "ASME Y14.26M-1989", "5.0"
  };
  static char *standard[] = {
    "None", "ISO", "AFNOR", "ANSI", "BSI", "CSA", "DIN", "JIS"
  };

  printf(" Parameter delimeter: \"%c\"\n", global->parmDelim);
  printf(" Record delimeter:    \"%c\"\n", global->recDelim);
  printf(" Sending system ID:   \"%s\"\n", global->sendID);
  printf(" File name:           \"%s\"\n", global->fileName);
  printf(" System ID:           \"%s\"\n", global->systemID);
  printf(" Preprocessor:        \"%s\"\n", global->preprocVers);
  printf(" Bits per integer:    %d\n", global->bits);
  printf(" Single magnitude:    %d\n", global->singleMagnitude);
  printf(" Single precision:    %d\n", global->singlePrecision);
  printf(" Double magnitude:    %d\n", global->doubleMagnitude);
  printf(" Double precision:    %d\n", global->doublePrecision);
  printf(" Receiving system ID: \"%s\"\n", global->receiveID);
  printf(" Model space scale:   %g\n", global->scale);
  printf(" Unit flag:           \"%s\"\n", unitFlag[global->unit-1]);
  printf(" Units:               \"%s\"\n", global->units);
  printf(" Maximum gradations:  %d\n", global->gradations);
  printf(" Maximum line weight: %g\n", global->weight);
  if (global->fileDate[0])
    printf(" File creation date:  \"%c%c/%c%c/%c%c %c%c:%c%c.%c%c\"\n",
	   global->fileDate[2], global->fileDate[3],
	   global->fileDate[4], global->fileDate[5],
	   global->fileDate[0], global->fileDate[1],
	   global->fileDate[7], global->fileDate[8],
	   global->fileDate[9], global->fileDate[10],
	   global->fileDate[11], global->fileDate[12]);
  else
    printf(" File modified date:  \"Unknown\"\n");
  printf(" Minimum resolution:  %g\n", global->resolution);
  printf(" Maximum coordinate:  %g\n", global->coordinate);
  printf(" Author's name:       \"%s\"\n", global->author);
  printf(" Organization name:   \"%s\"\n", global->org);
  printf(" IGES version:        \"%s\"\n", version[global->version-1]);
  printf(" Drafting standard:   \"%s\"\n", standard[global->standard]);
  if (global->modelDate[0])
    printf(" Model creation date: \"%c%c/%c%c/%c%c %c%c:%c%c.%c%c\"\n",
	   global->modelDate[2], global->modelDate[3],
	   global->modelDate[4], global->modelDate[5],
	   global->modelDate[0], global->modelDate[1],
	   global->modelDate[7], global->modelDate[8],
	   global->modelDate[9], global->modelDate[10],
	   global->modelDate[11], global->modelDate[12]);
  else
    printf(" Model modified date: \"Unknown\"\n");
}

void PrintIgesParameter(DirStruct *dir, FILE *fp, char *record, int *ip,
			char parmDelim, char recDelim)
{
  Type100 *circle;             /* Type 100 (circular arc) */
  Type102 *composite;          /* Type 102 (composite curve) */
  Type106 *copious;            /* Type 106 (copious data) */
  Type108 *plane;              /* Type 108 (plane) */
  Type110 *line;               /* Type 110 (line) */
  PowCurv *powCurv;            /* Type 112 (parametric spline curve) */
  PowSurf *powSurf;            /* Type 114 (parametric spline surface) */
  Type116 *point;              /* Type 116 (point) */
  Type118 *ruled;              /* Type 118 (ruled surface) */
  Type120 *surfRev;            /* Type 120 (surface of revolution) */
  Type122 *tabCyl;             /* Type 122 (tabulated cylinder) */
  Type124 *transf;             /* Type 124 (transformation matrix) */
  ParCurv *curv;               /* Type 126 (NURBS curve) */
  ParSurf *surf;               /* Type 128 (NURBS surface) */
  Type142 *cos;                /* Type 142 (curve on surface) */
  Type144 *trim;               /* Type 144 (trimmed surface) */
  int i;

  printf(" PD %d, Entity type: %d (DE %d)\n", dir->parameter, dir->type,
	 dir->DE);

  switch (dir->type) {
  case 100:                    /* circular arc */

    circle = ReadIges100(fp, record, ip, parmDelim, recDelim);
    printf("  Center point: %f %f\n", circle->x1, circle->y1);
    printf("  Start point:  %f %f\n", circle->x2, circle->y2);
    printf("  End point:    %f %f\n", circle->x3, circle->y3);
    free_garray1((char *)circle);
    break;

  case 102:                    /* composite curve */

    composite = ReadIges102(fp, record, ip, parmDelim, recDelim);
    printf("  Number of entities: %d\n", composite->n);
    printf("  DE:                 %d\n", composite->de[0]);
    for (i=1; i<composite->n; i++)
      printf("                      %d\n", composite->de[i]);
    free_iarray1(composite->de);
    free_garray1((char *)composite);
    break;

  case 106:                    /* copious data */
	    
    copious = ReadIges106(fp, record, ip, parmDelim, recDelim);
    printf("  Interpretation flag: %d\n", copious->ip);
    printf("  Number of tuples:    %d\n", copious->n);
    free_copious(copious);
    break;

  case 108:                    /* plane */

    plane = ReadIges108(fp, record, ip, parmDelim, recDelim);
    printf("  Coefficients A,B,C,D: %f %f %f", plane->a, plane->b, plane->c,
	   plane->d);
    free_garray1((char *)plane);
    break;

  case 110:                    /* line */

    line = ReadIges110(fp, record, ip, parmDelim, recDelim);
    printf("  Start point:     %f %f %f\n", line->x1, line->y1, line->z1);
    printf("  Terminate point: %f %f %f\n", line->x2, line->y2, line->z2);
    free_garray1((char *)line);
    break;

  case 112:                    /* parametric spline curve */

    powCurv = ReadIges112(fp, record, ip, parmDelim, recDelim);
    printf(" Number of segments: %d\n", powCurv->nsegmts);
    free_pgeom(powCurv);
    break;

  case 114:                    /* parametric spline surface */

    powSurf = ReadIges114(fp, record, ip, parmDelim, recDelim);
    printf("  Number of U segments: %d\n", powSurf->usegmts);
    printf("  Number of V segments: %d\n", powSurf->vsegmts);
    free_sgeom(powSurf);
    break;

  case 116:                    /* point */

    point = ReadIges116(fp, record, ip, parmDelim, recDelim);
    printf("  Coordinates: %f %f %f", point->x, point->y, point->z);
    free_garray1((char *)point);
    break;

  case 118:                    /* ruled surface */

    ruled = ReadIges118(fp, record, ip, parmDelim, recDelim);
    printf("  Pointer to DE of first curve:  %d\n", ruled->de1);
    printf("  Pointer to DE of second curve: %d\n", ruled->de2);
    printf("  Direction flag:                %d\n", ruled->dirflg);
    printf("  Developable surface flag:      %d\n", ruled->devflg);
    free_garray1((char *)ruled);
    break;

  case 120:                    /* surface of revolution */

    surfRev = ReadIges120(fp, record, ip, parmDelim, recDelim);
    printf("  Pointer to DE of axis:       %d\n", surfRev->de1);
    printf("  Pointer to DE of generatrix: %d\n", surfRev->de2);
    printf("  Start angle (radians):       %f\n", surfRev->start);
    printf("  Terminate angle (radians):   %f\n", surfRev->term);
    free_garray1((char *)surfRev);
    break;

  case 122:                    /* tabulated cylinder */

    tabCyl = ReadIges122(fp, record, ip, parmDelim,
			 recDelim);
    printf("  Pointer to DE of directrix:    %d\n", tabCyl->de);
    printf("  Terminate point of generatrix: %f %f %f\n", tabCyl->lx,
	   tabCyl->ly, tabCyl->lz);
    free_garray1((char *)tabCyl);
    break;

  case 124:                    /* transformation matrix */

    transf = ReadIges124(fp, record, ip, parmDelim,recDelim);
    printf("  Matrix: %f %f %f %f\n", transf->r[0][0], transf->r[0][1],
	   transf->r[0][2], transf->r[0][3], transf->t[0]);
    printf("          %f %f %f %f\n", transf->r[1][0], transf->r[1][1],
	   transf->r[1][2], transf->r[1][3], transf->t[1]);
    printf("          %f %f %f %f\n", transf->r[1][0], transf->r[2][1],
	   transf->r[2][2], transf->r[2][3], transf->t[2]);
    printf("          %f %f %f %f\n", transf->r[1][0], transf->r[3][1],
	   transf->r[3][2], transf->r[3][3], transf->t[3]);
    free_garray1((char *)transf);
    break;

  case 126:                    /* NURBS curve */

    curv = ReadIges126(fp, record, ip, parmDelim, recDelim);
    printf("  Upper index: %d\n", curv->ncontpts-1);
    printf("  Degree:      %d\n", curv->order-1);
    free_egeom(curv);
    break;

  case 128:                    /* NURBS surface */
	      
    surf = ReadIges128(fp, record, ip, parmDelim, recDelim);
    printf("  Upper index of U: %d\n", surf->ucontpts-1);
    printf("  Upper index of V: %d\n", surf->vcontpts-1);
    printf("  Degree of U:      %d\n", surf->uorder-1);
    printf("  Degree of V:      %d\n", surf->vorder-1);
    free_fgeom(surf);
    break;

  case 142:                    /* curve on parametric surface */
	      
    cos = ReadIges142(fp, record, ip, parmDelim, recDelim);
    printf("  Creation flag:             %d\n", cos->crtn);
    printf("  Pointer to DE of surface:  %d\n", cos->sptr);
    printf("  Pointer to DE of 2D curve: %d\n", cos->bptr);
    printf("  Pointer to DE of 3D curve: %d\n", cos->cptr);
    printf("  Representation preference: %d\n", cos->pref);
    free_garray1((char *)cos);
    break;

  case 144:                    /* trimmed (parametric) surface */

    trim = ReadIges144(fp, record, ip, parmDelim, recDelim);
    printf("  Pointer to DE of surface:           %d\n", trim->pts);
    printf("  Outer boundary flag:                %d\n", trim->n1);
    printf("  Number of inner boundaries:         %d\n", trim->n2);
    printf("  Pointer to DE of outer boundary:    %d\n", trim->pt0);
    if (trim->n2) {
      printf("  Pointer to DEs of inner boundaries: %d\n", trim->pti[0]);
      for (i=1; i<trim->n2; i++)
	printf("                                      %d\n", trim->pti[i]);
    }
    free_garray1((char *)trim->pti);
    free_garray1((char *)trim);
    break;

  default:                     /* unsupported type */

    printf("  Unknown type; don't know how to read!\n", dir->type);
    break;
  }
}

/********* PrintIgesTerminate() *********
* 1     Purpose
* 
*       Print IGES file Terminate structure.
* 
* 2     Specification
* 
*       #include "iges.h"
*       void PrintIgesTerminate(TermStruct *term);
* 
* 3     Description
* 
*       This function prints to the standard output unit the contents of a
*       structure containing the information from an IGES file Terminate
*       section.
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
*           1.  TermStruct * term
*               On entry:  the address of the Terminate structure.
*/

void PrintIgesTerminate(TermStruct *term)
{
  printf(" Start:     %d\n", term->start);
  printf(" Global:    %d\n", term->global);
  printf(" Directory: %d\n", term->directory);
  printf(" Parameter: %d\n", term->parameter);
}

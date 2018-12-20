/* Copyright (C) 1992 Massachusetts Institute of Technology, Cambridge, MA
   All rights reserved
*/

/* encode.c */

/* AddNextInteger()
 * AddNextReal()
 * AddNextString()
 * WriteIgesDirectory()
 * WriteIgesGlobal()
 * WriteIgesRecord()
 * WriteIgesStart()
 * WriteIgesTerminate()
*/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "iges.h"

/********* AddNextInteger() *********
* 1     Purpose
* 
*       Write integer value to IGES file.
* 
* 2     Specification
* 
*       #include "iges.h"
*       void AddNextInteger(FILE *fp, int value, char delim, char *record,
*                           int iDir, char cSection, int *iGlobal, int flush);
* 
* 3     Description
* 
*       This function writes an integer value starting at the current character
*       position in the current record.  If insufficient space remains on the
*       current record, then the record is flushed to the file and a new record
*       is started.
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
*               On entry:  the file structure returned by fopen() for the open
* 	      IGES file.
* 
*           2.  int value
*               On entry:  the integer value to be written to the file.
* 
*           3.  char delim
*               On entry:  the parameter delimiter.
* 
*           4.  char * record
*               On  entry:  the address of the character string containing
* 	      the current record,  before writing the integer.
* 
*               On exit: the address of the character string containing the
* 	      current record, after writing the integer.
* 
*           5.  int iDir
*               On entry:  the directory index of the current entity.
* 
*           6.  char cSection
*               On entry:  flag indicating the current section type:
* 
*                 cSection     File_section
*                 --------     -----------------
*                    'S'       Start section
*                    'G'       Global section
*                    'D'       Directory section
*                    'P'       Parameter section
*                    'T'       Terminate section
* 
*           7.  int * iGlobal
*               On entry:  the current line number.
* 
*               On exit:  the updated line number.
* 
*           8.  int flush
*               On entry:  flag controlling the flushing of the current record.
* 	      If flag = 1 then flush the current record after writing the
* 	      integer value.  Otherwise do not flush the record.
* 
* Functions referenced by AddNextInteger() are:
*  WriteIgesRecord()
*
* Functions that reference AddNextInteger() are:
*  WriteIges106()
*  WriteIges110()
*  WriteIges112()
*  WriteIges114()
*  WriteIges118()
*  WriteIges120()
*  WriteIges122()
*  WriteIges126()
*  WriteIges128()
*  WriteIges142()
*  WriteIges144()
*  WriteIgesGlobal()
*/

void AddNextInteger(FILE *fp, int value, char delim, char *record,
		    int iDir, char cSection, int *iGlobal, int flush)
{
  char s[73];
  int maxCol;

  if (cSection != PRAX_IGES_PARAMETER)
    maxCol = 72;
  else
    maxCol = 63;

  if (strlen(record) + sprintf(s, "%d%c", value, delim) > maxCol) {
    WriteIgesRecord(fp, record, iDir, cSection, ++(*iGlobal));
    record[0] = '\0';
  }
  strcat(record, s);
  if (flush)
    WriteIgesRecord(fp, record, iDir, cSection, ++(*iGlobal));
}

/********* AddNextReal() *********
* 1     Purpose
* 
*       Write real value to IGES file.
* 
* 2     Specification
* 
*       #include "iges.h"
*       void AddNextReal(FILE *fp, double value, char delim, char *record,
*                        int iDir, char cSection, int *iGlobal, int flush);
* 
* 
* 3     Description
* 
*       This function writes a real value starting at the current character
*       position in the current record.  If insufficient space remains on the
*       current record, then the record is flushed to the file and a new
*       record is started.
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
*               On entry:  the file structure returned by fopen() for the open
* 	      IGES file.
* 
*           2.  double value
*               On entry:  the real value to be written to the file.
* 
*           3.  char delim
*               On entry:  the parameter delimiter.
* 
*           4.  char * record
*               On  entry:  the address of the character string containing the
* 	      current record,  before writing the real.
* 
*               On exit: the address of the character string containing the
* 	      current record, after writing the real.
* 
*           5.  int iDir
*               On entry:  the directory index of the current entity.
* 
*           6.  char cSection
*               On entry:  flag indicating the current section type:
* 
*                 cSection     File_section
*                 --------     -----------------
*                    'S'       Start section
*                    'G'       Global section
*                    'D'       Directory section
*                    'P'       Parameter section
*                    'T'       Terminate section
* 
*           7.  int * iGlobal
*               On entry:  the current line number.
* 
*               On exit:  the updated line number.
* 
*           8.  int flush
*               On entry:  flag controlling the flushing of the current record.
* 	      If flag = 1 then flush the current record after writing the
* 	      integer value.  Otherwise do not flush the record.
* 
* Functions referenced by AddNextReal() are:
*  WriteIgesRecord()
*
* Functions that reference AddNextReal() are:
*  WriteIges106()
*  WriteIges110()
*  WriteIges112()
*  WriteIges114()
*  WriteIges120()
*  WriteIges122()
*  WriteIges126()
*  WriteIges128()
*  WriteIgesGlobal()
*/

void AddNextReal(FILE *fp, double value, char delim, char *record,
		 int iDir, char cSection, int *iGlobal, int flush)
{
  char str[73], *s = NULL;
  int i, j = 0, maxCol, n;

  if (cSection != PRAX_IGES_PARAMETER)
    maxCol = 72;
  else
    maxCol = 63;

  for (n=0; n<sprintf(str, "%#.15G", value); n++)
    if (str[n] == 'E') {
      j = n;
      break;
    }
  for (i=n-1; i>2; i--)
    if (str[i] != '0')
      break;
    else if (str[i-1] == '.')
      break;
  i++;
  if (j) {
    for (; str[j] != '\0'; i++,j++)
      str[i] = str[j];
  }
  str[i] = delim;
  str[i+1] = '\0';
  if (strlen(record) + i+1 > maxCol) {
    WriteIgesRecord(fp, record, iDir, cSection, ++(*iGlobal));
    record[0] = '\0';
  }
  strcat(record, str);
  if (flush)
    WriteIgesRecord(fp, record, iDir, cSection, ++(*iGlobal));
}

/********* AddNextString() *********
* 1     Purpose
* 
*       Write character string to IGES file.
* 
* 2     Specification
* 
*       #include "iges.h"
*       void AddNextString(FILE *fp, char *value, char delim, char *record,
*                          int iDir, char cSection, int *iGlobal, int flush);
* 
* 3     Description
* 
*       This function writes an string value starting at the current character
*       position in the current record.  If the current record will not contain
*       the entire string, only the appropriate portion is written and the
*       record is flushed to the file.  The remainder of the string is placed
*       at the beginning of a new record.
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
*               On entry:  the file structure returned by fopen() for the open
* 	      IGES file.
* 
*           2.  char * value
*               On entry:  the address of the character string to be written to
* 	      the file.
* 
*           3.  char delim
*               On entry:  the parameter delimiter.
* 
*           4.  char * record
*               On  entry:  the address of the character string containing the
* 	      current record,  before writing the character string.
* 
*               On exit: the address of the character string containing the
* 	      current record, after writing the character string.
* 
*           5.  int iDir
*               On entry:  the directory index of the current entity.
* 
*           6.  char cSection
*               On entry:  flag indicating the current section type:
* 
*                 cSection     File_section
*                 --------     -----------------
*                    'S'       Start section
*                    'G'       Global section
*                    'D'       Directory section
*                    'P'       Parameter section
*                    'T'       Terminate section
* 
*           7.  int * iGlobal
*               On entry:  the current line number.
* 
*               On exit:  the updated line number.
* 
* 
*           8.  int flush
*               On entry:  flag controlling the flushing of the current record.
* 	      If flag = 1 then flush the current record after writing the
* 	      integer value.  Otherwise do not flush the record.
* 
* Functions referenced by AddNextString() are:
*  WriteIgesRecord()
*
* Functions that reference AddNextString() are:
*  WriteIgesGlobal()
*/

void AddNextString(FILE *fp, char *value, char delim, char *record,
		   int iDir, char cSection, int *iGlobal, int flush)
{
  char holl[8], s[73];
  int h, i, ir, is, maxCol;

  if (cSection != PRAX_IGES_PARAMETER)
    maxCol = 72;
  else
    maxCol = 63;

  ir = strlen(record);
  is = strlen(value);
  h = (int)log10((double)is) + 1;
  sprintf(holl, "%dH", is);
  if (ir > maxCol-1 - h) {  /* can't break string during Hollerith */
    WriteIgesRecord(fp, record, iDir, cSection, ++(*iGlobal));
    record[0] = '\0';
  }
  strcat(record, holl);     /* add Hollerith count */
  i = 0;
  while (is > 0) {
    if ((ir = strlen(record)) > maxCol-1) {
      WriteIgesRecord(fp, record, iDir, cSection, ++(*iGlobal));
      record[0] = '\0';
      ir = 0;
    }
    strncat(record, &value[i], maxCol-ir);
    i += (maxCol - ir);
    is -= (maxCol - ir);
  }
  if ((ir = strlen(record)) > maxCol-1) {   /* add delimeter */
    WriteIgesRecord(fp, record, iDir, cSection, ++(*iGlobal));
    record[0] = '\0';
    ir = 0;
  }
  record[ir++] = delim;
  record[ir] = '\0';
  if (flush)
    WriteIgesRecord(fp, record, iDir, cSection, ++(*iGlobal));
}

/********* WriteIgesDirectory() **********
* 1     Purpose
* 
*       Write IGES file Directory entry.
* 
* 2     Specification
* 
*       #include "iges.h"
*       void WriteIgesDirectory(FILE *fp, DirStruct *dir, int iDirectory);
* 
* 3     Description
* 
*       This function writes an IGES file Directory entry.
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
*           2.  DirStruct * dir
*               On entry:  the address of the structure containing the Directory
* 	      entry.
* 
*           3.  int iDirectory
*               On entry:  the line number of the first line of the Directory
* 	      entry.
* 
* Functions referenced by WriteIgesDirectory() are:
*  WriteIgesRecord()
*
* Functions that reference WriteIgesDirectory() are:
*  SaveIgesCos()
*  SaveIgesCurv()
*  SaveIgesGrid()
*  SaveIgesList()
*  SaveIgesMinDist()
*  SaveIgesSurf()
*  SaveIgesTrim()
*  SaveIgesUv()
*/

void WriteIgesDirectory(FILE *fp, DirStruct *dir, int iDirectory)
{
  char record[73];

  sprintf(record, "%8d%8d%8d%8d%8d%8d%8d%8d%02.2d%02.2d%02.2d%02.2d",
	  dir->type, dir->parameter, dir->structure, dir->pattern,
	  dir->level, dir->view, dir->matrix, dir->associativity,
	  dir->blank, dir->subordinate, dir->use, dir->hierarchy);
  WriteIgesRecord(fp, record, 0, PRAX_IGES_DIRECTORY, iDirectory);
  sprintf(record, "%8d%8d%8d%8d%8d                %8.8s%8d",
	  dir->type, dir->weight, dir->color, dir->count,
	  dir->form, dir->label, dir->subscript);
  WriteIgesRecord(fp, record, 0, PRAX_IGES_DIRECTORY, ++iDirectory);
}

/********* WriteIgesGlobal() *********
* 1     Purpose
* 
*       Write IGES file Global section.
* 
* 2     Specification
* 
*       #include "iges.h"
*       void WriteIgesGlobal(FILE *fp, GlobalStruct *global);
* 
* 3     Description
* 
*       This function writes an IGES file Global section.
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
*           2.  GlobalStruct * global
*               On entry:  the address of the structure containing the Global
* 	      section.
* 
* Functions referenced by WriteIgesGlobal() are:
*  AddNextInteger()
*  AddNextReal()
*  AddNextString()
*
* Functions that reference WriteIgesGlobal() are:
*  SaveIgesCos()
*  SaveIgesCurv()
*  SaveIgesGrid()
*  SaveIgesList()
*  SaveIgesMinDist()
*  SaveIgesSurf()
*  SaveIgesTrim()
*  SaveIgesUv()
*/

int WriteIgesGlobal(FILE *fp, GlobalStruct *global)
{
  char record[73];
  int iGlobal = 0;

  sprintf(record, "1H%c%c1H%c%c", global->parmDelim, global->parmDelim,
	  global->recDelim, global->parmDelim);
  AddNextString(fp, global->sendID, global->parmDelim, record, 0,
		PRAX_IGES_GLOBAL, &iGlobal, 0);
  AddNextString(fp, global->fileName, global->parmDelim, record, 0,
		PRAX_IGES_GLOBAL, &iGlobal, 0);
  AddNextString(fp, global->systemID, global->parmDelim, record, 0,
		PRAX_IGES_GLOBAL, &iGlobal, 0);
  AddNextString(fp, global->preprocVers, global->parmDelim,record, 0,
		PRAX_IGES_GLOBAL, &iGlobal, 0);
  AddNextInteger(fp, global->bits, global->parmDelim, record, 0,
		 PRAX_IGES_GLOBAL, &iGlobal, 0);
  AddNextInteger(fp, global->singleMagnitude, global->parmDelim, record, 0,
		 PRAX_IGES_GLOBAL, &iGlobal, 0);
  AddNextInteger(fp, global->singlePrecision, global->parmDelim, record, 0,
		 PRAX_IGES_GLOBAL, &iGlobal, 0);
  AddNextInteger(fp, global->doubleMagnitude, global->parmDelim, record, 0,
		 PRAX_IGES_GLOBAL, &iGlobal, 0);
  AddNextInteger(fp, global->doublePrecision, global->parmDelim, record, 0,
		 PRAX_IGES_GLOBAL, &iGlobal, 0);
  AddNextString(fp, global->receiveID, global->parmDelim, record, 0,
		PRAX_IGES_GLOBAL, &iGlobal, 0);
  AddNextReal(fp, global->scale, global->parmDelim, record, 0,
	      PRAX_IGES_GLOBAL, &iGlobal, 0);
  AddNextInteger(fp, global->unit, global->parmDelim, record, 0,
		 PRAX_IGES_GLOBAL, &iGlobal, 0);
  AddNextString(fp, global->units, global->parmDelim, record, 0,
		PRAX_IGES_GLOBAL, &iGlobal,0);
  AddNextInteger(fp, global->gradations, global->parmDelim, record, 0,
		 PRAX_IGES_GLOBAL, &iGlobal, 0);
  AddNextReal(fp, global->weight, global->parmDelim, record, 0,
	      PRAX_IGES_GLOBAL, &iGlobal, 0);
  AddNextString(fp, global->fileDate, global->parmDelim,record, 0,
		PRAX_IGES_GLOBAL, &iGlobal, 0);
  AddNextReal(fp, global->resolution, global->parmDelim, record, 0,
	      PRAX_IGES_GLOBAL, &iGlobal, 0);
  AddNextReal(fp, global->coordinate, global->parmDelim, record, 0,
	      PRAX_IGES_GLOBAL, &iGlobal, 0);
  AddNextString(fp, global->author, global->parmDelim,record, 0,
		PRAX_IGES_GLOBAL, &iGlobal, 0);
  AddNextString(fp, global->org, global->parmDelim,record, 0,
		PRAX_IGES_GLOBAL, &iGlobal, 0);
  AddNextInteger(fp, global->version, global->parmDelim, record, 0,
		 PRAX_IGES_GLOBAL, &iGlobal, 0);
  AddNextInteger(fp, global->standard, global->parmDelim, record, 0,
		 PRAX_IGES_GLOBAL, &iGlobal, 0);
  if (global->modelDate[0])
    AddNextString(fp, global->modelDate, global->recDelim, record, 0,
		  PRAX_IGES_GLOBAL, &iGlobal, 1);
  else
    AddNextString(fp, global->fileDate, global->recDelim, record, 0,
		  PRAX_IGES_GLOBAL, &iGlobal, 1);

  return (iGlobal);
}

/********** WriteIgesRecord() *********
* 1     Purpose
* 
*       Write IGES file record.
* 
* 2     Specification
* 
*       #include "iges.h"
*       void WriteIgesRecord(FILE *fp, char *record, int iDir, char cSection,
*                            int iStart);
* 
* 3     Description
* 
*       This function writes a line of an IGES file.
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
*               On entry:  the address of a character string containing the
* 	      line to be written to the file.
* 
*           3.  int iDir
*               On entry:  if currently in the Parameter section, the line
* 	      number of the first line of the directory entry with which
* 	      this record is associated; otherwise, this value is ignored.
* 
*           4.  char cSection
*               On entry:  flag specifying the current file section.
* 
*                 cSection     File_Section
*                 --------     -----------------
*                    'S'       Start section
*                    'G'       Global section
*                    'D'       Directory section
*                    'P'       Parameter section
*                    'T'       Terminate section
* 
* 
*           5.  int iStart
*               On entry:  the current line number.
* 
* Functions that reference WriteIgesRecord() are:
*  AddNextInteger()
*  AddNextReal()
*  AddNextString()
*  WriteIgesDirectory()
*  WriteIgesStart()
*  WriteIgesTerminate()
*/

void WriteIgesRecord(FILE *fp, char *record, int iDir, char cSection,
		     int iStart)
{
  if (cSection != PRAX_IGES_PARAMETER)
    fprintf(fp, "%-72.72s%c%7d\n", record, cSection, iStart);
  else
    fprintf(fp, "%-63.63s  %7d%c%7d\n", record, iDir, cSection, iStart);
}

/********* WriteIgesTerminate() *********
* 1     Purpose
* 
*       Write IGES file Terminate section.
* 
* 2     Specification
* 
*       #include "iges.h"
*       void WriteIgesTerminate(FILE *fp, Tint iStart, int iGlobal,
*                               int iDirectory, int iParameter);
* 
* 3     Description
* 
*       This function writes an IGES file Terminate section.
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
*           2.  int iStart
*               On entry:  the number of lines in the Start section.
* 
*           3.  int iGlobal
*               On entry:  the number of lines in the Global section.
* 
*           4.  int iDirectory
*               On entry:  the number of lines in the Directory section.
* 
*           5.  int iParameter
*               On entry:  the number of lines in the Parameter section.
* 
* Functions referenced by WriteIgesTerminate() are:
*  WriteIgesRecord()
*
* Functions that reference WriteIgesTerminate() are:
*  SaveIgesCos()
*  SaveIgesCurv()
*  SaveIgesGrid()
*  SaveIgesList()
*  SaveIgesMinDist()
*  SaveIgesSurf()
*  SaveIgesTrim()
*  SaveIgesUv()
*/

void WriteIgesTerminate(FILE *fp, int iStart, int iGlobal,
			int iDirectory, int iParameter)
{
  char record[33];
  sprintf(record, "S%7dG%7dD%7dP%7d", iStart, iGlobal, iDirectory,
	  iParameter);
  WriteIgesRecord(fp, record, 0, PRAX_IGES_TERMINATE, 1);
}

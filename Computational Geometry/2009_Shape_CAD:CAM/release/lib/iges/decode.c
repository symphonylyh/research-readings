/* Copyright (C) 1992 Massachusetts Institute of Technology, Cambridge, MA
   All rights reserved
*/

/* decode.c */

/* DecodeIgesDirectory()
 * DecodeIgesGlobal()
 * DecodeIgesTerminate()
 * GetNextInteger()
 * GetNextReal()
 * GetNextString()
 * ParseIgesRecord()
 * ReadNextIgesRecord()
 * SkipToNext()
 * SkipToThis()
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h>
#include <math.h>
#include "iges.h"

/********* DecodeIgesDirectory() *********
* 1     Purpose
* 
*       Read IGES file Directory entry.
* 
* 2     Specification
* 
*       #include "iges.h"
*       void DecodeIgesDirectory(char *record, int *first, DirStruct **dir);
* 
* 3     Description
* 
*       This function reads a Directory entry from an IGES file and
*       allocates and initializes a Directory structure.
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
*           1.  char * record
*               On  entry:  the address of a character string containing the
* 	      first or second line of a Directory entry.
* 
*           2.  int * first
*               On  entry:  flag indicating whether record is the first or
* 	      second line of the Directory entry.
* 
*               On exit:  if first is equal to 1 on entry, it is set to 0 on
* 	      exit.  If first is equal to 0 on entry, it is set to 1 on exit.
* 
*           3.  DirStruct ** dir
*               On exit: the address of a variable that holds the address of a
* 	      newly allocated DirStruct.  If first is equal to 1 on entry,
* 	      the structure is allocated and the fields associated with
*               the first line of the Directory entry are set.  If first is
* 	      equal to 0 on entry, the fields associated with the second
* 	      line of the Directory entry are set.
* 
* Functions referenced by DecodeIgesDirectory() are:
*  gen_array1()
*
* Functions that reference DecodeIgesDirectory() are:
*  GetIgesDirList()
*/

void DecodeIgesDirectory(char *record, int *first, DirStruct **dir)
{
  char token[9];
  int status;

  token[8] = '\0';   /* see IGES manual for directory format */
  if (*first) {      /* first line of directory entry */
    *dir = (DirStruct *)gen_array1(1, sizeof(DirStruct));
    sscanf(record, "%8c", token);
    (*dir)->type = atoi(token);
    sscanf(&record[8], "%8c", token);  /* directory records are fixed format */
    (*dir)->parameter = atoi(token);
    sscanf(&record[16], "%8c", token);
    (*dir)->structure = atoi(token);
    sscanf(&record[24], "%8c", token);
    (*dir)->pattern = atoi(token);
    sscanf(&record[32], "%8c", token);
    (*dir)->level = atoi(token);
    sscanf(&record[40], "%8c", token);
    (*dir)->view = atoi(token);
    sscanf(&record[48], "%8c", token);
    (*dir)->matrix = atoi(token);
    sscanf(&record[56], "%8c", token);
    (*dir)->associativity = atoi(token);
    sscanf(&record[64], "%8c", token);
    status = atoi(token);
    (*dir)->blank = status / 1000000;
    (*dir)->subordinate = (status % 1000000) / 10000;
    (*dir)->use = (status % 10000) / 100;
    (*dir)->hierarchy = status % 100;
    *first = 0;
  }
  else {             /* second line of directory entry */
    sscanf(&record[8], "%8c", token);
    (*dir)->weight = atoi(token);
    sscanf(&record[16], "%8c", token);
    (*dir)->color = atoi(token);
    sscanf(&record[24], "%8c", token);
    (*dir)->count = atoi(token);
    sscanf(&record[32], "%8c", token);
    (*dir)->form = atoi(token);
    sscanf(&record[56], "%8c", (*dir)->label);
    (*dir)->label[8] = '\0';
    sscanf(&record[64], "%8c", token);
    (*dir)->subscript = atoi(token);
    *first = 1;
  }
}

/********* DecodeIgesGlobal() *********
* 1     Purpose
* 
*       Read IGES file Global section.
* 
* 2     Specification
* 
*       #include "iges.h"
*       void DecodeIgesGlobal(FILE *fp, char *record, GlobalStruct **global);
* 
* 3     Description
* 
*       This function reads the Global section from an IGES file and
*       allocates and initializes a Global structure.
* 
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
*               On entry:  file pointer of the IGES file.
* 
*           2.  char * record
*               On  entry:  the  address  of  a  character  string  containing
* 	      the  first  line  of  the  Global section.
* 
*           3.  GlobalStruct ** global
*               On  exit:  the address of a variable that holds the address of
* 	      a newly allocated and initialized GlobalStruct.
* 
* Functions referenced by DecodeIgesGlobal() are:
*  gen_array1()
*  GetNextInteger()
*  GetNextReal()
*  GetNextString()
*
* Functions that reference DecodeIgesGlobal() are:
*  OpenIgesFileCB()
*/

void DecodeIgesGlobal(FILE *fp, char *record, GlobalStruct **global)
{
  int ip = 0;
  char *token;

  *global = (GlobalStruct *)gen_array1(1, sizeof(GlobalStruct));

/* if      ,,        then parameter  ',' ';' */
/* record  1Haa,     and record       a  ';' */
/* starts  ,1Hb,     delimeters      ','  b  */
/* with    1Haa1Hba  are              a   b  */

  while (record[ip] == ' ')
    ip++;
  if (record[ip] == '1') {                       /* '1Hax', see IGES 2.2.3.1 */
    (*global)->parmDelim = record[ip+2];         /* set delim to 'a' */
    ip += 3;                                     /* set pointer to 'x' */
    while (record[ip] != (*global)->parmDelim)   /* find next delimeter */
      ip++;
  }
  else
    (*global)->parmDelim = ',';
  ip += 1;

  while (record[ip] == ' ')
    ip++;
  if (record[ip] == '1') {                       /* '1Hbx' */
    (*global)->recDelim = record[ip+2];          /* set delim to 'b' */
    ip += 3;                                     /* set pointer to 'x' */
    while (record[ip] != (*global)->parmDelim)   /* find next delimeter */
      ip++;
  }
  else
    (*global)->recDelim = ';';
  ip += 1;

  (*global)->unit = 3;          /* units specified by parameter 15 */
  (*global)->version = 3;       /* IGES Version 2.0 */
  (*global)->standard = 0;      /* no drafting standard specified */
  (*global)->gradations = 1;    /* default line weight */
  if (!GetNextString(fp, record, &ip, (*global)->parmDelim,
		       (*global)->recDelim, 72, &(*global)->sendID))
    return;
  if (!GetNextString(fp, record, &ip, (*global)->parmDelim,
		       (*global)->recDelim, 72, &(*global)->fileName))
    return;
  if (!GetNextString(fp, record, &ip, (*global)->parmDelim,
		       (*global)->recDelim, 72, &(*global)->systemID))
    return;
  if (!GetNextString(fp, record, &ip, (*global)->parmDelim,
		       (*global)->recDelim, 72, &(*global)->preprocVers))
    return;
  if (!GetNextInteger(fp, record, &ip, (*global)->parmDelim,
		      (*global)->recDelim, 72, &(*global)->bits))
    return;
  if (!GetNextInteger(fp, record, &ip, (*global)->parmDelim,
		      (*global)->recDelim, 72, &(*global)->singleMagnitude))
    return;
  if (!GetNextInteger(fp, record, &ip, (*global)->parmDelim,
		      (*global)->recDelim, 72, &(*global)->singlePrecision))
    return;
  if (!GetNextInteger(fp, record, &ip, (*global)->parmDelim,
		      (*global)->recDelim, 72, &(*global)->doubleMagnitude))
    return;
  if (!GetNextInteger(fp, record, &ip, (*global)->parmDelim,
		      (*global)->recDelim, 72, &(*global)->doublePrecision))
    return;
  if (!GetNextString(fp, record, &ip, (*global)->parmDelim,
		       (*global)->recDelim, 72, &(*global)->receiveID))
    return;
  if (!GetNextReal(fp, record, &ip, (*global)->parmDelim, (*global)->recDelim,
		   72, &(*global)->scale))
    return;
  if (!GetNextInteger(fp, record, &ip, (*global)->parmDelim,
		      (*global)->recDelim, 72, &(*global)->unit))
    return;
  if ((*global)->unit < 1 || (*global)->unit > 11)
    (*global)->unit = 3;
  if (!GetNextString(fp, record, &ip, (*global)->parmDelim,
		       (*global)->recDelim, 72, &(*global)->units))
    return;
  if (!GetNextInteger(fp, record, &ip, (*global)->parmDelim,
		      (*global)->recDelim, 72, &(*global)->gradations))
    return;
  if ((*global)->gradations < 1)
    (*global)->gradations = 1;
  if (!GetNextReal(fp, record, &ip, (*global)->parmDelim, (*global)->recDelim,
		   72, &(*global)->weight))
    return;
  if (!GetNextString(fp, record, &ip, (*global)->parmDelim,
		     (*global)->recDelim, 72, &token))
    return;
  strcpy((*global)->fileDate, token);
  free(token);
  if (!GetNextReal(fp, record, &ip, (*global)->parmDelim, (*global)->recDelim,
		   72, &(*global)->resolution))
    return;
  if (!GetNextReal(fp, record, &ip, (*global)->parmDelim, (*global)->recDelim,
		   72, &(*global)->coordinate))
    return;
  if (!GetNextString(fp, record, &ip, (*global)->parmDelim,
		       (*global)->recDelim, 72, &(*global)->author))
    return;
  if (!GetNextString(fp, record, &ip, (*global)->parmDelim,
		     (*global)->recDelim, 72, &(*global)->org))
    return;
  if (!GetNextInteger(fp, record, &ip, (*global)->parmDelim,
		      (*global)->recDelim, 72, &(*global)->version))
    return;
  if ((*global)->version < 1 || (*global)->version > 10)
    (*global)->version = 3;
  if (!GetNextInteger(fp, record, &ip, (*global)->parmDelim,
		      (*global)->recDelim, 72, &(*global)->standard))
    return;
  if ((*global)->standard < 0 || (*global)->standard > 7)
    (*global)->standard = 0;
  if (!GetNextString(fp, record, &ip, (*global)->parmDelim,
		     (*global)->recDelim, 72, &token)) {
    return;
  }
  strcpy((*global)->modelDate, token);
  free(token);
}

/********* DecodeIgesTerminate() *********
* 1     Purpose
* 
*       Read IGES file Terminate section.
* 
* 2     Specification
* 
*       #include "iges.h"
*       TermStruct *DecodeIgesTerminate(char *record);
* 
* 3     Description
* 
*       This function reads the Terminate section from an IGES file and
*       allocates and initializes a Terminate structure.
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
*           1.  char * record
*               On entry:  the address of a character string containing the
* 	      first line of the Terminate section.
* 
* 6     Return Values, Error Indicators and Warnings
* 
*       The address of a newly allocated and initialized TermStruct is
*       returned.   If  insufficient dynamic memory is available, the
*       function returns the value NULL.
* 
* Functions referenced by DecodeIgesTerminate() are:
*  gen_array1()
*/

TermStruct *DecodeIgesTerminate(char *record)
{
  TermStruct *term;

  term = (TermStruct *)gen_array1(1, sizeof(TermStruct));
  sscanf(record, "S%7dG%7dD%7dP%7d", &term->start, &term->global,
	 &term->directory, &term->parameter);

  return (term);
}

/********* GetNextInteger() *********
* 1     Purpose
* 
*       Read integer value from IGES file.
* 
* 2     Specification
* 
*       #include "iges.h"
*       int  GetNextInteger(FILE *fp, char *record, int *ip, char eop,
*                           char eor, int maxCol, int *value);
* 
* 3     Description
* 
*       This function reads an integer value starting at the current character
*       position in the current record of the IGES file.
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
*           2.  char * record
*               On  entry:  the address of the character string containing the
* 	      current record,  before reading the integer.
* 
*               On  exit:  the  address  of  the  character  string  containing
* 	      the  current  record,  before reading the integer.
* 
*           3.  int * ip
*               On entry:  the current character position, before reading the
* 	      integer.
* 
*               On exit:  the current character position, after reading the
* 	      integer.
* 
*           4.  char eop
*               On entry:  the parameter delimiter.
* 
*           5.  char eor
*               On entry:  the record delimiter.
* 
*           6.  int maxCol
*               On entry:  the maximum character position of the line.
* 
*           7.  int * value
*               On entry:  the address of a variable to contain the integer
* 	      value.
* 
*               On exit:  the address of a variable containing the integer
* 	      value read from the file.
* 
* 6     Return Values, Error Indicators and Warnings
* 
*       If the integer read operation was successful, the value 1 is returned;
*       otherwise, the value 0 is returned.
* 
* Functions referenced by GetNextInteger() are:
*  free_garray1()
*  ParseIgesRecord()
*
* Functions that reference GetNextInteger() are:
*  DecodeIgesGlobal()
*  ReadIges106()
*  ReadIges110()
*  ReadIges112()
*  ReadIges114()
*  ReadIges118()
*  ReadIges120()
*  ReadIges122()
*  ReadIges124()
*  ReadIges126()
*  ReadIges128()
*  ReadIges142()
*  ReadIges144()
*/

int GetNextInteger(FILE *fp, char *record, int *ip, char eop, char eor,
		     int maxCol, int *value)
{
  char *token;

  if (token = ParseIgesRecord(fp, record, ip, eop, eor, PRAX_IGES_INTEGER,
			      maxCol)) {
    *value = atoi(token);   /* decode the current token on the record */
    free_garray1(token);
    return (1);   /* return success */
  } 
  else {
    *value = 0;
    return (0);   /* return failure */
  }
}

/********* GetNextReal() *********
* 1     Purpose
* 
*       Read real value from IGES file.
* 
* 2     Specification
* 
*       #include "iges.h"
*       int  GetNextReal(FILE *fp, char *record, int *ip, char eop, char eor,
*                        int  maxCol,  double *value);
* 
* 3     Description
* 
*       This function reads an real value starting at the current character
*       position in the current record of the IGES file.
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
*           2.  char * record
*               On  entry:  the address of the character string containing
* 	      the current record,  before reading the real.
* 
*               On  exit:  the  address  of  the  character  string  containing
* 	      the  current  record,  before reading the real.
* 
*           3.  int * ip
*               On entry:  the current character position, before reading the
* 	      real.
* 
*               On exit:  the current character position, after reading the real.
* 
*           4.  char eop
*               On entry:  the parameter delimiter.
* 
*           5.  char eor
*               On entry:  the record delimiter.
* 
*           6.  int maxCol
*               On entry:  the maximum character position of the line.
* 
*           7.  double * value
*               On entry:  the address of a variable to contain the real value.
* 
*               On exit:  the address of a variable containing the real value
* 	      read from the file.
* 
* 6     Return Values, Error Indicators and Warnings
* 
*       If the real read operation was successful, the value 1 is returned;
*       otherwise, the value 0 is returned.
* 
* Functions referenced by GetNextReal() are:
*  free_garray1()
*  ParseIgesRecord()
*
* Functions that reference GetNextReal() are:
*  DecodeIgesGlobal()
*  ReadIges106()
*  ReadIges110()
*  ReadIges112()
*  ReadIges114()
*  ReadIges120()
*  ReadIges122()
*  ReadIges124()
*  ReadIges126()
*  ReadIges128()
*/

int GetNextReal(FILE *fp, char *record, int *ip, char eop, char eor,
		  int maxCol, double *value)
{
  char *token;

  if (token = ParseIgesRecord(fp, record, ip, eop, eor, PRAX_IGES_REAL,
			      maxCol)) {
    *value = atof(token);   /* decode the current token on the record */
    free_garray1(token);
    return (1);   /* return success */
  }
  else {
    *value = 0.0;
    return (0);   /* return failure */ 
  }
}

/********* GetNextString() *********
* 1     Purpose
* 
*       Read character string value from IGES file.
* 
* 2     Specification
* 
*       #include "iges.h"
*       int GetNextString(FILE *fp, char *record, int *ip, char eop, char eor,
*                         int  maxCol,  char **value);
* 
* 3     Description
* 
*       This function reads an character string value starting at the current
*       character position in the current record of the IGES file.
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
*           2.  char * record
*               On  entry:  the address of the character string containing the
* 	      current record,  before reading the string.
* 
*               On  exit:  the  address  of  the  character  string  containing
* 	      the  current  record,  before reading the string.
* 
*           3.  int * ip
* 	      On entry:  the current character position, before reading the
* 	      string.
* 
*               On exit:  the current character position, after reading the
* 	      string.
* 
*           4.  char eop
*               On entry:  the parameter delimiter.
* 
*           5.  char eor
*               On entry:  the record delimiter.
* 
*           6.  int maxCol
*               On entry:  the maximum character position of the line.
* 
*           7.  char ** value
*               On entry:  the address of a variable to contain the address
* 	      of the character string.
* 
*               On exit:  the address of a variable containing the address
* 	      of the character string read from the file.
* 
* 6     Return Values, Error Indicators and Warnings
* 
*       If the string read operation was successful, the value 1 is returned;
*       otherwise, the value 0 is returned.
* 
* Functions referenced by GetNextString() are:
*  ParseIgesRecord()
*
* Functions that reference GetNextString() are:
*  DecodeIgesGlobal()
*/

int GetNextString(FILE *fp, char *record, int *ip, char eop, char eor,
		    int maxCol, char **value)
{
  char *token;

  if (token = ParseIgesRecord(fp, record, ip, eop, eor, PRAX_IGES_STRING,
			      maxCol)) {
    *value = token;   /* return the current token on the record */
    return (1);   /* return success */
  }
  else {
    *value = NULL;
    return (0);   /* return failure */
  }
}

/********* ParseIgesRecord() *********
* 1     Purpose
* 
*       Return next token from IGES file record.
* 
* 2     Specification
* 
*       #include "iges.h"
*       char  *ParseIgesRecord(FILE *fp, char *record, int *ip, char eop,
*                              char  eor,  int  type,  int maxCol);
* 
* 3     Description
* 
*       This function returns the next token from the current character
*       position of the IGES file.
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
*           1.  FILE *fp
*               On entry:  the file pointer of the IGES file.
* 
*           2.  char * record
*               On entry:  the address of a character string containing the
* 	      current line of the IGES file.
* 
*           3.  int * ip
*               On entry:  the address of a variable containing the character
* 	      position in the current record, before reading the next token.
* 
*               On  exit:  the  address  of  a  variable  that  will  contain
* 	      the  character  position  in  the current record, after reading
* 	      the token.
* 
*           4.  char eop
*               On entry:  the parameter delimiter character.
* 
*           5.  char eor
*               On entry:  the record delimiter character.
* 
*           6.  int type
*               On entry: flag specifying the data type of the token to be read.
* 
*               type  Memonic            Data_Type
*               ----  -----------------  ----------------
*                 0   PRAX_IGES_INTEGER  Integer
*                 1   PRAX_IGES_REAL     Real
*                 2   PRAX_IGES_STRING   Character string
* 
*           7.  int maxCol
*               On entry:  the maximum line length of the IGES file record.
* 
* 6     Return Values, Error Indicators and Warnings
* 
*       The address of a character string containing the token.  This may be a
*       zero-length string if no further tokens existed in the file record.
* 
* 8     Further Comments
* 
*       The returned character string containing the token is newly allocated
*       from dynamic memory and must be explicitly deallocated using free().
* 
* Functions referenced by ParseIgesRecord() are:
*  gen_array1()
*  ReadNextIgesRecord()
*
* Functions that reference ParseIgesRecord() are:
*  GetNextInteger()
*  GetNextReal()
*  GetNextString()
*/

char *ParseIgesRecord(FILE *fp, char *record, int *ip, char eop, char eor,
		      int type, int maxCol)
{
  int h = 0, hh = 0, haveH, i = 0, n = 0, numeric = 0, string = 0;
  char *token = NULL;

  if (record[*ip] == eor)      /* 'eor' is the record delimeter */
    return (token);            /* 'ip' is the character position in the */
                               /* current record */
  while (1) {
    if (*ip >= maxCol) {       /* we have read entire record, get another */
      ReadNextIgesRecord(fp, record, 0);
      *ip = 0;
    }
    if (string && i < h) {     /* if current token is a string, apend the */
      token[i++] = record[(*ip)++];   /* current character */
      continue;
    }
    if (record[*ip] == eop) {  /* eop is the parameter delimeter */
      if (!token)
	token = (char *)gen_array1(1, 1);
      token[i] = '\0';         /* end the current token */
      numeric = 0;
      string = 0;
      (*ip)++;
      break;
    }
    if (record[*ip] == eor) {
      if (!token)
	token = (char *)gen_array1(1, 1);
      token[i] = '\0';         /* if end-of-record, end the current token */
      break;
    }
    if (record[*ip] == ' ') {  /* between tokens, white space is ignored */
      (*ip)++;
      continue;
    }
    if (numeric) {             /* current token is numeric */
      if (type == PRAX_IGES_REAL) {
	if (record[*ip] == 'd') record[*ip] = 'e'; /* convert exponential D */
	if (record[*ip] == 'D') record[*ip] = 'E'; /* to E for atof() */
      }
      token[i++] = record[(*ip)++];
      continue;
    }
    if (type == PRAX_IGES_STRING) {   /* current token is a string */
      string = 1;
      haveH = 0;                      /* catch Hollerith strings that */
      for (i = *ip+1; i<maxCol; i++)  /* are broken before the 'H' */
	if (record[i] == 'H') {       /* 'broken' meaning that the string */
	  haveH = 1;                  /* is split between two records */
	  break;
	}
      sscanf(&record[*ip], "%d%n", &h, &n);
      if (!haveH) {
	ReadNextIgesRecord(fp, record, 0);
	*ip = 0;
	if (record[0] == 'H')
	  n = 0;
	else {
	  sscanf(&record[*ip], "%d%n", &hh, &n);
	  h = h*pow(10, (int)(log10(hh)+1)) + hh;
	}
      }
      token = (char *)gen_array1(1, h+1);
      i = 0;
      *ip += (n+1);
    }
    else {
      numeric = 1;
      token = (char *)gen_array1(1, maxCol);
      i = 0;
      token[i++] = record[(*ip)++];
    }
  }
  return (token);
}

/********* ReadNextIgesRecord() *********
* 1     Purpose
* 
*       Read the next line from an IGES file.
* 
* 2     Specification
* 
*       #include "iges.h"
*       char ReadNextIgesRecord(FILE *fp, char *record, int reread);
* 
* 3     Description
* 
*       This function reads the next line from an IGES file.
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
*               On entry:  file pointer of the IGES file.
* 
*           2.  char * record
*               On entry:  the address of a character string that contains
* 	      the current line.
* 
*               On exit:  the address of a character string that contains
* 	      the next line.
* 
*           3.  int reread
*               On entry:  flag specifying whether or not to physically read
* 	      the next line (flag = 1) or to reread the current line
* 	      (flag = 0).
* 
* 6     Return Values, Error Indicators and Warnings
* 
*       The section flag of the line read from the IGES file is returned.
* 
*         value     File_Section
*         -----     -----------------
*          'S'      Start section
*          'G'      Global section
*          'D'      Directory section
*          'P'      Parameter section
*          'T'      Terminate section
* 
*       If an error occurred during the read operation or if end-of-file was
*       reached the value EOF is returned.
* 
* Functions that reference ReadNextIgesRecord() are:
*  CheckIgesFile()
*  GetIgesDirList()
*  OpenIgesFileCB()
*  ParseIgesRecord()
*  SkipToNext()
*  SkipToThis()
*/

char ReadNextIgesRecord(FILE *fp, char *record, int reread)
{
  if (!reread) {
    if (fscanf(fp, "%80c", record) == EOF)   /* read data record */
      return (EOF);
    if (fgetc(fp) == EOF)                    /* read newline */
      return (EOF);
  }
  return (record[72]);
}

/********* SkipToThis() *********
* 1     Purpose
* 
*       Set current position in IGES file.
* 
* 2     Specification
* 
*       #include "iges.h"
*       void SkipToNext(FILE *fp, char *record, int *ip, int iNumber,
*                       char cSection);
* 
* 3     Description
* 
*       This function sets the current position in the IGES file to either
*       the first line of a specified section (or the next line, if already
*       in that section), or a specified line number in the current section
*       (or the next section, if already past that line number in the current
*       section).
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
*               On entry:  the address of the file structure of the IGES file.
* 
*           2.  char * record
*               On entry:  the address of a character string containing the
* 	      current line, before repositioning.
* 
*               On  exit:  the address of a character string containing the
* 	      current line,  after repositioning.
* 
*           3.  int * ip
*               On exit:  the address of a variable containing the character
* 	      position in the current line set to 0.
* 
*           4.  int iNumber
*               On entry:  the desired line number.
* 
*           5.  int cSection
*               On entry:  the desired file section.
* 
*                 cSection     File_Section
*                 --------     -----------------
*                    'S'       Start section
*                    'G'       Global section
*                    'D'       Directory section
*                    'P'       Parameter section
*                    'T'       Terminate section
* 
* 8     Further Comments
* 
*       This function will only read forward from the current position in the
*       file.  If unsure of the relative direction of the desired position,
*       use the function SkipToThis().
* 
* Functions referenced by SkipToThis() are:
*  ReadNextIgesRecord()
*
* Functions that reference SkipToThis() are:
*  OpenIgesFileCB()
*  ReadIgesCos()
*  ReadIgesRuledSurf()
*  ReadIgesSurfRev()
*  ReadIgesTabCyl()
*  ReadIgesTrim()
*/

int SkipToThis(FILE *fp, char *record, int *ip, int iNumber,
		 char cSection)
{
  int n;
  char c;

  *ip = 0;
  sscanf(&record[73], "%7d", &n);
  if (record[72] != cSection || n > iNumber) {   /* wrong section or */
    rewind(fp);                                  /* past record */
    while (1) {
      if ((c = ReadNextIgesRecord(fp, record, 0)) == EOF)
	return (-1); /* end-of-file found before requested record */
      if (c == cSection) {
	sscanf(&record[73], "%7d", &n);
	break;                                   /* now at first record */
      }                                          /* of correct section */
    }
  }
  while (n < iNumber) {
    if ((c = ReadNextIgesRecord(fp, record, 0)) == EOF)
      return (-2);   /* end-of-file found before requested record */
    sscanf(&record[73], "%7d", &n);
  }
  return (0);        /* signal success */
}

/********* SkipToNext() *********
* 1     Purpose
* 
*       Set current position in IGES file.
* 
* 2     Specification
* 
*       #include "iges.h"
*       int SkipToThis(FILE *fp, char *record, int *ip, int iNumber,
*                      char cSection);
* 
* 3     Description
* 
*       This  function  sets  the  current  position  in  the  IGES  file  to
*       a  specified  line  in  a  specified section.
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
*               On entry:  the address of the file structure of the IGES file.
* 
* 
*           2.  char * record
*               On entry:  the address of a character string containing the
* 	      current line, before repositioning.
* 
*               On  exit:  the address of a character string containing the
* 	      current line,  after repositioning.
* 
*           3.  int * ip
*               On exit:  the address of a variable containing the character
* 	      position in the current line set to 0.
* 
*           4.  int iNumber
*               On entry:  the desired line number.
* 
*           5.  int cSection
*               On entry:  the desired file section.
* 
*                 cSection     File_Section
*                 --------     -----------------
*                    'S'       Start section
*                    'G'       Global section
*                    'D'       Directory section
*                    'P'       Parameter section
*                    'T'       Terminate section
* 
* Functions referenced by SkipToNext() are:
*  ReadNextIgesRecord()
*
* Functions that reference SkipToNext() are:
*  OpenIgesFileCB()
*/

void SkipToNext(FILE *fp, char *record, int *ip, int iNumber,
		char cSection)
{
  int n;

  while (1) {
    ReadNextIgesRecord(fp, record, 0);   /* read forward until we find */
    if (record[72] == cSection)          /* the requested section */
      break;
    sscanf(&record[73], "%7d", &n);
    if (n == iNumber)                    /* read until we find the */
      break;                             /* requested line number */
  }
  *ip = 0;
}

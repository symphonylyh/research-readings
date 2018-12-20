/* Copyright (C) 1992 Massachusetts Institute of Technology, Cambridge, MA
   All rights reserved
*/

/* utils.c */

/* AllocIgesDirectory()
 * AllocIgesGlobal()
 * FileCat()
 * FreeIgesDirectory()
 * FreeIgesDirList()
 * FreeIgesGlobal()
 * GetIgesDirEntry()
 * GetIgesDirList()
 * GetIgesEntity()
*/

#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <malloc.h>
#include "iges.h"
#include "editor.h"

#define nTypes 16

char *strdup(const char *);
/*int gethostname(char*, int); */

/********* AllocIgesDirectory() *********
* 1     Purpose
* 
*       Allocate IGES file Directory structure.
* 
* 2     Specification
* 
*       #include "iges.h"
*       DirStruct *AllocIgesDirectory(int type, int parameter, int structure,
*                                     int pattern, int level, int view,
* 				    int matrix, int associativity, int blank,
* 				    int subordinate, int use, int hierarchy,
* 				    int weight, int color, int count,
* 				    int form, char *label, int subscript,
* 				    int DE);
* 3     Description
* 
*       This function allocates a data structure that contains the information
*       in an entity Directory section entry.
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
*       Parameters that are undefined in the following list are unsupported by
*       the libiges.a library and should be specified as 0.
* 
*           1.  int type
*               On entry:  the entity type number.
* 
*           2.  int parameter
*               On entry:  the sequence number of the first line of parameter
* 	      data for this entity.
* 
*           3.  int structure
*               Note:  unsupported; specify as 0.
* 
*           4.  int pattern
*               Note:  unsupported; specify as 0.
* 
*           5.  int level
*               Note:  unsupported; specify as 0.
* 
*           6.  int view
*               Note:  unsupported; specify as 0.
* 
*           7.  int matrix
*               Note:  unsupported; specify as 0.
* 
*           8.  int associativity
*               Note:  unsupported; specify as 0.
* 
*           9.  int blank
*               Note:  unsupported; specify as 0.
* 
*          10.  int subordinate
*               Note:  unsupported; specify as 0.
* 
*          11.  int use
*               Note:  unsupported; specify as 0.
* 
*          12.  int hierarchy
*               Note:  unsupported; specify as 0.
* 
*          13.  int weight
*               Note:  unsupported; specify as 0.
* 
*          14.  int color
*               Note:  unsupported; specify as 0.
* 
*          15.  int count
*               On entry:  the count of the number of lines in the parameter
* 	      section for this entity.
* 
*          16.  int form
*               On entry:  format number to select entity-specific formatting
* 	      options.
* 
*          17.  char * label
*               On entry:  an eight character descriptive label for the entity.
* 
*          18.  int subscript
*               Note:  unsupported; specify as 0.
* 
*          19.  int DE
*               On entry:  the sequence number of the entry.
* 
* 6     Return Values, Error Indicators and Warnings
* 
*       If insufficient dynamic memory is available, the function returns the
*       value NULL.
* 
* 8     Further Comments
* 
*       Free the data structure with FreeIgesDirectory().
* 
* Functions referenced by AllocIgesDirectory() are:
*  gen_array1()
*
* Functions that reference AllocIgesDirectory() are:
*  SaveIgesCos()
*  SaveIgesCurv()
*  SaveIgesGrid()
*  SaveIgesList()
*  SaveIgesMinDist()
*  SaveIgesSurf()
*  SaveIgesTrim()
*  SaveIgesUv()
*/

DirStruct *AllocIgesDirectory(int type, int parameter, int structure,
			      int pattern, int level, int view, int matrix,
			      int associativity, int blank,
			      int subordinate, int use, int hierarchy,
			      int weight, int color, int count, int form,
			      char *label, int subscript, int DE)
{
  DirStruct *dir;

  if (dir = (DirStruct *)gen_array1(1, sizeof(DirStruct))) {
    dir->type = type;
    dir->parameter = parameter;
    dir->structure = structure;
    dir->pattern = pattern;
    dir->level = level;
    dir->view = view;
    dir->matrix = matrix;
    dir->associativity = associativity;
    dir->blank = blank;
    dir->subordinate = subordinate;
    dir->use = use;
    dir->hierarchy = hierarchy;
    dir->weight = weight;
    dir->color = color;
    dir->count = count;
    dir->form = form;
    strncpy(dir->label, label, 9);
    dir->subscript = subscript;
    dir->DE = DE;
  }
  return (dir);
}

/********** AllocateIgesGlobal() *********
* 1     Purpose
* 
*       Allocate IGES file Global data structure.
* 
* 2     Specification
* 
*       #include "iges.h"
*       GlobalStruct *AllocIgesGlobal(char parmDelim, char recDelim,
*                                     char *sendID, char *file Name,
* 				    char *systemID, char *preprocVers,
* 				    int bits, int singleMagnitude,
* 				    int single Precision, int doubleMagnitude,
* 				    int doublePrecision,  char *rceiveID,
* 				    double scale, int  unit,  char  *units,
* 				    int  gradations,  double  weight,
* 				    char  *fileDate,  double  resolution,
* 				    double  coordinate,  char  *author,
* 				    char  *org,  int  version,  int  standard,
* 				    char*modelDate);
* 3     Description
* 
*       This function allocates a data structure that contains the information
*       in a Global section.
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
*           1.  char parmDelim
*               On entry:  the parameter delimiter.
* 
*           2.  char recDelim
*               On entry:  the record delimiter.
* 
*           3.  char * sendID
*               On entry: the address of a character string containing the name
* 	      of the model specified by the file's entities.
* 
*           4.  char * fileName
*               On entry:  the address of a character string containing the file
* 	      name.
* 
*           5.  char * systemID
*               On  entry:  the address of a character string containing the
* 	      system that created the file.
* 
*           6.  char * preprocVers
*               On  entry:  the  address  of  a  character  string  containing
* 	      the  version  identifier  of  the system that created the file.
* 
*           7.  int bits
*               On entry:  the number of bits for integer representation on the
* 	      system that created the file.
* 
*           8.  int singleMagnitude
*               On entry:  the maximum power of ten representable by a single
* 	      precision number on the system that created the file.
* 
*           9.  int singlePrecision
*               On entry:  the number of significant digits in a single
* 	      precision number on the system that created the file.
* 
*          10.  int doubleMagnitude
*               On entry:  the maximum power of ten representable by a double
* 	      precision number on the system that created the file.
* 
*          11.  int doublePrecision
*               On entry:  the number of significant digits in a double
* 	      precision number on the system that created the file.
* 
*          12.  char * receiveID
*               On  entry:  the  address  of  a  character  string  containing
* 	      the  name  of  the  receiving system.
* 
*          13.  double scale
*               On entry:  the model space scale given as the ration of model
* 	      space units to real world units.
* 
*          14.  int unit
*               On entry:  flag specifying the unit system.
* 
*                 flag   Measuring system
*                 ----   -------------------------------------
*                   1    inches
*                   2    millimeters
*                   3    units named by units character string
*                   4    feet
*                   5    miles
*                   6    meters
*                   7    kilometers
*                   8    mils
*                   9    microns
*                  10    centimeters
*                  11    microinches
* 
*          15.  char * units
*               On entry:  the address of a character string containing the
* 	      name of the units.
* 
*                 units       Units
*                 -----       -----------
*                 "IN"        inches
*                 "INCH"      inches
*                 "MM"        millimeters
*                 "FT"        feet
*                 "MI"        miles
*                 "M"         meters
*                 "KM"        kilometers
*                 "MIL"       mils
*                 "UM"        microns
*                 "CM"        centimeters
*                 "UIN"       microinches
* 
*         16.   int gradations
*               On entry:  the number of distinct line thicknesses.
* 
*         17.   double weight
*               On entry:  the actual width of the thickness possible line.
* 
*         18.   char * fileDate
*               On entry:  the address of a character string containing the
* 	      file creation date.
* 
*         19.   double resolution
*               On entry:  the smallest discernible distance in model space.
* 
*         20.   double coordinate
*               On entry:  the absolute value of the maximum coordinate value.
* 
*         21.   char * author
*               On entry:  the address of a character string containing the
* 	      name of the author of the model.
* 
*         22.   char * org
*               On  entry:   the  address  of  a  character  string  containing
* 	      the  name  of  the  author's organization.
* 
*         23.   int version
*               On entry:  flag specifying the IGES version.
* 
*                 version     Version
*                 -------     ------------------------
*                     1       1.0
*                     2       ANSI Y14.26M - 1981
*                     3       2.0
*                     4       3.0
*                     5       ASME/ANSI Y14.26M - 1987
*                     6       4.0
*                     7       ASME Y14.26M - 1989
*                     8       5.0
* 
*         24.   int standard
*               On entry:  flag specifying the drafting standard.
* 
*               standard  Standard
*               --------  ------------------------------------------------------
*                    0     None
*                    1     ISO    International Organization for Standardization
*                    2     AFNOR  French Association for Standardization
*                    3     ANSI   American National Standards Institute
*                    4     BSI    British Standards Institute
*                    5     CSA    Canadian Standards Association
*                    6     DIN    German Institute for Standardization
*                    7     JIS    Japanese Institute for Standardization
* 
*          25.  char * modelDate
*               On entry:  the address of a character string containing the
* 	      model creation date.
* 
* 6     Return Values, Error Indicators and Warnings
* 
*       If insufficient dynamic memory is available, the function returns the
*       value NULL.
* 
* 8     Further Comments
* 
*       Free the data structure with FreeIgesDirectory().
* 
* Functions referenced by AllocIgesGlobal() are:
*  gen_array1()
*
* Functions that reference AllocIgesGlobal() are:
*  SaveIgesCos()
*  SaveIgesCurv()
*  SaveIgesGrid()
*  SaveIgesList()
*  SaveIgesMinDist()
*  SaveIgesSurf()
*  SaveIgesTrim()
*  SaveIgesUv()
*/

GlobalStruct *AllocIgesGlobal(char parmDelim, char recDelim, char *sendID,
			      char *fileName, char *systemID,
			      char *preprocVers, int bits,
			      int singleMagnitude, int singlePrecision,
			      int doubleMagnitude, int doublePrecision,
			      char *receiveID, double scale, int unit,
			      char *units, int gradations, double weight,
			      char *fileDate, double resolution,
			      double coordinate, char *author, char *org,
			      int version, int standard, char *modelDate)
{
  GlobalStruct *global;
  char *cwd, line[256];

  if (global = (GlobalStruct *)gen_array1(1, sizeof(GlobalStruct))) {
    global->parmDelim   = parmDelim;
    global->recDelim    = recDelim;
    global->sendID      = strdup(sendID);
    gethostname(line, 64);
    strcat(line, ":");
    if (fileName[0] != '/') {
      cwd = getcwd(NULL, 128);
      strcat(line, cwd);
      strcat(line, "/");
      free(cwd);
    }
    strcat(line, fileName);
    global->fileName    = strdup(line);
    global->systemID    = strdup(systemID);
    global->preprocVers = strdup(preprocVers);
    global->bits        = bits;
    global->singleMagnitude = singleMagnitude;
    global->singlePrecision = singlePrecision;
    global->doubleMagnitude = doubleMagnitude;
    global->doublePrecision = doublePrecision;
    global->receiveID   = strdup(receiveID);
    global->scale       = scale;
    global->unit        = unit;
    global->units       = strdup(units);
    global->gradations  = gradations;
    global->weight      = weight;
    strcpy(global->fileDate, fileDate);
    global->resolution  = resolution;
    global->coordinate  = coordinate;
    global->author      = strdup(author);
    global->org         = strdup(org);
    global->version     = version;
    global->standard    = standard;
    strcpy(global->modelDate, modelDate);
  }
  return (global);
}

/********* FileCat() *********
* 1     Purpose
* 
*       Concatenate two files.
* 
* 2     Specification
* 
*       #include "iges.h"
*       void FileCat(FILE *fp1, FILE *fp2);
* 
* 3     Description
* 
*       This function concatenates the file pointed to by fp2 onto the end of
*       the file pointed to by fp1.  File fp2 must be open for reading and
*       positioned at the beginning of the file.  File fp1 must be open for
*       writing and positioned at the end of the file.
* 
* 5     Parameters
* 
*           1.  FILE * fp1
*               On entry:  file pointer of the file to which will be
* 	      concatenated file fp2.
* 
*           2.  FILE * fp2
*               On entry:  file pointer of the file to be concatenated onto
* 	      the end of file fp1.
* 
* Functions that reference FileCat() are:
*  SaveIgesCos()
*  SaveIgesCurv()
*  SaveIgesGrid()
*  SaveIgesList()
*  SaveIgesMinDist()
*  SaveIgesSurf()
*  SaveIgesTrim()
*  SaveIgesUv()
*/

void FileCat(FILE *fp1, FILE *fp2)
{
  int c;

  while ((c = fgetc(fp2)) != EOF)
    fputc(c, fp1);
}

/********* FreeIgesDirectory() *********
 * Functions referenced by FreeIgesDirectory() are:
 *  free_garray1()
 *
 * Functions that reference FreeIgesDirectory() are:
 *  SaveIgesCos()
 *  SaveIgesCurv()
 *  SaveIgesGrid()
 *  SaveIgesList()
 *  SaveIgesMinDist()
 *  SaveIgesSurf()
 *  SaveIgesTrim()
 *  SaveIgesUv()
 */

void FreeIgesDirectory(DirStruct *dir)
{
  if (dir)
    free_garray1((char *)dir);
}

/********* FreeIgesDirList() *********
* 1     Purpose
* 
*       Deallocate IGES file Directory data structure.
* 
* 2     Specification
* 
*       #include "iges.h"
*       void FreeIgesDirectory(DirStruct *dir);
* 
* 3     Description
* 
*       This function deallocates an IGES file Directory structure.
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
*               On entry:  the address of the structure.
* 
* 8     Further Comments
* 
*       The structure should be allocated with AllocIgesDirectory
* 
* Functions referenced by FreeIgesDirList() are:
*  free_copious()
*  free_egeom()
*  free_fgeom()
*
* Functions that reference FreeIgesDirList() are:
*  OpenIgesFileCB()
*/

void FreeIgesDirList(EntryStruct *list)
{
  EntryStruct *prev = NULL, *this = NULL;

  this = list;
  while (this) {
    if (this->tmpfile) {        /* temp file name */
      unlink(this->tmpfile);
      free(this->tmpfile);
    }
    switch (this->dir->type) {  /* geometry structure */
    case 100:

      if (this->igeom)
	free_garray1(this->igeom);
      break;

    case 102:

      if (this->igeom) {
	free_iarray1(((Type102 *)(this->igeom))->de);
	free_garray1(this->igeom);
      }
      break;

    case 106:

      if (this->igeom)
	free_copious(this->igeom);
      break;

    case 108:

      if (this->igeom)
	free_garray1(this->igeom);
      break;

    case 110:

      if (this->igeom)
	free_garray1(this->igeom);
      break;

    case 112:

      if (this->igeom)
	free_egeom(this->igeom);
      break;

    case 114:

      if (this->igeom)
	free_egeom(this->igeom);
      break;

    case 116:

      if (this->igeom)
	free_garray1((char *)this->igeom);
      break;

    case 118:

      if (this->igeom)
	free_fgeom(this->igeom);
      break;

    case 120:

      if (this->igeom)
	free_fgeom(this->igeom);
      break;

    case 122:

      if (this->igeom)
	free_fgeom(this->igeom);
      break;

    case 126:

      if (this->igeom)
	free_egeom(this->igeom);
      break;

    case 128:

      if (this->igeom)
	free_fgeom(this->igeom);
      break;
    }
    if (this->dir)
      free(this->dir);          /* directory structure */
    prev = this;
    this = this->next;
    if (prev)
      free(prev);               /* directory entry */
  }
}

/********* FreeIgesGlobal() *********
* 1     Purpose
* 
*       Deallocate IGES file Global data structure.
* 
* 2     Specification
* 
*       #include "iges.h"
*       void FreeIgesGlobal(GlobalStruct *global);
* 
* 3     Description
* 
*       This function deallocates an IGES file Global structure.
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
*           1.  GlobalStruct * global
*               On entry:  the address of the structure.
* 
* 8     Further Comments
* 
*       The structure should be allocated with AllocIgesGlobal
* 
* Functions referenced by FreeIgesGlobal() are:
*  free_garray1()
*
* Functions that reference FreeIgesGlobal() are:
*  OpenIgesFileCB()
*  SaveIgesCos()
*  SaveIgesCurv()
*  SaveIgesGrid()
*  SaveIgesList()
*  SaveIgesMinDist()
*  SaveIgesSurf()
*  SaveIgesTrim()
*  SaveIgesUv()
*/

void FreeIgesGlobal(GlobalStruct *global)
{
  if (global->sendID)
    free (global->sendID);
  if (global->fileName)
    free (global->fileName);
  if (global->systemID)
    free (global->systemID);
  if (global->preprocVers)
    free (global->preprocVers);
  if (global->receiveID)
    free (global->receiveID);
  if (global->units)
    free (global->units);
  if (global->author)
    free (global->author);
  if (global->org)
    free (global->org);
  free_garray1((char *)global);
}

/********* GetIgesDirEntry() *********
* 1     Purpose
* 
*       Find entry in the IGES file Directory list.
* 
* 2     Specification
* 
*       #include "iges.h"
*       EntryStruct *GetIgesDirEntry(EntryStruct *list, int *iEntry);
* 
* 3     Description
* 
*       This  function  finds  an  entry  in  the  linked  list  of  Directory
*       entries  of  an  IGES  file.  The sequence numbers of the list entries,
*       stored in list->dir->DE, are 1,3,5,...
* 
* 5     Parameters
* 
*           1.  EntryStruct * list
*               On entry:  the address of the first entry in the Directory list.
* 
*           2.  int iEntry
*               On  entry:  the sequence number of the desired list entry.
* 	      This  number  will  be  an positive odd integer 2i - 1,  where
* 	      i is the positional index of the entry in the list.  Thus, for
* 	      the first entry, iEntry = 1; for the second, iEntry = 3; and
* 	      so on.
* 
* 6     Return Values, Error Indicators and Warnings
* 
*       The address of the entry is returned.  If the entry is not found the
*       value NULL is returned.
*/

EntryStruct *GetIgesDirEntry(EntryStruct *list, int iEntry)
{
  EntryStruct *this;
  int DE = -1;

  this = list;
  while (this) {
    DE += 2;
    if (DE == iEntry)
      return (this);
    this = this->next;
  }
  return (NULL);
}

/********* GetIgesDirList() *********
* 1     Purpose
* 
*       Read an IGES file Directory entry list.
* 
* 2     Specification
* 
*       #include "iges.h"
*       EntryStruct *GetIgesDirList(FILE *fp, char *record);
* 
* 3     Description
* 
*       This  function  reads  the  Directory  section  of  an  IGES  file  and
*       builds  a  linked  list  of  its entries.  The list is composed of a
*       sequence of EntryStruct structures, each of which points to a DirStruct,
*       which contains the entry for a single entity.
* 
* 5     Parameters
* 
*           1.  FILE * fp
*               On entry:  the file point of the IGES file.
* 
*           2.  char * record
*               On entry:  the address of a character string containing the
* 	      first line of the Directory section.
* 
* 6     Return Values, Error Indicators and Warnings
* 
*       The address of the linked list is returned.  If insufficient dynamic
*       memory is available, the value NULL is returned.
* 
* 8     Further Comments
* 
*       The entry list should be deallocated by FreeIgesDirList().
* 
* Functions referenced by GetIgesDirList() are:
*  DecodeIgesDirectory()
*  gen_array1()
*  ReadNextIgesRecord()
*
* Functions that reference GetIgesDirList() are:
*  OpenIgesFileCB()
*/

EntryStruct *GetIgesDirList(FILE *fp, char *record)
{
  EntryStruct *list = NULL, *prev = NULL, *this = NULL;
  DirStruct *dir;
  int DE = -1, first = 1, ip = 0;

  while (ReadNextIgesRecord(fp, record, 0) == PRAX_IGES_DIRECTORY) {
    DecodeIgesDirectory(record, &first, &dir);
    if (first) {
      dir->DE = (DE += 2);
      this = (EntryStruct *)gen_array1(1, sizeof(EntryStruct));
      if (!list)
	list = this;
      if (prev)
	prev->next = this;
      this->dir = dir;
      prev = this;
    }
  }
  return (list);
}

/********* GetIgesEntity() *********
* 1     Purpose
* 
*       Return name of IGES entity type.
* 
* 2     Specification
* 
*       #include "iges.h"
*       char *GetIgesEntity(int type);
* 
* 3     Description
* 
*       This function returns the address of a character string that contains
*       the name of an IGES entity type.
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
*           1.  int type
*               On entry:  the entity type number
* 
* 6     Return Values, Error Indicators and Warnings
* 
*       The address of a character string containing one of the following
*       labels is returned:
* 
*         type    Label
*         ----    -----------------
*         100     "Circular Arc"
*         102     "Composite Curve"
*         106     "Copious Data"
*         108     "Plane"
*         110     "Line"
*         112     "Param Curve"
*         114     "Param Surface"
*         116     "Point"
*         118     "Ruled Surface"
*         120     "Surf of Rev"
*         122     "Tabulated Cyl"
*         126     "NURBS Curve"
*         128     "NURBS Surface"
*         142     "Curve on Surface"
*         144     "Trimmed Surface"
*         else    "Unsupported!"
* 
* Functions that reference GetIgesEntity() are:
*  OpenIgesFileCB()
*/

char *GetIgesEntity(int type)
{
  static char *types[nTypes] = {
    "Circular Arc",
    "Composite Curve",
    "Copious Data",
    "Plane",
    "Line",
    "Param Curve",
    "Param Surface",
    "Point",
    "Ruled Surface",
    "Surf of Rev",
    "Tabulated Cyl",
    "NURBS Curve",
    "NURBS Surface",
    "Curve on Surface",
    "Trimmed Surface",
    "Unsupported!"
  };

  switch (type) {
  case 100:   /* circualr arc */

    return types[0];
    break;

  case 102:   /* composite curve */

    return types[1];
    break;

  case 106:   /* copious data */

    return (types[2]);
    break;

  case 108:   /* plane */

    return types[3];
    break;

  case 110:   /* copious data */

    return (types[4]);
    break;

  case 112:   /* parametric curve */

    return (types[5]);
    break;

  case 114:   /* parametric surface */

    return (types[6]);
    break;

  case 116:   /* point */

    return types[7];
    break;

  case 118:   /* ruled surface */

    return (types[8]);
    break;

  case 120:   /* surface of revolution */

    return (types[9]);
    break;

  case 122:   /* tabulated cylinder */

    return (types[10]);
    break;

  case 126:   /* NURBS curve */

    return (types[11]);
    break;

  case 128:   /* NURBS surface */

    return (types[12]);
    break;

  case 142:   /* curve on surface */

    return (types[13]);
    break;

  case 144:   /* trimmed surface */

    return (types[14]);
    break;
  }

  return (types[15]);   /* unsupported */
}

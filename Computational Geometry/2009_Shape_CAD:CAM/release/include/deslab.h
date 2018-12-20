/* Copyright (C) 1995 Massachusetts Institute of Technology, Cambridge, MA */
/* All rights reserved */

/* deslab.h */

#ifndef DESLAB_H
#define DESLAB_H

#include <stdio.h>
#include <string.h>

#ifdef __cplusplus
extern "C" {
#endif

/* undefined type */

#define DESLAB_NONE      -1

/* application types */

#define DESLAB_DESLAB     0
#define DESLAB_PRAXITELES 1
#define DESLAB_CLINE      2
#define DESLAB_CSURF      3
#define DESLAB_GEOD       4
#define DESLAB_ISO        5
#define DESLAB_LOC        6
#define DESLAB_MDIST      7
#define DESLAB_MINDIST    8
#define DESLAB_REFL       9

/* application names */

#define DESLAB_DESLAB_APP "deslab"
#define DESLAB_PRAX_APP   "praxiteles"
#define DESLAB_CLINE_APP  "cline"
#define DESLAB_CSURF_APP  "csurf"
#define DESLAB_GEOD_APP   "geod"
#define DESLAB_ISO_APP    "iso"
#define DESLAB_LOC_APP    "loc"
#define DESLAB_MDIST_APP  "mdist"
#define DESLAB_MINDIST_APP "mindist"
#define DESLAB_REFL_APP   "refl"

/* file types */

#define DESLAB_LIST        0
#define DESLAB_GRID        1
#define DESLAB_CURV        2
#define DESLAB_SURF        3
#define DESLAB_COS         4
#define DESLAB_IGES        5
#define DESLAB_MULT        6
#define DESLAB_LOGO        7
#define DESLAB_DIST        8
#define DESLAB_VECT        9
#define DESLAB_FUNC       10
#define DESLAB_UV         11
#define DESLAB_LOCAL      12
#define DESLAB_FACET      13
#define DESLAB_TRIM       14
#define DESLAB_GROUP      15
#define DESLAB_HULL       16
#define DESLAB_DIST_INSP  17
#define DESLAB_LOCAL_INSP 18
#define DESLAB_SCRIPT     19
#define DESLAB_LOG        20
#define DESLAB_PS         21
#define DESLAB_INSP       22

/* curve periodicity */

#define DESLAB_NON_PERIODIC 0
#define DESLAB_PERIODIC     1

/* message types */

#define DESLAB_INFO     0
#define DESLAB_WARNING	1
#define DESLAB_ERROR	2

#define DESLAB_FILE       0
#define DESLAB_EVALUATE   1
#define DESLAB_PIVOT      2
#define DESLAB_CONTINUITY 3

/* function prototypes */

/* add prototypes not in strict ANSI library headers */

char *cuserid(char *);
  /* int   gethostname(char *, int); */ 
int   gethostname(char *, size_t);   
int   strcasecmp(const char *, const char *);
int   strncasecmp(const char *, const char *, size_t);

/* deslab functions */

int   CheckDeslabFile(FILE *, int);
char *GetEntityName(int);
char *GetFormatName(int);
void  WriteDeslabHeader(FILE *, int, int, int, char *, char *, int);

#ifdef __cplusplus
}
#endif

#endif

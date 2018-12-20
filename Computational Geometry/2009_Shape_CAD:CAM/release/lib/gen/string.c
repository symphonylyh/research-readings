/************************************************************************
 *									*
			    Copyright (C) 1990
     Massachusetts Institute of Technology, Cambridge, Massachusetts.
			   All rights reserved.

     This portion of source code was prepared in the Ocean Engineering
     Design Laboratory by Bradley A. Moran.
 *									*
 ************************************************************************/

#include <ctype.h>
#include <stdio.h>
#include "gen.h"

#define LINESIZE 81

/* Function: stripwhite()
 * Purpose: Strip white space from a file
 * Method: Increment the file pointer forward by reading and
 *         discarding whitepspace one character at a time
 * Arguments:
 *  iop - file point of the open file
 */

void stripwhite (FILE *iop)	  /* strip white space from iop */
{
     register int c;

     /* read one character at a time until End-Of-File */

     while ((c = getc(iop)) != EOF)
	  if (!isspace(c)) {      /* if it's white space, ignore it */
	       ungetc (c, iop);   /* if not, move the pointer back one */
	       break;             /* and quit */
	  }
}

/* Function: fgetstring()
 * Purpose: Read string from file
 * Arguments:
 *  string -
 *  n - size of string
 *  ios - file pointer of open file
 */
/* Functions that reference fgetstring() are:
 *  fileopen()
 *  ReadParSurf()
 */

char *fgetstring (char *string, int n, FILE *ios)

/* Original declaration follows: */
/* char *string;		   target string */
/* int n;			   size of input string */
/* FILE *ios;			   input stream */

{
     char line[LINESIZE];	/* # defined to be 81, can be larger */
     char format[6];		/* for telling sscanf how big to go */
     char *rc;

     stripwhite (ios);          /* ignore any white space */

     if ((rc = fgets(line, sizeof(line), ios)) != (char *) NULL) {
	  sprintf (format, "%%%ds", n);
	  sscanf (line, format, string);
     }
     return (rc);
}

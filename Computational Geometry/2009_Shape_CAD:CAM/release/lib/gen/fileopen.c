/************************************************************************
 *									*
		Copyright (C) 1989 by Bradley A. Moran
			 All rights reserved.

    The  author hereby grants to  MIT  permission to reproduce and  to
    distribute  copies of  this source  code  document in whole or  in
    part.
 *									*
 ************************************************************************/

#include <stdio.h>
#include "gen.h"

/* Function: fileopen()
 * Purpose: Open a file
 * Method: Open the named file.  If an error occurs, repeatedly
 *         prompt for a new file name to entered via the standard
 *         input unit.  If an exclamtion mark "!" is entered, then
 *         no file is opened.
 * Arguments:
 *  filename - the name of the file
 *  n - the length of the character string filename
 *  type - the file open type (for a complete list do "man fopen")
 * Return: The file pointer of the opened file
 */
/* Functions referenced by fileopen() are:
 *  fgetstring()
 */
/* Functions that reference fileopen() are:
 *  ReadParCurv()
 *  ReadParSurf()
 */

FILE *fileopen (char *filename, int n, char *type)

/* Open a file. filename may be relative or absolute, n is the size
   of memory allocated for filename, type is standard i/o type string
   (e.g. "r", "rb", "wb", etc.) */

{
     FILE *fp;

     while ((fp = fopen (filename, type)) == (FILE *) NULL) {
				/* open failure */
	  perror (filename);
	  fputs ("Enter another file (! to abort) : ", stderr);
	  fflush (stderr);
	  fgetstring (filename, n, stdin);
	  if (*filename == (unsigned)'!')
	       break;
     }
     return (fp);
}

# ifdef DEBUG

main (int argc, char **argv)

{
     FILE *fp, *fileopen ();
     char str[26];

     fputs ("File name : ", stdout);
     fgetstring (str, sizeof(str), stdin);
     fp = fileopen (str, sizeof(str), "r");
}

# endif

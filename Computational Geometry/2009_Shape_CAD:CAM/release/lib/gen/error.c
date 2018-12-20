/************************************************************************
 *									*
			Copyright (C) 1989 by
	Massachusetts Institute of Technology, Cambridge, MA
			 All rights reserved
 *									*
 ************************************************************************/

#include <stdio.h>
#include <stdlib.h>

/* Function: errormsg()
 * Purpose: Display error message
 * Method: Print the error message on the standard output unit and
 *         then use the kernal function exit() to terminate the
 *         execution with an exit code.
 * Arguments:
 *  code - the exit code returned by the terminating program.
 *         If the code is 1, then only print a warning message and
 *         continue, do not terminate execution.
 *  msg - textual error message
 */
/* Functions that reference errormsg() are:
 *  Aij()
 *  Aijqs()
 *  Aijqs_per()
 *  Aijqs_z()
 *  build_offset()
 *  curvature()
 *  curvature1()
 *  dbl_array1()
 *  dbl_array2()
 *  dbl_array3()
 *  egeomalloc1()
 *  egeomalloc2()
 *  evalbsp()
 *  evalbsp_per()
 *  fgeomalloc1()
 *  fgeomalloc2()
 *  flt_array1()
 *  flt_array2()
 *  flt_array3()
 *  int_array1()
 *  int_array2()
 *  int_sample_gencyl()
 *  lf2comp()
 *  loft_rational()
 *  nbasisd()
 *  partan_dir()
 *  pgeomalloc()
 *  pretransmat()
 *  ptr_array1()
 *  ptr_array2()
 *  ReadParCurv_Per()
 *  ReadParSurf_Peru()
 *  sample_offset()
 *  sgeomalloc()
 *  sht_array1()
 *  sht_array2()
 *  transmat()
 *  vectalloc()
 *  vec_array1()
 *  vec_array2()
 *  WriteParCurv()
 *  WriteParCurv_Per()
 *  WriteParSurf()
 *  WriteParSurf_Peru()
 */

int errormsg (int code, char *msg)	/* print to stderr */
{
     (void) fprintf (stderr, "\n%s\n", msg);

     if (code != 1) {		/* decide whether to exit */
	  (void) fprintf (stderr, "Error code is %d.  Exiting...\n\r", code);
	  exit (code);
     }
     return 0;
}

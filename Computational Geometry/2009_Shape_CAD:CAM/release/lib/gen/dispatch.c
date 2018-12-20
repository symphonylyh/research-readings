/************************************************************************
 *									*
      Copyright (C) 1989 Massachusetts Institute of Technology,
			    Cambridge, MA.
			 All rights reserved.

	  Department of Ocean Engineering Design Laboratory
 *									*
 ************************************************************************/

#include <errno.h>
#include <string.h>
#include <gen.h>

/* Function: dispatch()
 * Purpose: Set the geometry type
 * Method: Use the specified input string to identity a particular
 *         geometry type.
 * Arguments:
 *  string - the name of the type
 * Return: The numerical flag of the type
 */

int dispatch (char *string)

{
     extern int errno;
     int type;
     
     if (!strcmp(string, "ACurve"))              /* Curve */
	  type = ACurve;
     else if (!strcmp(string, "ASurface"))       /* Surface */
	  type = ASurface;
     else if (!strcmp(string, "PCurveOpen"))     /* Non-periodic curve */
	  type = PCurveOpen;
     else if (!strcmp(string, "PCurvePer"))      /* Periodic curve */
	  type = PCurvePer;
     else if (!strcmp(string, "PSurfaceOpen"))   /* Non-periodic surface */
	  type = PSurfaceOpen;
     else if (!strcmp(string, "PSurfacePerU"))   /* surface periodic i u */
	  type = PSurfacePerU;
     else if (!strcmp(string, "PSurfacePerV"))   /* surface periodic in v */
	  type = PSurfacePerV;
     else if (!strcmp(string, "PSurfacePerUV"))  /* Surface Periodic in u,v */
	  type = PSurfacePerUV;
     else {
	  errno = ENODEV;
	  type = ENOGEOM;
     }
     return(type);
}

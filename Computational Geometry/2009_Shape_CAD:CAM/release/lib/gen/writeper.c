/************************************************************************
 *									*
			Copyright (C) 1989 by
	Massachusetts Institute of Technology, Cambridge, MA
			 All rights reserved
 *									*
 ************************************************************************/
/************************************************************************
 *									*
			 Modification History

     12 Apr 89 - B.A. Moran, fixed bug in number knots
     12 Dec 88 - S.T. Tuohy, original author
 *									*
 ************************************************************************/

#include <stdio.h>
#include "gen.h"

/* Function: WriteParCurv_Per()
 * Purpose: Write a periodic NURBS curve to a file
 * Method: The NURBS curve structure is defined in the header file
 *         gen.h
 * Arguments:
 *  fp - file pointer of the open file to contain the curve
 *  egm - address of the curve structure
 */
/* Functions referenced by WriteParCurv_Per() are:
 *  errormsg()
 */

void WriteParCurv_Per(FILE *fp, ParCurv *egm)
{
  int fpg = FALSE,i;
  char fname[60],msg[100];
 
  if (fp == NULL) {   /* if open file not specified, prompt for name */
    fpg = TRUE;
    printf("\n\rWriteParCurv_Per() output file name:  "), scanf("%s",fname);
    if ((fp = fopen(fname,"w")) == NULL) {
      sprintf(msg,"in WriteParCurv_Peru(), could not open:  %s",fname);
      errormsg(0,msg);
    }
  }

  /* write order and number of control points */
  fprintf(fp,"%d %d\n",egm->order,egm->ncontpts);

  for (i = 0; i < egm->ncontpts + 1; i++)   /* knot vector */
    fprintf(fp,"%+.16le\n",egm->knots[i]);

  for (i = 0; i < egm->ncontpts; i++)       /* control points */
    fprintf(fp,"%+.16le %+.16le %+.16le \n\t%+.16le\n",egm->contpts[i]->x,
	   egm->contpts[i]->y,egm->contpts[i]->z,egm->contpts[i]->w);

  if (fpg == TRUE) fclose(fp);
}

/* Function: WriteParSurf_Peru()
 * Purpose: Write a NURBS surface periodic in u to a file
 * Method: The NURBS surface structure is defined in the header file
 *         gen.h
 * Arguments:
 *  fp - file pointer of the open file to contain the surface
 *  fgm - address of the surface structure
 */
/* Functions referenced by WriteParSurf_Peru() are:
 *  errormsg()
 */

void WriteParSurf_Peru (FILE *fp, ParSurf *fgm)
{
     int fpg = FALSE, i, j;
     char fname[61], msg[100];
 
     if (fp == NULL) {   /* if open file not specified, prompt for name */
	  fpg = TRUE;
	  fputs ("WriteParSurf_Peru() output file name : ", stdout);
	  scanf("%60s%*c", fname);
	  if ((fp = fopen (fname, "w")) == NULL) {
	       sprintf (msg, "in WriteParSurf_Peru(), could not open:  %s",
			fname);
	       errormsg (0, msg);
	  }
     }

     /* write u,v order and number of u,v control points */
     fprintf (fp, "%d %d %d %d\n", fgm->uorder, fgm->vorder,
		     fgm->ucontpts,fgm->vcontpts);

     for (i = 0; i < fgm->ucontpts + 1; i++)             /* u knot vector */
       fprintf(fp,"%+.16le\n",fgm->uknots[i]);
     for (i = 0; i < fgm->vorder + fgm->vcontpts; i++)   /* v knot vector */
       fprintf(fp,"%+.16le\n",fgm->vknots[i]);
     
     for (i = 0; i < fgm->ucontpts; i++)                 /* control points */
       for (j = 0; j < fgm->vcontpts; j++) 
	 fprintf(fp,"%+.16le %+.16le %+.16le\n\t%+.16le\n",
		 fgm->contpts[i][j]->x, fgm->contpts[i][j]->y,
		 fgm->contpts[i][j]->z, fgm->contpts[i][j]->w);
     
     if (fpg == TRUE)
	  fclose(fp);
}

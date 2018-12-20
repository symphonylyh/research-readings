/************************************************************************
 *									*
 *   Copyright  (C)  1988  Massachusetts  Institute  of  Technology,    *
 *   Cambridge, MA.  All rights reserved.   Ocean Engineering Design    *
 *   Laboratory.  Last modified  by  Bradley A. Moran  on 08 Dec 88.    *
 *                                                                      *
 ************************************************************************/

#include <stdio.h>
#include "gen.h"

/* Function: WriteParSurf()
 * Purpose: Write a non-periodic NURBS surface to a file
 * Method: The NURBS surface structure is defined in the header file
 *         gen.h
 * Arguments:
 *  fp - file pointer of the open file to contain the surface
 *  fgm - address of the surface structure
 */
/* Functions referenced by WriteParSurf() are:
 *  errormsg()
 */
/* Functions that reference WriteParSurf() are:
 *  main()
 *  sample_csurf()
 */

void WriteParSurf(FILE *fp, ParSurf *fgm)

{
  int fpg = FALSE,i,j;
  char fname[60],msg[100];
 
  if (fp == NULL) {   /* if no open file is specified, prompt for name */
    fpg = TRUE;
    printf("\n\rWriteParSurf() output file name:  "), scanf("%s",fname);
    if ((fp = fopen(fname,"w")) == NULL) {
      sprintf(msg,"in WriteParSurf(), could not open:  %s",fname);
      errormsg(0,msg);
    }
  }

   /* write u,v order and number of u,v control points */
  fprintf(fp,"%d %d %d %d\n", fgm->uorder, fgm->vorder,
	  fgm->ucontpts,fgm->vcontpts);

  for (i = 0; i < fgm->uorder + fgm->ucontpts; i++)   /* u knot vector */
    fprintf(fp,"%+.16le\n",fgm->uknots[i]);
  for (i = 0; i < fgm->vorder + fgm->vcontpts; i++)   /* v knot vector */
    fprintf(fp,"%+.16le\n",fgm->vknots[i]);

  for (i = 0; i < fgm->ucontpts; i++)                 /* control points */
    for (j = 0; j < fgm->vcontpts; j++) 
      fprintf(fp,"%+.16le %+.16le %+.16le \n\t%+.16le\n",
	      fgm->contpts[i][j]->x, fgm->contpts[i][j]->y,
	      fgm->contpts[i][j]->z,fgm->contpts[i][j]->w);

  if (fpg == TRUE) fclose(fp);
}

/* Function: WriteParCurv()
 * Purpose: Write a non-periodic NURBS curve to a file
 * Method: The NURBS curve structure is defined in the header file
 *         gen.h
 * Arguments:
 *  fp - file pointer of the open file to contain the curve
 *  egm - address of the curve structure
 */
/* Functions referenced by WriteParCurv() are:
 *  errormsg()
 */
/* Functions that reference WriteParCurv() are:
 *  main()
 *  sample_csurf()
 */

void WriteParCurv(FILE *fp, ParCurv *egm)
{
  int fpg = FALSE,i;
  char fname[60],msg[100];
 
  if (fp == NULL) {   /* if an open file is not specified, prompt for name */
    fpg = TRUE;
    printf("\n\rWriteParCurv() output file name:  "), scanf("%s",fname);
    if ((fp = fopen(fname,"w")) == NULL) {
      sprintf(msg,"in WriteParCurv(), could not open:  %s",fname);
      errormsg(0,msg);
    }
  }

  /* write order and number of control points */
  fprintf(fp,"%d %d\n", egm->order, egm->ncontpts);

  for (i = 0; i < egm->order + egm->ncontpts; i++)   /* knot vector */
    fprintf(fp,"%+.16le\n",egm->knots[i]);

  for (i = 0; i < egm->ncontpts; i++)                /* control points */
    fprintf(fp,"%+.16le %+.16le %+.16le \n\t%+.16le\n",egm->contpts[i]->x,
	   egm->contpts[i]->y,egm->contpts[i]->z,egm->contpts[i]->w);

  if (fpg == TRUE) fclose(fp);
}

/* Function: Writevector()
 * Purpose: Write a vector to a file
 * Method: The file is not tested to see if it is valid or open
 * Arguments:
 *  fp - file pointer of an open file to contain the vector
 *  v - address of the vector structure
 */

void Writevector(FILE *fp, vector *v)
{
  fprintf(fp,"%+.16le %+.16le %+.16le %+.16le\n", v->x, v->y, v->z, v->w);
}

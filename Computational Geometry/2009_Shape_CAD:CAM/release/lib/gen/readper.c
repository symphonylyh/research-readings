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

     12 Apr 89 - B.A. Moran fixed bug in number of knots
     10 Apr 89 - B.A. Moran fixed bug in reading order
     01 Mar 89 - S.T. Tuohy author
 *									*
 ************************************************************************/

#include <malloc.h>
#include <stdio.h>
#include "gen.h"

/* Because the file may contain comments (introduced by a pound sign "#")
 * we use GetNextToken() to get the next non-comment token.  This is
 * returned as a character string.  We then use the standard library
 * function sscanf() to convert from a character string to a numeric
 * value
 */

/* Function: ReadParCurv_Per()
 * Purpose: Read a periodic NURBS curve from a file
 * Method: The NURBS curve structure is defined in the header file "gen.h"
 * Arguments:
 *  fp - file pointer of the open file containing the curve
 *  egm - address of a NURBS curve structure to hold the curve
 *        If egm is NULL, a new structure is allocated
 * Return: the address of the structure holding the curve
 */
/* Functions referenced by ReadParCurv_Per() are:
 *  copyvector()
 *  dbl_array1()
 *  errormsg()
 *  free_darray1()
 *  free_varray1()
 *  gen_array1()
 *  GetNextToken()
 *  vec_array1()
 */
/* Functions that reference ReadParCurv_Per() are:
 *  ReadDeslabCos()
 *  ReadDeslabCurv()
 */

ParCurv *ReadParCurv_Per (FILE *fp, ParCurv *egm)
{
  vector in;
  int i, fpg = FALSE;
  char *cp, fname[60], msg[100], token[80];
  
  if (fp == NULL) {
    fpg = TRUE;
    printf("\n\rReadParCurv_Per() input file name:  ");
    scanf("%60s", fname);
    if ((fp = fopen(fname, "r")) == NULL) {
      sprintf(msg, "in ReadParCurv_Per(), could not find:  %s", fname);
      errormsg(0, msg);
    }
  }
  if (!egm) {			/* allocate new egeom */
    egm = (ParCurv *) gen_array1(1, sizeof(ParCurv));
    if (!egm)
      errormsg(2, "allocation failure in ReadParCurv_Per()");
    egm->kmem = egm->pmem = 0;
  }

  GetNextToken(fp, token);
  sscanf(token, "%d", &egm->order);
  GetNextToken(fp, token);
  sscanf(token, "%d", &egm->ncontpts);

  if (egm->kmem == 0) {		/* allocate knots if no mem */
    egm->knots = dbl_array1((unsigned) egm->ncontpts + 1);
  }
  else if (egm->kmem < egm->ncontpts + 1) { /* free and allocate if */
					    /* not enough mem */
    free_darray1 (egm->knots);
    egm->knots = dbl_array1((unsigned) egm->ncontpts + 1);
  }
  egm->kmem = egm->ncontpts + 1;
  for (i = 0; i < egm->ncontpts + 1; i++) {   /* knot vector */
    GetNextToken(fp, token);
    sscanf (token, "%lf", &egm->knots[i]);
  }

  if (egm->pmem == 0)		/* allocate points if none */
    egm->contpts = vec_array1 ((unsigned) egm->ncontpts);
  else if (egm->pmem < egm->ncontpts) {
    free_varray1 (egm->contpts, (unsigned) egm->ncontpts);
    egm->contpts = vec_array1 ((unsigned) egm->ncontpts);
  }
  egm->pmem = egm->ncontpts;

  for (i = 0; i < egm->ncontpts; i++) {       /* control points */
    GetNextToken(fp, token);
    sscanf (token, "%lf", &in.x);
    GetNextToken(fp, token);
    sscanf (token, "%lf", &in.y);
    GetNextToken(fp, token);
    sscanf (token, "%lf", &in.z);
    GetNextToken(fp, token);
    sscanf (token, "%lf", &in.w);
    egm->contpts[i] = copyvector(&in, egm->contpts[i]);
  }
  if (fpg == TRUE)
    fclose (fp);
  egm->type = PCurvePer;
  return (egm);
}

/* Function: ReadParSurf()
 * Purpose: Read a non-periodic NURBS surface from a file
 * Method: The NURBS surface structure is defined in the header file "gen.h"
 * Arguments:
 *  fp - file pointer of the open file containing the surface
 *  fgm - address of a NURBS curve structure to hold the surface
 *        If egm is NULL, a new structure is allocated
 * Return: the address of the structure holding the surface
 */
/* Functions referenced by ReadParSurf_Peru() are:
 *  copyvector()
 *  dbl_array1()
 *  errormsg()
 *  free_darray1()
 *  free_varray2()
 *  gen_array1()
 *  GetNextToken()
 *  vec_array2()
 */

ParSurf *ReadParSurf_Peru (FILE *fp, ParSurf *fgm)
{
  vector in;
  int i, j, fpg = FALSE;
  char *cp, fname[60], msg[100], token[80];

  if (fp == NULL) {
    fpg = TRUE;
    printf("\n\rReadParSurf_Peru() input file name:  ");
    scanf("%s", fname);
    if ((fp = fopen(fname, "r")) == NULL) {
      sprintf(msg, "in ReadParSurf_Peru(), could not find:  %s", fname);
      errormsg(0, msg);
    }
  }
  if (!fgm) {
    fgm = (ParSurf *) gen_array1(1, sizeof(ParSurf));
    if (!fgm)
      errormsg(2, "allocation failure in ReadParSurf_Peru()");
    fgm->ukmem = fgm->vkmem = fgm->upmem = fgm->vpmem = 0;
  }
  GetNextToken(fp, token);
  sscanf(token, "%d", &fgm->uorder);
  GetNextToken(fp, token);
  sscanf(token, "%d", &fgm->vorder);
  GetNextToken(fp, token);
  sscanf(token, "%d", &fgm->ucontpts);
  GetNextToken(fp, token);
  sscanf(token, "%d", &fgm->vcontpts);

  if (fgm->ukmem == 0) {
    fgm->ukmem = fgm->ucontpts + 1;
    fgm->uknots = dbl_array1 ((unsigned) fgm->ukmem);
  }

  if (fgm->ukmem < fgm->ucontpts + 1) {
    fgm->ukmem = fgm->ucontpts + 1;
    free_darray1 (fgm->uknots);
    fgm->uknots = dbl_array1 ((unsigned) fgm->ukmem);
  }
  
  if (fgm->vkmem == 0) {
    fgm->vkmem = fgm->vcontpts + fgm->vorder;
    fgm->vknots = dbl_array1 ((unsigned) fgm->vkmem);
  }

  if (fgm->vkmem < fgm->vcontpts + fgm->vorder) {
    fgm->vkmem = fgm->vcontpts + fgm->vorder;
    free_darray1 (fgm->vknots);
    fgm->vknots = dbl_array1 ((unsigned) fgm->vkmem);
  }
  
  if (fgm->upmem == 0 || fgm->vpmem == 0) {
    fgm->contpts = vec_array2((unsigned) fgm->ucontpts, (unsigned)
			      fgm->vcontpts);
    fgm->upmem = fgm->ucontpts;
    fgm->vpmem = fgm->vcontpts;
  }

  if (fgm->upmem < fgm->ucontpts || fgm->vpmem < fgm->vcontpts) {
    free_varray2(fgm->contpts, (unsigned) fgm->ucontpts, (unsigned)
		 fgm->vcontpts);
    fgm->contpts = vec_array2((unsigned) fgm->ucontpts, (unsigned)
			      fgm->vcontpts);
    fgm->upmem = fgm->ucontpts;
    fgm->vpmem = fgm->vcontpts;
  }

  for (i = 0; i < fgm->ucontpts + 1; i++) {             /* u knot vector */
    GetNextToken(fp, token);
    sscanf(token, "%lf", &fgm->uknots[i]);
  }
  for (i = 0; i < fgm->vorder + fgm->vcontpts; i++) {   /* v knot vector */
    GetNextToken(fp, token);
    sscanf(token, "%lf", &fgm->vknots[i]);
  }
  for (i = 0; i < fgm->ucontpts; i++)                   /* control points */
    for (j = 0; j < fgm->vcontpts; j++) {
      GetNextToken(fp, token);
      sscanf(token, "%lf", &in.x);
      GetNextToken(fp, token);
      sscanf(token, "%lf", &in.y);
      GetNextToken(fp, token);
      sscanf(token, "%lf", &in.z);
      GetNextToken(fp, token);
      sscanf(token, "%lf", &in.w);
      fgm->contpts[i][j] = copyvector(&in, fgm->contpts[i][j]);
    }
  if (fpg == TRUE)
    fclose(fp);
  fgm->type = PSurfacePerU;
  return(fgm);
}

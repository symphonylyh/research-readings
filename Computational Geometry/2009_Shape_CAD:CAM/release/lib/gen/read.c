/************************************************************************
 *									*
		Copyright (C) 1989 by Bradley A. Moran
			 All rights reserved.

    The  author hereby grants to  MIT  permission to reproduce and  to
    distribute  copies of  this source  code  document in whole or  in
    part.
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

static char fname[81];

/* Function: ReadParCurv()
 * Purpose: Read a non-periodic NURBS curve from a file
 * Method: The NURBS curve structure is defined in the header file "gen.h"
 * Arguments:
 *  fp - file pointer of the open file containing the curve
 *  egm - address of a NURBS curve structure to hold the curve
 *        If egm is NULL, a new structure is allocated
 * Return: the address of the structure holding the curve
 */
/* Functions referenced by ReadParCurv() are:
 *  copyvector()
 *  dbl_array1()
 *  egeomalloc1()
 *  fileopen()
 *  free_darray1()
 *  free_varray1()
 *  GetNextToken()
 *  vec_array1()
 */
/* Functions that reference ReadParCurv() are:
 *  main()
 *  ReadDeslabCos()
 *  ReadDeslabCurv()
 *  ReadDeslabLocal()
 *  ReadDeslabMinDist()
 *  ReadTrimSurf()
 *  sample_csurf()
 */

ParCurv *ReadParCurv (FILE *fp, ParCurv *egm)
{
  vector in;
  int i, fpg = FALSE;
  int ord, n_pts;
  char *cp, token[80];
  
  if (fp == (FILE *) NULL)
  {				/* open input stream if not done */
    fputs ("Aperiodic B-Spline curve file name : ", stdout);
    scanf ("%60s%*c", fname);
    fp = fileopen (fname, sizeof(fname), "r");
    fpg = TRUE;			/* set file pointer flag */
  }

  GetNextToken(fp, token);
  sscanf (token, "%d", &ord);
  GetNextToken(fp, token);
  sscanf (token, "%d", &n_pts);

  if (!egm)			/* allocate egm if not done */
    egm = egeomalloc1 (ord, n_pts);

  if (egm->kmem < ord + n_pts)
  {				/* knot vector insufficient */
    if (egm->kmem)
      free_darray1 (egm->knots);
    egm->knots = dbl_array1 ((unsigned) ord + n_pts);
    egm->kmem = ord + n_pts;
  }
  if (egm->pmem < n_pts)
  {				/* vector pointers insufficient */
    if (egm->pmem)
      free_varray1 (egm->contpts, (unsigned) egm->ncontpts);
    egm->contpts = vec_array1 ((unsigned) egm->ncontpts);
    egm->pmem = n_pts;
  }
  egm->type = PCurveOpen;
  egm->order = ord;        /* order */
  egm->ncontpts = n_pts;   /* number of control points */
  for (i = 0; i < egm->order + egm->ncontpts; i++) {   /* knot vector */
    GetNextToken(fp, token);
    sscanf (token, "%lf", &egm->knots[i]);
  }
  for (i = 0; i < egm->ncontpts; i++) {   /* control points */
    GetNextToken(fp, token);
    sscanf (token, "%lf", &in.x);
    GetNextToken(fp, token);
    sscanf (token, "%lf", &in.y);
    GetNextToken(fp, token);
    sscanf (token, "%lf", &in.z);
    GetNextToken(fp, token);
    sscanf (token, "%lf", &in.w);
    egm->contpts[i] = copyvector (&in, egm->contpts[i]);
  }

  if (fpg == TRUE)
    fclose (fp);		/* file opened locally, close */
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
/* Functions referenced by ReadParSurf() are:
 *  copyvector()
 *  dbl_array1()
 *  fgeomalloc1()
 *  fgetstring()
 *  fileopen()
 *  free_darray1()
 *  free_varray2()
 *  GetNextToken()
 *  vec_array2()
 */
/* Functions that reference ReadParSurf() are:
 *  main()
 *  ReadDeslabLocal()
 *  ReadDeslabMinDist()
 *  ReadDeslabSurf()
 *  ReadTrimSurf()
 *  sample_csurf()
 */

ParSurf *ReadParSurf (FILE *fp, ParSurf *fgm)

{
  vector in;
  int i, j, fpg = FALSE, rc = NULL;
  int uord, vord, u_p, v_p;
  char *cp, token[80];
  
  if (fp == NULL)
  {				/* open file if not specified */
    fputs ("Aperiodic B-Spline patch file name : ", stdout);
    fgetstring (fname, sizeof(fname), stdin);
    fp = fileopen (fname, sizeof(fname), "r");
    fpg = TRUE;
  }

  GetNextToken(fp, token);
  sscanf (token, "%d", &uord);
  GetNextToken(fp, token);
  sscanf (token, "%d", &vord);
  GetNextToken(fp, token);
  sscanf (token, "%d", &u_p);
  GetNextToken(fp, token);
  sscanf (token, "%d", &v_p);

  if (!fgm)			/* allocate fgm if not specified */
    fgm = fgeomalloc1 (uord, vord, u_p, v_p);

  fgm->uorder = uord;    /* u order */
  fgm->ucontpts = u_p;   /* number of u control points */
  fgm->vorder = vord;    /* v order */
  fgm->vcontpts = v_p;   /* number of v control points */
  if (fgm->ukmem < uord + u_p) {
    if (fgm->ukmem)
      free_darray1 (fgm->uknots);
    fgm->uknots = dbl_array1 ((unsigned) uord + u_p);
    fgm->ukmem = uord + u_p;
  }
  if (fgm->vkmem < vord + v_p) {
    if (fgm->vkmem)
      free_darray1 (fgm->vknots);
    fgm->vknots = dbl_array1 ((unsigned) vord + v_p);
    fgm->vkmem = vord + v_p;
  }
  if (fgm->upmem < u_p || fgm->vpmem < v_p)
  {
    if (fgm->contpts)
      free_varray2 (fgm->contpts, (unsigned) fgm->ucontpts, (unsigned)
		    fgm->vcontpts);
    fgm->contpts = vec_array2 ((unsigned) u_p, (unsigned) v_p);
    fgm->upmem = u_p, fgm->vpmem = v_p;
  }
  fgm->type = PSurfaceOpen;
  for (i = 0; i < fgm->uorder + fgm->ucontpts; i++) {   /* u knot vector */
    GetNextToken(fp, token);
    sscanf (token, "%lf", &fgm->uknots[i]);
  }
  for (i = 0; i < fgm->vorder + fgm->vcontpts; i++) {   /* v knot vector */
    GetNextToken(fp, token);
    sscanf (token, "%lf", &fgm->vknots[i]);
  }
  for (i = 0; i < fgm->ucontpts; i++)                   /* control points */
    for (j = 0; j < fgm->vcontpts; j++) {
      GetNextToken(fp, token);
      sscanf (token, "%lf", &in.x);
      GetNextToken(fp, token);
      sscanf (token, "%lf", &in.y);
      GetNextToken(fp, token);
      sscanf (token, "%lf", &in.z);
      GetNextToken(fp, token);
      sscanf (token, "%lf", &in.w);
      fgm->contpts[i][j] = copyvector (&in, fgm->contpts[i][j]);
    }

  if (fpg == TRUE)		/* close file if opened in this routine */
    fclose (fp);

  return (fgm);
}

/* Function: GetNextToken()
 * Purpose: Get the next token from a file
 * Method: Files may contain comments, which start with a pound sign "#"
 *         and continue till the end of the line.
 * Arguments:
 *  fp - file pointer of the open file
 *  token - address of a character string to hold the next token
 * Return: The number of characters in the token.  If End-Of-File is
 *         reached, return EOF, which is usually -1
 */
/* Functions that reference GetNextToken() are:
 *  main()
 *  ReadDeslabFacet()
 *  ReadDeslabHull()
 *  ReadDeslabLocal()
 *  ReadDeslabMinDist()
 *  ReadDeslabMult()
 *  ReadDeslabVect()
 *  ReadGridSurf()
 *  ReadInitFile()
 *  ReadListCurv()
 *  ReadParCurv()
 *  ReadParCurv_Per()
 *  ReadParSurf()
 *  ReadParSurf_Peru()
 *  ReadParUv()
 *  ReadPraxDirectory()
 *  ReadTrimSurf()
 *  sample_csurf()
 */

int GetNextToken(FILE *fp, char *token)
{
  int i, iret, n;
  char c;

  /* recall that fscanf(fp, "%s", token) reads the next sequence of
   * non-whitespace characters.  Also, it returns the number of characters
   * read, unless End-Of-File is reached, when it returns EOF
   */

  iret = fscanf(fp, "%s", token); /* read the next token */

  while (token[0] == (unsigned)'#' && iret != EOF) { /* start of comment */
    c = token[1];                              /* keep reading 1 character */
    while (c != (unsigned)'\n' && iret != EOF) /* at a time until newline */
      iret = fscanf(fp, "%c", &c);
    iret = fscanf(fp, "%s", token);   /* next token after comment */
    /* make sure to check if it is a comment */
  }

  if (token[0] == '"') {          /* string token, may have embdedded blanks */
    for (i=1; token[i]; i++)      /* shift by one to remove initial " */
      token[i-1] = token[i];

    i -= 2;
    while (token[i] != '"' && iret != EOF)   /* there are embedded blanks */
      iret = fscanf(fp, "%c", &token[++i]);  /* continue until closing " */
    token[i] = '\0';
  }

  return iret;
}

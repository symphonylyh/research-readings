/* Copyright (C) 1992 Massachusetts Institute of Technology, Cambridge, MA
   All rights reserved

   Last edited: April 11, 2003 for v.11.1

   v.11.1: (1) WriteTrimSurf(...) changed.

*/

/* write.c */

/* WriteFacet()
   WriteFacetSurface()
   WriteGridSurf()
   WriteInsp()
   WriteListCurv()
   WriteParCurv2()
   WriteParCurv2_Per()
   WriteParSurf2()
   WriteParUv()
   WriteTrimSurf()
*/

#include <stdio.h>
#include "gen.h"
#include "editor.h"

void WriteListCurv(FILE *fp, ListCurv *list, short e_format)
{
  int i;

  fprintf(fp, "%d\n", list->npoints);
  for (i=0; i<list->npoints; i++)
    if (e_format)
      fprintf(fp, "%+.15e %+.15e %+.15e\n", list->pts[i][0], list->pts[i][1],
	      list->pts[i][2]);
    else
      fprintf(fp, "%9f %9f %9f\n", list->pts[i][0], list->pts[i][1],
	      list->pts[i][2]);
}

void WriteGridSurf(FILE *fp, GridSurf *grid, short e_format)
{
  int i, j;

  fprintf(fp, "%d %d\n", grid->ncol, grid->nrow);
  for (i=0; i<grid->ncol; i++)
    for (j=0; j<grid->nrow; j++)
      if (e_format)
	fprintf(fp, "%+.15e %+.15e %+.15e\n", grid->pts[i][j][0],
		grid->pts[i][j][1], grid->pts[i][j][2]);
      else
	fprintf(fp, "%9f %9f %9f\n", grid->pts[i][j][0],
		grid->pts[i][j][1], grid->pts[i][j][2]);
}

void WriteParCurv2(FILE *fp, ParCurv *egeom, short e_format)
{
  int i;
 
  fprintf(fp, "%d %d\n", egeom->order, egeom->ncontpts);
  for (i=0; i < egeom->order + egeom->ncontpts; i++)
    if (e_format)
      fprintf(fp, "%+.15e\n", egeom->knots[i]);
    else
      fprintf(fp, "%9f\n", egeom->knots[i]);
  for (i=0; i < egeom->ncontpts; i++)
    if (e_format)
      fprintf(fp,"%+.15e %+.15e %+.15e\n\t%+.15e\n", egeom->contpts[i]->x,
	      egeom->contpts[i]->y, egeom->contpts[i]->z,
	      egeom->contpts[i]->w);
    else
      fprintf(fp, "%9f %9f %9f %9f\n", egeom->contpts[i]->x,
	      egeom->contpts[i]->y, egeom->contpts[i]->z,
	      egeom->contpts[i]->w);
}

void WriteParCurv2_Per(FILE *fp, ParCurv *egeom, short e_format)
{
  int i;

  fprintf(fp,"%d %d\n", egeom->order, egeom->ncontpts);
  for (i=0; i < egeom->ncontpts + 1; i++)
    if (e_format)
      fprintf(fp,"%+.15e\n",egeom->knots[i]);
    else
      fprintf(fp,"%9f\n",egeom->knots[i]);
  for (i=0; i < egeom->ncontpts; i++)
    if (e_format)
      fprintf(fp,"%+.15e %+.15e %+.15e\n\t%+.15e\n", egeom->contpts[i]->x,
	      egeom->contpts[i]->y, egeom->contpts[i]->z,
	      egeom->contpts[i]->w);
    else
      fprintf(fp,"%9f %9f %9f %9f\n", egeom->contpts[i]->x,
	      egeom->contpts[i]->y, egeom->contpts[i]->z,
	      egeom->contpts[i]->w);
}

void WriteParSurf2(FILE *fp, ParSurf *fgeom, short e_format)
{
  int i, j;
 
  fprintf(fp,"%d %d %d %d\n", fgeom->uorder, fgeom->vorder,
	  fgeom->ucontpts, fgeom->vcontpts);
  for (i=0; i < fgeom->uorder + fgeom->ucontpts; i++)
    if (e_format)
      fprintf(fp, "%+.15e\n", fgeom->uknots[i]);
    else
      fprintf(fp, "%9f\n", fgeom->uknots[i]);
  for (i=0; i < fgeom->vorder + fgeom->vcontpts; i++) 
    if (e_format)
      fprintf(fp, "%+.15e\n", fgeom->vknots[i]);
    else
      fprintf(fp, "%9f\n", fgeom->vknots[i]);
  for (i=0; i < fgeom->ucontpts; i++)
    for (j = 0; j < fgeom->vcontpts; j++)
      if (e_format)
	fprintf(fp, "%+.15e %+.15e %+.15e\n\t%+.15e\n",
		fgeom->contpts[i][j]->x, fgeom->contpts[i][j]->y,
		fgeom->contpts[i][j]->z, fgeom->contpts[i][j]->w);
      else
	fprintf(fp, "%9f %9f %9f %9f\n",
		fgeom->contpts[i][j]->x, fgeom->contpts[i][j]->y,
		fgeom->contpts[i][j]->z, fgeom->contpts[i][j]->w);
}

void WriteParUv(FILE *fp, ParUv *ugeom, short e_format)
{
  int i;

  if (e_format) {
    fprintf(fp, "%+.15e %+.15e\n", ugeom->umin, ugeom->umax);
    fprintf(fp, "%+.15e %+.15e\n", ugeom->vmin, ugeom->vmax);
  }
  else {
    fprintf(fp, "%9f %9f\n", ugeom->umin, ugeom->umax);
    fprintf(fp, "%9f %9f\n", ugeom->vmin, ugeom->vmax);
  }
  fprintf(fp, "%d\n", ugeom->npoints);
  for (i=0; i<ugeom->npoints; i++)
    if (e_format)
      fprintf(fp, "%+.15e %+.15e\n", ugeom->pts[i][0], ugeom->pts[i][1]);
    else
      fprintf(fp, "%9f %9f\n", ugeom->pts[i][0], ugeom->pts[i][1]);
  fprintf(fp, "%d\n", ugeom->is_open);
}

void WriteFacet(FILE *fp, int nPts, double **pts, int *iEnd,
		int nAdj, int *iAdj, short e_format)
{
  int i, j, n = 0;

  fprintf(fp, "%d\n", nPts);
  for (i=0; i<nPts; i++)
    if (e_format)
      fprintf(fp, "%+.15e %+.15e %+.15e\n", pts[i][0], pts[i][1], pts[i][2]);
    else
      fprintf(fp, "%9f %9f %9f\n", pts[i][0], pts[i][1], pts[i][2]);
  for (i=0; i<nPts; i++)
    fprintf(fp, "%d\n", iEnd[i]+1);
  for (i=0; i<nPts; i++) {
    fprintf(fp, "%d", iAdj[n]+1);
    for (j=n+1; j<=iEnd[i]; j++)
      fprintf(fp, " %d", iAdj[j]+1);
    fprintf(fp, "\n");
    n = j;

  }
}

void WriteFacetSurf(FILE *fp, FacetSurf *tgeom, short e_format)
{
  int i, j, n = 0;

  fprintf(fp, "%d\n", tgeom->nPts);
  for (i=0; i<tgeom->nPts; i++)
    if (e_format)
      fprintf(fp, "%+.15e %+.15e %+.15e\n", tgeom->pts[i][0],
	      tgeom->pts[i][1], tgeom->pts[i][2]);
    else
      fprintf(fp, "%9f %9f %9f\n", tgeom->pts[i][0], tgeom->pts[i][1],
	      tgeom->pts[i][2]);
  for (i=0; i<tgeom->nPts; i++)
    fprintf(fp, "%d\n", tgeom->iEnd[i]+1);
  for (i=0; i<tgeom->nPts; i++) {
    fprintf(fp, "%d", tgeom->iAdj[n]+1);
    for (j=n+1; j<=tgeom->iEnd[i]; j++)
      fprintf(fp, " %d", tgeom->iAdj[j]+1);
    fprintf(fp, "\n");
    n = j;
  }
}

void WriteInsp(FILE *fp, int nPts, double **pts, double **norms,
	       short e_format)
{
  int i;

  fprintf(fp, "%d\n", nPts);
  for (i=0; i<nPts; i++)
    if (e_format)
      fprintf(fp, "%+.15e %+.15e %+.15e\n\t%.15e %+.15e %+.15e\n",
	      pts[i][0], pts[i][1], pts[i][2],
	      norms[i][0], norms[i][1], norms[i][2]);
    else
      fprintf(fp, "%9f %9f %9f  %9f %9f %9f\n",
	      pts[i][0], pts[i][1], pts[i][2],
	      norms[i][0], norms[i][1], norms[i][2]);
}



/***********************************************
  "WriteTrimSurf(...)" is modified in v.11.1
**********************************************/
void WriteTrimSurf(FILE *fp, TrimSurf *trim, short e_format)
{
  int i, j;

  WriteParSurf2(fp, trim->fgeom, e_format);
  fprintf(fp, "%d\n%d\n", trim->bLoop, trim->nLoops);
  for (i=0; i<trim->nLoops; i++) {
    fprintf(fp, "%d\n", trim->loops[i].nCurves);
    for (j=0; j<trim->loops[i].nCurves; j++)
      WriteParCurv2(fp, trim->loops[i].egeoms[j], e_format);
  }

  /*************** 11.1 begin ***************/
  if(trim->c_loops) { /* C-curve info. available */
    PrintLog("--> \"C-curves are provided explicitly.\"\n");
    for (i=0; i<trim->nLoops; i++) {
      fprintf(fp, "%d\n", trim->c_loops[i].nCurves);
      for (j=0; j<trim->c_loops[i].nCurves; j++)
	WriteParCurv2(fp, trim->c_loops[i].egeoms[j], e_format);
    }
  }

  else {
    PrintLog("--> \"C-curves are not provided explicitly.\"\n");
  }
  /*************** 11.1 end *****************/
}


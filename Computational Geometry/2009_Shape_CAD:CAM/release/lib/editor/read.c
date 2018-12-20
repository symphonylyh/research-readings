/* Copyright (C) 1992 Massachusetts Institute of Technology, Cambridge, MA
   All rights reserved

   Last edited: April 11, 2003 for v.11.1

   v.11.1: (1) ReadTrimSurf(...) changed.
*/

/* read.c */

/* ReadFacetSurface()
   ReadGridSurf()
   ReadListCurv()
   ReadParUv()
   ReadTrimSurf()
*/

#include <stdio.h>
#include <malloc.h>
#include "gen.h"
#include "editor.h"

FacetSurf *ReadFacetSurf(FILE *fp, FacetSurf *facet)
{
  FacetSurf *tgeom;
  char token[80];
  int i;

  if (facet)
    tgeom = facet;
  else
    tgeom = (FacetSurf *)gen_array1(1, sizeof(FacetSurf));

  GetNextToken(fp, token);
  sscanf(token, "%d", &tgeom->nPts);
  tgeom->pts = dbl_array2(tgeom->nPts, 3);
  for (i=0; i<tgeom->nPts; i++) {
    GetNextToken(fp, token);
    sscanf(token, "%lf", &tgeom->pts[i][0]);
    GetNextToken(fp, token);
    sscanf(token, "%lf", &tgeom->pts[i][1]);
    GetNextToken(fp, token);
    sscanf(token, "%lf", &tgeom->pts[i][2]);
  }
  tgeom->iEnd = int_array1(tgeom->nPts);
  for (i=0; i<tgeom->nPts; i++) {
    GetNextToken(fp, token);
    sscanf(token, "%d", &tgeom->iEnd[i]);
    tgeom->iEnd[i] -= 1;
  }
  tgeom->nAdj = tgeom->iEnd[tgeom->nPts-1] + 1;
  tgeom->iAdj = int_array1(tgeom->nAdj);
  for (i=0; i<tgeom->nAdj; i++) {
    GetNextToken(fp, token);
    sscanf(token, "%d", &tgeom->iAdj[i]);
    tgeom->iAdj[i] -= 1;
  }

  return tgeom;
}

ListCurv *ReadListCurv(FILE *fp, ListCurv *list)
{
  ListCurv *lgeom;
  int i, j;
  char token[80];

  if (list)
    lgeom = list;
  else
    lgeom = (ListCurv *)gen_array1(1, sizeof(ListCurv));

  lgeom->knottype = 0;
  lgeom->isotype = 0;
  GetNextToken(fp, token);
  sscanf(token, "%d", &lgeom->npoints);
  lgeom->knots = lgeom->npoints;
  lgeom->pts = dbl_array2(lgeom->npoints, 3);
  for (i=0; i<lgeom->npoints; i++)
    for (j=0; j<3; j++) {
      GetNextToken(fp, token);
      sscanf(token, "%lf", &lgeom->pts[i][j]);
    }

  return (lgeom);
}

GridSurf *ReadGridSurf(FILE *fp, GridSurf *grid)
{
  GridSurf *ggeom;
  int i, j, k;
  char token[80];

  if (grid)
    ggeom = grid;
  else
    ggeom = (GridSurf *)gen_array1(1, sizeof(GridSurf));

  ggeom->knottype = 0;
  ggeom->isotype = 0;
  GetNextToken(fp, token);
  sscanf(token, "%d", &ggeom->ncol);
  GetNextToken(fp, token);
  sscanf(token, "%d", &ggeom->nrow);
  ggeom->uknots = ggeom->nrow;
  ggeom->vknots = ggeom->ncol;
  ggeom->pts = dbl_array3(ggeom->ncol, ggeom->nrow, 3);
  for (i=0; i<ggeom->ncol; i++)
    for (j=0; j<ggeom->nrow; j++)
      for (k=0; k<3; k++) {
	GetNextToken(fp, token);
	sscanf(token, "%lf", &ggeom->pts[i][j][k]);
      }

  return (ggeom);
}

ParUv *ReadParUv(FILE *fp, ParUv *uv)
{
  ParUv *ugeom;
  int i;
  char token[80];

  if (uv)
    ugeom = uv;
  else
    ugeom = (ParUv *)gen_array1(1, sizeof(ParUv));

  GetNextToken(fp, token);
  sscanf(token, "%lf", &ugeom->umin);
  GetNextToken(fp, token);
  sscanf(token, "%lf", &ugeom->umax);
  GetNextToken(fp, token);
  sscanf(token, "%lf", &ugeom->vmin);
  GetNextToken(fp, token);
  sscanf(token, "%lf", &ugeom->vmax);
  GetNextToken(fp, token);
  sscanf(token, "%d", &ugeom->npoints);
  ugeom->pts = dbl_array2(ugeom->npoints, 2);
  for (i=0; i<ugeom->npoints; i++) {
    GetNextToken(fp, token);
    sscanf(token, "%lf", &ugeom->pts[i][0]);
    GetNextToken(fp, token);
    sscanf(token, "%lf", &ugeom->pts[i][1]);
  }
  GetNextToken(fp, token);
  sscanf(token, "%hd", &ugeom->is_open);

  return (ugeom);
}



/***********************************************
  "ReadTrimSurf(...)" is modified in v.11.1
**********************************************/
TrimSurf *ReadTrimSurf(FILE *fp, TrimSurf *trim)
{
  int i, j;
  char token[80];

  if (!trim)
    trim = (TrimSurf *)gen_array1(1, sizeof(TrimSurf));

  trim->fgeom = ReadParSurf(fp, NULL);
  GetNextToken(fp, token);
  sscanf(token, "%hd", &trim->bLoop);
  GetNextToken(fp, token);
  sscanf(token, "%d", &trim->nLoops);
  trim->loops = (TrimLoop *)gen_array1(trim->nLoops, sizeof(TrimLoop));
  for (i=0; i<trim->nLoops; i++) {
    GetNextToken(fp, token);
    sscanf(token, "%d", &trim->loops[i].nCurves);
    trim->loops[i].egeoms = (ParCurv **)gen_array1(trim->loops[i].nCurves,
						   sizeof(ParCurv *));
    for (j=0; j<trim->loops[i].nCurves; j++)
      trim->loops[i].egeoms[j] = ReadParCurv(fp, NULL);
  }

  /*************** 11.1 begin ***************/
  /* GetNextToken(...) is defined in "lib/gen/read.c" */
  /* and it reads all the consecutive #-comment lines */
  /* If it encounters EOF, it returns EOF(=-1).       */
  if( GetNextToken(fp, token) == EOF ) { /* No C-curves provided */
    trim->c_loops = NULL; /* set it NULL */
    PrintLog("--> \"C-curves are not provided explicitly.\"\n");
    /*--> PrintLog throws a message to both stdout and .log file. */
  }

  else {
    PrintLog("--> \"C-curves are provided explicitly.\"\n");
    trim->c_loops = (TrimLoop *)gen_array1(trim->nLoops, sizeof(TrimLoop));
    for (i=0; i<trim->nLoops; i++) {
    /*GetNextToken(fp, token);*/ /*---> should be commented out as
                                        GetNextToken was already called in the
                                        above EOF check - it can't be called
                                        twice in a row. Input function must be
                                        called in between. So, GetNextToken
                                        is placed at the last of this for-loop
                                        - see below. */
      sscanf(token, "%d", &trim->c_loops[i].nCurves);
      trim->c_loops[i].egeoms = (ParCurv**)gen_array1(trim->c_loops[i].nCurves,
						      sizeof(ParCurv *));
      for (j=0; j<trim->c_loops[i].nCurves; j++)
	trim->c_loops[i].egeoms[j] = ReadParCurv(fp, NULL);

      GetNextToken(fp, token); /* Yes, here the GetNextToken is */
    }

  }
  /*************** 11.1 end *****************/

  return (trim);
}

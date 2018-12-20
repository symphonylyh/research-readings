/* Copyright (C) 1992 Massachusetts Institute of Technology, Cambridge, MA
   All rights reserved
*/

/* rev_param.c */

/* RevCurvParam()
   RevSurfUParam()
   RevSurfVParam()
   SwitchFacetSurf()
   SwitchSurfParam()
*/

#include "gen.h"
#include "editor.h"

void RevCurvParam(ParCurv *egeom)
{
  ParCurv *egm;
  int i;

  egm = copyegeom(egeom, 0);

  for(i=0; i < egm->order + egm->ncontpts;i ++)
    egeom->knots[i] =  1.0 - egm->knots[egm->order + egm->ncontpts - 1 - i];

  for(i=0; i < egm->ncontpts; i++)
    copyvector(egm->contpts[egm->ncontpts-1-i], egeom->contpts[i]);

  free_egeom(egm);
}

void RevSurfUParam(ParSurf *fgeom)
{
  ParSurf *fgm;
  int i, j;

  fgm = copyfgeom(fgeom, 0);

  for(i=0; i<fgm->uorder + fgm->ucontpts; i++)
    fgeom->uknots[i] = 1.0 - fgm->uknots[fgm->uorder + fgm->ucontpts - 1 - i];

  for(i=0; i< fgm->ucontpts; i++)
    for(j=0; j< fgm->vcontpts; j++)
      copyvector(fgm->contpts[fgm->ucontpts-1-i][j], fgeom->contpts[i][j]);

  free_fgeom(fgm);
}

void RevSurfVParam(ParSurf *fgeom)
{
  ParSurf *fgm;
  int i, j;

  fgm = copyfgeom(fgeom, 0);

  for(i=0; i<fgm->vorder + fgm->vcontpts; i++)
    fgeom->vknots[i] = 1.0 - fgm->vknots[fgm->vorder + fgm->vcontpts - 1 - i];

  for(i=0; i< fgm->ucontpts; i++)
    for(j=0; j< fgm->vcontpts; j++)
      copyvector(fgm->contpts[i][fgm->vcontpts-1-j], fgeom->contpts[i][j]);

  free_fgeom(fgm);
}

void SwitchFacetSurf(FacetSurf *tgeom)
{
  SwitchFacetSurf2(tgeom->nPts, tgeom->pts, tgeom->iEnd, tgeom->iAdj);
}

void SwitchFacetSurf2(int nPts, double **pts, int *iEnd, int *iAdj)
{
  double t;
  int *adj, i, j, last, n, n1, n2;

  for (i=0; i<nPts; i++) {   /* switch x and y coordinates */
    t = pts[i][0];
    pts[i][0] = pts[i][1];
    pts[i][1] = t;
  }

  /* reverse the order of the adjacency lists to preserve counter-
   * clockwise order
   * for point i, the array iAdj[j], j from iEnd[i-1] to iEnd[i], contains
   * the indices of the adjacent points 
   * recall that a -1 at the end of the list indicates that the
   * point is on the boundary
   */

  last = -1;                        /* ending index of last point */
  for (i=0; i<nPts; i++) {
    n1 = last + 1;                                    /* start index of list */
    n2 = ((iAdj[iEnd[i]] < 0) ? iEnd[i]-1 : iEnd[i]); /* end index of list */
    n = n2 - last;                  /* size of adjacency list */
    adj = int_array1(n);            /* temp copy of list */

    for (j=0; j<n; j++)             /* copy adjacency list */
      adj[j] = iAdj[j+n1];
    for (j=0; j<n; j++)             /* reverse the order of the adjacencies */
      iAdj[j+n1] = adj[n-1-j];

    free_iarray1(adj);              /* free temp storage */
    last = iEnd[i];                 /* update end index for next iteration */
  }
}

void SwitchSurfParam(ParSurf **fgeom)
{
  ParSurf *fgm;
  int i, j;

  fgm = fgeomalloc1((*fgeom)->vorder, (*fgeom)->uorder, (*fgeom)->vcontpts, 
		    (*fgeom)->ucontpts);
  for (i=0; i < (*fgeom)->uorder + (*fgeom)->ucontpts; i++)
    fgm->vknots[i] = (*fgeom)->uknots[i];

  for (i=0; i < (*fgeom)->vorder + (*fgeom)->vcontpts; i++)
    fgm->uknots[i] = (*fgeom)->vknots[i];

  for (i=0; i<(*fgeom)->ucontpts; i++)
    for (j=0; j<(*fgeom)->vcontpts; j++)
      copyvector((*fgeom)->contpts[i][j], fgm->contpts[j][i]);

  free_fgeom(*fgeom);
  (*fgeom) = copyfgeom(fgm, 0);
}

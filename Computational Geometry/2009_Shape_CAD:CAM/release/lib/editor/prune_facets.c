/* Copyright (C) 1994 Massachusetts Institute of Technology, Cambridge, MA
   All rights reserved
*/

/* prune_facets.c */

/* PruneFacets()
*/

#include <stdio.h>
#include <unistd.h>
#include <malloc.h>
#include "gen.h"
#include "editor.h"

extern char *tempnam(const char *, const char *);

/* iOut[i][0] is 1 if point i is on alpha-hull */
/* iOut[i][1] is previous point (before i) on alpha-hull */
/* iOut[i][2] is next point (after i) on alpha-hull */

void PruneFacets(int nPts, double **pts, int **iEnd, int **iAdj, int nAdj,
		 int nOut, int **iOut)
{
  FILE *fp;
  int end, i, *iEnd2, j, k1, k2, kk, start;
  char *tmpfil;

  iEnd2 = int_array1(nPts);      /* allocate working copies */

  tmpfil = tempnam("/usr/tmp", "prune");
  fp = fopen(tmpfil, "w");

  for (i=kk=0; i<nPts; i++) {    /* for each point ... */
    start = (i ? (*iEnd)[i-1]+1 : 0);
    end = (*iEnd)[i];
    if ((*iAdj)[end] > -1 && !iOut[i][0]) {

      /* Case A. Point not on boundary and not on alpha-hull */
      /*         Copy the entire adjacency list */

      for (j=start; j<=end; j++,kk++)
	fprintf(fp, "%d\n", (*iAdj)[j]);
    }
    else if ((*iAdj)[end] < 0 && iOut[i][0]) {

      /* Point on boundary and on alpha-hull */

      if ((*iAdj)[(*iEnd)[iOut[i][1]]] < 0 &&
	  (*iAdj)[(*iEnd)[iOut[i][2]]] < 0) {

	/* Case B. Locally convex, i.e. alpha-hull predecessor and successor */
	/*         points on also on boundary */
	/*         Copy the entire adjacency list */

	for (j=start; j<=end; j++,kk++)
	  fprintf(fp, "%d\n", (*iAdj)[j]);
      }
      else {

	/* Case C. Locally concave, i.e. alpha-hull predecessor or successor */
	/*         points not both on boundary */
	/*         Copy adjacency list from successor to predecessor */

	k1 = k2 = -1;
	for (j=start; j<end; j++) {
	  if ((*iAdj)[j] == iOut[i][2]) k1 = j;
	  if ((*iAdj)[j] == iOut[i][1]) k2 = j;
	}

	if (k1 > -1 && k2 > -1) {
	  for (j=k1; j<=k2; j++,kk++)
	    fprintf(fp, "%d\n", (*iAdj)[j]);
	  fprintf(fp, "%d\n", -1);
	  kk++;
	}
	else
	  printf("PruneFacets(1): i = %d k1,k2 = %d %d iOut = %d %d\n", i+1,
		 k1+1, k2+1, iOut[i][1]+1, iOut[i][2]+1);
      }
    }
    else { 

      /* Case D. Point inside boundary but on alpha-hull. */
      /*         Copy adjacency list from successor to predecessor */

      k1 = k2 = -1;
      for (j=start; j<=end; j++) {
	if ((*iAdj)[j] == iOut[i][2]) k1 = j;
	if ((*iAdj)[j] == iOut[i][1]) k2 = j;
      }

      if (k1 > -1 && k2 > -1) {
	if (k1 <= k2) {
	  for (j=k1; j<=k2; j++,kk++)
	    fprintf(fp, "%d\n", (*iAdj)[j]);
	}
	else {
	  for (j=k1; j<=end; j++,kk++)
	    fprintf(fp, "%d\n", (*iAdj)[j]);
	  for (j=start; j<=k2; j++,kk++)
	    fprintf(fp, "%d\n", (*iAdj)[j]);
	}
	fprintf(fp, "%d\n", -1);
	kk++;
      }
      else
	printf("PruneFacets(2): i = %d k1,k2 = %d %d iOut = %d %d\n", i+1,
	       k1+1, k2+1, iOut[i][1]+1, iOut[i][2]+1);
    }
    iEnd2[i] = kk-1;
  }
  free_iarray1(*iEnd);           /* delete old arrays */
  free_iarray1(*iAdj);
  fclose(fp);

  *iEnd = iEnd2;                 /* substitute new iEnd array */

  nAdj = (*iEnd)[nPts-1] + 1;    /* build new adjacency matrix */
  *iAdj = int_array1(nAdj);

  fp = fopen(tmpfil, "r");
  for (i=0; i<nAdj; i++)
    fscanf(fp, "%d", &(*iAdj)[i]);
  fclose(fp);

  unlink(tmpfil);
  free(tmpfil);
}

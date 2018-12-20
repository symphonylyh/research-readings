/* Copyright (c) 1996 by Massachusetts Institute of Technology
 * All rights reserved
 */

/* Constrained Delaunay triangulation
 * See A. K. Cline and R. J. Renka, "A Constrained Two-Dimensional
 * Triangulation and the Solution of Closest Node Problems in the
 * Presence of Barriers," SIAM Journal of Numerical Analysis 27(5):
 * 1305-1321, October 1990.
 */

/* This program uses a version of "tripack.f" modified by deleting
 * several functions and subroutines that are not needed.
 * From the TRIPACK package, as implemented by R. J. Renka,
 * Department of Computer Science, University of North Texas.
 */

#include <stdio.h>
#include <unistd.h>
#include <malloc.h>
#include "gen.h"

/* define the SVR4/XOPEN function "tempnam" */
extern char *tempnam(const char *, const char *);

/* create Deslab Facet structure
 */

void toFacet(int ncc, int *lcc, int n, int *list, int *lptr, int *lend,
	     int **iEnd, int **iAdj)
{
  FILE *fp;
  int c, i, j, k,  /* loop counters */
      lmn, lmx,
      lp,          /* index of next node in adjacency list */
      lpl,         /* index of last node in adjacency list */
      m,
     *nabor,       /* list of adjacent nodes */
      nd,          /* node index */
     *temp;        /* temporary copy of adjacency list */
  int i1, i2, i3, imax, imin, nEnd;
  char quit, *tmpfil;

  *iEnd = int_array1(n);
  nEnd = -1;

  tmpfil = tempnam("/usr/tmp", "cdt");
  fp = fopen(tmpfil, "w");
  for (i=0; i< (ncc ? lcc[0]-1 : n); i++) {
    lp = lpl = lend[i];   /* index of last node in adjacency list */

    k = 0;                /* size of this adjacency list */
    do {                  /* build adjacency list */
      lp = lptr[lp-1];    /* pointer to index of next node in list */
      nd = list[lp-1];    /* index of next node in list */
      k++;
    } while (lp != lpl);

    if (nd <= 0)          /* this flags the boundary nodes */
      k++;

    nabor = int_array1(k+1);

    k = 0;                /* size of this adjacency list */
    do {                  /* build adjacency list */
      lp = lptr[lp-1];    /* pointer to index of next node in list */
      nd = list[lp-1];    /* index of next node in list */
      nabor[k++] = nd;
    } while (lp != lpl);

    if (nd <= 0) {        /* this flags the boundary nodes */
      nabor[k-1] = -nd;
      nabor[k++] = 0;
    }

    /* print adjacency list indices */

    for (j=0; j<k; j++)   /* print non-constraint curve nodes */
      fprintf(fp, " %2d", nabor[j]);
    fprintf(fp, "\n");

    free_iarray1(nabor);

    nEnd += k;
    (*iEnd)[i] = nEnd;    /* keep track of size of adjacency list array */
  }

  for (c=0; c<ncc; c++) {
    lmn = lcc[c] - 1;
    lmx = ((c < ncc-1) ? lcc[c+1]-2 : n-1);

    for (i=lmn; i<=lmx; i++) {
      lp = lpl = lend[i];   /* index of last node in adjacency list */

      k = 0;                /* size of this adjacency list */
      do {                  /* build adjacency list */
	lp = lptr[lp-1];    /* pointer to index of next node in list */
	nd = list[lp-1];    /* index of next node in list */
	k++;
      } while (lp != lpl);

      if (nd <= 0)          /* this flags the boundary nodes */
	k++;

      nabor = int_array1(k+1);
      temp  = int_array1(k+1);

      k = 0;                /* size of this adjacency list */
      do {                  /* build adjacency list */
	lp = lptr[lp-1];    /* pointer to index of next node in list */
	nd = list[lp-1];    /* index of next node in list */
	nabor[k++] = nd;
      } while (lp != lpl);

      if (nd <= 0) {        /* this flags the boundary nodes */
	nabor[k-1] = -nd;
	nabor[k++] = 0;
      }

      if (i == lmn)         /* index of the preceding node in list */ 
	nd = lmx+1;
      else
	nd = i;

      for (j=0; j<k; j++)
	if (nd == nabor[j]) {   /* find the preceding node in list */
	  lp = j;
	  break;
	}

      m = 0;                    /* reorder list starting with prev node */
      for (j=lp; j<k; j++,m++)
	temp[m] = nabor[j];
      for (j=0; j<lp; j++,m++)
	temp[m] = nabor[j];

      if (i == lmx)             /* nd is index of next node on curve */
	nd = lmn + 1;           /* if i is last of list, nd is the first */
      else
	nd = i + 2;             /* else nd is the next in the list */

      /* recall that the node index i is offset-0 */
      /* while the adjacency list nabor (and temp) is offset-1 */

      nabor[0] = temp[0];
      fprintf(fp, " %2d", nabor[0]);     /* print first adjacent node */

      quit = 0;
      for (j=1; j<k && !quit; j++) {
	if (temp[j] == nd) {
	  nabor[j] = temp[j];
	  quit = 1;
	}
	else if (temp[j] == 0) {
	  nabor[j] = temp[j];
	  quit = 1;
	}
	else if (lmn < temp[j] && temp[j] <= lmx+1) {
	  i1 = i+1;
	  i2 = temp[j-1];
	  i3 = temp[j];

	  imin = imax = i1;
	  if (i2 < imin) imin = i2;
	  if (i3 < imin) imin = i3;
	  if (i2 > imax) imax = i2;
	  if (i3 > imax) imax = i3;

	  if (imin > lmn && imax <= lmx+1 &&
	      ((imin == i1 && imax == i3) ||
	       (imin == i2 && imax == i1) ||
	       (imin == i3 && imax == i2))) {
	    nabor[j] = 0;
	    quit = 1;
	  }
	  else
	    nabor[j] = temp[j];
	}
	else
	  nabor[j] = temp[j];
	fprintf(fp, " %2d", nabor[j]);     /* print adjacent node */
      }
      if (nabor[j-1] != 0) {   /* this is a constraint curve node */
	nabor[j++] = 0;        /* make sure that its list ends with a 0 */
	fprintf(fp, " %2d", nabor[j-1]);   /* print adjacent node */
      }
      fprintf(fp, "\n");

      free_iarray1(nabor);
      free_iarray1(temp);

      nEnd += j;
      (*iEnd)[i] = nEnd;    /* keep track of size of adjacency list array */
    }
  }

  *iAdj = int_array1(nEnd+1);
  fclose(fp);

  fp = fopen(tmpfil, "r");                 /* save adjacency lists */
  for (i=0; i<=nEnd; i++) {
    fscanf(fp, "%d", &(*iAdj)[i]);
    (*iAdj)[i] -= 1;
  }
  fclose(fp);
  unlink(tmpfil);
  free(tmpfil);
}

/* C version of the Fortran subroutine TRPRNT implemented by Cline and
 * Renka in their TRIPACK package
 */

void trprnt(int ncc, int *lcc, int n, float *x, float *y, int *list,
	    int *lptr, int *lend)
{
  int i, j, k,     /* loop counters */
      lp,          /* index of next node in adjacency list */
      lpl,         /* index of last node in adjacency list */
      na,          /* number of arcs */
      nabor[30],   /* list of adjacent nodes */
      nb = 0,      /* number of boundary nodes */
      nd,          /* node index */
      nt;          /* number of triangles */

  printf("                           ADJACENCY SETS,    N = %5d\n", n);
  printf(" NODE     X(NODE)        Y(NODE)");
  printf("                    NEIGHBORS OF NODE\n");

  for (i=0; i<n; i++) {
    lp = lpl = lend[i];   /* index of last node in adjacency list */

    k = 0;                /* size of this adjacency list */
    do {                  /* build adjacency list */
      lp = lptr[lp-1];    /* pointer to index of next node in list */
      nd = list[lp-1];    /* index of next node in list */
      nabor[k++] = nd;
    } while (lp != lpl);

    if (nd <= 0) {        /* this flags the boundary nodes */
      nabor[k-1] = -nd;
      nabor[k++] = 0;
      nb++;
    }

    /* print x, y, and adjacency list indices */

    printf(" %4d  %13.6e  %13.6e     ", i+1, x[i], y[i]);
    for (j=0; j<k; j++)
      printf("%5d", nabor[j]);
    printf("\n");
  }

  nt = 2*n - nb - 2;
  na = nt + n - 1;
  printf(" NB = %4d BOUNARY NODES     NA = %5d ARCS     NT = %5d TRIANGLES\n",
	 nb, na, nt);
  printf(" NCC = %3d CONSTRAINT CURVES\n", ncc);
  printf("          ");
  for (i=0; i< (ncc < 14 ? ncc : 14); i++)
    printf("%5d", lcc[i]);
  printf("\n");
}

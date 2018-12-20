/* Copyright (C) 2000 Massachusetts Institute of Technology, Cambridge, MA
   All rights reserved
*/

/* delaunay_to_edge.c */

/* DelaunayToEdgeU()
   DelaunayToEdgeV()
*/

#include <stdlib.h>
#include "gen.h"
#include "editor.h"

void DelaunayToEdgeU(int nPts_o, double **pts_o, int *nPts, double ***pts,
                    int **iEnd, int *nAdj, int **iAdj)
{
  
/*******************************************************************
  The left half of the measured data points (0<u<0.5) are wrapped
  around to (1.0<u<1.5). Then a new facet surface is generated
  by using the data with 0.5<=u<1.5
 *******************************************************************/
 
  int nPts_c,*iEnd_c,nAdj_c,*iAdj_c; 
  double **pts_c;  /* The coordinates of points after wrapped around*/ 
  double **pts_a;  /* The coordinates of points to be added*/
  int ntedge,ne; 
  int i,j,k,m;
  int **edge;  /*  the array for the end points of the edges  */
  double u1,v1,d1,u2,v2,d2;
  int n,n1,n2; 

  nPts_c = nPts_o;
  iEnd_c = int_array1(nPts_c);
  pts_c = dbl_array2(nPts_c, 3);   
  for (j=0;j<nPts_c;j++) {
    if (pts_o[j][0]<0.5) pts_c[j][0]=pts_o[j][0]+1.0;
    else 
      pts_c[j][0]=pts_o[j][0];
    pts_c[j][1]=pts_o[j][1];
    pts_c[j][2]=pts_o[j][2];
  }
  
  k = Delaunay(nPts_c, pts_c, iEnd_c, &nAdj_c, &iAdj_c, 0);  

/*******************************************************************
 Compute all the edges which intersects with the plane u=1.0
*******************************************************************/
  ntedge=0;
  for (i=0;i<nPts_c;i++) {   /* search for all node points */
    k=iEnd_c[i];
    if (i==0) j=0;
    else
      j=iEnd_c[i-1]+1; 
    for (m=j;m<=k;m++) {
      n = iAdj_c[m];
      if ((n!=-1) && ((pts_c[i][0]-1.0)*(pts_c[n][0]-1.0)<0.0)) {
	ntedge++;
      }
    }
  }

  ntedge=ntedge/2;     /** Each edge is repeated once  **/
  
  edge = int_array2(ntedge,2);

  ne = 0;
  for (i=0;i<nPts_c;i++) {   /* search for all node points  */
    k=iEnd_c[i];
    if (i==0) j=0;
    else
      j=iEnd_c[i-1]+1;
    for (m=j;m<=k;m++) {
      n = iAdj_c[m];
      if ((n!=-1) && ((pts_c[i][0]-1.0)*(pts_c[n][0]-1.0)<0.0)) {
	if (InBdryEdge(i,n,edge,ntedge)==0) {
	  edge[ne][0]=i;
	  edge[ne][1]=n;
          ne++;
	}
      }
    }
  }

/**************************************************************************
 Compute the intersection points of the edges with the plane u=1.0 
**************************************************************************/
  pts_a = dbl_array2(ntedge, 3); 
  for (i=0;i<ntedge;i++) {
    n1 = edge[i][0];
    n2 = edge[i][1];
    u1 = pts_c[n1][0]; v1 = pts_c[n1][1]; d1 = pts_c[n1][2];
    u2 = pts_c[n2][0]; v2 = pts_c[n2][1]; d2 = pts_c[n2][2]; 
    pts_a[i][0] = 1.0;
    pts_a[i][1] = v2 - (v2-v1)*(u2-1.0)/(u2-u1);
    pts_a[i][2] = d2 - (d2-d1)*(u2-1.0)/(u2-u1);
  }
 
/**************************************************************************
  Add the intersection points to original point (before u shift) lists 
  and redo the triangulation by using all the points
**************************************************************************/
  *nPts = nPts_o + 2*ntedge;
  *iEnd = int_array1(*nPts);
  *pts = dbl_array2(*nPts, 3);  
  for (j=0;j<nPts_o;j++) {
    for (k=0;k<3;k++) 
      (*pts)[j+ntedge][k]=pts_o[j][k];
  }

  for (j=0;j<ntedge;j++) {
      (*pts)[j][0]=pts_a[j][0]-1.0;
      (*pts)[j][1]=pts_a[j][1];
      (*pts)[j][2]=pts_a[j][2];
      for (k=0;k<3;k++) 
	(*pts)[j+ntedge+nPts_o][k]=pts_a[j][k];
  }

  k = Delaunay(*nPts, *pts, *iEnd, nAdj, iAdj, 0);  
}


void DelaunayToEdgeV(int nPts_o, double **pts_o, int *nPts, double ***pts,
                    int **iEnd, int *nAdj, int **iAdj)
{
  
/*******************************************************************
  The lower half of the measured data points (0<v<0.5) are wrapped
  around to (1.0<v<1.5). Then a new facet surface is generated
  by using the data with 0.5<=v<1.5
 *******************************************************************/
 
  int nPts_c,*iEnd_c,nAdj_c,*iAdj_c; 
  double **pts_c;  /* The coordinates of points after wrapped around*/ 
  double **pts_a;  /* The coordinates of points to be added*/
  int ntedge,ne; 
  int i,j,k,m;
  int **edge;  /*  the array for the end points of the edges  */
  double u1,v1,d1,u2,v2,d2;
  int n,n1,n2; 

  nPts_c = nPts_o;
  iEnd_c = int_array1(nPts_c);
  pts_c = dbl_array2(nPts_c, 3);   
  for (j=0;j<nPts_c;j++) {
    if (pts_o[j][1]<0.5) pts_c[j][1]=pts_o[j][1]+1.0;
    else 
      pts_c[j][1]=pts_o[j][1];
    pts_c[j][0]=pts_o[j][0];
    pts_c[j][2]=pts_o[j][2];
  }
  
  k = Delaunay(nPts_c, pts_c, iEnd_c, &nAdj_c, &iAdj_c, 0);  

/*******************************************************************
 Compute all the edges which intersects with the plane u=1.0
*******************************************************************/
  ntedge=0;
  for (i=0;i<nPts_c;i++) {   /* search for all node points */
    k=iEnd_c[i];
    if (i==0) j=0;
    else
      j=iEnd_c[i-1]+1; 
    for (m=j;m<=k;m++) {
      n = iAdj_c[m];
      if ((n!=-1) && ((pts_c[i][1]-1.0)*(pts_c[n][1]-1.0)<0.0)) {
	ntedge++;
      }
    }
  }

  ntedge=ntedge/2;     /** Each edge is repeated once  **/
  
  edge = int_array2(ntedge,2);

  ne = 0;
  for (i=0;i<nPts_c;i++) {   /* search for all node points  */
    k=iEnd_c[i];
    if (i==0) j=0;
    else
      j=iEnd_c[i-1]+1;
    for (m=j;m<=k;m++) {
      n = iAdj_c[m];
      if ((n!=-1) && ((pts_c[i][1]-1.0)*(pts_c[n][1]-1.0)<0.0)) {
	if (InBdryEdge(i,n,edge,ntedge)==0) {
	  edge[ne][0]=i;
	  edge[ne][1]=n;
          ne++;
	}
      }
    }
  }

/**************************************************************************
 Compute the intersection points of the edges with the plane v=1.0 
**************************************************************************/
  pts_a = dbl_array2(ntedge, 3); 
  for (i=0;i<ntedge;i++) {
    n1 = edge[i][0];
    n2 = edge[i][1];
    u1 = pts_c[n1][0]; v1 = pts_c[n1][1]; d1 = pts_c[n1][2];
    u2 = pts_c[n2][0]; v2 = pts_c[n2][1]; d2 = pts_c[n2][2]; 
    pts_a[i][1] = 1.0;
    pts_a[i][0] = u2 - (u2-u1)*(v2-1.0)/(v2-v1);
    pts_a[i][2] = d2 - (d2-d1)*(v2-1.0)/(v2-v1);

  }
 
/**************************************************************************
  Add the intersection points to original point (before u shift) lists 
  and redo the triangulation by using all the points
**************************************************************************/
  *nPts = nPts_o + 2*ntedge;
  *iEnd = int_array1(*nPts);
  *pts = dbl_array2(*nPts, 3);  
  for (j=0;j<nPts_o;j++) {
    for (k=0;k<3;k++) 
      (*pts)[j+ntedge][k]=pts_o[j][k];
  }

  for (j=0;j<ntedge;j++) {
      (*pts)[j][0]=pts_a[j][0];
      (*pts)[j][1]=pts_a[j][1]-1.0;
      (*pts)[j][2]=pts_a[j][2];
      for (k=0;k<3;k++) 
	(*pts)[j+ntedge+nPts_o][k]=pts_a[j][k];
  }

  k = Delaunay(*nPts, *pts, *iEnd, nAdj, iAdj, 0);  
}

int InBdryEdge(int i, int j, int **edge, int ntedge)
{
  int k,flag;
  k = 0;
  flag = 0;
  
  while ((flag==0) && (k<ntedge)) {
    if (((i==edge[k][0])&&(j==edge[k][1]))||((i==edge[k][1])&&(j==edge[k][0]))) {
      flag = 1;
    }
    else 
      k++;
  }

  return flag;

}

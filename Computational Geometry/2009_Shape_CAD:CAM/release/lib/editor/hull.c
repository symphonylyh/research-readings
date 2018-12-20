/* Copyright (C) 1994 Massachusetts Institute of Technology, Cambridge, MA
   All rights reserved
*/

/* hull.c */

/* AddAdj()
   ComparAdj()
   NonConvexEdge()
*/

#include <math.h>
#include "editor.h"

/* adjs[i][0] and adjs[i][1] are the vertices (v1,v2) of edge i */
/* adjs[i][2] is the third vertex (v3) for triangle (v1,v2,v3)  */
/* adjs[i][3] is the third vertex (V3) for triangle (v2,v1,V3)  */

void AddAdj(int **adjs, int v1, int v2, int v3, int *nadj)
{
  int found, i;

/* although each edge is specified twice, as (v1,v2) and (v2,v1), we only  */
/* want each edge in the list once, so if edge (v1,v2) is already in list */
/* as (v2,v1), add v3 */

  found = 0;
  for (i=0; i<*nadj; i++)
    if (v1 == adjs[i][1] && v2 == adjs[i][0]) {
      adjs[i][3] = v3;
      found = 1;
    }

/* else add edge (v1,v2) to list */

  if (!found) {
    adjs[*nadj][0] = v1;
    adjs[*nadj][1] = v2;
    adjs[*nadj][2] = v3;
    *nadj += 1;
  }
}

/* sort adjacency list */

int ComparAdj(const void *c1, const void *c2)
{
  int *adj1, *adj2;
  int i;

  adj1 = (int *)c1;
  adj2 = (int *)c2;
  for (i=0; i<4; i++)
    if (adj1[i] < adj2[i])
      return -1;
    else if (adj1[i] > adj2[i])
      return  1;
  return 0;
}

/* determine if edge is non-convex */

int NonConvexEdge(double **pts, int *adjs)
{
  double com, cph, cth, e1[3], e2[3], e3[3], len, n[3], som, sph, sth;

/* the triangles adjacent to edge (v1,v2) = (adjs[0],adjs[1]) are: */
/* (v1,v2,v3) = (adjs[0],adjs[1],adjs[2]) and */
/* (v2,v1,v4) = (adjs[1],adjs[0],adjs[3])     */

  e1[0] = pts[adjs[0]][0] - pts[adjs[1]][0];   /* edge (v2,v1)  */
  e1[1] = pts[adjs[0]][1] - pts[adjs[1]][1];   /* of (v2,v1,v4) */
  e1[2] = pts[adjs[0]][2] - pts[adjs[1]][2];

  e2[0] = pts[adjs[3]][0] - pts[adjs[1]][0];   /* edge (v2,v4)  */
  e2[1] = pts[adjs[3]][1] - pts[adjs[1]][1];   /* of (v2,v1,v4) */
  e2[2] = pts[adjs[3]][2] - pts[adjs[1]][2];

  n[0] = e1[1]*e2[2] - e1[2]*e2[1];   /* normal of (v2,v1,v4 */
  n[1] = e1[2]*e2[0] - e1[0]*e2[2];   /* is cross product */
  n[2] = e1[0]*e2[1] - e1[1]*e2[0];

  len = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
  if (len) {
    n[0] /= len;                      /* normalize to unit vector */
    n[1] /= len;
    n[2] /= len;
  }

/* transform edge (v1,v2) so that v1 = (0,0,0) and v2 is (0,0,z) */
/*                                and v3 is on xy-plane */
/* 1. translate v1 = 0 and v2 to v2-v1 and v3 to v3-v1 */

  e1[0] = pts[adjs[1]][0] - pts[adjs[0]][0];   /* e1 = v2-v1 */
  e1[1] = pts[adjs[1]][1] - pts[adjs[0]][1];
  e1[2] = pts[adjs[1]][2] - pts[adjs[0]][2];

  e3[0] = pts[adjs[2]][0] - pts[adjs[0]][0];   /* e3 = v3-v1 */
  e3[1] = pts[adjs[2]][1] - pts[adjs[0]][1];
  e3[2] = pts[adjs[2]][2] - pts[adjs[0]][2];

/* 2. rotate v2 about y-axis so that v2 is in yz-plane */

  len = sqrt(e1[0]*e1[0] + e1[2]*e1[2]);
  if (len) {
    cth = e1[2]/len;   /* cos(theta), where theta is rotation about y-axis */
    sth = e1[0]/len;   /* sin(theta), where theta is rotation about y-axis */
  }
  else {
    cth = 1.0;
    sth = 0.0;
  }

  e2[0] = 0.0;       /* e2 is e1 rotated about y-axis by -theta */
  e2[1] = e1[1];
  e2[2] = len;

  e1[0] = e3[0]*cth - e3[2]*sth;   /* e1 is e3 rotated about y-axis */
  e1[1] = e3[1];
  e1[2] = e3[0]*sth + e3[2]*cth;

/* 3. rotate v2 about x-axis so that v2 is in xz-plane */

  len = sqrt(e2[1]*e2[1] + e2[2]*e2[2]);
  if (len) {
    cph =  e2[2]/len;  /* cos(phi), where phi is rotation about x-axis */
    sph = -e2[1]/len;  /* sin(phi), where phi is rotation about x-axis */
  }
  else {
    cph = 1.0;
    sph = 0.0;
  }
  e3[0] =  e1[0];                  /* e3 is e1 rotated about x-axis */
  e3[1] =  e1[1]*cph + e1[2]*sph;
  e3[2] = -e1[1]*sph + e1[2]*cth;

/* 4. rotate e3 about z-axis to put e3 on xz-plane */

  len = sqrt(e3[0]*e3[0] + e3[1]*e3[1]);
  if (len) {
    com = e3[0]/len;  /* cos(omega), where omega is rotation about z-axis */
    som = e3[1]/len;  /* sin(omega), where omega is rotation about z-axis */
  }
  else {
    com = 1.0;
    som = 0.0;
  }

/* triangle (v1,v2,v3) is now in xz-plane with normal parallel to y-axis */
/* rotate normal about y by -theta and about x by -phi */
/* recall that cos(-theta) = cos(theta) and sin(-theta) = -sin(theta) */

  e1[0] =  n[0]*cth - n[2]*sth;
  e1[1] =  n[1];
  e1[2] =  n[0]*sth + n[2]*cth;

  e2[0] =  e1[0];
  e2[1] =  e1[1]*cph + e1[2]*sph;
  e2[2] = -e1[1]*sph + e1[2]*cph;

  n[0]  =  e2[0]*com + e2[1]*som;
  n[1]  = -e2[0]*som + e2[1]*com;
  n[2]  =  e2[2];

  if (n[0] >= 0.0)   /* if normal is parallel to or intersects  +y-axis */
    return 0;        /* edge is non-convex */

  return 1;          /* else edge is convex */
}

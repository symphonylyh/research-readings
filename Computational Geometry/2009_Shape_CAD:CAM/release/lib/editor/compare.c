/* Copyright (C) 1997 Massachusetts Institute of Technology, Cambridge, MA
   All rights reserved
*/

/* compare.c */

/* check_cv()
   check_one_dir()
   CheckOneDir()
   check_patch_u()
   CheckPatchContinuity()
   CheckPatchU()
   com_knot_vector()
   compare_curves()
   compare_patches()
   CompatePatches()
   GetCurvTang()
   patch_with_com_knots()
   RaiseAndMerge()
   SampleContinuity()
*/

#include <math.h>
#include "gen.h"
#include "bspl.h"
#include "editor.h"

void compare_curves(ParCurv **egeom, int np, int *POS_FLAG, double *c0err,
		    int *g1, int *g2)
{
  ParCurv *egeom1, *egeom2;
  int or, nc;
  short i;

  for (i=0; i<np-1; i++)
    if(egeom[i]->order != egeom[i]->ncontpts){
      or = egeom[i]->order;
      nc = egeom[i]->ncontpts+egeom[i]->order-1;
      egeom1 = egeomalloc1(or,  nc);
      or = egeom[i+1]->order;
      nc = egeom[i+1]->ncontpts+egeom[i+1]->order-1;
      egeom2 = egeomalloc1(or, nc);
      subdiv_cv(egeom[i], egeom1, 2);
      subdiv_cv(egeom[i+1], egeom2, 1);
      check_cv(egeom1, egeom2, POS_FLAG, &c0err[i], &g1[i], &g2[i]);
      free_egeom(egeom1);
      free_egeom(egeom2);
    }
    else
      check_cv(egeom[i], egeom[i+1], POS_FLAG, &c0err[i], &g1[i], &g2[i]);
}

void compare_patches(ParSurf ***fgeom, int nps, int npt, int *UPOS_FLAG,
		     int *VPOS_FLAG, double *c0u, double *c0v, int *g1u,
		     int *g1v, int *g2u, int *g2v)
{
  ParSurf ***nf;
  short i, j;

  check_one_dir(fgeom, nps, npt, UPOS_FLAG, VPOS_FLAG, c0u, g1u, g2u);

  nf = (ParSurf ***)gen_array2(nps, npt, sizeof(ParSurf *));
  for(i=0; i<npt; i++)
    for(j=0; j<nps; j++)
      nf[j][i] = swapsurf(fgeom[i][j]);

  check_one_dir(nf, npt, nps, UPOS_FLAG, VPOS_FLAG, c0v, g1v, g2v);
  
  free_fgeom_array2(nf, npt, nps);
}

void check_one_dir(ParSurf ***fgeom, int nps, int npt, int *UPOS_FLAG,
		   int *VPOS_FLAG, double *c0, int *g1, int *g2)
{
  ParSurf *fgeom1, *fgeom2, *fg1, *fg2;
  double *cknots;
  short i, j, k, nk;
  
  for (j=k=0; j<npt; j++)
    for (i=0; i<nps-1; i++,k++) {
      if (fgeom[j][i]->ucontpts != fgeom[j][i]->uorder) {
	/* make sure that cknots s large enough to hold the largest */
	/* possible merged v knot vector */
	cknots = dbl_array1(fgeom[j][ i ]->vorder + fgeom[j][ i ]->vcontpts +
			    fgeom[j][i+1]->vorder + fgeom[j][i+1]->vcontpts);
	nk = com_knot_vector(fgeom[j][i], fgeom[j][i+1], cknots);
	fg1 = patch_with_com_knots(fgeom[j][i], cknots, nk);
	fg2 = patch_with_com_knots(fgeom[j][i+1], cknots, nk);
	free_darray1(cknots);
	fgeom1 = fgeomalloc1(fg1->uorder, fg1->vorder, fg1->ucontpts+
			     fg1->uorder-1, fg1->vcontpts);
	fgeom2 = fgeomalloc1(fg2->uorder, fg2->vorder, fg2->ucontpts+
			     fg2->uorder-1, fg2->vcontpts);
	
	subdiv_bez(fg1, fgeom1, 2, 0);
	subdiv_bez(fg2, fgeom2, 1, 0);
	check_patch_u(fgeom1, fgeom2, UPOS_FLAG, VPOS_FLAG, &c0[k], &g1[k],
		      &g2[k]);
	free_fgeom(fgeom1);
	free_fgeom(fgeom2);
	free_fgeom(fg1);
	free_fgeom(fg2);
      }
      else
	check_patch_u(fgeom[j][i], fgeom[j][i+1], UPOS_FLAG, VPOS_FLAG,
		      &c0[k], &g1[k], &g2[k]);
    }
}

void check_patch_u(ParSurf *fgeom1, ParSurf *fgeom2, int *UPOS_FLAG,
		   int *VPOS_FLAG, double *c0, int *g1, int *g2)
{
  vector *T1, *T2, *K1, *K2, *tmp, *norm1, *norm2;
  double maxErr=0.0, dis;
  int vc, uc;
  short i, j;

  T1=vectalloc();
  T2=vectalloc();
  K1=vectalloc();
  K2=vectalloc();
  tmp = vectalloc();

  vc = fgeom1->vcontpts;
  uc = fgeom1->ucontpts;
  *UPOS_FLAG = 1;
  for (i=0; i<vc; i++) {
    dis = mag(sub_vect(fgeom1->contpts[uc-1][i], fgeom2->contpts[0][i]));
    if (dis > maxErr)
      maxErr = dis;
    if (mag(sub_vect(fgeom1->contpts[uc-1][i], fgeom2->contpts[0][i])) >
       EDITOR_JBSPL_TOL) {
      *UPOS_FLAG = 0;
      break;
    }
  }
  *c0 = maxErr;
  
  T1 = revalderivsurf(fgeom1, 1.0, 0.0, 1, 0);
  T2 = revalderivsurf(fgeom2, 0.0, 0.0, 1, 0);
  vect_colinear(T1, T2);
  *g1 = 1;
  for (i=0; i<vc; i++) {
    sub_vect1(fgeom1->contpts[uc-2][i], fgeom1->contpts[uc-1][i], T1);
    sub_vect1(fgeom2->contpts[1][i], fgeom2->contpts[0][i], T2);
    sub_vect1(fgeom1->contpts[uc-2][i],fgeom2->contpts[1][i], tmp);
    vect_colinear(T1,  tmp);
    if (vect_colinear(T1, T2) == 0) {
      *g1 = 0;
      break;
    }
  }

  for (i=0; i<11; i++) {
    norm1 = normalsurf(fgeom1, 1.0, 0.1*i);
    norm2 = normalsurf(fgeom2, 0.0, 0.1*i);
    vect_colinear(norm1, norm2);
    vectfree(norm1);
    vectfree(norm2);
  }
  T1 = revalderivsurf(fgeom1, 1.0, 0.0, 2, 0);
  T2 = revalderivsurf(fgeom2, 0.0, 0.0, 2, 0);
  sub_vect1(T2, T1, K1);
  sub_vect1(fgeom2->contpts[1][i], fgeom2->contpts[0][i], K2);
  vect_colinear(K1, K2);
  *g2 = 1;
  for (i=0; i<vc; i++) {
     sub_vect1(fgeom2->contpts[2][i], fgeom2->contpts[1][i], K2);
     sub_vect1(K2, fgeom2->contpts[1][i], K2);
     add_vect1(K2, fgeom2->contpts[0][i], K2);

     sub_vect1(fgeom1->contpts[uc-3][i], fgeom1->contpts[uc-2][i], K1);
     sub_vect1(K1, fgeom1->contpts[uc-2][i], K1);
     add_vect1(K1, fgeom1->contpts[uc-1][i], K1);

     sub_vect1(K2, K1, T1);
     sub_vect1(fgeom1->contpts[uc-2][i], fgeom1->contpts[uc-1][i], T2);

     if (vect_colinear(T1, T2) == 0) {
       *g2 = 0;
       break;
     }
  }

  vectfree(T1);
  vectfree(T2);
  vectfree(K1);
  vectfree(K2);
}

int com_knot_vector(ParSurf *fgeom1, ParSurf *fgeom2, double cknots[])
{
  short i, n, k1, k2;

  for (i=0; i<fgeom1->vorder; i++)
    cknots[i] = 0.0;
  n=fgeom1->vorder;
  k1=fgeom1->vorder;
  k2=fgeom2->vorder;
  while (1) {
    if (k1 == fgeom1->vcontpts && k2 == fgeom2->vcontpts)
      break;
    if (fgeom1->vknots[k1] < fgeom2->vknots[k2]){
      cknots[n] = fgeom1->vknots[k1];
      k1++; n++;
    }
    else if (fgeom1->vknots[k1] == fgeom2->vknots[k2]) {
      cknots[n] = fgeom1->vknots[k1];
      k1++; k2++; n++;
    }
    else{
      cknots[n] = fgeom2->vknots[k2];
      k2++; n++;
    }
  }
  for (i=n; i<fgeom1->vorder+n; i++)
    cknots[i] = 1.0;

  return n+fgeom1->vorder;
}

ParSurf *patch_with_com_knots(ParSurf *fgeom, double cknots[], int nk)
{
  ParSurf *nf;
  short i;

  nf = fgeomalloc1(fgeom->uorder, fgeom->vorder, fgeom->ucontpts, 
		   nk-fgeom->vorder);
		   
  for (i=0; i<nf->uorder+nf->ucontpts; i++)
    nf->uknots[i]=fgeom->uknots[i];
  for (i=0; i<nk; i++)
    nf->vknots[i] = cknots[i];

  surfoslo3(fgeom, nf, 4);

  return nf;
}

void check_cv(ParCurv *egeom1, ParCurv *egeom2, int *POS_FLAG, double *c0err,
	      int *g1, int *g2)
{
  vector *T1, *T2, *K1, *K2;
  double maxErr = 0.0, dis;
  int uc;

  T1 = vectalloc();
  T2 = vectalloc();
  K1 = vectalloc();
  K2 = vectalloc();

  uc = egeom1->ncontpts;
  *POS_FLAG = 1;
  
  dis = mag(sub_vect(egeom1->contpts[uc-1], egeom2->contpts[0]));
  if (dis > maxErr)
    maxErr = dis;
  if (mag(sub_vect(egeom1->contpts[uc-1], egeom2->contpts[0])) >
      EDITOR_JBSPL_TOL)
    *POS_FLAG = 0;

  *c0err = maxErr;
  
  T1 = GetCurvTang(egeom1, 1.0, 1);
  T2 = GetCurvTang(egeom2, 0.0, 1);

  if (vect_colinear(T1, T2) == 0)
    *g1 = 0;
  else
    *g1 = 1;

  sub_vect1(egeom2->contpts[2], egeom2->contpts[1], K2);
  sub_vect1(K2, egeom2->contpts[1], K2);
  add_vect1(K2, egeom2->contpts[0], K2);
  
  sub_vect1(egeom1->contpts[uc-3], egeom1->contpts[uc-2], K1);
  sub_vect1(K1, egeom1->contpts[uc-2], K1);
  add_vect1(K1, egeom1->contpts[uc-1], K1);
  
  sub_vect1(K2, K1, T1);
  sub_vect1(egeom1->contpts[uc-2], egeom1->contpts[uc-1], T2);

  if (vect_colinear(T1, T2) == 0)
    *g2 = 0;
  else
    *g2 = 1;
  
  vectfree(T1);
  vectfree(T2);
  vectfree(K1);
  vectfree(K2);
}

vector *GetCurvTang(ParCurv *egeom, float t, int is_open)
{
  vector *tang;

  if (is_open)
    tang = rbspeval(egeom, t, 1);
  else
    tang = rbspeval_per(egeom, t, 1);

  return (tang);
}

void ComparePatches(ParSurf ***fgeoms, int nps, int npt, int nsegs,
		    double *c0u, double *c0ut, double *c0v, double *c0vt,
		    double *c1u, double *c1ut, double *c1v, double *c1vt,
		    double *c2u, double *c2ut, double *c2v, double *c2vt)
{
  ParSurf ***temps;
  int i, j;

 /* check continuity in U */
  CheckOneDir(fgeoms, nps, npt, nsegs, c0u, c0ut, c1u, c1ut, c2u, c2ut);

  temps = (ParSurf ***)gen_array2(nps, npt, sizeof(ParSurf *));
  for (i=0; i<npt; i++)                /* temporarily swap U and V */
    for (j=0; j<nps; j++)
      temps[j][i] = swapsurf(fgeoms[i][j]);

  /* check continuity in V */
  CheckOneDir(temps, npt, nps, nsegs, c0v, c0vt, c1v, c1vt, c2v, c2vt);

  free_fgeom_array2(temps, npt, nps);   /* free temporary copy */
}

void CheckOneDir(ParSurf ***fgeoms, int nps, int npt, int nsegs, double *c0u,
		 double *c0t, double *c1u, double *c1t, double *c2u,
		 double *c2t)
{
  ParSurf *fgeom1, *fgeom2;
  int i, j;

  /* check continuity of patches adjacent in the U direction */

  for (j=0; j<npt; j++)
    for (i=1; i<nps; i++) {       
      RaiseAndMerge(fgeoms[j][i-1], fgeoms[j][i], &fgeom1, &fgeom2);
      CheckPatchU(fgeom1, fgeom2, nsegs, &c0u[i-1], &c0t[i-1], &c1u[i-1],
		  &c1t[i-1], &c2u[i-1], &c2t[i-1]);
      free_fgeom(fgeom1);
      free_fgeom(fgeom2);
    }
}

void RaiseAndMerge(ParSurf *fgeom1, ParSurf *fgeom2, ParSurf **fgm1,
                   ParSurf **fgm2)
{
  ParSurf *tmp1, *tmp2;
  double eps = 1.0e-12, *knots;
  int k, n;

  if (fgeom1->vorder != fgeom2->vorder) {  /* V orders don't match */
    if (fgeom1->vorder > fgeom2->vorder) { /* raise V order of 2nd surface */
      tmp1 = copyfgeom(fgeom1, NULL);      /* 1st surface is ok */

      n = fgeom1->vorder   - fgeom2->vorder;
      k = fgeom2->vcontpts - fgeom1->vorder + 1;
      if (k < 1) k = 1;
      tmp2 = fgeomalloc1(fgeom2->uorder, fgeom2->vorder+n, fgeom2->ucontpts,
                          fgeom2->vcontpts + fgeom1->vorder*k);
      copyfgeom(fgeom2, tmp2);
      raise_surf(tmp2, eps, 0, n);
    }
    else {                                 /* raise V order of 1st surface */
      n = fgeom2->vorder   - fgeom1->vorder;
      k = fgeom1->vcontpts - fgeom2->vorder + 1;
      if (k < 1) k = 1;
      tmp1 = fgeomalloc1(fgeom1->uorder, fgeom1->vorder+n, fgeom1->ucontpts,
                          fgeom1->vcontpts + fgeom2->vorder*k);
      copyfgeom(fgeom1, tmp1);
      raise_surf(tmp1, eps, 0, n);

      tmp2 = copyfgeom(fgeom2, NULL);      /* 2nd surface is ok */
    }
  }
  else {
    tmp1 = fgeom1;
    tmp2 = fgeom2;
  }

  /* if both surfaces are not Bezier patches, merge their knot vectors */
  if (tmp1->vorder != tmp1->vcontpts || tmp2->vorder != tmp2->vcontpts) {
    knots = dbl_array1(2*tmp1->vorder + tmp1->vcontpts + tmp2->vcontpts);
    n = com_knot_vector(tmp1, tmp2, knots);
    *fgm1 = patch_with_com_knots(tmp1, knots, n);
    *fgm2 = patch_with_com_knots(tmp2, knots, n);
    free_darray1(knots);
  }
  else {
    *fgm1 = tmp1;
    *fgm2 = tmp2;
  }
}

void CheckPatchU(ParSurf *fgeom1, ParSurf *fgeom2, int nsegs, double *c0u,
		 double *c0t, double *c1u, double *c1t, double *c2u,
		 double *c2t)
{
  double curv, d, dv, r1[7], r2[7], rad, v, v0;
  int i, j, nspans;

  /* note: fgeom1 and fgeom2 have same order and knot vector in V */

  nspans = fgeom1->vcontpts - fgeom1->vorder + 1;

  *c0u = *c0t = 0.0;
  *c1u = *c1t = 0.0;
  *c2u = *c2t = 0.0;
  for (i=0; i<nspans; i++) {
    v0 = fgeom1->vknots[fgeom1->vorder + i - 1];
    dv = fgeom1->vknots[fgeom1->vorder + i] - v0;
    for (j = (i ? 1 : 0); j<nsegs; j++) {
      v = v0 + dv*(double)j/(double)(nsegs-1);

      CheckPatchContinuity(fgeom1, fgeom2, v, &d, &rad, &curv, r1, r2);

      if (fabs(d) > fabs(*c0u)) {        /* find and save extrema */
	*c0u = d;
	*c0t = v;
      }
      if (fabs(rad) > fabs(*c1u)) {
	*c1u = rad;
	*c1t = v;
      }
      if (fabs(curv) > fabs(*c2u)) {
	*c2u = curv;
	*c2t = v;
      }
    }
  }
}

void CheckPatchContinuity(ParSurf *fgeom1, ParSurf *fgeom2, double v,
			  double *d, double *rad, double *curv,
			  double *d1, double *d2)
{
  vector *norm1, *norm2, *r1[5], *r2[5];
  double cth, e, l, normu1, normu2;
  int i;

  evalrsurf(fgeom1, 1.0, v, 4, &r1[0]);
  evalrsurf(fgeom2, 0.0, v, 4, &r2[0]);

  *d = distance(r1[0], r2[0]);       /* positional differance */
  d1[0] = r1[0]->x/r1[0]->w;
  d1[1] = r1[0]->y/r1[0]->w;
  d1[2] = r1[0]->z/r1[0]->w;
  d2[0] = r2[0]->x/r2[0]->w;
  d2[1] = r2[0]->y/r2[0]->w;
  d2[2] = r2[0]->z/r2[0]->w;

  norm1 = normalsurf(fgeom1, 1.0, v);
  norm2 = normalsurf(fgeom2, 0.0, v);

  cth = dot(norm1, norm2);           /* cosine of angle between vectors */
  if (cth > 9.99999e-1)
    *rad = 0.0;
  else
    *rad = acos(cth);                /* angle (radians) between vectors */

  d1[3] = norm1->x/norm1->w;
  d1[4] = norm1->y/norm1->w;
  d1[5] = norm1->z/norm1->w;
  d2[3] = norm2->x/norm2->w;
  d2[4] = norm2->y/norm2->w;
  d2[5] = norm2->z/norm2->w;

  e = dot(r1[1], r1[1]);             /* fundamental coefficients */
  l = dot(norm1, r1[4]);
  normu1 = -l/e;                     /* normal u curvature */
  e = dot(r2[1], r2[1]);
  l = dot(norm2, r2[4]);
  normu2 = -l/e;
  *curv = fabs(normu2-normu1);       /* curvature differance */
  d1[6] = normu1;
  d2[6] = normu2;

  vectfree(norm1);
  vectfree(norm2);

  for (i=0; i<5; i++) {
    vectfree(r1[i]);
    vectfree(r2[i]);
  }
}

void SampleContinuity(ParSurf *fgm1, ParSurf *fgm2, int nsegs, double ***pts,
		      int *n, double **knots)
{
  ParSurf *fgeom1, *fgeom2;
  double curv, d, dv, r1[7], r2[7], rad, v, v0;
  int i, j, k, nspans;

  RaiseAndMerge(fgm1, fgm2, &fgeom1, &fgeom2);

  /* note: fgeom1 and fgeom2 have same order and knot vector in V */

  nspans = fgeom1->vcontpts - fgeom1->vorder + 1;

  *n = nspans*(nsegs-1) + 1;
  *pts = dbl_array2(*n, 3);
  *knots = dbl_array1(*n);

  for (i=k=0; i<nspans; i++) {
    v0 = fgeom1->vknots[fgeom1->vorder + i - 1];
    dv = fgeom1->vknots[fgeom1->vorder + i] - v0;
    for (j = (i ? 1 : 0); j<nsegs; j++,k++) {
      v = v0 + dv*(double)j/(double)(nsegs-1);
      (*knots)[k] = v;

      CheckPatchContinuity(fgeom1, fgeom2, v, &d, &rad, &curv, r1, r2);
      (*pts)[k][0] = d;
      (*pts)[k][1] = rad*RAD_TO_DEG;
      (*pts)[k][2] = curv;
    }
  }
  free_fgeom(fgeom1);
  free_fgeom(fgeom2);
}

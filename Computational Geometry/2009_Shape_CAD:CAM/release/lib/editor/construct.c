/* Copyright (C) 1998 Massachusetts Institute of Technology, Cambridge, MA
   All rights reserved
*/

/* construct.c */

/* add_curv()
   add_surf()
   construct_bspl()
   construct_jbspl()
   renew_knots()
   unify_knots()
   unify_knots_cv()
*/

#include "bspl.h"
#include "gen.h"
#include "editor.h"
#include <math.h>

ParCurv *construct_bspl(ParCurv **egeoms, int np, int POS_FLAG, int join_key)
{
  ParCurv *bspls;
  vector *tmp;
  int ncontpts=0, order, uk, uc, cpos;
  short i, j;

  order = egeoms[0]->order;

  for (i=0; i<np; i++)
    ncontpts += egeoms[i]->ncontpts;

  bspls = egeomalloc1(order, ncontpts);
 
  for (i=uc=uk=0; i<np; i++) {
    add_curv(egeoms, np, bspls, i, uk, uc);
    
    uk = uk + egeoms[i]->ncontpts + egeoms[i]->order;
    uc = uc + egeoms[i]->ncontpts;
  }

  if (POS_FLAG == 0 && join_key) {
    tmp = vectalloc();
    cpos = 0-1;  
    for (i=0; i<np-1; i++) {
      cpos += egeoms[i]->ncontpts;
      add_vect1(bspls->contpts[cpos], bspls->contpts[cpos+1], tmp);
      scale_vect1(0.5, tmp, tmp);
      bspls->contpts[cpos] = copyvector(tmp, bspls->contpts[cpos]); 
      bspls->contpts[cpos+1] = copyvector(tmp, bspls->contpts[cpos+1]); 
    }
    vectfree(tmp);
  }

  unify_knots_cv(egeoms, np, bspls);

  return bspls;
}

void add_curv(ParCurv **egeoms, int np, ParCurv *bspls, int mu, int uk, int uc)
{
  short i, j;
  
  for (i=uc; i<uc+egeoms[mu]->ncontpts; i++)
    bspls->contpts[i] = copyvector(egeoms[mu]->contpts[i-uc],
				   bspls->contpts[i]);
}

void unify_knots_cv(ParCurv **egeoms, int np, ParCurv *bspls)
{
  ParCurv *tmp_curv, *tmp;
  double s, uintval, vintval, *initial, *length;
  short i, j, k;

  length = dbl_array1(np);
  initial = dbl_array1(np);

  for (i=0; i<np; i++)
    length[i] = arc_length(egeoms[i]);
  
  uintval = 0.0;
  for (i=0; i<np; i++)
    uintval += length[i];
  
  s = 0.0;
  for (i=0; i<np; i++) {
    initial[i] = s/uintval;
    s += length[i];
  }

  for (i=k=0; i<np; i++)
    for (j=0; j<egeoms[i]->ncontpts; j++,k++)
      bspls->knots[k] = egeoms[i]->knots[j]*length[i]/uintval + initial[i];

  free_darray1(length);
  free_darray1(initial);

  for (i=0; i<egeoms[np-1]->order; i++,k++)
    bspls->knots[k] = 1.0;
}

ParSurf *construct_jbspl(ParSurf ***fgeom, int nps, int npt, int UPOS_FLAG,
			 int VPOS_FLAG, int join_key)
{
  ParSurf *bspl_surf, *surf;
  vector *tmp;
  int ucontpts=0, vcontpts=0, uorder, vorder, uc, vc, cpos;
  short i, j;

  uorder = fgeom[0][0]->uorder;
  vorder = fgeom[0][0]->vorder;

  for (i=0; i<nps; i++)
    ucontpts += fgeom[0][i]->ucontpts;
  for (i=0; i<npt; i++)
    vcontpts += fgeom[i][0]->vcontpts;

  bspl_surf = fgeomalloc1(uorder, vorder, ucontpts, vcontpts);

  for (j=vc=0; j<npt; j++) {
    for (i=uc=0; i<nps; i++) {
      add_surf(fgeom, i, j, uc, vc, bspl_surf);
      uc += fgeom[j][i]->ucontpts;
    }
    vc += fgeom[j][0]->vcontpts;
  }

  if (UPOS_FLAG == 0 && join_key) {
    tmp = vectalloc();
    cpos = 0-1;  
    for (i=0; i<nps-1; i++) {
      cpos += fgeom[0][i]->ucontpts;
      for (j=0; j<bspl_surf->vcontpts; j++) {
	add_vect1(bspl_surf->contpts[cpos][j],bspl_surf->contpts[cpos+1][j],
		  tmp);
	scale_vect1(0.5, tmp, tmp);
	bspl_surf->contpts[cpos][j] = 
	  copyvector(tmp, bspl_surf->contpts[cpos][j]); 
	bspl_surf->contpts[cpos+1][j] = 
	  copyvector(tmp, bspl_surf->contpts[cpos+1][j]); 
      }
    }
    vectfree(tmp);
  }
  if (VPOS_FLAG == 0) {
    tmp = vectalloc();
    cpos = 0-1;
    for (i=0; i<npt-1; i++) {
      cpos += fgeom[i][0]->vcontpts;
      for (j=0; j<bspl_surf->ucontpts; j++) {
	add_vect1(bspl_surf->contpts[j][cpos],bspl_surf->contpts[j][cpos+1],
		  tmp);
	scale_vect1(0.5, tmp, tmp);
	bspl_surf->contpts[j][cpos] = 
	  copyvector(tmp, bspl_surf->contpts[j][cpos]); 
	bspl_surf->contpts[j][cpos+1] = 
	  copyvector(tmp, bspl_surf->contpts[j][cpos+1]); 
      }
    }
    vectfree(tmp);
  }

  return bspl_surf;
}

void add_surf(ParSurf ***fgeom, int mu, int nv, int uc, int vc, ParSurf *bspls)
{
  short i, j;
  
  for (i=uc; i<uc+fgeom[nv][mu]->ucontpts; i++)
    for (j=vc; j<vc+fgeom[nv][mu]->vcontpts; j++)
      bspls->contpts[i][j] = copyvector(fgeom[nv][mu]->contpts[i-uc][j-vc],
					bspls->contpts[i][j]);
}

void unify_knots(ParSurf ***fgeom, int nps, int npt, ParSurf *bspls)
{
  ParSurf *tmp_surf;
  double *initial, *length, s, uintval, vintval;
  short i, j, k;

  /* normalize the knot vectors to the interval [0,1] */

  length = dbl_array1(nps);
  for(i=0; i<nps; i++) {
    length[i] = 0.0;
    for(j=0; j<npt; j++)
      length[i] += avg_arc_length(fgeom[j][i]);
  }
  
  uintval = 0.0;
  for(i=0; i<nps; i++)
    uintval += length[i];
  
  initial = dbl_array1(nps);
  s = 0.0;
  for(i=0; i<nps; i++){
    initial[i] = s/uintval;
    s += length[i];
  }

  for (i=k=0; i<nps; i++)
    for (j=0; j<fgeom[0][i]->ucontpts; j++,k++)
      bspls->uknots[k] = fgeom[0][i]->uknots[j]*length[i]/uintval + initial[i];

  free_darray1(length);
  free_darray1(initial);

  for (i=0; i<fgeom[0][nps-1]->uorder; i++,k++)
    bspls->uknots[k] = 1.0;

  length = dbl_array1(npt);
  for (i=0; i<npt; i++) {
    length[i] = 0.0;
    for (j=0; j<nps; j++) {
      tmp_surf = swapsurf(fgeom[i][j]);
      length[i] += avg_arc_length(tmp_surf);
      free_fgeom(tmp_surf);
    }
  }
  
  vintval = 0.0;
  for (i=0; i<npt; i++)
    vintval += length[i];

  initial = dbl_array1(npt);
  s = 0.0;
  for (i=0; i<npt; i++) {
    initial[i] = s/vintval;
    s += length[i];
  }

  for (i=k=0; i<npt; i++)
    for (j=0; j<fgeom[i][0]->vcontpts; j++,k++)
      bspls->vknots[k] = fgeom[i][0]->vknots[j]*length[i]/vintval + initial[i];

  free_darray1(length);
  free_darray1(initial);

  for (i=0; i<fgeom[0][0]->vorder; i++,k++)
    bspls->vknots[k] = 1.0;
}

/* To make the patches in the same direction have the same knot vector */

void renew_knots(ParSurf ***fgeom, int nps, int npt)
{
  ParSurf *temp;
  double *knots, mn;
  int i, *index, j, k, m, n,flag,jk;

  /* the surfaces to be merged form an array of surfaces fgeom[j][i] */
  /* the rows of the array are in the V parametric direction */
  /* the columns of the array are in the U parametric direction */

  for (i=0; i<nps; i++) {        /* for this u column of the surfs matrix */
    if (npt > 1) {               /* merge the u knots of the surfs in row */
    for (j=n=0; j<npt; j++)      /* if more than 1 surfs in v direction */
      n += (fgeom[j][i]->ucontpts - fgeom[j][i]->uorder);
    /* n is the maximum possible size of the merged knot vector */
    
    knots = dbl_array1(2*fgeom[0][i]->uorder + n);
    index = int_array1(npt);     /* index[i] is pointer to current */
    for (j=0; j<npt; j++)        /* knot in each knot vector */
      index[j] = fgeom[0][i]->uorder;  /* start with first internal knot */

    for (j=k=0; j<fgeom[0][i]->uorder; j++)
      knots[k++] = 0.0;

    m = 1;
    while (m >= 0) {             /* while there are still unused knots... */
      mn = 1.0;
      m = -1;
      for (j=0; j<npt; j++)      /* find next smallest knot value */
	if (index[j] < fgeom[j][i]->ucontpts) /* if knots are left in vector */
	  if (fgeom[j][i]->uknots[index[j]] < mn) { /* if value is smallest */
	    mn = fgeom[j][i]->uknots[index[j]];     /* save value */
	    m = j;                                  /* save index of vector */
	  }
      if (m >= 0) {
	index[m]++;
        jk=index[m]-2;
	if ((mn > knots[k-1] + 1.0e-6)||(fabs(knots[k-1]-fgeom[m][i]->uknots[jk])<1.0e-8)){
                      /* if this value is far enough away */
	                               /* from the previous value or they
                                          are from the same surface*/
	  /*
          printf("&&&&&&&&&&&&*******\n");
          printf("m=%2d  mn=%10.8f\n",m,mn);
	  printf("k=%2d, knots[k-1]=%10.8f\n",k,knots[k-1]);
	  printf("jk=%2d\n",jk);
	  printf("fgeom[m][i]->uknots[index[m]-2]=%10.8f\n",
		 fgeom[m][i]->uknots[jk]);
	  */
	  knots[k++] = mn; 
	}
      }
    }
    for (j=0; j<fgeom[0][i]->uorder; j++)
      knots[k++] = 1.0;

    free_iarray1(index);

    for (j=0; j<npt; j++) {           /* use Oslo to update knot vector */
      temp = copyfgeom(fgeom[j][i], NULL);
      free_fgeom(fgeom[j][i]);
      fgeom[j][i] = fgeomalloc1(temp->uorder, temp->vorder, k - temp->uorder,
				temp->vcontpts);
      copyfgeom(temp, fgeom[j][i]);
      fgeom[j][i]->ucontpts = k - temp->uorder;
      for (m=0; m<k; m++)
	fgeom[j][i]->uknots[m] = knots[m];
      surfoslo3(temp, fgeom[j][i], 4);
      free_fgeom(temp);
    }

    free_darray1(knots);
    }
  }

  for (j=0; j<npt; j++) {        /* for this v row of the surfs matrix */
    if (nps > 1) {               /* merge the v knots of the surfs column */
    for (i=n=0; i<nps; i++)      /* if more than 1 surfs in u direction */
      n += (fgeom[j][i]->vcontpts - fgeom[j][i]->vorder);
    /* n is the maximum possible size of the merged knot vector */
    
    knots = dbl_array1(2*fgeom[j][0]->uorder + n);
    index = int_array1(nps);     /* index[i] is pointer to current */
    for (i=0; i<nps; i++)        /* knot in each knot vector */
      index[i] = fgeom[j][0]->vorder;  /* start with first internal knot */

    for (i=k=0; i<fgeom[j][0]->vorder; i++)
      knots[k++] = 0.0;

    m = 1;
    while (m >= 0) {             /* while there are still unused knots... */
      mn = 1.0;
      m = -1;
      for (i=0; i<nps; i++)      /* find next smallest knot value */
	if (index[i] < fgeom[j][i]->vcontpts) /* if knots are left in vector */
	  if (fgeom[j][i]->vknots[index[i]] < mn) { /* if value is smallest */
	    mn = fgeom[j][i]->vknots[index[i]];     /* save value */
	    m = i;                                  /* save index of vector */
	  }
      if (m >= 0) {
	index[m]++;
	if ((mn > knots[k-1] + 1.0e-6) 
	    ||(fabs(knots[k-1]-fgeom[j][m]->vknots[index[m]-2])<1.0e-8)) {
	  /*
             printf("k=%2d, knots[k-1]=%10.8f\n",k,knots[k-1]);
             printf("index[m]-2=%2d\n",index[m]-2);
	     printf("fgeom[j][m]->vknots[index[m]-2]=%10.8f\n",
		    fgeom[j][m]->vknots[index[m]-2]);
	  */
             knots[k++] = mn;     /* if this value is far enough away */
	                          /* from the previous value or are from
                                     the same surface*/
	}
      }
    }
    for (i=0; i<fgeom[j][0]->vorder; i++)
      knots[k++] = 1.0;

    free_iarray1(index);

    for (i=0; i<nps; i++) {           /* use Oslo to update knot vector */
      temp = copyfgeom(fgeom[j][i], NULL);
      free_fgeom(fgeom[j][i]);
      fgeom[j][i] = fgeomalloc1(temp->uorder, temp->vorder, temp->ucontpts,
				k - temp->vorder);
      copyfgeom(temp, fgeom[j][i]);
      fgeom[j][i]->vcontpts = k - temp->vorder;
      for (m=0; m<k; m++)
	fgeom[j][i]->vknots[m] = knots[m];
      surfoslo3(temp, fgeom[j][i], 4);
      free_fgeom(temp);
    }

    free_darray1(knots);
    }
  }
}

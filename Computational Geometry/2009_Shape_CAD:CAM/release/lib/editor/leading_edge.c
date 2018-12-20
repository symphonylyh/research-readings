/* Copyright (C) 1992 Massachusetts Institute of Technology, Cambridge, MA
   All rights reserved
*/

/* leading_edge.c */

/* eval_LE2D()
   eval_LE3D()
   ExtractLeadingEdge2D()
   ExtractLeadingEdge3D()
   lefunct1()
*/

#include "gen.h"
#include "bspl.h"
#include "appr.h"
#include "editor.h"

/*********************** obsolete *****************************************
void e04jaf_(int *, int *, double *, double *, double *, double *, double *,
	     int *, double *, int *, int *);
***************************************************************************/
/*********************** for NAG Marks 19,20 ******************************/
void e04jyf_(int *, int *, void(int*,double*,double*,int*,double*),
             double *, double *, double *, double *, int *, int *, double *, 
             int *, int *, double *, int *);
/**************************************************************************/


static ParSurf *fgeom;
static double g_u = 1.0;

ParCurv *ExtractLeadingEdge2D(ParSurf *fgm, double *eps, short nbeg, 
			      short nmax)
{
  ParCurv *egeom;

  fgeom = fgm;

  /************* obsolete ****************
  set_fun1_ptr(lefunct1);
  **************************************/

  egeom = approx_fnbc(eval_LE2D, eps, 4, 0, nbeg, nmax);

  return (egeom);
}

ParCurv *ExtractLeadingEdge3D(ParSurf *fgm, double *eps, short nbeg, 
			      short nmax)
{
  ParCurv *egeom;

  fgeom = fgm;

  /************* obsolete ****************
  set_fun1_ptr(lefunct1);
  **************************************/

  egeom = approx_fnbc(eval_LE3D, eps, 4, 0, nbeg, nmax);

  return (egeom);
} 

vector *eval_LE2D(double u, int ideriv)
{
  int n, ibound, liw, lw, ifail;
  /******* obsolete *****
  double bl[2], bu[2], x[2], iw[4], w[25], f;
  **********************/
  /******** for NAG Marks 19,20 **********/
  double bl[2], bu[2], x[2], w[25], f;
  int iw[4];
  /***************************************/

  vector *v;
  /********* for NAG Marks 19,20 ********/
  int iuser[1];   double user[1];
  /**************************************/
  
  g_u = u;
  v = vectalloc();
  
  if (ideriv == 0) {
    n = 1;
    ibound = 0;
    bl[0] = 0.5;
    bu[0] = 0.9;
    x[0] = 0.5;
    liw = 3;
    lw = 13;
    ifail = 1;

    /************* obsolete ****************
    e04jaf_(&n, &ibound, bl, bu, x, &f, iw, &liw, w, &lw, &ifail);
    **************************************/
    /********* for NAG Marks 19,20 ********/
    e04jyf_(&n, &ibound, lefunct1, bl, bu, x, &f, iw, &liw, w, &lw, 
            iuser, user, &ifail);
    /**************************************/

    v->x = u;
    v->y = x[0];
    v->z = 0.0;
    v->w = 1.0;
  }
  else
    v->x = v->y = v->z = v->w = 1.0;

  return (v);
}

vector *eval_LE3D(double u, int ideriv)
{
  int n, ibound, liw, lw, ifail;
  /***** obsolete *****
  double bl[2], bu[2], x[2], iw[4], w[25], f, t;
  ********************/
  /***** for NAG Marks 19,20 *******/
  double bl[2], bu[2], x[2], w[25], f, t;
  int iw[4];
  /*********************************/

  vector *v;
  /********* for NAG Marks 19,20 ********/
  int iuser[1];   double user[1];
  /**************************************/
  
  t = u;
  if (ideriv == 0) {
    n=2;
    ibound = 0;
    bl[0] = 0.2;
    bu[0] = 0.9;
    bl[1] = 0.5;
    bu[1] = 0.9;
    /*  We use bl[1] and bu[1] to fix v  */
    x[0] = 0.5;
    x[1] = 0.3;
    liw = 4;
    lw = 25;
    ifail = 1;
    
    bl[0] = bu[0] = t;

    /************* obsolete ****************
    e04jaf_(&n, &ibound, bl, bu, x, &f, iw, &liw, w, &lw, &ifail);
    **************************************/
    /********* for NAG Marks 19,20 ********/
    e04jyf_(&n, &ibound, lefunct1, bl, bu, x, &f, iw, &liw, w, &lw, 
            iuser, user, &ifail);
    /**************************************/

    v = revalderivsurf(fgeom, x[0], x[1], 0, 0);
  }
  else {
    v = vectalloc();
    v->x = v->y = v->z = v->w = 1.0;
  }
  return(v);
}

void lefunct1(int *n, double *xc, double *fc, int *iuser, double *user)
{
  double kmin, kmax;
  
  curvature1(fgeom, g_u, xc[0], &kmax, &kmin);
  
  *fc = -kmax;
}

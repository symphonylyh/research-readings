/* ***************************************************************************
 Copyright (C) 1996 Massachusetts Institute of Technology all rights reserved 
	Programmer: George A. Kriezis
**************************************************************************** */
#include <stdio.h>
#include <math.h>
#include "gen.h"
#include "bspl.h"

/*****************************************************************************
*                                  merge_knotv()
******************************************************************************
* 
* 1   Purpose
*     This function merges the knot vectors of n NURBS curves and updates 
*     their control points to reflect the new augmented common knot vector.
* 
* 2   Specification
*     #include "bspl.h"
*     void merge_knotv(ParCurv **eg, int is, int n, double *knot, int *index, 
*                      double eps)
* 
* 3   Description
*     This routine merges the knot vectors of n NURBS curves (stored in eg[is] 
*     through eg[is + n - 1]) and forms a common knot vector. This routine 
*     returns the number of knots in the common knot vector.
* 
* 5   Parameters
*       1.ParCurv ** eg
*         On entry: an array of NURBS curve data structures of length  n which 
* 	          contain the geometry for the NURBS curves whose knots are to 
* 		  be merged.
*       2.int is
*         On entry: index of starting element in array eg[].
*       3.int n
*         On entry: number of curves to consider in knot merging.
*       4.double * knots
*         On entry: array of length k, where k = the total number of knots in 
* 	          the merged knot vector. The value of nu is such that: 
* 		  max (eg[i] -> uorder + eg[i] -> ncontpts)
* 		   i
* 		  <= nu <=
* 		  sum{i=is..(is+n-1)}(eg[i] -> uorder+eg[i] -> ncontpts)
*         On exit: array containing common knot vector.
*       5.int * index
*         On entry: array of length n used as workspace.
*       6.double eps
*         On entry: numerical tolerance for determining when double precision 
* 	          numbers are equal.
* 
* 6   Return Values, Error Indicators and Warnings
*     Not applicable.
* 
* 7   Accuracy
*     Not applicable.
* 
* 8   Further Comments
*     Knot addition is performed using surfolso3 (bspl _ oslo3.c). Knot 
*     merging is performed as in routine merge_knotv (bspl _ merge.c).
*     This routine is usually followed by routine addpoints (bspl-addpoints.c) 
*     which uses the oslo algorithm to update the control points of the above 
*     curves based on the new common knot vector.
* 
* 9   Functions referenced by merge_knotv() are:
*     find_nextknot()
*  
* 10  Functions that reference merge_knotv() are:
*     BlendSurfCB()
*     camber_surface()
*     ConvexHullTest()
*     convexhull_test()
*     global_pos_error()
*     LoftSurfCB()
*     MergeKnots()
*     merge_can_off()
*     merge_fgeom()
*     ParSurf_approx()
* 
******************************************************************************/

int merge_knotv(ParCurv *eg[], int is, int n, double kn[], int index[],
		double eps)
/* Merge the knot vectors of n B-spline curves (stored in eg[is]
 * through eg[is+n-1]) and form a common knot vector. index is a working 
 * space of length n */
{
  int k,i ;    /* k is the current length of the common knot vector. */

  k = 0 ; 
  /* set the first knot in each knot vector as the current knot in each 
   * knot vector */
  for (i=0;i<n;i++) index[i] = 0 ;  
  /* merge the next available knot into the common knot vector.
   * if the first knot vector reaches its end, then stop. */
  while (find_nextknot(eg,is,n,index,kn,k,eps))
         ++k;

  return(++k) ;
}

/*****************************************************************************
*                              find_nextknot()
******************************************************************************
* 
* 1   Purpose
*     This function is a subroutine of merge_knotv(), which merges the knot
*     vectors of a set of B-spline curves. It finds the next available knot.
* 
* 2   Specification
*     #include "bspl.h"
*     int find_nextknot(ParCurv *eg[], int is, int n, int index[], double kn[],
* 		      int k, double eps)
* 
* 3   Description
*     This function finds the current minimum knot among the n knot vectors 
*     and compile this knot into the common knot vector. 
* 
* 4   References
*     Not applicable
* 
* 5   Parameters
*       1.ParCurv *eg[]
*         On entry: an array of NURBS curve data structures of length  n which 
*                   contain the geometry for the NURBS curves whose knots are 
*                   to be merged.
*       2.int is
*         On entry: index of starting element in array eg[].
*       3.int n
*         On entry: number of curves to consider in knot merging.
*       4.double kn[]
*         On entry: array of length k, where k = the total number of knots in 
*                   the merged knot vector. The value of nu is such that: 
*                   max (eg[i] -> uorder + eg[i] -> ncontpts)
*                    i
*                   <= nu <=
*                   sum{i=is..(is+n)}eg[i] -> uorder+eg[i] -> ncontpts
*         On exit: array containing common knot vector.
*       5.int index[]
*         On entry: array whose elements indicate the current knot in
* 	            corresponding knot vector.
*         On exit:  updated after compiling a knot into the common knot vector.
*       6.double eps
*         On entry: numerical tolerance for determining when double precision 
*                   numbers are equal.
* 
* 6   Return Values, Error Indicators and Warnings
*     If the current knot in the first knot vector reaches the end, return 0;
*     otherwise return 1.
* 
* 7   Accuracy
*     Not applicable.
* 
* 8   Further Comments
*     see also merge_knotv() and find_minknot().
* 
* 9   Functions referenced by find_nextknot() are:
*     find_minknot()
* 
* 10  Functions that reference find_nextknot() are:
*     merge_knotv()
* 
******************************************************************************/

int find_nextknot(ParCurv *eg[], int is, int n, int index[], double kn[],
		  int k, double eps)
/* Compile the next knot (k) to the common knot vector (kn) out of n knot
 * vectors belonging to eg[i]. Consider two knots as identical if their
 * values differ by an amount less than eps */
{
  int i,j ;  /* i is the index of the curve whose current knot is to be 
              * merged into the common knot vector next.
              * j is the index of the knot, which is to be merged next,
	      * in the knot vector of i-th curve. */

  /* find the index of the curve whose current knot is to be merged next. */
  i=find_minknot(eg,is,n,index,eps) ;
  /* get the index of the knot to be merged next. */
  j = index[i]-1 ;
  kn[k] = eg[i+is]->knots[j] ;    /* put the knot into the common knot vector */
  /* if the current knot in 1st curve reaches its end, return 0;
   * otherwise, return 1 */
  if (index[0]==eg[is]->ncontpts+eg[is]->order) return(0) ; else return(1) ;
}

/*****************************************************************************
                             find_minknot()
******************************************************************************
*
* 1   Purpose
*     This function is a subroutine of find_nextknot(). It finds the minimum
*     knot.
* 
* 2   Specification
*     #include "bspl.h"
*     int find_minknot(ParCurv *eg[], int is, int n, int index[], double eps)
* 
* 3   Description
*     This function finds the index of the B-spline curve which has the 
*     current minimum knot among the n curves whose knot vectors are to be
*     merged. 
* 
* 4   References
*     Not applicable
* 
* 5   Parameters
*       1.ParCurv *eg[]
*         On entry: an array of NURBS curve data structures of length  n which 
*                   contain the geometry for the NURBS curves whose knots are 
*                   to be merged.
*       2.int is
*         On entry: the index of starting element in the above array.
*       3.int n
*         On entry: the number of curves to be merged
*       4.int index[]
*         On entry: the array whose elements indicate the current knot in 
* 	            corresponding knot vector.
*         On exit:  updated after finding the minimum knot which is to be
*                   merged.
*       5.double eps
*         On entry: tolerance for comparison of the equality of two double
*                   values.
* 
* 6   Return Values, Error Indicators and Warnings
*     Not applicable.
* 
* 7   Accuracy
*     Not applicable.
* 
* 8   Further Comments
*     see also find_nextknot() and merge_knotv().
* 
* 9   Functions referenced by find_minknot() are:
* 
* 10  Functions that reference find_minknot() are:
*     find_nextknot()
*     size_nextknot()
* 
******************************************************************************/

int find_minknot(ParCurv *eg[], int is, int n, int index[], double eps)
/* Find the minimum knot out of n knot vectors belonging to eg[0] through
 * eg[n-1] starting at the position index[0] through index[n-1],
 * respectively */
{
  int i,j ;   /* j will be the index of the curve which has the current
	       * minimum knot. */
  double t,minknot ;   /* minknot is the current minimum value of knot. */

  /* determine the minimum knot */
  minknot = 1.0/eps ;    /* set the current minimum value as a very large
			  * number */
  for (i=0;i<n;i++) {
      t = eg[i+is]->knots[index[i]] ;       /* get the current knot from i-th
					     * curve */
      if (t<minknot) { minknot = t ; j = i ; }  /* compare with the currnet
						 * minimum */
      }

  /* update index */
  for (i=0;i<n;i++) { 
      /* if the current knot in i-th curve is equal to the minimum, 
       * skip to the next knot. */
      t = eg[i+is]->knots[index[i]]-minknot ; 
      if (t<eps) index[i]++ ;
    }

  return(j) ;
}

/*************************************************************************
*                            addpoints()
**************************************************************************
* 
* 1   Purpose
*     This function updates the knot vectors and control points of an array 
*     of b-spline curves.
* 
* 2   Specification
*     #include "bspl.h"
*     void addpoints(ParCurv **egeom, int is, int n, double *knot, int k)
* 
* 3   Description
*     This routine updates the knot vectors and the control points of n 
*     B-spline curves (stored in eg[is] through eg[is + n - 1]) using the 
*     Oslo algorithm for subdividing the curves in the new knot vector.
* 
* 4   References
*     [1]W. Boehm, Inserting New Knots Into B-Spline Curves, CAD, 
*     (12)4:199-201, 1980.
* 
* 5   Parameters
*       1.ParCurv ** egeom
*         On entry: an array of pointers to NURBS curve.
*       2.int is
*         On entry: integer defining starting element of NURBS curve array 
* 	to update.
*       3.int n
*         On entry: number of NURBS B-spline curves to update.
*       4.double * knot
*         On entry: array containing common knot vector.
*       5.int k
*         On entry: the number of knots in common knot vector.
* 
* 8   Further Comments
*     The routine assumes that the array of NURBS curves have adequate 
*     space to fit all the data after the knot insertion. 
*     See also curve_oslo1 (bspl _ olso1.c).
* 
* 9   Functions referenced by addpoints() are:
*     curve_oslo1()
*     dbl_array2()
*     free_darray2()
*     vectalloc()
* 
* 10  Functions that reference addpoints() are:
*     camber_surface()
*     ConvexHullTest()
*     convexhull_test()
*     global_pos_error()
*     LoftSurfCB()
*     MergeKnots()
*     ParSurf_approx()
* 
***************************************************************************/

void addpoints(ParCurv *eg[], int is, int n, double kn[], int k)
/* update the knot vectors of the curves eg[is] through eg[is+n-1] (re-
 * presenting the border lines and the derivatives thereof). 
 * kn is the common knot vector and k the number of the knots. */
{
  int i,j,l ;
  double **inpoly;   /* control points of a curve before update */
  double **outpoly;  /* control points of a curve after update */

  inpoly = dbl_array2((unsigned)k,4);
  outpoly = dbl_array2((unsigned)k,4);

  for (i=0;i<n;i++) {
      l = i+is ;    /* index of the current curve to be updated */
      for (j=0;j<eg[l]->ncontpts;j++) {
	  /* get the control points of the curve to be updated */
          inpoly[j][0] = eg[l]->contpts[j]->x ;
          inpoly[j][1] = eg[l]->contpts[j]->y ;
          inpoly[j][2] = eg[l]->contpts[j]->z ;
          inpoly[j][3] = eg[l]->contpts[j]->w ;
          }
      /* update the control points by calling knot insertion routine 
       * curve_oslo1() */
      curve_oslo1(eg[l]->order,eg[l]->ncontpts,k,eg[l]->knots,kn,inpoly,
		  outpoly);
      eg[l]->ncontpts = k-eg[l]->order ;          /* number of control points
					           * after update */
      for (j=0;j<k;j++) eg[l]->knots[j] = kn[j] ; /* new knot vector */
      for (j=0;j<eg[l]->ncontpts;j++) {   /* update the control points of the
					   * curve */
          if (eg[l]->contpts[j]==NULL) eg[l]->contpts[j] = vectalloc() ;
          eg[l]->contpts[j]->x = outpoly[j][0] ;
          eg[l]->contpts[j]->y = outpoly[j][1] ;
          eg[l]->contpts[j]->z = outpoly[j][2] ;
          eg[l]->contpts[j]->w = outpoly[j][3] ;
          }
      }
 
  free_darray2(inpoly);
  free_darray2(outpoly);
}

/***************************************************************************
*                               merge_fgeom()
****************************************************************************
* 
* 1   Purpose
*     This function merges the knot vectors of two NURBS surfaces and updates 
*     their control points to reflect the new augmented knot vectors.
* 
* 2   Specification
*     #include "bspl.h"
*     void merge_fgeom(ParSurf *fg1, ParSurf *fg2, ParCurv **eg, double *kn, 
*                      double *kl, int *index, double eps)
* 
* 3   Description
*     This routine merges the knot vectors of the NURBS surfaces and uses the 
*     surface knot addition algorithm to update the control points of the 
*     above surfaces based on the new augmented knot vectors.
* 
* 4   References
*     [1] W. Boehm, Inserting New Knots Into B-Spline Curves, CAD,
*         (12)4:199-201, 1980.
* 
* 5   Parameters
*       1.ParSurf * fg1
*         On entry: NURBS surface data structure containing the geometry of 
* 	          the first surface.
*         On exit: NURBS surface data structure containing the geometry of the 
* 	          first surface with control points added to reflect the 
* 		  augmented knot vector.
*       2.ParSurf * fg2
*         On entry: NURBS surface data structure containing the geometry of the
* 	          second surface.
*         On exit: NURBS surface data structure containing the geometry of the 
* 	          second surface with control points added to reflect the 
* 		  augmented knot vector.
*       3.ParCurv ** eg
*         On entry: an array of pointers to NURBS curve data structures of 
* 	          length >= 2 used as workspace. Each curve data structure 
* 		  must be allocated for a NURBS of 
* 		      order >= max ( fg1 -> uorder, fg1 -> vorder)
* 		  and 
* 		      number of control points >= max ( fg1->ucontpts, 
* 		      fg2 -> ucontpts, fg1 -> vcontpts, fg2 -> vcontpts).
*       4.double * kn
*         On entry: array of length nu, where nu = the total number of knots 
* 	          in the merged knot vector in the u parameter. The value of 
* 		  nu is such that:
* 		  fg1 -> uorder + max (fg1 -> ucontpts, fg2 -> ucontpts)
* 		  <= nu <= 
* 		  fg1 -> uorder+(fg1 -> ucontpts + fg2 -> ucontpts)
*         On exit: array containing common knot vector for the u parameter.
*       5.double * kl
*         On entry: array of length nv, where nv = the total number of knots in
* 	          the merged knot vector in the v parameter. The value of nv 
* 		  is such that:
* 		  fg1 -> vorder + max (fg1 -> vcontpts, fg2 -> vcontpts)
* 		  <= nv <=
* 		  fg1 -> vorder+(fg1 -> vcontpts + fg2 -> vcontpts)
*         On exit: array containing common knot vector for the v parameter.
*       6.int * index
*          On entry: array of length  max(nu; nv) used as workspace.
*       7.double eps
*          On entry: numerical tolerance for determining when double precision 
* 	           numbers are equal.
* 
* 6   Return Values, Error Indicators and Warnings
*     Not applicable.
* 
* 7   Accuracy
*     Not applicable.
* 
* 8   Further Comments
*     Knot addition is performed using surfolso3 (bspl _ oslo3.c). Knot 
*     merging is performed as in routine merge_knotv (bspl _ merge.c).
* 
* 9   Functions referenced by merge_fgeom() are:
*     addpoints_surf()
*     extract_edge()
*     merge_knotv()
* 
* 10  Functions that reference merge_fgeom() are:
*     ConvexHullSurf()
*     convexhull_surf()
* 
******************************************************************************/

void merge_fgeom(ParSurf *fg, ParSurf *fg0, ParCurv *eg[], double kn[],
		 double kl[], int index[], double eps)
/* merge the knot vectors of fg, fg0. Add knots as required */
{
  int i,j,n1,n2,uo,vo,up,vp,up0,vp0;  /* indices and some temporary
				       * variables */
  int ku,kv ;           /* number of knots in common knot vector in u,v */

  /* extract the edge v=0 of the 1st surface as eg[0], which is a B-spline 
   * curve in u with the same order as the surface in u. */
  extract_edge(eg[0],fg,1,0) ; 
  /* extract the edge v=0 of the 2nd surface as eg[1], which is a B-spline 
   * curve in u with the same order as the surface in u. */
  extract_edge(eg[1],fg0,1,0) ;
  /* merge the knot vectors in u of the two surfaces and get the common knot
   * vector in u. */
  ku = merge_knotv(eg,0,2,kn,index,eps) ;
  /* extract the edge u=0 of the 1st surface as eg[0], which is a B-spline 
   * curve in v with the same order as the surface in v. */
  extract_edge(eg[0],fg,0,0) ; 
  /* extract the edge u=0 of the 2nd surface as eg[1], which is a B-spline 
   * curve in v with the same order as the surface in v. */
  extract_edge(eg[1],fg0,0,0) ;
  /* merge the knot vectors in v of the two surfaces and get the common knot
   * vector in v. */
  kv = merge_knotv(eg,0,2,kl,index,eps) ;

  uo = fg0->uorder ; vo = fg0->vorder ; up = fg->ucontpts ; vp = fg->vcontpts ;
  up0 = fg0->ucontpts ; vp0 = fg0->vcontpts ;

  /* update the control points of the 1st surface */
  n1 = uo+up-1 ;   /* number of control points in u direction of the 1st 
		    * surface before update */
  i = 0 ; 
  while (++i<n1 && fg->uknots[i]==kn[i]) ;   /* i now points to the 1st knots
                      * which are unequal between the old knot vector and the
		      * common knot vector in u */
  n2 = vo+vp-1 ;     /* number of control points in v direction of the 1st
                      * surface before update */
  j = 0 ; 
  while (++j<n2 && fg->vknots[j]==kl[j]) ;   /* same as i but in v direction */
  /* if either of the knot vectors in u,v directions is not the same as the 
   * corresponding common knot vector, update the control points */
  if (i<n1 || j<n2) addpoints_surf(fg,kn,kl,ku,kv) ;

  /* update the control points of the 2nd surface similarly */
  n1 = uo+up0-1 ; i = 0 ; while (++i<n1 && fg0->uknots[i]==kn[i]) ;
  n2 = vo+vp0-1 ; j = 0 ; while (++j<n2 && fg0->vknots[j]==kl[j]) ;
  if (i<n1 || j<n2) addpoints_surf(fg0,kn,kl,ku,kv) ;
}

/***************************************************************************
*                          addpoints_surf()
****************************************************************************
* 
* 1   Purpose
*     This function updates the knot vectors and control points of a NURBS 
*     surface.
* 
* 2   Specification
*     #include "bspl.h"
*     void addpoints_surf(ParSurf *fgeom, double *kn, double *kl, int ku,
*                         int kv)
* 
* 3   Description
*     This routine updates the knot vectors and the control points of a NURBS 
*     surface using the Oslo algorithm for adding knots to a NURBS surface 
*     given new knot vectors.
* 
* 
* 4   References
*     [1] W. Boehm, Inserting New Knots Into B-Spline Curves, CAD, 
*     (12)4:199-201, 1980.
* 
* 5   Parameters
*        1.ParSurf * fgeom
*          On entry: pointer to NURBS surface
*        2.double * kn
*          On entry: array containing common knot vector in the u direction
*        3.double * kl
*          On entry: array containing common knot vector in the v direction
*        4.int ku
*          On entry: number of knots in common knot vector in the u direction
*        5.int kv
*          On entry: number of knots in common knot vector in the v direction
* 
* 6   Return Values, Error Indicators and Warnings
*     Not applicable.
* 
* 7   Accuracy
*     Not applicable.
* 
* 8   Further Comments
*     This routine calls routine surfoslo3 (bspl _ oslo3.c). It assumes that 
*     the NURBS surface data structure has adequate space reserved to 
*     accommodate the addition of points.
* 
* 9   Functions referenced by addpoints_surf() are:
*     copyfgeom()
*     fgeomalloc1()
*     free_fgeom()
*     surfoslo3()
* 
* 10  Functions that reference addpoints_surf() are:
*     merge_fgeom()
* 
***************************************************************************/

void addpoints_surf(ParSurf *fg, double kn[], double kl[], int ku, int kv) 
/* compute the vertices of fg with the refined knot vectors kn, kl having
 * ku, kv knots */
{
  int i,j ; 
  struct fgeom *fgold;   /* temporarily to hold the information of the old
			  * surface */
  /* assume fg adequate initial space */
  fgold = fgeomalloc1(fg->uorder,fg->vorder,fg->ucontpts,fg->vcontpts);

  copyfgeom(fg,fgold);   /* copy the old surface to fgold */

  /* update the number of control points in u,v */
  fg->ucontpts = ku-fg->uorder;
  fg->vcontpts = kv-fg->vorder;
  /* update the knot vectors in u,v */
  for (i=0;i<ku;i++) fg->uknots[i] = kn[i] ;
  for (i=0;i<kv;i++) fg->vknots[i] = kl[i] ;
  /* updata the control points by calling knot insertion routine surfoslo3() */
  surfoslo3(fgold,fg,4);

  free_fgeom(fgold);
}

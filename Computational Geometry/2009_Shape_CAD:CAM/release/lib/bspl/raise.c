/* ***************************************************************************
 Copyright (C) 1996 Massachusetts Institute of Technology all rights reserved 
	Programmer: George A. Kriezis
**************************************************************************** */
#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include "gen.h"
#include "bspl.h"

/*****************************************************************************
*                               compute_lim()
******************************************************************************
* 
* 1   Purpose
*     This function calculates the array lambda and a new knot vector useful
*     for degree elevation (by one) of a NURBS curve.
* 
* 2   Specification
*     #include "bspl.h"
*     int compute_lim(ParCurv *egeom, double **lim, double *knot, double eps)
* 
* 3   Description
*     This routine computes the lambda matrix and the new knot vector for 
*     elevating the degree of a NURBS curve by one. This routine is used in 
*     conjunction with routine raise_byone() to degree elevate a NURBS curve 
*     by one. The routine returns the number of knots in the new knot vector. 
*     For explanation of the lambda matrix consult the reference list.
* 
* 4   References
*     [1]E. Cohen, T. Lyche, and L. Schumaker. Algorithms for degree raising 
*     of splines, ACM Transactions on Graphics 4(3):171-181, July 1985.
* 
* 5   Parameters
*       1.ParCurv * egeom
*         On entry: NURBS curve geometry data structure of curve that
* 	          eventually will be degree elevated
*       2.double ** lim
*         On exit: two dimensional array which contains values for the lambda 
* 	         matrix used in degree elevation
*       3.double * knot
*         On exit: array of new knot vector corresponding to higher order curve
*       4.double eps
*         On entry: tolerance used in comparing two knots values
* 
* 6   Return Values, Error Indicators and Warnings
*     Not applicable.
* 
* 7   Accuracy
*     Not applicable.
* 
* 8   Further Comments
*     Space for array lim should be provided from the calling subprogram. 
*     Array lim should be of length  [ncontpts + order - 1][ncontpts]. 
*     See also raise_byone (bspl _ raise.c).
* 
* 9   Functions referenced by compute_lim() are:
*     dbl_array2()
*     free_darray2()
* 
* 10  Functions that reference compute_lim() are:
*     ApproxCurvCB()
*     ConvexHullTest()
*     convexhull_test()
*     LoftSurfCB()
*     MergeKnots()
*     raise_surf()
* 
****************************************************************************/

int compute_lim(ParCurv *egeom, double **lim, double kv[], double eps)
/* raise the order of egeom by one */
{
  int i,j,m;             /* loop indices */
  int k,l,n,o,n1,n2 ;    /* k -- number of knots inserted after degree
                          *      elevation.
			  * l -- index for counting know knot vector
			  * n -- number of control points of the old curve
			  * o -- order of the old curve
			  * n1,n2 -- results of some inequalities  */
  double **aim;          /* alpha matrix, which contains discrete B-spline. */
  double d,d1,d2,s ;     /* working variables */

  /* get the number of control points and the order of the curve. */
  n = egeom->ncontpts ; o = egeom->order ;
  /* allocate memory */
  aim = dbl_array2((unsigned)(2*n+2),(unsigned)(2*n-o+1));

  /* compute the new knot vector - store in kv */
  l = 0;      /* initialize the index */
  /* increase the multiplicity of each knot by 1 */      
  for (i=0;i<n+o;i++) {
      /* everytime the value of knot change, we insert a knot with value
       * knots[i-1] before it. Notice the condition in 'if' also deals with
       * the boundary at i=0. */
      if (i && egeom->knots[i]>egeom->knots[i-1]+eps) 
         kv[l++] = egeom->knots[i-1] ;
      kv[l++] = egeom->knots[i] ;
      }	
  kv[l] = egeom->knots[n+o-1] ; 
  k = l-o-n ;

  /* compute a sub(i) sup(1) (j) and l sub(i) sup(1) (j) */
  for (i=0;i<n+o;i++) 
      for (j=0;j<n+k;j++)
          if (egeom->knots[i]<=kv[j] && kv[j]<egeom->knots[i+1])
             aim[i][j] = lim[i][j] = 1.0 ; 
	  else aim[i][j] = lim[i][j] = 0.0 ;

  /* compute a sub(i) sum(m) (j) and l sub(i) sup(m) (j) iteratively */
  for (m=2;m<=o;m++) 
      for (i=0;i<n+o-m;i++) {
          d1 = egeom->knots[i+m-1]-egeom->knots[i] ; n1 = d1<eps ;
          d2 = egeom->knots[i+m]-egeom->knots[i+1] ; n2 = d2<eps ;
          for (j=0;j<n+k;j++) {
              if (aim[i][j]<eps) d = 0.0 ; 
	      else {
		    s = kv[j+m-1]-egeom->knots[i] ;
                    if (s<eps && n1) d = 0.0 ; 
		    else d = s*aim[i][j]/d1 ;
                   }
              if (aim[i+1][j]>eps) {
                 s = egeom->knots[i+m]-kv[j+m-1] ;
                 if (s>eps || !n2) d +=s*aim[i+1][j]/d2 ;
                 }
              aim[i][j] = d ;
              if (lim[i][j]>eps) {
                 s = kv[j+m]-egeom->knots[i] ;
                 if (s>eps || !n1) d += s*lim[i][j]/d1 ; 
                 }
              if (lim[i+1][j]>eps) {
                 s = egeom->knots[i+m]-kv[j+m] ;
                 if (s>eps || !n2) d += s*lim[i+1][j]/d2 ; 
                 }
              lim[i][j] = d ;
              }
          }

  /* free memory */
  free_darray2(aim);

  /* return the number of knots inserted. */
  return(k) ;
}

/****************************************************************************
*                                raise_byone()
*****************************************************************************
* 
* 1   Purpose
*     This function elevates the degree of a NURBS curve by one.
* 
* 2   Specification
*     #include "bspl.h"
*     void raise_byone(ParCurv *egeom, double **lim, double *knot, int k)
* 
* 3   Description
*     This routine elevates the degree of a NURBS curve by one, computing the
*     new control points after order raising and updating the data structure.
* 
* 4   References
*     [1] E. Cohen, T. Lyche, and L. Schumaker. Algorithms for degree raising
*         of splines, ACM Transactions on Graphics 4(3):171-181, July 1985.
* 
* 5   Parameters
*       1.ParCurv * egeom
*         On entry: NURBS curve data structure containing geometry of the 
* 	          original curve.
*         On exit: NURBS curve data structure containing geometry of the 
* 	         degree elevated curve.
*       2.double * lim
*         On entry: matrix containing the lambda matrix used in degree
*                 elevation.
* 	          For explanation of the lambda matrix consult the reference 
* 		  list.
*       3.double * knot
*         On entry: array of length k containing the knot vector of the degree 
* 	          elevated curve.
*       4.int k
*         On entry: the length of the array knot.
* 
* 6   Return Values, Error Indicators and Warnings
*     Not applicable.
* 
* 7   Accuracy
*     Not applicable.
* 
* 8   Further Comments
*     To compute the elements of the matrix lim and array knot, consult 
*     routine compute_lim (bspl _ raise.c).
* 
* 9   Functions referenced by raise_byone() are:
*     add_vect1()
*     copyegeom()
*     copyvector()
*     egeomalloc1()
*     free_egeom()
*     free_varray1()
*     scale_vect1()
*     vectalloc()
*     vec_array1()
* 
* 10  Functions that reference raise_byone() are:
*     ApproxCurvCB()
*     ConvexHullTest()
*     convexhull_test()
*     LoftSurfCB()
*     MergeKnots()
*     raise_surf()
* 
******************************************************************************/

void raise_byone(ParCurv *egeom, double **lim, double kv[], int k)
/* assign new edge data - compute new vertices after order raising */
{
  struct egeom *egeom1;      /* never used! (glshen) */
  vector **q,*v ;            /* v -- working array
			      * q -- control points of the new curve  */ 
  double d ;                 /* factor in the formula which calculates the new
			      * control points. */
  int i,j,n ;                /* number of control points of the old curve */
  q = vec_array1(egeom->ncontpts+k); 

  d = 1.0/(double)egeom->order ; 
  n = egeom->ncontpts ;

  /* construct another curve to keep the information of the old curve */
  egeom1=egeomalloc1(egeom->order,egeom->ncontpts);
  copyegeom(egeom,egeom1);

  /* set the order and the number of control points of the new curve */
  egeom->order += 1; egeom->ncontpts += k;

  v = vectalloc();
  /* set the knot vector of the new curve */
  for (i=0;i<egeom->ncontpts+egeom->order;i++) egeom->knots[i] = kv[i] ;
  /* calculate the control points of the new curve */
  for (j=0;j<egeom->ncontpts;j++) {
      /* initialize each control point */
      q[j]->x = q[j]->y = q[j]->z = 0.0 ; q[j]->w = 1.0 ;
      /* calculate the control points by the formula in reference */
      for (i=0;i<n;i++) {
          scale_vect1(lim[i][j],egeom->contpts[i],v) ;
          add_vect1(v,q[j],q[j]) ;
          }
      scale_vect1(d,q[j],q[j]) ;
      }
  /* set the control points of the new curve */
  for (j=0;j<egeom->ncontpts;j++) {
      if (egeom->contpts[j] == NULL)
         egeom->contpts[j] = vectalloc();
         copyvector(q[j],egeom->contpts[j]) ; 
      }

  /* free memory */
  free_varray1(q,egeom->ncontpts);
  free_egeom(egeom1);
  free ((char *)v) ;
}

/*****************************************************************************
*                                raise_surf()
******************************************************************************
*
* 1   Purpose
*     This function elevates the degree of a NURBS surface by one.
*
* 2   Specification
*     #include "bspl.h"
*     void raise_surf(ParSurf *fgeom, double eps, int u, int v)
* 
* 3   Description
*     This routine elevates the degree of a NURBS surface in either or both 
*     directions by a specified degree. It computes the new control points 
*     after order raising and updating the data structure.
* 
* 4   References
*     [1] E. Cohen, T. Lyche, and L. Schumaker. Algorithms for degree raising 
*         of splines, ACM Transactions on Graphics 4(3):171-181, July 1985.
* 
* 5   Parameters
*        1.ParSurf * fgeom
*          On entry: NURBS surface data structure containing geometry of the 
* 	           original surface.
*          On exit: NURBS surface data structure containing geometry of the 
* 	           degree elevated surface.
*        2.double * lim
*          On entry: matrix containing the lambda matrix used in degree
*                  elevation.
* 	           For explanation of the lambda matrix consult the reference l
* 		   ist.
*        3.double eps
*          On entry: tolerance used to determine the equivalence of two double 
* 	           precision numbers.
*        4.int u
*          On entry: number of degrees to raise surface in the u direction
*        5.int v
*          On entry: number of degrees to raise surface in the v direction
* 
* 6   Return Values, Error Indicators and Warnings
*     Not applicable.
* 
* 7   Accuracy
*     Not applicable.
* 
* 8   Further Comments
*     Not applicable.
* 
* 9   Functions referenced by raise_surf() are:
*     compute_lim()
*     dbl_array1()
*     dbl_array2()
*     egeomalloc2()
*     extract_edge()
*     free_darray1()
*     free_darray2()
*     free_egeom()
*     raise_byone()
*     store_vertices()
* 
* 10  Functions that reference raise_surf() are:
*     ApproxSurfCB()
*     ConvexHullSurf()
*     convexhull_surf()
*  
****************************************************************************/

void raise_surf(ParSurf *fgeom, double eps, int u, int v)
/* raise the orders of fgeom by u,v */
{
  double **lim,*kv;
  struct egeom *eg ; int i,j,k ;
  int uorder,vorder,ucontpts,vcontpts;

  /* maximum possible numbers of control points in two directions after
   * degree elevation */
  ucontpts= fgeom->ucontpts + u*(fgeom->ucontpts - fgeom->uorder +1);
  vcontpts= fgeom->vcontpts + v*(fgeom->vcontpts - fgeom->vorder +1);
  if (vcontpts > ucontpts)
     ucontpts = vcontpts;

  /* orders of the new surface */
  uorder = fgeom->uorder + u;
  vorder = fgeom->vorder + v;
  if (vorder > uorder)
     uorder = vorder;

  /* allocate memory which is large enough to contain a curve extracted from 
   * the surface in either direction */
  eg = egeomalloc2((short)uorder,(short)ucontpts) ;
  /* allocate memory for lamda matrix and new knot vector */
  lim = dbl_array2((unsigned)(uorder+ucontpts),(unsigned)ucontpts);
  kv = dbl_array1((unsigned)(uorder+ucontpts));

  /* if elevate degree in u */
  if (u) 
     for (i=0;i<u;i++) {   /* elevate degree u times */
         /* extract the edge in u direction */
         extract_edge(eg,fgeom,1,0) ; 
	 /* calculate the lamda matrix and new knot vector in u direction.
	  * k is the number of knots inserted. */ 
	 k = compute_lim(eg,lim,kv,eps) ;
         /* calculate the control points due to degree elevation in u
	  * direction. */
         for (j=0;j<fgeom->vcontpts;j++) {
	     extract_edge(eg,fgeom,1,j) ; 
	     raise_byone(eg,lim,kv,k) ;
	     /* set the control points */
	     store_vertices(eg,fgeom,1,j) ;
	     }
         /* increase the order and the number of control points in u */
         fgeom->uorder += 1 ; 
	 fgeom->ucontpts += k ;
	 /* set the knot vector in u */
         for (j=0;j<fgeom->uorder+fgeom->ucontpts;j++) 
	     fgeom->uknots[j] = kv[j] ;
         }

  /* if elevate degree in v*/
  if (v) 
     for (i=0;i<v;i++) {    /* elevate degree v times */
         /* extract the edge in v direction */
         extract_edge(eg,fgeom,0,0) ; 
         /* calculate the lamda matrix and new knot vector in v direction.
	  * k is the number of knots inserted. */
	 k = compute_lim(eg,lim,kv,eps) ;
	 /* calculate the control points due to degree elevation in v
	  * direction. */
	 for (j=0;j<fgeom->ucontpts;j++) {
	     extract_edge(eg,fgeom,0,j) ; 
	     raise_byone(eg,lim,kv,k) ;
	     /* set the control points */
	     store_vertices(eg,fgeom,0,j) ;
	     }
	 /* increase the order and the number of control points in v */
	 fgeom->vorder += 1 ; 
	 fgeom->vcontpts += k ;
	 /* set the knot vector v */
	 for (j=0;j<fgeom->vorder+fgeom->vcontpts;j++) 
	     fgeom->vknots[j] = kv[j] ;
         }

  /* free memory */
  free_egeom(eg) ;
  free_darray2(lim);
  free_darray1(kv);
}

/***************************************************************************
*                                 store_vertives()
****************************************************************************
* 
* 1   Purpose 
*     This function is a subroutine of degree elevation routine. It assigns the
*     control points of a curve to one row or column of the control points of a
*     surface
* 
* 2   Specification
*     #include "bspl.h"
*     void store_vertices(ParCurv *eg, ParSurf *fgeom, int n0, int k)
* 
* 3   Description
*     This function put the control points of a B-spline curve to one column of
*     the control net of a B-spline surface. The index of the column and the 
*     direction is specified.
* 
* 4   References
*     Not applicable
*  
* 5   Parameters
*     1.ParCurv *eg
*       On entry: the B-spline curve whose control points is to be used.
*     2.ParSurf *fgeom
*       On entry: the B-spline surface whose control points is to be assigned.
*     3.int n0
*       On entry: indicates in which direction the control points are to be 
*                 assigned.
*     4.int k
*       On entry: the index of row or column where the control points are to
*                 be assigned.
* 
* 6   Return Values, Error Indicators and Warnings
*     Not applicable
* 
* 7   Accuracy
*     Not applicable
* 
* 8   Further Comments
*     Not applicable
* 
* 9   Functions referenced by bernrs() are:
*     copyvector()
*     vectalloc()
*  
* 10  Functions that reference bernrs() are: 
*     raise_surf()
* 
***************************************************************************/

void store_vertices(ParCurv *eg, ParSurf *fgeom, int n0, int k)
{
  int i ;

  /* if the direction is u */
  if (n0) 
     for (i=0;i<eg->ncontpts;i++) {
         if (fgeom->contpts[i][k]==NULL) 
	    fgeom->contpts[i][k] = vectalloc() ;
         copyvector(eg->contpts[i],fgeom->contpts[i][k]) ;
         }
  /* if the direction is v */
  else  
     for (i=0;i<eg->ncontpts;i++) {
         if (fgeom->contpts[k][i]==NULL) 
	    fgeom->contpts[k][i] = vectalloc() ;
         copyvector(eg->contpts[i],fgeom->contpts[k][i]) ;
         }
}

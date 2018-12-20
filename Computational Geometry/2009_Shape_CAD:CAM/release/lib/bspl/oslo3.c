/* ***************************************************************************
 Copyright (C) 1996 Massachusetts Institute of Technology all rights reserved 
	Programmer: George A. Kriezis
**************************************************************************** */
# include <stdio.h>
# include "gen.h"
# include "bspl.h"

/**************************************************************************
*                             surfoslo3()
***************************************************************************
* 
* 1   Purpose
*     This function add knots to a NURBS surface according to a new augmented 
*     knot set.
* 
* 2   Specification
*     #include "bspl.h"
*     void surfoslo3(ParSurf *ofgeom, ParSurf *nfgeom, int ncoord)
* 
* 3   Description
*     This subroutine computes the set of control points for a new surface 
*     from the set of control points from a supplied surface and a new knot 
*     vector. The new knot vector is supplied in the allocated surface which 
*     on exit, contains the new surface. It uses the first algorithm given 
*     in Cohen(1980) and follows the basic method outlined in Lyche(1985).
* 
* 4   References
*     [1] E. Cohen, T. Lyche and R. F. Riesenfeld. Discrete B-Splines and 
*         Subdivision Techniques in Computer-Aided Geometric Design and 
* 	Computer Graphics, Computer Graphics and Image Processing, 
* 	14:87-111, 1980
*     [2] T. Lyche, E. Cohen, and K. Morken.  Knot Line Refinement Algorithms 
*         for Tensor Product B-Spline Surfaces, Computer Aided Geometric 
* 	Design, 2:133-139, 1985.
* 
* 5   Parameters
*        1.ParSurf * ofgeom
*          On entry: NURBS surface data structure containing geometry to which 
* 	           algorithm is to be applied.
*        2.ParSurf * nfgeom
*          On entry: NURBS surface data structure containing new augmented knot 
* 	           vector.
*          On exit: the geometry of the surface with knots and control points 
* 	           added.
*        3.int ncoord
*          On entry: the number of coordinates to which knot addition is to be 
* 	           applied.
* 
* 6   Return Values, Error Indicators and Warnings
*     At this time, it is suggested that ncoord = 4.
* 
* 7   Accuracy
*     Not applicable.
* 
* 8   Further Comments
*     Not applicable.
* 
* 9   Functions referenced by surfoslo3() are:
*     soslo()
*  
* 10  Functions that reference surfoslo3() are:
*     addpoints_surf()
*     BsplineToBezierSurf()
*     build_offset()
*     EditKnotsSurf()
*     merge_can_off()
*     patch_with_com_knots()
*     recover_patch()
*     renew_knots()
*     SplitSurf()
*     subbezier()
*     subdiv_bez()
*     SurfRev_to_ParSurf()
*  
*************************************************************************/

void surfoslo3(ParSurf *ofgeom, ParSurf *nfgeom, int ncoord)
{
  /* call soslo() with appropriate variables. */
  soslo(ofgeom->uorder,ofgeom->vorder,ofgeom->ucontpts,ofgeom->vcontpts,
        nfgeom->uorder+nfgeom->ucontpts,nfgeom->vorder+nfgeom->vcontpts,
	ofgeom->uknots,ofgeom->vknots,nfgeom->uknots,nfgeom->vknots,
	ofgeom->contpts,nfgeom->contpts,ncoord);

}

/******************************************************************************
*                             soslo()
*******************************************************************************
* 
* 1   Purpose
*     This function calculates the control points of a new B-spline surface
*     which is obtained by augmenting the knot vector of an old B-spline
*     surface.
* 
* 2   Specification
*     void soslo(int uk, int vk, int ucontpts, int vcontpts, int uq, int vq,
*                double utau[], double vtau[], double ut[], double vt[],
* 	       vector ***inpoly, vector ***outpoly, int ncoord)
* 
* 3   Description
*     This function calculates the control points of a new B-spline surface by
*     transformation in both u,v directions.
*                      P = transpose(Mu) * P0 * Mv
*     where P is the new control points,
*           P0 is the old control points,
* 	  Mu is the transformation matrix in u,
* 	  Mv is the transformation matrix in v.
* 
* 4   References
*     [1] E. Cohen, T. Lyche and R. F. Riesenfeld. Discrete B-Splines and 
*         Subdivision Techniques in Computer-Aided Geometric Design and 
*         Computer Graphics, Computer Graphics and Image Processing, 14:87-111,
*         1980
* 
* 5   Parameters
*     1.int uk
*       On entry: the order in u direction of the original B-spline surface.
*     2.int vk 
*       On entry: the order in v direction of the original B-spline surface.
*     3.int ucontpts
*       On entry: number of control points in u direction of the original
*               surface.
*     4.int vcontpts
*       On entry: number of control points in v direction of the original
*               surface.
*     5.int uq
*       On entry: the number of knots in u direction of the new surface.
*     6.int vq
*       On entry: the number of knots in v direction of the new surface.
*     7.double utau[]
*       On entry: the original knot vector in u direction
*     8.double vtau[]
*       On entry: the origianl knot vector in v direction
*     9.double ut[]
*       On entry: the augamented knot vector in u direction
*     10.double vt[]
*       On entry: the augamented knot vector in v direction
*     11.vector ***inpoly
*       On entry: the control points of the old surface.
*     12 vector ***outpoly
*       On exit : the control points of the new surface
*     13.int ncoord
*        On entry: the number of coordinates to which knot addition is to be 
* 	         applied.
* 		 ncoord = 1,  only x-coordinate.
* 		 ncoord = 2,  only x&y coordinates.
*                  ncoord = 3,  x,y,z coordinates.
*                  ncoord = 4,  all the coordinates.
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
* 9   Functions referenced by soslo() are:
*     dbl_array1()
*     dbl_array2()
*     discrete_bsp()
*     find()
*     free_darray1()
*     free_darray2()
*     free_iarray1()
*     int_array1()
*     vectalloc()
* 
* 10  Functions that reference soslo() are:
*     surfoslo3()
* 
**************************************************************************/

void soslo(int uk, int vk, int ucontpts, int vcontpts, int uq, int vq,
      double utau[], double vtau[], double ut[], double vt[],
      vector ***inpoly, vector ***outpoly, int ncoord)
{
  int k,i,j,l;               /* loop indices */
  int nucontpts,nvcontpts;   /* number of control points in u and v of the new
			      * surface, repectively. */
  int mu,jstart;     /* some indices */
  int ncontpts,nknots,order;  /* to be used as maximum size */
  int *jsarray;      /* array to be used as indices */
  vector *v;         /* used to temporarily contain a control point */    
  double *au;        /* a row of the transformation matrix in u, only
		      * containing nonzero elements */
  double **av;       /* transformation matrix in v, only containing nonzero
		      * elements */
  double **temp;     /* temporary working space */
  double **alpha;    /* to be used to contain a value of a dicrete B-spline */

  nucontpts = uq - uk;   /* number of ctrl pts in u direction of the new 
			  * surface */
  nvcontpts = vq - vk;   /* number of ctrl pts in v direction of the new
			  * surface */

  /* determine the maximum size that the working spaces could be. */
  if (nucontpts >= nvcontpts)
     ncontpts=nucontpts;
  else
     ncontpts=nvcontpts;
  if (uq >= vq)
     nknots=uq;
  else
     nknots=vq;
  if (uk >= vk)
     order=uk;
  else
     order=vk;

  /* allocate memory */
  au = dbl_array1((unsigned)order);
  jsarray=int_array1((unsigned)ncontpts);
  av= dbl_array2((unsigned)ncontpts,(unsigned)order);
  temp = dbl_array2(4,(unsigned)ncontpts);
  alpha = dbl_array2((unsigned)nknots,(unsigned)(order+1));

  /* calculate the transformation matrix in v direction, i.e., the discrete
   * B-spline in v direction. */
  for (i=0; i<nvcontpts; i++)
      {
       /* find index mu such that vtau[mu] <= vt[i] <= vtau[mu+1]. */
       mu = find(vcontpts, vtau, vt[i]);
       /* set all alpha to zero */
       for (k=0; k<nknots; k++)
	    for (j=0; j<order+1; j++)
	        alpha[k][j]=0.0;
       /* evaluate the discrete B-spline at index i. */
       discrete_bsp(mu,i,vk,vt,vtau,alpha);
       /* set the i-th row of the transformation matrix. */       
       jsarray[i] = mu - vk + 1;  /* the index of the first element in i-th
				   * row which is nonzero. */
       for (j=0; j<vk; j++)
           av[i][j] = alpha[j+jsarray[i]][vk];
     }

  /* calculate the transformation matrix in u direction, and finally
   * calculate the control points of the new surface. */
  for (i=0; i<nucontpts; i++)
      {
       /* find index mu such that utau[mu] <= ut[i] <= tau[mu+1]. */
       mu = find(ucontpts, utau, ut[i]);
       /* set all alpha to zero */
       for (k=0; k<nknots; k++)
	   for (j=0; j<order+1; j++)
	       alpha[k][j]=0.0;
       /* evaluate the discrete B-spline at index i. */
       discrete_bsp(mu, i, uk, ut, utau, alpha);
       /* set the i-th column of the transformation matrix in u. */
       jstart = mu - uk + 1;   /* the index of the first element in the
				* current row which is nonzero. */
       for (j=0; j<uk; j++)
           au[j] = alpha[j+jstart][uk];

       /* pre-multiply the control points by the transformation matrix in v and
        * put the result in matrix temp. */ 
       for (j=0; j<vcontpts; j++)
           {
	    for (l=0; l<4; l++)
	        temp[l][j] = 0.0;

            for (l=0; l<uk; l++)
	        {
	         if (ncoord >= 1)
	            temp[0][j] += au[l]*inpoly[l+jstart][j]->x;
	         if (ncoord >= 2)
	            temp[1][j] += au[l]*inpoly[l+jstart][j]->y;
	         if (ncoord >= 3)
	            temp[2][j] += au[l]*inpoly[l+jstart][j]->z;
	         if (ncoord >= 4)
	            temp[3][j] += au[l]*inpoly[l+jstart][j]->w;
	        }
           }
    
      /* multiply the above result by the transformation matrix in v, and get
       * the control points. */
      for (j=0; j<nvcontpts; j++)
          {
	   /* initialize */
	   if (outpoly[i][j] == NULL)
	      outpoly[i][j] = vectalloc();
	      v = outpoly[i][j];
	      v->x = 0.0;
	      v->y = 0.0;
	      v->z = 0.0;
	      v->w = 0.0;
       
           for (l=0; l<vk; l++)
	       {
	        if (ncoord >= 1)
	           v->x += temp[0][l+jsarray[j]]*av[j][l];
	        if (ncoord >= 2)
	           v->y += temp[1][l+jsarray[j]]*av[j][l];
	        if (ncoord >= 3)
	           v->z += temp[2][l+jsarray[j]]*av[j][l];
	        if (ncoord >= 4)
	           v->w += temp[3][l+jsarray[j]]*av[j][l];
	       }
          }
  }

  /* free memory */
  free_iarray1(jsarray);
  free_darray1(au);
  free_darray2(av);
  free_darray2(temp);
  free_darray2(alpha);
}

/*
 *main(){
 *surfaceoslotest();
 *}
 *surfaceoslotest()
 *{
 *struct face *face;
 *struct fgeom *ofgeom, *nfgeom;
 *
 *face = readface(NULL);
 *ofgeom = face->geom;
 *face = readface(NULL);
 *nfgeom = face->geom;
 *writefgeom(stdout, ofgeom);
 *writefgeom(stdout, nfgeom);
 *
 *hsoslo1(ofgeom,nfgeom,4);
 *
 *writefgeom(stdout, nfgeom);
 *}
 */
/*
 *double cinpoly[20][3] = {{0.0,0.0,0.0},
 *		      {100.0,100.0,100.0},
 *	           	{200.0,200.0,200.0},
 *			{300.0,300.0,300.0},
 *			{400.0,400.0,400.0},
 *			{500.0,500.0,500.0}};
 *
 *double tau[10]={0.0,0.0,0.0,0.0,1.0,2.0,3.0,3.0,3.0,3.0};
 *double t[11]={0.0,0.0,0.0,0.0,1.0,1.3,2.0,3.0,3.0,3.0,3.0};
 *
 *int k = 4;
 *int n = 6;
 *int q = 11;
 */
/*
 *curveoslotest()
 *{
 *int i;
 *double outpoly[30][3];
 *printf("enter k, n, q\n");
 *printf("order = %d ncontpts  = %d  n0: of newknots  = %d\n",k,n,q);
 *
 *printf("original knots\n");
 *for(i=0; i<n+k; i++)
 *printf("%lf ", tau[i]);
 *printf("\n");
 *
 *for(i=0; i<n; i++)
 *printf("original contpts = %d %lf %lf %lf\n", i+1,cinpoly[i][0],
 *       cinpoly[i][1],cinpoly[i][2]);
 *
 *curve_oslo1(k,n,q,tau,t,cinpoly,outpoly);
 *
 *printf("final knots\n");
 *for(i=0; i<q; i++)
 *printf("%lf ", t[i]);
 *printf("\n");
 *
 *for(i=0; i<q-k; i++)
 *printf("final contpts = %d %lf %lf %lf\n",i+1, outpoly[i][0],outpoly[i][1],
 *       outpoly[i][2]);
 *}
 *
*/









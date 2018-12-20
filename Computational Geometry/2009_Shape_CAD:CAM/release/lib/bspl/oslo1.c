/* ***************************************************************************
 Copyright (C) 1996 Massachusetts Institute of Technology all rights reserved 
	Programmer: George A. Kriezis
**************************************************************************** */
# include "gen.h"
# include "bspl.h"

/* needs matrix.o for linking */

/****************************************************************************
*                                curve_oslo1()
*****************************************************************************
* 
* 1   Purpose
*     This function converts a NURBS curve to another one with a new given 
*      knot vector using knot insertion.
* 
* 2   Specification
*     #include "bspl.h"
*     void curve_oslo1(int order, int ncontpts, int q, double *tau, double *t, 
*                      double **inpoly, double **outpoly)
* 
* 3   Description
*     This subroutine returns a new set of control points corresponding to an 
*     augmented knot vector for a NURBS curve that is geometrically identical 
*     to a NURBS defined by a set of control points and a corresponding knot 
*     vector. For a description of the algorithm, see the reference list.
* 
* 4   References
*     [1] E. Cohen, T. Lyche and R. F. Riesenfeld. Discrete B-Splines and 
*         Subdivision Techniques in Computer-Aided Geometric Design and 
* 	  Computer Graphics, Computer Graphics and Image Processing, 14:87-111,
*         1980
* 
* 5   Parameters
*        1.int order
*          On entry: order of the original NURBS curve.
*        2.int ncontpts
*          On entry: the number of control points in the original NURBS curve.
*        3.int q
*          On entry: the number of elements in the new augmented knot vector.
*        4.double * tau
*          On entry: array of length order + ncontpts containing the original 
* 	           knot vector. In this array, tau[0] is the first knot and 
* 		   tau[order + ncontpts - 1] is the last knot.
*        5.double * t
*          On entry: array of length q containing the new augmented knot
*                  vector.
* 	           In this array, tau[0] is the first knot and tau[q - 1] is 
* 		   the last knot.
*        6.double ** inpoly
*          On entry: two dimensional array length ncontpts x 4 containing the 
* 	           original control points.  In this matrix, the element 
* 		   inpoly[i][j] refers to the ith control point i = 0,...,
* 		   ncontpts - 1 and jth coordinate where 0 = x, 1 = y, 2 = z, 
* 		   and 3 defines the homogeneous coordinate.
*        7.double ** outpoly
*          On entry: two dimensional array length q - order x 4.
*          On exit: contains the new set of control points. In this matrix, 
* 	           the element outpoly[i][j] refers to the ith control point
* 		   i = 0, ...,q - order - 1 and jth coordinate where 0 = x,
* 		   1 = y,  2 = z, and 3 defines the homogeneous coordinate.
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
* 9   Functions referenced by curve_oslo1() are:
*     dbl_array2()
*     free_darray2()
*     matrixmu()
*     transform()
* 
* 10  Functions that reference curve_oslo1() are:
*     addpoints()
*     BsplineToBezier()
*     CurveOrientation()
*     EditKnotsCos()
*     EditKnotsCurv()
*     offset_curve()
*     recover_cv()
*     split_bspl()
*     SubdivideCos()
*     SubdivideCurv()
*     subdiv_cv()
* 
*****************************************************************************/

void curve_oslo1(int order, int ncontpts, int q, double tau[], double t[],
		 double **inpoly, double **outpoly)

/* int order;     	 order of the curve */
/* int ncontpts;  	 number of control points  */
/* int q;       	 number of knots in the new knot vector */
/* double **inpoly;   	 array of original control points, [ncontpts][3] */
/* double tau[];         vector of original knots, [ncontpts+order] */
/* double t[];         	 vector of new augmented knots, [q] */
/* double **outpoly;  	 array of new control points,  [q-order] */
{
  int i,j,k;
  double **au;          /* matrix containing the values of the corresponding
			 * discrete B-spline. */
  
  /* allocate memory */
  au = dbl_array2((unsigned)(q-order),(unsigned)q);

  /* evaluate the discrete B-spline at t[j], 0 <= j <= q. */
  transform(order,ncontpts,q, tau,t, au);

  /* calculate the control points corresponding to the augamented knot
   * vector. */
  matrixmu(q-order,ncontpts,4, q,4,4, &au[0][0],&inpoly[0][0],&outpoly[0][0]);

  /* free memory */
  free_darray2(au);
}

/**************************************************************************/
/* not needed anymore ****
 * void surface_oslo1(uorder,vorder,ucontpts,vcontpts,uq,vq,utau,vtau,ut,vt,
 * 		inpoly,d1,outpoly,d2)
 * 
 * int d1, d2;		 middle dimension of inpoly and outpoly 
 * int uorder,vorder;     	 order of the curve in u and v 
 * int ucontpts,vcontpts;   number of control points in u and v 
 * int uq,vq;       	 number of knots in the new knot vector 
 * double inpoly[];    	 array of original control points, dim [][d1][4] 
 * double utau[],vtau[];    vector of original knots,  
 * double ut[],vt[];        vector of new augmented knots, (uq) (vq)
 * double outpoly[];   	 array of new control points,  dim [][d2][4] 
 * 
 * {
 * int i,j,k;
 * double **au,**av,**output,**avt,**temp;
 * au = dbl_array2((unsigned)(uq-uorder),(unsigned)uq);
 * av = dbl_array2((unsigned)(vq-vorder),(unsigned)vq);
 * output = dbl_array2((unsigned)(uq-uorder),(unsigned)(vq-vorder));
 * avt = dbl_array2((unsigned)vq,(unsigned)(vq-vorder));
 * temp = dbl_array2((unsigned)(uq-uorder),(unsigned)(vq-vorder));
 * 
 * transform(uorder,ucontpts,uq, utau,ut, au);
 * transform(vorder,vcontpts,vq, vtau,vt, av);
 * 
 * for(k=0; k<4; k++)
 * 	{
 * 	for(i=0; i<ucontpts; i++)
 * 		for(j=0; j<vcontpts; j++)
 * 			output[i][j] = inpoly[i*d1*4 + j*4 + k];
 * 
 * matrixmu(uq-uorder,ucontpts,vcontpts,uq,vq-vorder,uq-uorder,au,output,temp);
 * matrix_transpose(vq-vorder,vcontpts,av,vq,avt,vq-vorder);
 * matrixmu(uq-uorder,vcontpts,vq-vorder,NCONTPTS,vq-vorder,,temp,avt,output);
 * 
 * 	for(i=0; i<uq - uorder; i++)
 * 		for(j=0; j<vq - vorder; j++)
 * 			outpoly[i*d2*4 + j*4 + k] = output[i][j];
 * 	}
 * free_darray2(au);
 * free_darray2(av);
 * free_darray2(avt);
 * free_darray2(output);
 * free_darray2(temp);
 *}
 */	

/***************************************************************************
*                                transform()
****************************************************************************
*
* 1   Purpose
*     This function calculates the matrix of discrete B-spline, which is to be
*     used for calculating the control points after the augamentation of knot
*     vector.
*
* 2   Specification
*     void transform(int k, int n, int q, double tau[], double t[], double **a)
* 
* 3   Description
*     This function first calculates the discrete B-spline values for each knot
*     in the augamented knot vector, and then constructs the matrix of discrete
*     B-spline.
* 
* 4   References 
*     [1] E. Cohen, T. Lyche and R. F. Riesenfeld. Discrete B-Splines and 
*         Subdivision Techniques in Computer-Aided Geometric Design and 
*         Computer Graphics, Computer Graphics and Image Processing, 14:87-111,
*         1980
* 
* 5   Parameters
*     1.int k
*       On entry: the order of the original B-spline curve.
*     2.int n
*       On entry: the length of the original knot vector
*     3.int q 
*       On entry: the length of the augamented knot vector
*     4.double tau[]
*       On entry: the array containing the original knot vector
*     5.double t[]
*       On entry: the array containing the augamented knot vector
*     6.double **a
*       On exit:  the matrix containing the discrete B-spline
* 
* 6   Return Values, Error Indicators and Warnings
*     Not applicable.
* 
* 7   Accuracy
*     Not applicable.
* 
* 8   Further Comments
*     This is a subroutine of function curve_oslo1(), see also curve_oslo1().
* 
* 9   Functions referenced by transform() are:
*     dbl_array2()
*     discrete_bsp()
*     find()
*     free_darray2()
* 
* 10  Functions that reference transform() are:
*     curve_oslo1()
* 
****************************************************************************/

void transform(int k, int n, int q, double tau[], double t[], double **a)
{
  int i,j,l,mu;         /* indices */
  double **alpha;       /* the matrix of discrete B-spline for t[j], each 
			 * column corresponds to an order <= order k. */
				      
  /* allocate memory */
  alpha = dbl_array2((unsigned)q,(unsigned)(k+1));

  /* initialize the discrete B-spline matrix. */
  for (i=0; i<q-k; i++)
      for (j=0; j<q; j++)
	  a[i][j] = 0.0;

  /* the discrete B-spline matrix has (q-k) rows, which is the number of
   * control points of the new curve. */
  for (j=0; j<q-k; j++)
      {
       /* find the index mu such that tau[mu] <= t[j] <= tau[mu+1]. */
       mu = find(n, tau, t[j]); 
       /* set all alpha to zero */
       for (i=0; i<q; i++)
	   for (l=0; l<k+1; l++)
	       alpha[i][l]=0.0;
       /* calculate the alpha matrix for t[j] */
       discrete_bsp(mu, j, k, t, tau, alpha);
       /* assign the j-th row of the a matrix, in which only (mu-k+1)-th to
	* (mu)-th elements are nonzero. */
       for (i=mu-k+1; i<=mu; i++)
	   a[j][i] = alpha[i][k];
       }
  
  /* free memory */
  free_darray2(alpha);
}

/*******************************************************************/
/*
 * static find ( kn, tau, t1)
 * int  kn;
 * double tau[],t1;
 * {
 * int i,mu;
 *
 * for(i=0; i < kn; i++)
 * if(t1>=tau[i])
 * 	mu = i;
 * 
 * return(mu);
 * }
 */

/****************************************************************************
*                                 discrete_bsp()
*****************************************************************************
* 
* 1   Purpose
*     This function evaluates a discrete B-spline.
* 
* 2   Specification 
*     void discrete_bsp(int mu, int j, int k, double t[], double tau[],
* 		      double **alpha)
* 
* 3   Description
*     This function uses the iteration algorithm proposed in the following 
*     reference to evaluate a discrete B-spline at a fixed index j.
* 
* 4   References
*     [1] E. Cohen, T. Lyche and R. F. Riesenfeld. Discrete B-Splines and 
*         Subdivision Techniques in Computer-Aided Geometric Design and 
*         Computer Graphics, Computer Graphics and Image Processing, 14:87-111,
*         1980
* 
* 5   Parameters
*     1.int mu
*       On entry: the index of the knot in the original knot vector, such that
*                 for a given index j, tau[mu] <= t[j] <= tau[mu+1].
*     2.int j
*       On entry: the index where the value of the discrete B-spline is to be
*                 evaluated.
*     3.int k
*       On entry: the order of the original B-spline curve
*     4.double t[]
*       On entry: the array containing the augamented knot vector.
*     5.double tau[]
*       On entry: the array containing the original knot vector.
*     6.double **alpha
*       On exit:  the matrix containing the values of the discrete B-spline at
*                 index j.
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
* 9   Functions referenced by discrete_bsp() are: None
* 
* 10  Functions that reference discrete_bsp() are:
*     soslo()
*     transform()
* 
*****************************************************************************/

void discrete_bsp(int mu, int j, int k, double t[], double tau[],
		  double **alpha)
{
  int r, i, mu2;      /* r -- order of the discrete B-spline
                       * i -- loop index
		       * mu2 -- index */
  double beta_1, beta, d1, d2, tj;       /* working variables */

  /* initialization */
  alpha[mu][1] = 1.0;
  mu2 = mu;

  /* each column of matrix alpha corresonding to the value fo the discrete
   * B-spline with order r+1, which varies from 1 to k. */ 
  for (r=1; r<k; r++)
      {
       /* boundary */
       beta_1 = 0.0;
       /* evaluate by iteration algorithm, see the reference. */
       tj = t[j+r];
       for (i=mu2; i<=mu; i++)
	   {
	    d1 = tj - tau[i];
	    d2 = tau[i+r] - tj;
	    beta = alpha[i][r]/ (d1 + d2);
	    alpha[i-1][r+1] = d2*beta + beta_1;
	    beta_1 = d1*beta;
	   }
       alpha[mu][r+1] = beta_1;
       mu2 = mu2 - 1;
      }
}	

/*   surface test 1   
 * 
 * double inpoly[2][2][4] = { {0.0,0.0,0.0,1.0}, {0.0,100.0,0.0,1.0},
 * 		      {100.0,0.0,0.0,1.0},{100.0, 100.0, 0.0,1.0}};
 * double  utau[4]= {0.0,0.0,1.0,1.0};
 * double vtau[4]= {0.0,0.0,1.0,1.0};
 * double ut[4]= {0.0,0.0,1.0,1.0};
 * double vt[5]= {0.0,0.0,0.5,1.0,1.0};
 * 
 * int uk = 2;
 * int vk = 2;
 * int un = 2;
 * int vn = 2;
 * int uq = 4;
 * int vq = 5;
 */

/*  surface test 2  
 * 
 * double inpoly[4][4][4] = { { {0.,0.,0.},{0.,25.,50.},{0.,75.,50.},
 *                              {0.,100.,0.} },
 * 		 { {25.,0.,0.},{25.,25.,50.},{25.,75.,50.},{25.,100.,0.} },
 * 		 { {75.,0.,0.},{75.,25.,50.},{75.,75.,50.},{75.,100.,0.} },
 * 		 { {100.,0.,0.},{100.,25.,50.},{100.,75.,50.},{100.,100.,0.}}};
 * double  utau[8]={0.0,0.0,0.0,0.0,3.0,3.0,3.0,3.0};
 * double vtau[8]={0.0,0.0,0.0,0.0,3.0,3.0,3.0,3.0};
 * double  ut[9]={0.0,0.0,0.0,0.0,1.5,3.0,3.0,3.0,3.0};
 * double vt[9]={0.0,0.0,0.0,0.0,1.5,3.0,3.0,3.0,3.0};
 * 
 * int uk = 4;
 * int vk = 4;
 * int un = 4;
 * int vn = 4;
 * int uq = 9;
 * int vq = 9;
 */
/*
 * surfaceoslotest()
 * {
 * double outpoly[5][5][4];
 * 
 * surface_oslo1(uk,vk,un,vn,uq,vq,utau,vtau,ut,vt,inpoly,2,outpoly,5);
 * 
 * }
 */
/*
 * double cinpoly[20][3] = {{0.0,0.0,0.0},
 * 		      {100.0,100.0,100.0},
 * 	           	{200.0,200.0,200.0},
 * 			{300.0,300.0,300.0},
 * 			{400.0,400.0,400.0},
 * 			{500.0,500.0,500.0}};
 * 
 * double tau[10]={0.0,0.0,0.0,0.0,1.0,2.0,3.0,3.0,3.0,3.0};
 * double t[11]={0.0,0.0,0.0,0.0,1.0,1.3,2.0,3.0,3.0,3.0,3.0};
 * 
 * int k = 4;
 * int n = 6;
 * int q = 11;
 */
/*
 * curveoslotest(void)
 * {
 * int i;
 * double outpoly[30][3];
 * printf("enter k, n, q\n");
 * printf("order = %d ncontpts  = %d  n0: of newknots  = %d\n",k,n,q);
 * 
 * printf("original knots\n");
 * for(i=0; i<n+k; i++)
 * printf("%lf ", tau[i]);
 * printf("\n");
 * 
 * for(i=0; i<n; i++)
 * printf("original contpts = %d %lf %lf %lf\n", i+1,cinpoly[i][0],
 *        cinpoly[i][1],cinpoly[i][2]);
 * 
 * curve_oslo1(k,n,q,tau,t,cinpoly,outpoly);
 * 
 * printf("final knots\n");
 * for(i=0; i<q; i++)
 * printf("%lf ", t[i]);
 * printf("\n");
 * 
 * for(i=0; i<q-k; i++)
 * printf("final contpts = %d %lf %lf %lf\n",i+1, outpoly[i][0],outpoly[i][1],
 *        outpoly[i][2]);
 * }
 * 
*/

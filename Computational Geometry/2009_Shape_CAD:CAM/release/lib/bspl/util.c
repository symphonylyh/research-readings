/* ***************************************************************************
 Copyright (C) 1996 Massachusetts Institute of Technology all rights reserved 
	Programmer: George A. Kriezis
**************************************************************************** */
#include "gen.h"
#include "bspl.h"

/*****************************************************************************
*                               nodes_bsp()
******************************************************************************
* 
* 1   Purpose
*     This function calculates the node values of a given knot vector.
* 
* 2   Specification
*     #include "bspl.h"
*     void nodes_bsp(double *knot, int number, int order, double *node)
* 
* 3   Description
*     This routine evaluates the nodes corresponding to the given knot vector 
*     using the relation:
*                                   1    
*                     node[i]= ----------- sum(k=i..(i+order-1)) knot[k]
*                               order - 1
* 
* 4   References
*     [1] P. G. Alourdas. Shape Creation, Interrogation and Fairing Using 
*         B-Splines, Engineer's Thesis, Massachusetts Institute of Technology, 
* 	Department of Ocean Engineering, Cambridge, Massachusetts, 1989.
* 
* 5   Parameters
*        1.double * knot
*          On entry: array of length >= number + order containing the knot
*                    vector.
*        2.int number
*          On entry: number of control points or weights corresponding to knot 
* 	             vector.
*        3.int order
*          On entry: order of B-spline corresponding to the knot vector.
*        4.double * node
*          On entry: array of length  n.
*          On exit: elemental values are the computed nodes.
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
* 9   Functions referenced by nodes_bsp() are: None
* 
* 10  Functions that reference nodes_bsp() are:
*     approx_fnbc()
*     approx_fnbc_knot()
*     avg_arc_length()
*     BlendSurfCB()
*     build_offset()
*     build_offsets_int()
*     calc_knots_hartley()
*     camber_surface()
*     convexhull_test()
*     fit_curve()
*     integral_build_offset()
*     opt_param()
*     ParSurf_approx()
*     PowCurv_to_ParCurv()
*     PowSurf_to_ParSurf_int()
*     PowSurf_to_ParSurf_loft()
*     sample_blend()
*     sample_surf_test()
* 
******************************************************************************/
void nodes_bsp(double knot[], int number, int order, double node[])
{
  int i,k;
  for (i=0; i<number; i++){
      node[i]=0.0;
      for (k=1; k<order; k++)
          node[i] += knot[i+k];
      node[i] /= (double)(order-1);
      }
}

/******************************************************************************
*                           bsp_basis()
*******************************************************************************
* 
* 1   Purpose
*     This function computes the value of a B-spline basis function.
* 
* 2   Specification
*     #include ""
*     double2 bsp_basis(int i, int k, double u_val, double *knots)
* 
* 3   Description
*     This routine evaluates the ith B-spline basis function of order k at 
*     t = u_val based on the knot vector {t[i]} defined as:
* 
*       N(i,1) = 1 if t[i]<=t<t[i+1], 
*              = 0 otherwise.
* 
*                      t - t[i]                     t[i+k]- t
*       N[i,k](t) =----------------N[i,k-1](t) + ----------------N[i+1,k-1](t)
*                   t[i+k-1]- t[i]                t[i+k]- t[i+1]
* 
* 	   if k>1
* 
*     The value of the basis function is returned from the routine.
* 
* 4   References
*     [1]C. De Boor. A Practical Guide to Splines, Springer, New York, 1978.
* 
* 5   Parameters
*       1.int i
*         On entry: identifies the B-spline basis function to evaluate
*       2.int k
*         On entry: order of B-spline
*       3.double u_val
*         On entry: parameter value at which the B-spline basis is evaluated
*       4.double * knots
*         On entry: array containing the B-spline knot vector elements
* 
* 6   Return Values, Error Indicators and Warnings
*     Not applicable.
* 
* 7   Accuracy
*     Not applicable.
* 
* 8   Further Comments
*     This routine checks for u_val between knots[i] and knots[i + k - 1]. 
*     If u_val not in this range, then it returns zero.
* 
* 9   Functions referenced by bsp_basis() are: None
* 
* 10  Functions that reference bsp_basis() are:
*     bsp_basis()
* 
******************************************************************************/

double bsp_basis(int i, int k, double u_val, double knots[])
{
  int i1, k1;
  double v, v1, d_n;

  if (k==1)  {
     v = 0.0;
     if ((knots[i]-1.e-12) <= u_val && u_val <= knots[i+1])
        v = 1.0;
     }
  else  {
     v = 0.0;
     d_n = knots[i+k-1] - knots[i];
     if (d_n != 0.0)  {
        v1 = (u_val - knots[i]) / d_n;
        k1 = k - 1;
        v = v1 * bsp_basis(i,k1,u_val,knots);
        }
     d_n = knots[i+k] - knots[i+1];
     if(d_n != 0.0)  {
        v1 = (knots[i+k] - u_val) / d_n;
        i1 = i + 1;
        k1 = k - 1;
        v = v + v1 * bsp_basis(i1,k1,u_val,knots);
       }
    }

  return(v);
}  /*  end of bsp_basis  */

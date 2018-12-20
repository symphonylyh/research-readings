/***************************************************************************
 *                                                                         *
                           Copyright (C) 1992 by
            Massachusetts Institute of Technology, Cambridge, MA
                             All Rights Reserved
 *                                                                         *
 **************************************************************************/
/*  Filename:  blnd_curv.c

    Written by:  Peter C. Filkins
    Date:        21 January 1991
===========================================================================
    File Modification History:
           13 February 1991 - changed to calculate first and second
                              fundamental forms separately
	   12 March 1991    - added computation of principal curvatures.
	                      Does not recompute surface information like
			      curvature1() (Note: not rqd for blending)
===========================================================================
    Description:  computes the normal and principal curvatures of 
                  rational surfaces.
===========================================================================
    Subroutines called:       dot()
    Library dependencies:     libgen.a
===========================================================================
    Arguments:  vector *r[]   an array of vectors through 2nd derivative
                              on the given surface
		vector *n     surface normal
		double *dir   tangent direction of normal section expressed
		              in terms of the surface parameters (du,dv)
===========================================================================
    Variable list:  E, F, G   Elements of First Fundamental Form of Surface
                    L, M, N   Elements of Second Fundamental Form
===========================================================================*/
#include <stdio.h>
#include <math.h>
#include "gen.h"
#include "editor.h"

double *fundamentals(vector *n, vector **r)
{
     double *f;

     f = dbl_array1(6);

     f[0] = dot( r[1], r[1] );           /* E */
     f[1] = dot( r[1], r[2] );           /* F */
     f[2] = dot( r[2], r[2] );           /* G */

     f[3] = dot( n, r[4] );              /* L */
     f[4] = dot( n, r[3] );              /* M */
     f[5] = dot( n, r[5] );              /* N */

     return (f);
   }


double normal_curvature(double *f, double *dir)
{
     double du2, dv2, dudv;
     double fund1, fund2;

     du2 = dir[0] * dir[0];
     dudv = dir[0] * dir[1];
     dv2 = dir[1] * dir[1];

     fund1 = f[0]*du2 + 2.0*f[1]*dudv + f[2]*dv2;
     fund2 = f[3]*du2 + 2.0*f[4]*dudv + f[5]*dv2;

     return (fund2/fund1);
   }


void principal_curv(double *f, double *k)
{
     double Gdet, Ddet, H, K, tempmod, root;

     Gdet = f[0]*f[2] - f[1]*f[1];
     Ddet = f[3]*f[5] - f[4]*f[4];

     H = ((2.0*f[1]*f[4])-(f[0]*f[5] + f[2]*f[3])) / (2.0*Gdet);
     K = Ddet / Gdet;

     tempmod = (H * H) - K;
     tempmod = fabs( tempmod );

     if ( tempmod >= 0.0 && tempmod < ZERO )  {
       k[0] = H;
       k[1] = H;
     }
     else if ( tempmod >= ZERO )  {
       root = sqrt( tempmod );
       k[0] = H + root;
       k[1] = H - root;
     }
     else if (tempmod < 0.0)
       printf("H*H < K in principal_curv()\n");
     return;
   }

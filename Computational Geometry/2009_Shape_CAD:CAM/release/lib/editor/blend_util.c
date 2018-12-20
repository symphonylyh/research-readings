/*************************************************************************
 *                                                                       *
                         Copyright (C) 1992 by
         Massachusetts Institute of Technology, Cambridge, MA
                          All Rights Reserved
 *                                                                       *
 ************************************************************************/
/*  util_blnd.c */
/*  A file for utilities to support blend operations

    Written by:  Peter C. Filkins
    Date:        05 November 1990

    Modified:    30 January 1991 - added new functions and rewrote 
                                   existing function
===========================================================================
    Functions:
	 surf_data()   evaluates the position, first and second parametric
	               derivatives on the surface at a specified point on
		       the linkage curve lying on the surface
	 surf_normal() computes the normal vector to the surface.  Function
	               is an alternative to the one in the bspl lib in that
		       it does not recompute the first parametric deriv but
		       takes them as arguments.  Handles degenerate patch
		       per Faux and Pratt, p. 236.
**************************************************************************/

#include <math.h>
#include <malloc.h>
#include "gen.h"
#include "bspl.h"
#include "editor.h"

double vec_angle(vector *v1, vector *v2)
{
  double angle, u, v;

  u = mag(v1);
  v = mag(v2);
  if ((u<1.0e-09) && (v<1.0e-09)) 
    angle = 0.0;
  else
    angle = acos(dot(v1,v2)/(u*v));

  return (angle);
}

vector *surf_data(ParCurv *egm, ParSurf *fgm, double u, vector **v)
{
  vector *t, *n;

  t = rbspeval(egm, u, 0);
  evalrsurf(fgm, t->x, t->y, 5, v);
  vectfree(t);
 
  n = surf_normal(v[1], v[2], v[3]);

  return (n);
}

vector *surf_normal(vector *ru, vector *rv, vector *ruv)
{
  vector *n;

  if (mag(ru) < ZERO)
    n = cross(ruv, rv);
  else if (mag(rv) < ZERO)
    n = cross(ru, ruv);
  else
    n = cross(ru, rv);
  
  unitvector1(n, n);

  return(n);
}

ParCurv *link_alloc(short bc1, short bc2, short *m)
{
  ParCurv *egm; 
  short j, order, ncontpts, test;
  
  test = bc1 + bc2;
  spline_par(test, &order, &ncontpts);
     
  /* allocate memory for cross-link curve  */
  egm = egeomalloc1(order, ncontpts);

  /* endpoint knot multiplicity */
  for (j=0; j<order; j++)  {
    egm->knots[j] = 0.0;
    egm->knots[j+ncontpts] = 1.0;
  }
  /* initialize remaining knots */
  /* tangent at surf 1, curvature at surf 2 */
  if (bc1 == 1 && bc2 == 2)  egm->knots[4] = 2.0/3.0;
                                 /* tangent at surf 2, curvature at surf 1 */
  if (bc1 == 2 && bc2 == 1)  egm->knots[4] = 1.0/3.0;

  if (bc1 == 2 && bc2 == 2)  {
    egm->knots[4] = 1.0/3.0;
    egm->knots[5] = 2.0/3.0;
  }
  *m = ncontpts;
  return ( egm );
}

void spline_par(short test, short *ord, short *n)
{
  if (test == 4) {
    *ord = 4;
    *n = 6;
  }
  else if (test == 3) {
    *ord = 4;
    *n = 5;
  }
  else if (test == 2) {
    *ord = 4;
    *n = 4;
  }
  else if (test == 1) {
    *ord = 3;
    *n = 3;
  }
  else {
    *ord = 2;
    *n = 2;
  }
}

double bias_factor1(short f, double t, double surf1, double surf2 )
{       /* distance function */
  double b;

  if (f == 0)
    b = t / surf1;
  else 
    b = t / surf2;
  
  return (b);
}

double bias_factor(short f, double t, double surf1, double surf2)
{       /*  constant value bias */
  double b;

  if (f == 0)
    b = surf1;
  else
    b = surf2;
  
  return (b);
}

vector *bias( vector *t, vector *v1, vector *v2, double scale)
{
/* vector *t        tangent */
/* vector *v1, *v2  position on two linkage curves */
/* double scale     bias factor */
  vector *p, *t_minus;
  double angle1, angle2;

  p = sub_vect(v2, v1);
  t_minus = scale_vect(-1.0, t);
  angle1 = vec_angle(t, p );
  angle2 = vec_angle(t_minus, p);
  
  if (angle2 < angle1)
    scale = -scale;
  
  scale_vect1(scale, t, p);
  vectfree(t_minus);
  
  return( p );
}

/***************************************************************************
 *                                                                         *
                           Copyright (C) 1991 by
            Massachusetts Institute of Technology, Cambridge, MA
                             All Rights Reserved
 *                                                                         *
 **************************************************************************/
/*  Filename:  curv_cont.c

    Written by:  Seamus T. Tuohy (Peter C. Filkins)
    Date:        23 January 1991
===========================================================================
    File Modification History:
           13 February 1991 - modified curvature functions
	   08 March 1991    - revised to return normal curvature to
	                      to calling function for use in sampling
			      blend surface
	   12 March 1991    - revised to return principal curvatures 
	                      for sampling instead of normal
--------------------- Seamus takes over ----------------------------------
           24 July  1991    - made a correction to the calc of x2l and y2l
                              also made correction so that local coords
			      are rotated into different planes
===========================================================================
    Description: computes the value of the second derivative at the 
                 endpoints of a rational B-spline based on the normal
		 curvature of the surface at that point.  The second
		 derivative is computed in a local frame and then 
		 transformed to the global coordinates of the object.          
===========================================================================
    Subroutines called:    partan_dir(), fundamentals(), normal_curvature()
                           mag(), vector_dircos()
    Library dependencies:  libgen.a, libeditor.a
===========================================================================
    Arguments:  vector *v[]     surface data through second derivative
                vector *n       normal to the surface
		vector *tan_v[] vector pair defining the tangent plane
		vector *qtan    biased tangent vector in normal section
=========================================================================*/

vector *curv_cont(vector **v, vector *n, vector **tan_v, vector *qtan)
{
  vector *sec_deriv;
  double *dudv, kn, qt, qt3;
  double x2l, y2l;
  double *f;

  /*  find the direction of the normal section in terms of surf parameters */
  dudv = partan_dir(tan_v[0], v[1], v[2]);

  /*  compute elements of first and second fundamental forms  */
  f = fundamentals(n, v);

  /*  compute the normal curvature of the section */  
  kn = normal_curvature(f, dudv);
  free_darray1(f);
  free_darray1(dudv);

  qt = mag(qtan);
  qt3 = qt * qt * qt;

  /*  find the components of the 2nd deriv in local coord system */
  x2l = -qt * kn * qtan->y;
  y2l = qt * kn * qtan->x; 

  /*  transform to global coord system using direction cosines  */
  sec_deriv = vec_rot_fix(x2l, y2l, tan_v[0], n); 

  return (sec_deriv);
}

vector *vec_rot(double x2, double y2, vector **vec)
/* double x2, y2  second derivative values in local coord */
/* vector *vec[]  vectors defining tangent plane */
{
  vector *d1, *ans;

  ans = vector_dircos(vec[0], 0);
  d1 = copyvector(ans, 0);
  vector_dircos(vec[1], ans);

  ans->x = (d1->x * x2) + (ans->x * y2);
  ans->y = (d1->y * x2) + (ans->y * y2);     
  ans->z = (d1->z * x2) + (ans->z * y2);
  ans->w = 1.0;

  vectfree(d1);
  return (ans);
}

vector *vec_rot_fix(double x2, double y2, vector *v0, vector *v1)
/* double x2, y2   second derivative values in local coord */
/* vector *v0,*v1  vectors defining tangent plane */
{
  vector *d1, *ans;

  ans = vector_dircos(v0, 0);
  d1 = copyvector(ans, 0);
  vector_dircos(v1, ans);

  ans->x = (d1->x * x2) + (ans->x * y2);
  ans->y = (d1->y * x2) + (ans->y * y2);     
  ans->z = (d1->z * x2) + (ans->z * y2);
  ans->w = 1.0;

  vectfree(d1);
  return (ans);
}

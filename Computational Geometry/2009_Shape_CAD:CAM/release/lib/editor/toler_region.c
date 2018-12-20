/************************************************************************
 *									*
			Copyright (C) 1992 by
	Massachusetts Institute of Technology, Cambridge, MA
			 All rights reserved
 *									*
 ************************************************************************/

#include <stdio.h>
#include "gen.h"
#include "editor.h"

#define MAXSECT 500

short tolerance_region(ParSurf *plate, ParSurf **offset, ParSurf **canals,
		       ParSurf **corners, double tz,  int *nsec, int *nofs,
		       double eps, double fc, int *ncan, int *noff, int *nm)

/* Build the tolerance region around a deformed plate, represented by
   the rational B-spline surface patch plate.
   tz is the value of the tolerance.
   offset are the approximating B-splines to plate's offsets.
   offset[0] = "upper" offset (in the direction to plate's normal).
   offset[1] = "lower" offset (in a direction opposite to that
   of plate's normal).
   canals are rational B-splines approximating the canal surfaces, i.e.
   the offsets of the four plate edges.  canals[0] = offset to v = 0 .
   canals[1] = offset to v = 1 . canals[2] = offset to u = 0 .
   canals[3] = offset to u = 1 .
   corners are sphere segments (offsets to the four plate corners) .
   corners[0] = offset to u = v = 0 . corners[1] = offset to u = 1 , v = 0 .
   corners[2] = offset to u = 0 , v = 1 . corners[3] = offset to u = v = 1 .
   nsec is an integer array containing the initial number of
   generatrices used to build the lateral canal surfaces along u and v,
   respectively.
   nofs is an integer array containing the initial number of internal
   knots along u,v during the first iteration step of the offset building
   process.
   eps is the maximum allowable position difference between
   fc is a factor by which to divide eps if the edge merging operation
   disturbs convergence between approximating and exact surfaces.
   approximating and true surfaces.
   ncan[0] is the number of sampled canal surface sections between
   consecutive interpolated sections.
   ncan[1] is the number of sampled points on each section.
   noff[0], noff[1] are the number of sampling points on the
   approximating and the true offset surfaces of plate in the u and v
   direction, respectively.
   nm is an integer array containing the number of sampling points after
   the modification of the tolerance region edges.
   nm[0] is the number of sampling points in the direction along the edge.
   nm[1] is the number of sampling points in the first nontrivial span
   of the knot vector "normal" to the edge.
   The routine returns the achieved tolerance. */

{
  ParCurv ***gen;
  double **pars,**par,t0,u,v,ep;
  int *m,*m0,l,s,uv,cm;

/* Meaning of the most important local variables (in alphabetical order)

cm       =  2*s + uv (see below).
            cm = 0 corresponds to the edge v = 0 .
            cm = 1 corresponds to the edge v = 1 .
            cm = 2 corresponds to the edge u = 0 .
            cm = 3 corresponds to the edge u = 1 .
ep       =  current maximum position error between exact and
            approximating surfaces, returned by this routine.
gen      =  two dimensional array of ParCurv (4 x MAXSECT)
            gen[i][.] contains the interpolated sections of canals[i].
l        =  logical variable. If false iterate to build approximating
            canal and offset surfaces with a tighter tolerance.
m        =  integer array of dimension 2.
            m[0] contains the number of interpolated sections of
                 canals[0], canals[1].
            m[1] contains the number of interpolated sections of
                 canals[2], canals[3].
            m is formed after the merging operation applied to pairs of
            opposite canal surfaces.
m0       =  integer array of dimension 4.
            m[i] contains the number of interpolated sections
            (generatrices) used to build canals[i] during the first
            approximation step.
par      =  2-dimensional double array (2 x MAXSECT) .
            par[0][.] contains the locations of the interpolated
            sections (u-parameter values of plate), which are used to
            build canals[0], canals[1] .
            par[1][.] contains the locations of the interpolated
            sections (v-parameter values of plate), which are used to
            build canals[2], canals[3] .
pars     =  2-dimensional double array (4 x MAXSECT).
            pars[i][.] contains the locations (u or v parameter values
            of plate) of the interpolated sections, which have been used to
            build canals[i] during the first approximation step.
s        =  loop control variable.
            s = 0 characterizes a v = ct. edge.
            s = 1 characterizes a u = ct. edge.
uv       =  loop control variable.
            uv = 0  corresponds to u = 0  or  v = 0 .
            uv = 1  corresponds to u = 1  or  v = 1 .

*/

/* Allocate memory */

  gen = (ParCurv ***)gen_array2(4, MAXSECT, sizeof(ParCurv **));
  pars = dbl_array2(4, MAXSECT);
  par = dbl_array2(2, MAXSECT);
  m0 = int_array1(4);
  m = int_array1(2);
  l = 1;
  ep = eps;

  do {
    if (!l) ep /= fc ;

/* Build approximating canal surfaces */

    printf(" Approximate canal surfaces\n");
    for (s=0; s<2; s++)
      for (uv=0;uv<2;uv++) {
	cm = 2*s+uv;
	t0 = (double)uv;
	if (!l)
	  free_fgeom(canals[cm]);
	if ((canals[cm] = canal_to_rat(plate, t0, s, ncan[0], ncan[1], ep,
				       pars[cm], gen[cm], &m0[cm], nsec[s],
				       tz)) == NULL)
	  return (1);
      }

/* Merge opposite pairs of canal surfaces */

    if (merge_canals(plate, canals, gen, pars, par, m0, m, tz))
      return (1);

/* Compute approximating offset surfaces. */

    if (!l)
      for (s=0; s<2; s++)
	free_fgeom(offset[s]);
    if (integral_offset(plate, offset, tz, nofs, ep, noff))
      return (1);
    
/* Merge knot vectors of canal surfaces and offsets */

    printf(" Merge knot vectors of canal surfaces and offsets\n");
    if (merge_can_off(canals, offset))
      return (1);

/* Merge corresponding edges of canals, offsets */

    printf(" Merge edges of canal surfaces and offsets\n");
    merge_tol_edges(plate, canals, offset, tz);

/* Sample first span of canals, offset knot vector to check if position
difference is still within ep */

    printf(" Sample span of knot vector\n");
    l = sample_mod_surf(plate, canals, offset, nm, eps, tz);
  } while (!l);
  
/* Compute corner offsets */

  printf(" Compute corner offsets\n");
  for (cm=0; cm<4; cm++)
    corners[cm] = corner_offset(canals, cm, tz);

/* Free auxiliary space */

  free_darray2(pars);
  free_darray2(par);
  free_earray2(gen, 4, MAXSECT);
  free_iarray1(m0);
  free_iarray1(m);

  return(0);
}

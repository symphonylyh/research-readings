/***********************************************************************************
 Copyright (C) 2008 Massachusetts Institute of Technology, Cambridge, MA
 All rights reserved

						geod.c

  Given an open NURBS surface, the geodesics are computed at a specified 
  number of (nsegu & nsegv) points along the surface boundary (in v and u) 
  w/ initial directions orthogonal to the boundary curves.

	Make:
	make geod

	Run:
	geod [-s step_size] [-u number_of_segments_per_u-knot_span] [-v number_of_segments_per_v-knot_span] [-o output_file_for_geod]   input_surface_file 

    Output: 
	Resulting geodesics in .VECT format 
 
    Example:
           prompt> geod -s 0.01  -u 20 -v 20  -o geod.VECT s.SURF
  *********************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h>
#include <math.h>
#include <getopt.h>
#include "gen.h"
#include "bspl.h"
#include "geod.h"


int main(unsigned argc, char *argv[])
{
  FILE *fp;
  ParSurf *fgeom;
  double h = 0.1;
  int c;
  short errflg = 0, nsegu = 20, nsegv = 20;
  char *input = NULL, *output = NULL, version = 0;

  if (argc > 1) {
    while ((c = getopt(argc, argv, "Vs:u:v:o:")) != -1)
      switch (c) {
      case 's':
	h = atof(optarg);
	if (h < 0.01) h = 0.01;
	if (h > 1.00) h = 1.0;
	break;
      case 'u':
	nsegu = atoi(optarg);
	if (nsegu <   1) nsegu =   1;
	if (nsegu > 100) nsegu = 100;
	break;
      case 'v':
	nsegv = atoi(optarg);
	if (nsegv <  1)  nsegv =   1;
	if (nsegv > 100) nsegv = 100;
	break;
      case 'o':
	output = optarg;
	break;
      case 'V':       
	version = 1;
	break;
      default:
	errflg = 1;
      }
    input = argv[optind];
  }
  else
    errflg = 1;

  if (errflg) {
    printf("usage: geod [-s step_size] [-u number_of_segments_per_u-knot_span] [-v number_of_segments_per_v-knot_span] [-o output_file_for_geod]   input_surface_file\n");
    return (-1);
  }


  if (fp = fopen(input,"r")) {
    fgeom = ReadParSurf(fp, 0);
    fclose(fp);
  }
  else {
    printf("geod: can't open input file \"%s\"!\n", input);
    return(-2);
  }

  if (output) {
    if (fp = fopen(output, "w")) {
      DrawGeodesics(fgeom, nsegu, nsegv, h, fp);
      fclose(fp);
    }
    else
      printf("geod: can't open output file \"%s\"!\n", output);
  }
  else {
    DrawGeodesics(fgeom, nsegu, nsegv, h, stdout);
  }
  return (0);
}

void DrawGeodesics(ParSurf *fgeom, short nsegu, short nsegv, double h,
		   FILE *fp)
{
  vector *vect;
  double x,u,v,du,dv,y[4],yout[4],dy[4];
  short i,j,n;

  n = 4; x = 0.0;

  u = (double)fgeom->uknots[fgeom->ucontpts];
  v = (double)fgeom->vknots[fgeom->vcontpts];
  du = u/(double)(nsegu+1);
  dv = v/(double)(nsegv+1);

  for (i=1; i<=nsegu; i++) {
    y[1] = yout[1] = 0.0;
    y[0] = yout[0] = i*du;
    find_init(fgeom, yout);
    vect = revalderivsurf(fgeom, yout[0], yout[1], 0, 0);
    fprintf(fp, "M %f %f %f ", vect->x, vect->y, vect->z);
    free(vect);
    while (1) {
      for (j=0; j<4; j++)
	y[j] = yout[j];
      fn(fgeom, x, y, dy);
      rk4(fgeom, y, dy, n, x, h, yout, fn);
      if ((yout[0] > 1.0 || yout[0] < 0.0) || (yout[1] > 1.0 || yout[1] < 0.0))
	break;
      vect = revalderivsurf(fgeom, yout[0], yout[1], 0, 0);
      fprintf(fp, "D %f %f %f ", vect->x, vect->y, vect->z);
      free(vect);
    }
    trn(yout, y);
    vect = revalderivsurf(fgeom, yout[0], yout[1], 0, 0);
    fprintf(fp, "D %f %f %f ", vect->x, vect->y, vect->z);
    free(vect);
  }

  for(i=1; i<=nsegv; i++){
    y[0] = yout[0] = 0.0; y[1] = yout[1] = i*dv;
    find_init(fgeom, yout);
    vect = revalderivsurf(fgeom, yout[0], yout[1], 0, 0);
    fprintf(fp, "M %f %f %f ", vect->x, vect->y, vect->z);
    free(vect);
    while (1) {
      for (j=0; j<4; j++)
	y[j] = yout[j];
      fn(fgeom, x, y, dy);
      rk4(fgeom, y, dy, n, x, h, yout, fn);
      if ((yout[0] > 1.0 || yout[0] < 0.0) || (yout[1] > 1.0 || yout[1] < 0.0))
        break;
      vect = revalderivsurf(fgeom, yout[0], yout[1], 0, 0);
      fprintf(fp, "D %f %f %f ", vect->x, vect->y, vect->z);
      free(vect);
    }
    trn(yout,y);
    vect = revalderivsurf(fgeom, yout[0], yout[1], 0, 0);
    fprintf(fp, "D %f %f %f ", vect->x, vect->y, vect->z);
    free(vect);
  }
  fprintf(fp, "E");
}

void trn(double *y, double *yp)
{
  if (y[0] < 0.0) {
    y[1] = yp[1] + ((y[1]-yp[1])*(0.0-yp[0]))/(y[0]-yp[0]);
    y[0] = 0.0;
 }
  else if (y[0] > 1.0) {
    y[1] = yp[1] + ((y[1]-yp[1])*(1.0-yp[0]))/(y[0]-yp[0]);
    y[0] = 1.0;
  }
  else if (y[1] < 0.0) {
    y[0] = yp[0] + ((y[0]-yp[0])*(0.0-yp[1]))/(y[1]-yp[1]);
    y[1] = 0.0;
  }
  else if (y[1] > 1.0) {
    y[0] = yp[0] + ((y[0]-yp[0])*(1.0-yp[1]))/(y[1]-yp[1]);
    y[1] = 1.0;
  }
}

void rk4(ParSurf *fgeom, double *y, double *dy, short n, double x, double h,
	 double *yout, void (*derivs)(ParSurf *, double, double *, double *))
{
  double xh,hh,h6,*dym,*dyt,*yt;
  short i;
  
  dym = dbl_array1(n);
  dyt = dbl_array1(n);
  yt = dbl_array1(n);

  hh = 0.5*h;
  h6 = h/6.0;
  xh = x + hh;

  for (i=0; i<n; i++)
    yt[i] = y[i] + hh*dy[i];

  if ((yt[0] > 1.0 || yt[0] < 0.0) || (yt[1] > 1.0 || yt[1] < 0.0)) {
    yout[0] = yt[0];
    yout[1] = yt[1];
    return;
  }

  (*derivs)(fgeom, xh, yt, dyt);
  
  for (i=0; i<n; i++) 
    yt[i] = y[i] + hh*dyt[i];

  if ((yt[0] > 1.0 || yt[0] < 0.0) || (yt[1] > 1.0 || yt[1] < 0.0)) {
    yout[0] = yt[0];
    yout[1] = yt[1];
    return;
  }

  (*derivs)(fgeom, xh, yt, dym);

  for (i=0; i<n; i++) {
    yt[i] = y[i] + h*dym[i];
    dym[i] += dyt[i];
  }

  if ((yt[0] > 1.0 || yt[0] < 0.0) || (yt[1] > 1.0 || yt[1] < 0.0)) {
    yout[0] = yt[0];
    yout[1] = yt[1];
    return;
  }

  (*derivs)(fgeom, x+h, yt, dyt);

  for (i=0; i<n; i++)
    yout[i] = y[i]+h6*(dy[i]+dyt[i]+2.0*dym[i]);

  free_darray1(yt);
  free_darray1(dyt);
  free_darray1(dym);
}

void find_init(ParSurf *fgeom, double *y)
{
  vector *r[3],*ru,*rv;

  evalrsurf(fgeom, y[0], y[1], 2, r);

  ru = r[1];
  rv = r[2];

  if (y[0] == 0.0) {
    y[2] = 1.0 / sqrt(dot(ru, ru));
    y[3] = 0.0;
  }
  else if (y[1] == 0.0) {
    y[2] = 0.0;
    y[3] = 1.0/sqrt(dot(rv,rv));
  }
  free(r[0]);
  free(r[1]);
  free(r[2]);
}

void fn(ParSurf *fgeom, double x, double *y, double *dy)
{
  vector *r[6],*ru,*rv,*ruu,*rvv,*ruv;
  vector *r1,*r2,*r3,*r4,*r5,*r6,*n;
  double dd,d11,d12,d22,g11,g12,g22;
  short i;

  evalrsurf(fgeom, y[0], y[1], 5, r);

  ru  = r[1];
  rv  = r[2];
  ruv = r[3];
  ruu = r[4];
  rvv = r[5];

  n = cross(ru, rv);
  dd = mag(n);
  unitvector1(n, n);
  
  r1 = cross(ru, ruu);
  r2 = cross(ru, ruv);
  r3 = cross(ru, rvv);

  r4 = cross(ruu, rv);
  r5 = cross(ruv, rv);
  r6 = cross(rvv, rv);

  g11 = dot(n, r1) / dd;
  g12 = dot(n, r2) / dd;
  g22 = dot(n, r3) / dd;

  d11 = dot(n, r4) / dd;
  d12 = dot(n, r5) / dd;
  d22 = dot(n, r6) / dd;

  dy[0] = y[2];
  dy[1] = y[3];
  dy[2] = -(d11*y[2]*y[2]+2.0*d12*y[2]*y[3]+d22*y[3]*y[3]);
  dy[3] = -(g11*y[2]*y[2]+2.0*g12*y[2]*y[3]+g22*y[3]*y[3]);

  free(n);
  free(r1);
  free(r2);
  free(r3);
  free(r4);
  free(r5);
  free(r6);
  for(i=0;i<6;i++)
    free(r[i]);
}

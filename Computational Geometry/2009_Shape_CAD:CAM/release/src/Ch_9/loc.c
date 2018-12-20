/***********************************************************************************
 Copyright (C) 2008 Massachusetts Institute of Technology, Cambridge, MA
 All rights reserved

						loc.c

  Given an open NURBS surface, the lines of principal curvatures are 
  computed at a specified (nsegu & nsegv) points

	Make:
	make loc

	Run:
	loc [-s step_size] [-u number_of_segments_per_u-knot_span] [-v number_of_segments_per_v-knot_span]
	    [-x output_file_for_max_loc] [-n output_file_for_min_loc]   input_surface_file 

    Output: 
	Resulting lines of curvatures in .VECT format
 
    Example:
           prompt> loc -s 0.01  -u 20 -v 20  -x loc-max.VECT -n loc-min.VECT  s.SURF
  *********************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h>
#include <math.h>
#include <getopt.h>
#include "gen.h"
#include "bspl.h"
#include "loc.h"

static vector *vpu,*vpv;
static double dyp[2];
static short iflag, jflag;

int main(unsigned argc, char *argv[])
{
  FILE *fp;
  ParSurf *fgeom;
  double h = 0.1;
  int c;
  short errflg = 0, nsegu = 20, nsegv = 20;
  char *input = NULL, *out1 = NULL, *out2 = NULL, version = 0;

  vpu = vectalloc();
  vpv = vectalloc();
  clear_vector(vpu);
  clear_vector(vpv);

  if (argc > 1) {
    while ((c = getopt(argc, argv, "Vs:u:v:x:n:")) != -1)
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
      case 'x':
	out1 = optarg;
	break;
      case 'n':
	out2 = optarg;
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

  if (errflg || input == NULL || out1 == NULL || out2 == NULL) {
    printf("usage: loc [-s step_size] [-u number_of_segments_per_u-knot_span] [-v number_of_segments_per_v-knot_span] [-x output_file_for_max_loc] [-n output_file_for_min_loc]   input_surface_file\n");
    return(-1);
  }

  if (fp = fopen(input, "r")) {
    fgeom = ReadParSurf(fp, 0);
    fclose(fp);
  }
  else {
    printf("loc: can't open input file \"%s\"!\n", input);
    return(-2);
  }

  fp = fopen(out1, "w");
  iflag = 1;
  
  DrawLinesOfCurvature(fgeom, nsegu, nsegv, h, fp);
  fclose(fp);

  fp = fopen(out2, "w");
  iflag = 2;
  
  DrawLinesOfCurvature(fgeom, nsegu, nsegv, h, fp);
  fclose(fp);

  return (0);
}

void DrawLinesOfCurvature(ParSurf *fgeom, short nsegu, short nsegv, double h,
			  FILE *fp)
{
  int i,n;
  double u0,u1,v0,v1;
  double x,u,v,du,dv,y[2],yout[2],dy[2];
  vector *vect[1];
  
  n = 2;
  x = 0.0;

  u0 = fgeom->uknots[fgeom->uorder-1];
  v0 = fgeom->vknots[fgeom->vorder-1];
  u1 = fgeom->uknots[fgeom->ucontpts];
  v1 = fgeom->vknots[fgeom->vcontpts];

  u = u1 - u0;
  v = v1 - v0;
  du = u/(double)(nsegu+1);
  dv = v/(double)(nsegv+1);

  for (i=1; i<=nsegu; i++){
    y[1] = yout[1] = v0;
    y[0] = yout[0] = u0 + i*du;
    clear_vector(vpu);
    clear_vector(vpv);
    evalrsurf(fgeom, yout[0], yout[1], 0, vect);
    fprintf(fp, "M %f %f %f ", vect[0]->x, vect[0]->y, vect[0]->z);
    free(vect[0]);
    while (1) {
      y[0] = yout[0];
      y[1] = yout[1];
      fn(fgeom, x,y,dy);
      rk4(fgeom, y, dy, n, x, h, yout, fn);
      if ((yout[0] > u1 || yout[0] < u0) || (yout[1] > v1 || yout[1] < v0))
        break;
      evalrsurf(fgeom,yout[0],yout[1],0,vect);
      fprintf(fp, "D %f %f %f ", vect[0]->x, vect[0]->y, vect[0]->z);
      free(vect[0]);
    }
    trn(fgeom, yout, y);
    evalrsurf(fgeom, yout[0], yout[1], 0, vect);
    fprintf(fp, "D %f %f %f ", vect[0]->x, vect[0]->y, vect[0]->z);
    free(vect[0]);
  }

  for (i=1; i<=nsegv; i++){
    y[0] = yout[0] = u0; y[1] = yout[1] = v0 + i*dv;
    clear_vector(vpu);
    clear_vector(vpv);
    evalrsurf(fgeom, yout[0], yout[1], 0, vect);
    fprintf(fp, "M %f %f %f ", vect[0]->x, vect[0]->y, vect[0]->z);
    free(vect[0]);
    while (1) {
      y[0] = yout[0];
      y[1] = yout[1];
      fn(fgeom, x, y, dy);
      rk4(fgeom, y, dy, n, x, h, yout, fn);
      if ((yout[0] > u1 || yout[0] < u0) || (yout[1] > v1 || yout[1] < v0))
        break;
      evalrsurf(fgeom, yout[0], yout[1], 0, vect);
      fprintf(fp, "D %f %f %f ", vect[0]->x, vect[0]->y, vect[0]->z);
      free(vect[0]);
    }
    trn(fgeom, yout, y);
    evalrsurf(fgeom, yout[0], yout[1], 0, vect);
    fprintf(fp, "D %f %f %f ", vect[0]->x, vect[0]->y, vect[0]->z);
    free(vect[0]);
  }
  
  for (i=1; i<=nsegu; i++) {
    y[1] = yout[1] = v1;
    y[0] = yout[0] = u0 + i*du;
    clear_vector(vpu);
    clear_vector(vpv);
    evalrsurf(fgeom, yout[0], yout[1], 0, vect);
    fprintf(fp, "M %f %f %f ", vect[0]->x, vect[0]->y, vect[0]->z);
    free(vect[0]);
    while (1) {
      y[0] = yout[0];
      y[1] = yout[1];
      fn(fgeom, x, y, dy);
      rk4(fgeom, y, dy, n, x, h, yout, fn);
      if ((yout[0] > u1 || yout[0] < u0) || (yout[1] > v1 || yout[1] < v0))
        break;
      evalrsurf(fgeom, yout[0], yout[1], 0, vect);
      fprintf(fp, "D %f %f %f ", vect[0]->x, vect[0]->y, vect[0]->z);
      free(vect[0]);
    }
    trn(fgeom, yout, y);
    evalrsurf(fgeom, yout[0], yout[1], 0, vect);
    fprintf(fp, "D %f %f %f ", vect[0]->x, vect[0]->y, vect[0]->z);
    free(vect[0]);
  }

  for (i=1; i<=nsegv; i++) {
    y[0] = yout[0] = u1;
    y[1] = yout[1] = v0 + i*dv;
    clear_vector(vpu);
    clear_vector(vpv);
    evalrsurf(fgeom, yout[0], yout[1], 0, vect);
    fprintf(fp, "M %f %f %f ", vect[0]->x, vect[0]->y, vect[0]->z);
    free(vect[0]);
    while (1) {
      y[0] = yout[0];
      y[1] = yout[1];
      fn(fgeom, x, y, dy);
      rk4(fgeom, y, dy, n, x, h, yout, fn);
      if ((yout[0] > u1 || yout[0] < u0) || (yout[1] > v1 || yout[1] < v0))
        break;
      evalrsurf(fgeom, yout[0], yout[1], 0, vect);
      fprintf(fp, "D %f %f %f ", vect[0]->x, vect[0]->y, vect[0]->z);
      free(vect[0]);
    }
    trn(fgeom, yout, y);
    evalrsurf(fgeom, yout[0], yout[1], 0, vect);
    fprintf(fp, "D %f %f %f ", vect[0]->x, vect[0]->y, vect[0]->z);
    free(vect[0]);
  }
  fprintf(fp, "E");
}

void trn(ParSurf *fgeom, double *y, double *yp)
{
  double u0,u1,v0,v1;

  u0 = fgeom->uknots[fgeom->uorder-1];
  v0 = fgeom->vknots[fgeom->vorder-1];
  u1 = fgeom->uknots[fgeom->ucontpts];
  v1 = fgeom->vknots[fgeom->vcontpts];

  if ( y[0] < u0 ){
    y[1] = yp[1] + ((y[1]-yp[1])*(u0-yp[0]))/(y[0]-yp[0]);
    y[0] = u0;
  }
  else if ( y[0] > u1 ){
    y[1] = yp[1] + ((y[1]-yp[1])*(u1-yp[0]))/(y[0]-yp[0]);
    y[0] = u1;
  }
  else if ( y[1] < v0 ){
    y[0] = yp[0] + ((y[0]-yp[0])*(v0-yp[1]))/(y[1]-yp[1]);
    y[1] = v0;
  }
  else if ( y[1] > v1 ){
    y[0] = yp[0] + ((y[0]-yp[0])*(v1-yp[1]))/(y[1]-yp[1]);
    y[1] = v1;
  }
}

void rk4(ParSurf *fgeom, double *y, double *dy, short n, double x, double h,
	 double *yout, void (*derivs)(ParSurf *, double, double *, double *))
{
  double u0,u1,v0,v1;
  double xh,hh,h6,*dym,*dyt,*yt;
  short i;

  u0 = fgeom->uknots[fgeom->uorder-1];
  v0 = fgeom->vknots[fgeom->vorder-1];
  u1 = fgeom->uknots[fgeom->ucontpts];
  v1 = fgeom->vknots[fgeom->vcontpts];

  dym = dbl_array1(n);
  dyt = dbl_array1(n);
  yt = dbl_array1(n);
  
  hh = 0.5*h;
  h6 = h/6.0;
  xh = x + hh;

  for (i=0; i<n; i++) 
    yt[i] = y[i] + hh*dy[i];

  if ((yt[0] > u1 || yt[0] < u0)||(yt[1] > v1 || yt[1] < v0)){
    yout[0] = yt[0]; yout[1] = yt[1]; return; }

  (*derivs)(fgeom, xh, yt, dyt);
  
  for (i=0; i<n; i++) 
    yt[i] = y[i] + hh*dyt[i];
  
  if ((yt[0] > u1 || yt[0] < u0)||(yt[1] > v1 || yt[1] < v0)){
    yout[0] = yt[0]; yout[1] = yt[1]; return; }
  
  (*derivs)(fgeom, xh, yt, dym);

  for (i=0; i<n; i++){
    yt[i] = y[i] + h*dym[i];
    dym[i] += dyt[i];
  }

  if ((yt[0] > u1 || yt[0] < u0)||(yt[1] > v1 || yt[1] < v0)){
    yout[0] = yt[0]; yout[1] = yt[1]; return; }

  (*derivs)(fgeom, x+h, yt, dyt);

  for (i=0; i<n; i++) 
    yout[i] = y[i]+h6*(dy[i]+dyt[i]+2.0*dym[i]);

  free_darray1(yt);
  free_darray1(dyt);
  free_darray1(dym);
}

void fn(ParSurf *fgeom, double x, double *y, double *dy)
{
  vector *r[6],*ru,*rv,*ruu,*rvv,*ruv,*n;
  double d11,d12,d21,d22,g11,g12,g21,g22;
  double H,K,Gdet,Ddet,kmax,kmin,kk,bb,tempmod;
  short i;

  evalrsurf(fgeom, y[0], y[1], 5, r);

  ru  = r[1];
  rv  = r[2];
  ruv = r[3];
  ruu = r[4];
  rvv = r[5];

  n = cross(ru, rv);
  unitvector1(n, n);

  g11 = dot(ru, ru);
  g12 = g21 = dot(ru, rv);
  g22 = dot(rv, rv);
  
  d11 = dot(n, ruu);
  d12 = d21 = dot(n, ruv);
  d22 = dot(n, rvv);

  Gdet = g11*g22 - g12*g21;
  Ddet = d11*d22 - d12*d21;

  H = ((g12*d21 + g21*d12) - (g11*d22 + g22*d11))/(2.0*Gdet);

  K = Ddet/Gdet;

  tempmod = H*H - K;
  tempmod = fabs(tempmod);
  
  if(tempmod >= 0.0 && tempmod < ZERO) {
    kmax = H;
    kmin = H;
  }
  else if(tempmod >= ZERO) {
    kmax = H + sqrt(tempmod);
    kmin = H - sqrt(tempmod);
  }

  if (iflag == 1)
    kk = kmax;
  else if (iflag == 2)
    kk = kmin;

  bb = sqrt(g11*(d12+kk*g12)*(d12+kk*g12)
	    -2.0*g12*(d12+kk*g12)*(d11+kk*g11)
	    +g22*(d11+kk*g11)*(d11+kk*g11));

  bb = 1.0/bb;

  dy[0] =  bb*(d12+kk*g12);
  dy[1] = -bb*(d11+kk*g11);

  find_dir(fgeom, ru, rv, y, dy);

  if (jflag == 2) {
    dy[0] = -dy[0];
    dy[1] = -dy[1];
  }

  copyvector(ru, vpu);
  copyvector(rv, vpv);
  dyp[0] = dy[0];
  dyp[1] = dy[1];

  free(n);
  for(i=0;i<6;i++)
    free(r[i]);

}

void find_dir(ParSurf *fgeom, vector *ru, vector *rv, double *y, double *dy)
{
  vector *vv1,*vv2,*vv3,*vv4,*vv5,*vv6,*rhs,*lhs;
  double a1,a2;
  double u0,u1,v0,v1;

  u0 = fgeom->uknots[fgeom->uorder-1];
  v0 = fgeom->vknots[fgeom->vorder-1];
  u1 = fgeom->uknots[fgeom->ucontpts];
  v1 = fgeom->vknots[fgeom->vcontpts];

  if ((y[0] == u0 && dy[0] < 0.0) || (y[1] == v0 && dy[1] < 0.0)) {
    jflag = 2;
    return;
  }

  if ((y[0] == u1 && dy[0] > 0.0) || (y[1] == v1 && dy[1] > 0.0)) {
    jflag = 2;
    return;
  }

  vv1 = vectalloc();
  vv2 = vectalloc();
  vv3 = vectalloc();
  vv4 = vectalloc();
  vv5 = vectalloc();
  vv6 = vectalloc();
  rhs = vectalloc();
  lhs = vectalloc();
  
  scale_vect1(-dyp[0], vpu, vv1);
  scale_vect1(-dyp[1], vpv, vv2);
  scale_vect1(-dy[0], ru, vv3);
  scale_vect1(-dy[1], rv, vv4);
  scale_vect1(dyp[0], vpu, vv5);
  scale_vect1(dyp[1], vpv, vv6);

  add_vect1(vv1, vv2, lhs);
  add_vect1(lhs, vv3, lhs);
  add_vect1(lhs, vv4, lhs);

  add_vect1(vv5, vv6, rhs);
  add_vect1(rhs, vv3, rhs);
  add_vect1(rhs, vv4, rhs);

  a1 = mag(lhs);
  a2 = mag(rhs);

  if (a1 < a2)
    jflag = 2;
  else
    jflag = 1;

  free(vv1);
  free(vv2);
  free(vv3);
  free(vv4);
  free(vv5);
  free(vv6);
  free(lhs);
  free(rhs);
}

void clear_vector(vector *vect)
{
  vect->x = 0.0;
  vect->y = 0.0;
  vect->z = 0.0;
  vect->w = 1.0;
}

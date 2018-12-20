/************************************************************************
 *									*
			 Copyright (C) 1989, 1990, 1992, 1994
     Massachusetts Institute of Technology, Cambridge, Massachusetts.
			   All rights reserved.

 *									*
 ************************************************************************/

/* bspl.h header file - bspl library */

#ifndef BSPL_H_
#define BSPL_H_

#include "gen.h"

#ifdef __cplusplus
extern "C" {
#endif

#define PRECISION 1.0e-14
/* function declaration   */
void bernrs(double2[],int,double2,double2,double2,double2[],int*);
void bernrs1(double2[],int,double2,double2[],int*);
int rootisolate(double2[],double2[],int,double2,double2,double2**,int*);
void splitpoly2(double2[],int,double2,double2,double2,double2[],double2[]);
void splitpoly(double2[],int,double2,double2[],double2[]);
int findichange(double2[],int);
void berntomono(ParSurf*,double2**,int*,int*);
int signchange(double [], int);

int boehm_curve(ParCurv*,double2,int);
void boehm_surface(ParSurf *, double, int, char);
void extract_edge(ParCurv*,ParSurf*,int,int);
void writevect(FILE *, vector *);
void alloc_fgeompts(ParSurf *);
void alloc_egeompts(ParCurv *);

double2 abssumarray(double2[],int);
double2 sum(double2,double2);
double2 findroot(double2[],double[],int,double2,double2);

double2 zevalderivsurf(ParSurf*,double2,double2,int,int);
double2 zevaldsurf(ParSurf*,int,int,int,int,double**,double**);
double2 combine(int,int);
double2 factorial(int);
double2 absmaxarray(int,int,double2[],int);
double2 bsp_basis(int,int,double2,double2[]);

void convexbox(ParSurf*,double2*,double2*,double2*,double2*,double2*,double2*);
int curvature(ParSurf*,double2,double2,double2*,double2*);
int curvature1(ParSurf*,double2,double2,double2*,double2*);

void nbasisd(int,int,double2,double2[],double2**);
void Aijqs(int,int,int,double2[],int,int,int,double2[],vector***,double2***);
void Aijqs_z(int,int,int,double2[],int,int,int,double2[],vector***,double2***);
void evalrsurf(ParSurf*,double2,double2,int,vector*[]);
void Aij(int,int,int,double[],vector*[],double***);

void nodes_bsp_per(double2[],int,int,double2[]);
void nbasisd_per(int,int,double,double2[],double2**,int*);
void Aijqs_per(int,int,int,int*,double2[],int,int,int,double2[],vector***,
	       double2***);
void Aij_per(int,int,int,int[],double[],vector*[],double***);
void evalrsurf_per(ParSurf*,double2,double2,int,vector*[]);
void n_matrix_per(double **, double *, double *, short, short);

vector *evalbsp(ParCurv*,double);
vector *evalderivbsp(ParCurv*,double,int);
vector *rbspeval(ParCurv*,double,int);
vector *evalsurf(ParSurf*,double,double);
vector *normalcurv(ParCurv *, double, int);
vector *normalsurf(ParSurf*,double,double);
vector *revalderivsurf(ParSurf*,double,double,int,int);
vector *evalderivsurf(ParSurf*,double,double,int,int);
vector *normalsurf1(ParSurf*,double,double);
vector *evalsurf1(ParSurf*,int,int,double2**,double2**);
vector *evalsurf1_per(ParSurf*,int,double**,double**,int*);
vector *evaldersurf(ParSurf*,int,int,int,int,double**,double**);
void eval_surface_bounded(vector*[], ParSurf*, int, double[]);

void surface_translate(double *, double *, double *, double *,
		       double *, double *);
void surface_scale_trans(ParSurf *, double, double, double, double);

void boundbox_c(ParCurv *, double []);
void boundbox(ParSurf *, double []);
int compare_boxes(double [],double []);

int merge_knotv(ParCurv*[],int,int,double2[],int[],double2);
int find_nextknot(ParCurv*[],int,int,int[],double2[],int,double2);
int find_minknot(ParCurv*[],int,int,int[],double2);
void addpoints(ParCurv*[],int,int,double2[],int);

void merge_fgeom(ParSurf*,ParSurf*,ParCurv*[],double2[],double2[],int[],
		 double);
void addpoints_surf(ParSurf*,double2[],double2[],int,int);

void nodes_bsp(double2[],int,int,double2[]);
int lf2comp(double2,double2,double2);
void bernmat(int,double2[],int);
int find (int,double2[],double2);
void monotobern(double**,int,int,double2,double2,double2,double2,ParSurf*,int);
int linear_trans(int,double2**,double2,double2);
void curve_oslo1(int,int,int,double2[],double2[],double2**,double2**);
void discrete_bsp(int,int,int,double2[],double2[],double2**);
void surfoslo3(ParSurf*,ParSurf*,int);
void soslo(int,int,int,int,int,int,double2[],double2[],
	  double2[],double2[],vector***,vector***,int);

ParCurv *line_gen(vector**);
ParCurv *circle_gen(vector **, double);
ParCurv *arc_gen(vector **, double);
int      arc_gen_intersect(vector *, vector *, vector *, vector *, vector *,
			   double *, double *);

ParSurf *plane_gen(vector**);
ParSurf *cylinder_gen(vector**,double);
ParSurf *sphere_gen(vector*,double);
ParSurf *torus_gen(vector*[],double[]);
ParSurf *cone_gen(vector**,double);

void hodograph_surf(ParSurf *, ParSurf *, char);
ParCurv *ParCurv_hodograph(ParCurv *);

int compute_lim(ParCurv*,double2**,double2[],double2);
void raise_byone(ParCurv*,double2**,double2[],int);
void raise_surf(ParSurf*,double2,int,int);
void store_vertices(ParCurv*,ParSurf*,int,int);

void unionsrt(double2[],int,double2[],int,double2[],int*,double2);
void insertval(double2[],int,int,double2,double2[]);
int member(double2[],int,double2,double2);

void subbezier(ParSurf*,double2,double2,double2,double2,ParSurf*);
void subdivbs(ParCurv*,ParCurv*,double2);
void subdivids(ParSurf *, ParSurf *, double, char,
	       int);

vector *evalbsp_per(ParCurv*,double);
vector *rbspeval_per(ParCurv*,double,int);
vector *evalderivbsp_per(ParCurv*,double,int);
vector *evaldersurf_per(ParSurf*,int,int,int,double**,double**,int*);
vector *revalderivsurf_per(ParSurf*,double,double,int,int);
vector *evalderivsurf_per(ParSurf*,double,double,int,int);
vector *evalsurf_per(ParSurf*,double,double);

void  transform(int,int,int,double[],double[],double**);

void copymat4(double[][4], double[][4]);
void translate1(double, double, double, vector *, vector *);
void chgsys(double[][4], vector *);
void chgsyspre(double[][4], vector *);
void transmat( vector *, vector *, vector *, vector *, double[][4]);
void pretransmat(vector *, vector *, vector *, vector *,double[][4]);
void midpoint1(vector *v1, vector *v2, vector *n);
void ParSurf_rotrans (ParSurf *, vector *, vector *);
/*****
void translate(vector *, vector *);
void rotate(vector *, vector *, vector *, vector *);
*****/

#ifdef __cplusplus
}
#endif

#endif

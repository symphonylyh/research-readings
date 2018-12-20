/************************************************************************
 *									*
			 Copyright (C) 1989, 1990
     Massachusetts Institute of Technology, Cambridge, Massachusetts.
			   All rights reserved.

 *									*
 ************************************************************************/

/* appr.h header file - appr library */

#ifndef APPR_H_
#define APPR_H_

#include "editor.h"
#include "gen.h"

#ifdef __cplusplus
extern "C" {
#endif

ParCurv *approx_fnbc(vector*(*)(double,int),double2[],int,int,int,int);
void approx_fnbc2(int, int);
int  approx_fnbc3(int *);
ParCurv *approx_fnbc_knot(vector*(*)(double,int),double2[],int,int,int,int,double2*);
ParCurv *approx_fn_per(vector*(*)(double,int),double2[],int,int,int);
ParCurv *approx_fnbc_linear(vector *(*)(double, int), double [], int, int);
int sample_linear (vector *(*)(double, int), ParCurv *, double [], int *, int);

void approx_cubic_curv(ParCurv*,double2**,double2*,int,int,int,double2*);
void calc_par_chord(double2**,double2*,int);
void calc_par_uniform(double2**,double2*,int);
void calc_par_hartley(double2**,double2*,int);
void calc_par_foley(double2**,double2*,int);
void find_error_curv1(ParCurv*,double2**,double2*,int,double2*);

void approx_cubic_surf(ParSurf*,double2***,double2*,int,int,int,int,int,double2*);
void find_error_surf1(ParSurf*,double2***,double2*,double2*,int,int,double2*);
double approx_surf_error(double *, double *, double *);

ParSurf *ApproxOpenSurf    (ParSurf *, int, int, double *);
ParCurv *ApproxOpenCurvIso (ParSurf *, double, double *, int, double *);
ParCurv *ApproxOpenCurv    (ParCurv *, int, double *);

int checkacc(ParCurv*,ParCurv*,double*,int*,double,int,double);
int checkoff(double2,ParCurv*,ParCurv*,double2,double2);
double funccurve(double);

int sample_offset(ParSurf*,ParSurf*,int,int,double2,double2,int,double2*,int);
void calc_bounds_offs(ParSurf*,double2*,double2*,int*,int*);
void swap_vectors(double2*,double2*,int);
void swap_ivectors(int*,int*,int);
ParCurv *repeat_offset(ParCurv *, double);

ParCurv *ParCurv_approx(ParCurv*,int,double[],double[],int,int,int,int);
void set_last(double);
ParCurv *convexhull_test(ParCurv*,ParCurv*,double[],double[],int,int,int);

vector *geod(double, int);
vector *geod_par(double,int);
void fgeodpar(double *, double [], double []);
void fcn(double *, double [], double []);
vector *offnorm(double,int);
vector *offnorm2d(double,int);

ParCurv *cos_conv_nurbs(ParSurf*,ParCurv*,int,double[],int,int,int);
ParCurv *offset_geod(ParSurf*,ParCurv*,double,double[],int,int,int);
ParCurv *offset_geod_par(ParSurf*,ParCurv*,double,double[],int,int,int);
ParCurv *offset_normal(ParSurf*,ParCurv*,double,double[],int,int,int);
ParCurv *offset_curv2d(ParCurv*,double,double[],int,int,int);

void find_dist_curv(double[],double*,int,vector*,ParCurv*);
void find_estim_curv(double[],int,vector*,ParCurv*);
 
ParCurv *orthog_proj (ParSurf*,ParCurv*,double*,int,int,double,double);
ParCurv *orthog_proj_par (ParSurf*,ParCurv*,double*,int,int,double,double);
void init_proj (ParSurf*,ParCurv*);
void set_output_ptr (void(*)(double*,double*));
vector **project_curve_onsurf (double*,int,double);

void convexhull_surf(ParSurf*,ParSurf*,double*[],int*,double[],double[],int,int,int);
void output_2D (double*,double*);
void output_3D (double*,double*);
void compute_diffeq (double*,double*,double*);

void fair_cubic_surf(ParSurf*,ParSurf*,int,int,int,int,int,int,int,int,int,int,
		     double*,double*,double*);
void fair_knot(ParCurv*,int);
void calc_disc_surf(ParSurf*,double*);
void find_dev_surf(ParSurf*,ParSurf*,double*);
void calc_disc_curv(ParCurv*,double*);
void store_edge(ParCurv*,ParSurf*,int,int);

ParCurv *fit_curve (double*,int,int,vector**(*)(double*,int,double),double,double);
double funcsurf(double,double,double,double);
double guncsurf(double,double,double,double);
double funcsurfp(double,double,double,double);
double guncsurfp(double,double,double,double);

int knot_holzle(ParCurv*,double[],double,int);
int knot_holzle_fn(vector*(*)(double,int),double[],double,int);

void calc_hartley(vector**,double*,int,int);
void calc_hartley_per(vector**,double*,int,int);
void fitnodes (vector*[],vector*(*)(double,int),double[],int);
void solve_gram_vector(double**,vector**,int,int,int,int);
void solve_gram_pervector(double**,vector**,int,int);
vector *eval_fn(double, int);
void eval_fnbc(ParCurv*,vector*(*)(double,int),int);
void solve_bc(ParCurv*,double**,vector**,int,double**,double**);
void matr_m_vect(vector*[],double**,int,int,vector*[],int);

void inter_curv(ParCurv*,double**,int,double*);
void calc_knots_hartley(double**,double*,double*,int,int);
void calc_knots_chord(double**,double*,double*,int,int);
void calc_knots_uniform(double**,double*,double*,int,int);
void calc_knots_foley(double**,double*,double*,int,int);
void calc_gram_curv(double*,double*,int,double**,int*,int*,int);
void solve_gram(double**,double**,int,int,int);
void find_error_curv(ParCurv*,double**,double*,double*);
void transfer_matrix(double**,double**,int,int,int);
void transpose_matrix(double**,double**,int,int);

void inter_percurv(ParCurv*,double**,double*);
void calc_knots_per(double**,double*,int,int);
void calc_nodes_per(double*,double*,int,int);
void calc_gram_percurv(double*,double*,int,int,double**);
void solve_gram_percurv(double**,double**,int);
void find_error_percurv(ParCurv*,double**,double*,double*);

void inter_surf(ParSurf*,double***,int,double*);
void average_uknots(double*,double*,int,int,int);
void average_vknots(double*,double*,int,int,int);
void find_error_surf(ParSurf*,double***,double*,double*,double*);
double inter_surf_error(double *, double *, double *);
void calc_gram_surf(double*,double*,double*,double*,int,int,int,int,
		    double**,int*,int*);

ParSurf *interp_loft_surf(GridSurf *, int, int, int, double, double *,
			  double *);
void find_error_loft(ParSurf *, int, int, double ***, int, double *, double **,
		     double *);
double inter_loft_error(double *, double *, double *);
void inter_curv_nodes(ParCurv *, double **, int, double *, double *);

void hartley_judd (vector**,ParCurv*,double*,int);
void interp_ortho(ParCurv*,vector*[],double[],int);
void implstline(vector*,vector*,double*,double*,double*);

ParSurf *loft_rational(ParCurv**,double*,double*,int,int);
void calc_bandwidth(double*,double*,int,int,int*,int*);

double NRroot1(double,double(*)(double),double);
void NRroot2(double,double,double,double,double,double,double*,double*,
	     double(*)(double,double,double,double),
	     double(*)(double,double,double,double),double);

void putknots(double*,ParCurv*);
void putcontpts(double**,ParCurv*);
void getcontpts(ParCurv *,double**);
ParCurv *offset_curve(ParCurv*,double,double,int);
void copyknots(ParCurv*,ParCurv*);

ParSurf *offset_surf(ParSurf*,double,double,int,int,int);
ParSurf *rational_offset(ParSurf*,double,double**,int[],double,int[],int,double*,int);

void build_offset(ParSurf*,ParSurf*,double,int);
int check_off_rat(ParSurf*);

int sample (vector*(*)(double,int),ParCurv*,double[],double[],int*,int,int);
int compar(double*,double*);
int sample_per (vector*(*)(double,int),ParCurv*,double[],double[],int*,int);
void del_dbl_knots(double[],int,int*);
void sample_error(double *, double *, double *);

double2 deriv_dev (vector *(*)(double, int), ParCurv *, double);
double2 deriv_dev_per (vector *(*)(double, int), ParCurv *, double);
double2 curv_dev (vector *(*)(double, int), ParCurv *, double);
double curv_dev_per (vector*(*)(double,int),ParCurv*,double);

int sample_surf_test (ParSurf*,ParSurf*,double*[2],double[],int*,int,int,int);
ParSurf *check_surf_param(ParSurf*,ParSurf*);

int sample_proj(ParCurv*,double*,vector**,double,int*,int,
		vector**(*)(double*,int,double),double);
ParSurf *ParSurf_approx(ParSurf*,int[],double[],double[],int,int,int[],int[]);
vector *eval_iso(double,int);
vector *eval_ParCurv(double,int);

void find_surf_dist(double[],int,int,vector*,ParSurf*);
void find_estim(double[],int,int);
void fdist(int,double[],double*,double[]);
/********** obsolete ***********
void fdist_surf(int*,double*,double*,double*);
*******************************/
/********** for NAG Marks 19,20 *********/
void fdist_surf(int*,double*,double*,double*,int *,double *);
/****************************************/
/***** don't use 1st derivatives 8/8/95 *****
void fdist_curv(int *, double [], double *fc, double []);
*****/
/********** obsolete ***********
void fdist_curv(int *, double [], double *fc);
*******************************/
/********** for NAG Marks 19,20 *********/
void fdist_curv(int *, double [], double *fc, int *, double *);
/****************************************/
vector *inters_stlines(vector*,vector*,vector*,vector*);
vector *parallel(vector*,vector*,double);

int ptatinf(vector*);

ParSurf *scattered_fit(double*,double*,double*,double*,double*,double*,
		       double,int, int, int);
double scattered_surf_error(void);
ParCurv *curve_fit(double*,double*,double*,double*,double*,double,int,int);
ParCurv *ParCurv_iso(ParSurf *,double,int,ParCurv *);

#ifdef __cplusplus
}
#endif

#endif

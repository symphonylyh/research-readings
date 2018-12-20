/* Copyright (c) Massachusetts Insititute of Technology, 1996
 * All rights reserved
 */

#ifndef INSP_H
#define INSP_H

#include "editor.h"
#include "gen.h"

ParCurv *ApproxCamber(int, int, vector *, double2*);
ParCurv *ApproxThickness(int, vector *, double2 *);
double   AxialDisplacement (vector *, vector *);
double   BladeAngleDisplacement (vector *, vector *);
void 	 CalcHubTip (ParSurf *, int, int, int, int, double *, double *);
double   callback_GetBone(double*);
double   callback_GetCamberPoint(double*);
int      camberSurface(ParSurf *, FacetSurf *, int *, double **, vector ***,
		  double **, int, double, int, double, double, int, int, int,
		  int, ParSurf **, ParCurv ***, ParCurv ***, ParCurv ***,
		  ParCurv **, ParCurv **, int);
double   ChordLengthDeviation (double, double);
void	 compute_ortho_proj(double *, double *, double *);
void	 con_confun (int *, int *, int *, int *, int *,	double *, double *,
		double *, int *, int *,	double *);
void	 con_confun_blade (int *, int *, int *, int *, int *, double *,
		double *, double *, int *, int *, double *);
void	 con_confun_trim (int *, int *, int *, int *, int *, double *,
			  double *, double *, int *, int *, double *);
void	 con_localize (ParSurf **, vector **, double **, int *, double *,
	 	short *, double *, double *, double *, short, short *,
		int *, int *);
void	 con_localize_blade (ParSurf **, vector **, double **, double **,
		short *, double *, short *, double *, double *, double *,
		short *, int *, int *);
void	 con_localize_trim (TrimSurf **, vector **, double **, int *, double *,
	 	short *, double *, double *, double *, short, short *,
		int *, int *);
void	 con_objfun (int *, int *, double *, double *,
		double *, int *, int *, double *);
void	 con_objfun_trim (int *, int *, double *, double *,
		double *, int *, int *, double *);
vector  *ConeToPlane2d (double, double, vector *, vector *);
ParCurv *ConicalSection (ParSurf *, FacetSurf *, double, double, double, int,
			 int, int, double ***, int *, int);
vector  *ConicalSectionApprox (double, int);
void     ConicalSectionEval (double *, double *);
int	 conloc_error (int *, int *);
vector  *CylinderToPlane2d (double, vector *, vector *);
/********* obsolete **********
void	 EvalMaxCamber(int *, double *, double *);
*****************************/
/********* for NAG Marks 19,20 **********/
void	 EvalMaxCamber(int *, double *, double *, int *, double *);
/****************************************/
/********* obsolete **********
void	 EvalMaxThickness(int *, double *, double *);
*****************************/
/********* for NAG Marks 19,20 **********/
void	 EvalMaxThickness(int *, double *, double *, int *, double *);
/****************************************/
/********* obsolete **********
void	 fdist_surf1 (int *, double *, double *, double *);
*****************************/
/********* for NAG Marks 19,20 **********/
void	 fdist_surf1 (int *, double *, double *, double *, int *, double *);
/****************************************/
/********* obsolete **********
void     fdist_surf2 (int *, double *, double *, double *);
*****************************/
/********* for NAG Marks 19,20 **********/
void     fdist_surf2 (int *, double *, double *, double *, int *, double *);
/****************************************/
/********* obsolete **********
void     fdist_trim1 (int *, double *, double *);
*****************************/
/********* for NAG Marks 19,20 **********/
void     fdist_trim1 (int *, double *, double *, int *, double *);
/****************************************/
double   find_estim_trim (vector *, ParSurf *, ParCurv *);
void	 find_estim2 (double [], int, int, vector ***pts, int);
double   find_min_trim_curv (double *, vector *, ParCurv *, ParSurf *);
double	 find_surf_dist1 (double *, short, short, vector *, ParSurf *);
void	 find_surf_dist2 (double [], int, int, vector *, ParSurf *,
			  vector ***, int);
int      find_trim_dist (double *, int, int, vector *, TrimSurf *);
double	*footpoint (double *, vector *, ParSurf *);
double	 footpoint_curve (double, vector *, ParCurv *);
int	 footpoint_error (void);
void	 footpoint_fcn (int *, double *, double *, double *, int *, int *);
double	*footpoint_trim (double *, vector *, TrimSurf *);
ParCurv *GenerateCamberLine(ParCurv *, double, double, ParCurv **);
double2  GetBone(ParCurv *, vector*, vector*, double);
ParCurv *GetCamberLines    (ParCurv *, ParCurv *, double, ParCurv **,
			    ParCurv **, double *, int);
vector*  GetCamberPoint(ParCurv *, vector *, vector *, double *);
double   GetConicalSection (ParSurf *, FacetSurf *, double, double, double,
			    double, double, int, ParCurv **, ParCurv **, int,
			    double ***, int *, int);
void     GetFeatures (ParCurv *, ParCurv *, ParCurv *, ParCurv *, double,
		      double, double, int, int, double *);
double   GetPlanarSection  (ParSurf *, FacetSurf *tgeom, vector *, double,
			    double, double, double, int, ParCurv **,
			    ParCurv **, int, double ***, int *, int);
double   GetRadialSection  (ParSurf *, FacetSurf *, double, double, double,
			    double, int, ParCurv **, ParCurv **, int,
			    double ***, int *, int);
double	 global_pos_error(ParCurv *, ParCurv *);
int	 init_ortho_proj(ParSurf *, ParCurv *, int, int, double *);
void     inter_curv_leo(ParCurv *, double **, double *);
ParCurv	*InvertCurve (ParCurv *, int, int, int);
ParCurv	*InvertCurve2d (ParCurv *, int, int, int);
ParSurf	*InvertSurface (ParSurf *, int, int, int);
vector  *InvertVector (vector *, int, int, int);
double	 lead_chord(ParCurv *, double *);
double	 localize_blade2 (ParSurf ***, vector **, double **, double **,
		int *, int, double *, int **, double *, double *,
		int *, int, int *, int *);
double	 localize_sparse_2d (ParCurv **, double **, double *, int, int,
		double *, int *, double *, double *, int *, int *);
double	 localize_sparse2 (ParSurf **, vector **, double **, int *,
		int, double *, int *, double *, double *, int,
		int *, int, int *, int *);
/********* obsolete **********
void	 localize_sumsq (int *, double *, double *, double *);
*****************************/
/**************** for NAG Marks 19,20 *****************/
void	 localize_sumsq (int *, double *, double *, double *, int *, double *);
/******************************************************/
/********* obsolete **********
void	 localize_sumsq_2d (int *, double *, double *, double *);
*****************************/
/**************** for NAG Marks 19,20 *****************/
void	 localize_sumsq_2d (int *, double *, double *, double *, int *, double *);
/******************************************************/
/********* obsolete **********
void	 localize_sumsq_blade (int *, double *, double *, double *);
*****************************/
/**************** for NAG Marks 19,20 *****************/
void	 localize_sumsq_blade (int *, double *, double *, double *, int *, double *);
/******************************************************/
void	 localize_sumsq_trim (int *, double *, double *); /*** , double *);*/
double	 localize_trim (TrimSurf **, vector **, double **, int *,
		int, double *, int *, double *, double *, int,
		int *, int, int *, int *);
ParSurf *LoftCamberSurface (int, int, ParCurv **, double *, int);
ParSurf *LoftThickness     (int, int, ParCurv **, double *, double *, int);
/****************** obsolete ***************************
void	 max_chord(int *, double *, double *, double *);
*******************************************************/
/**************** for NAG Marks 19,20 *****************/
void	 max_chord(int *, double *, double *, double *, int *, double *);
/******************************************************/
/****************** obsolete ***************************
void	 max_rbspcur(int *, double *, double *);
*******************************************************/
/**************** for NAG Marks 19,20 *****************/
void	 max_rbspcur(int *, double *, double *, int *, double *);
/******************************************************/
double   MaxCamber(ParCurv *, double *);
double   MaxThickness(ParCurv *, double *);
double	 maxcur(ParCurv *, double *);
double	 maxcur_csurf(ParCurv *, double *, double *);
vector  *MidpointCone(double, double, vector *, vector *);
vector  *MidpointCyl(double, vector *, vector *);
void	 output_ortho_proj(double *, double *);
double   Pitch (double, double);
double   PitchAngle (vector *, vector*);
double   PitchAngleDeviation (double, double);
double   PitchDeviation (double, double);
vector	*point_eval(vector *, ParSurf *, double *);
double	**project_curve_on_surf(double *, int *, double);
double   RakeKerwin (vector *, int);
double   RakeNavseaCyl (vector *, int, double, double);
double   RakeNavseaPlan (vector *, vector *, int);
double	 rbspcur(ParCurv *, double);
ParCurv *PlanarSection (ParSurf *, FacetSurf *, vector *, vector *, double,
			 int, int, int, double ***, int *, int);
vector  *PlanarSectionApprox (double, int);
void     PlanarSectionEval (double *, double *);
vector  *Plane2dToCone (double, double, vector *, vector *);
vector  *Plane2dToCylinder (double, vector *, vector *);
vector  *Plane2dToPlane3d (vector *, vector *, vector *, vector *);
vector  *Plane3dToPlane2d (vector *, vector *, vector *, vector *);
ParCurv *RadialSection (ParSurf *, FacetSurf *, double, double, int, int, int,
			double ***, int *, int);
vector  *RadialSectionApprox (double, int);
void     RadialSectionEval (double *, double *);
double   RakeDeviation (vector *, vector *, vector *);
void	 set_conloc_blade_ptr (void (*fun_ptr)(int));
void	 set_conloc_fun_ptr (void (*fun_ptr)(int));
void	 set_conloc_trim_ptr (void (*fun_ptr)(int));
void	 set_unconloc_2d_ptr (void (*fun_ptr)(int));
void	 set_unconloc_blade_ptr (void (*fun_ptr)(int));
void	 set_unconloc_fun_ptr (void (*fun_ptr)(int));
void     set_unconloc_trim_ptr(void (*fun_ptr)(int));
void	 signed_distance (vector *, double *, ParSurf *, double *);
double   SkewCyl (double, double, double);
double   SkewAngle (vector *, vector *, vector *, vector *, double);
double   SkewAngleDeviationCyl (double, double);
double   SkewAngleDeviationPlan (vector *, vector *, vector *, vector *,
				 double);
double   SkewDeviation (vector *, vector *, vector *, vector *);
double   SkewPlan (vector *, vector *, vector *);
int	 SpanwiseIsU (ParSurf *, int);
int	 SpanwiseU (ParSurf *, int, int);
ParSurf *TransposeSurface (ParSurf *, int, int, int);
vector  *TransposeVector (vector *, int, int, int);
double	 vect_length(double *);

#endif





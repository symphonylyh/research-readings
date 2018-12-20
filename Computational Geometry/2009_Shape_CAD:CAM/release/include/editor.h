/* Copyright (C) 1992 Massachusetts Institute of Technology, Cambridge, MA
   All rights reserved

   Last edited: April 11, 2003 for v.11.1

   v.11.1: (1) C-curves added to "TrimSurf"
*/

/* editor.h */

#ifndef EDITOR_H
#define EDITOR_H

#include "gen.h"

#ifdef __cplusplus
extern "C" {
#endif

/* constants */

#define EDITOR_EMIN	1.0e-30
#define EDITOR_EMAX	1.0e+30
#define EDITOR_ZERO_6	1.0e-6
#define EDITOR_ZERO_10	1.0e-10

#define LOG_NO_PRINT	0
#define LOG_FILE	1
#define LOG_STDOUT	2
#define LOG_FILE_STDOUT	3

#define DEG_TO_RAD      0.017453292519943295
#define RAD_TO_DEG      57.295779513082323

#define EDITOR_2_PI     6.28318530717958647692
#define EDITOR_PI       3.14159265358979323846
#define EDITOR_PI_2     1.57079632679489661923

#define EDITOR_TRIM_MAX 1000
#define EDITOR_TRIM_U	0
#define EDITOR_TRIM_V	1
#define EDITOR_TRIM_EXTRA 1.0e-6

#define EDITOR_LOW 0
#define EDITOR_UPP 1
#define EDITOR_SURF_NT_INTRSCT 0
#define EDITOR_SURF_INTRSCT 1
#define EDITOR_U 0
#define EDITOR_V 1
#define EDITOR_NO_ROOT 0
#define EDITOR_WITHIN_TOL 1
#define EDITOR_SUBDIVISION 2
#define EDITOR_BEZIER_CLIP 3
#define EDITOR_CONVEX_HULL_INTERSECT 1
#define EDITOR_MAX_INTERSECT_POINT 100
#define EDITOR_MRVAL 0.80
#define EDITOR_TOL 1.0e-3
#define EDITOR_TOLKNOT 1.0e-10
#define EDITOR_ZERO_CVX 1.0e-10
#define EDITOR_ZERO_MIN 1.0e-8
#define EDITOR_ZERO_MIN_SOL 1.0e-12
#define EDITOR_ZEROD 1.0e-10
#define EDITOR_JBSPL_TOL 1.0e-8

#define EDITOR_XPOS 0
#define EDITOR_XNEG 1
#define EDITOR_YPOS 2
#define EDITOR_YNEG 3
#define EDITOR_ZPOS 4
#define EDITOR_ZNEG 5

#define EDITOR_B_PER_BLOCK GEN_B_PER_BLOCK
#define EDITOR_KB_PER_B GEN_KB_PER_B
#define EDITOR_MB_PER_B GEN_MB_PER_B

#define EDITOR_CONICAL 2
#define EDITOR_PLANAR  1
#define EDITOR_RADIAL  0

/* list structure */

typedef struct {
  int		npoints;	/* number of points */
  int		knottype;	/* not used */
  int		knots;		/* not used */
  int		isotype;	/* not used */
  double	**pts;		/* list of (x,y,z) points */
} ListCurv;

/* grid structure */

typedef struct {
  int		ncol;		/* number of columns */
  int		nrow;		/* number of rows */
  int		knottype;	/* not used */
  int		uknots;		/* not used */
  int		vknots;		/* not used */
  int		isotype;	/* not used */
  double	***pts;		/* grid of (x,y,z) points */
} GridSurf;

/* uv structure */

typedef struct {
  double	umin;		/* parametric space */
  double	umax;
  double	vmin;
  double	vmax;
  int		npoints;	/* number of points */
  short		is_open;	/* periodicity flag */
  double	**pts;		/* uv points */
} ParUv;

/* ruled surface structure */

typedef struct rgeom { 
  short		form;		/* form */
  ParCurv	*de1;		/* first boundary curve */
  ParCurv	*de2;		/* second boundary curve */
  short		dirflg;		/* direction flag */
  short		devflg;		/* developable flag */
} RulSurf;

/* surface of revolution */

typedef struct {
  double	axis[6];	/* axis of revolution */
  ParCurv	*de;		/* generatrix curve */
  double	start;		/* start angle (radians) */
  double	term;		/* terminate angle (radians) */
} SurfRev;

/* tabulated cylinder */

typedef struct {
  ParCurv	*de;		/* directrix curve */
  double	lx;		/* terminate point of generatrix */
  double	ly;
  double	lz;
} TabCyl;

/* trimmed surface */

typedef struct {
  int		nCurves;
  ParCurv	**egeoms;
} TrimLoop;

typedef struct {
  ParSurf	*fgeom;
  short		bLoop;
  int		nLoops;
  TrimLoop	*loops;   /* B-curves */
  /********** 11.1 begin **********/
  TrimLoop	*c_loops; /* C-curves */
  /********** 11.1 end ************/
} TrimSurf;

/* constrained Delaunay triangulation */

typedef struct adj_elem {
  int adj;
  struct adj_elem *prev;
  struct adj_elem *next;
} AdjElem;

typedef struct {
  int nAdjs;
  int bound;
  AdjElem *first;
  AdjElem *last;
} AdjType;

/* facet structure */

typedef struct {
  double **pts;   /* coordinates */
  int nAdj,       /* size of adjacency array iAdj */
      nPts;       /* number of points */
  int   *iAdj,    /* adjacency array */
        *iEnd;    /* pointers to adjacencies for a point i */
} FacetSurf;

/* function prototypes */

void	Add(short, float *, float *, float *);
void	add_curv(ParCurv **, int, ParCurv *, int, int, int);
void	add_surf(ParSurf ***, int, int, int, int, ParSurf *);
void	AddAdj(int **, int, int, int, int *);

void    addcst_(int *, int *, int *, float *, float *, int *, int *, int *,
		int *, int *, int *);
short	addsections_canal(ParSurf *, ParSurf *, ParCurv **, double *,
		 int *, double *, int, double, int, double);
void	AddU(FILE *, double, int *);
int	AlphaShape(int, double **, int *, int *, double, int *, int *);
void	approx_cubic_curv1(ParCurv *, double **, double *, double *,
		int, int, short, double *);
short	approx_user_knots(ParCurv *, double **, double *, int, int,
		int, double *, double *);
double	arc_length(ParCurv *);
void	arr_vec(int, double **, vector **);
double	avg_arc_length(ParSurf *);
double	axis_and_terminate(ParSurf *, double *);
int	bezier_clip_decision_egeom(ParCurv *, int, double *, double *);
void	bezier_knots(int order, double, double, double *);
vector	*****BezierContpts(vector ***, short, short, short, short);
double	**BezierKnots(double *, short, short);
vector	*bias( vector *, vector *, vector *, double);
double	bias_factor(short, double, double, double);
double	bias_factor1(short, double, double, double);
int	BsplineToBezier(ParCurv *, ParCurv ***);
double	*BsplineToBezierKnots( double *, short, short, short *);
ParSurf	***BsplineToBezierSurf(ParSurf *, int *, int *);
short	build_offsets_int(ParSurf **, ParSurf *, double);
double	**build_w_table(vector **, double, int, double, double);
int	calc_convex_hull(int, int, double *, double *, double *, double *);
void	calc_disc(ParCurv *, double *, int *);
int	calc_interval_egeom(ParCurv *, int, double *, double *T);
int	calc_intsct(int, double *, double *, int, double *, double *,
		double *, int *);
ParSurf *canal_to_int(ParCurv *, short, short, double, double *, double *,
		ParCurv **, int *, int *, short, short, double, double);
ParSurf	*canal_to_rat(ParSurf *, double, int, int, int, double, double *,
		ParCurv **, int *, int, double);
void	check_cv(ParCurv *, ParCurv *, int *, double *, int *, int *);
void	check_one_dir(ParSurf ***, int, int, int *, int *, double *,
		int *, int *);
void    CheckOneDir(ParSurf ***, int, int, int, double *, double *, double *,
		    double *, double *, double *);
void	check_patch_u(ParSurf *, ParSurf *, int *, int *, double *,
		int *, int *);
int     CheckCurvOpen(ParCurv *);
int	checknot(ParCurv *);
void    CheckPatchContinuity(ParSurf *, ParSurf *, double, double *, double *,
			     double *, double *, double *);
void    CheckPatchU(ParSurf *, ParSurf *, int, double *, double *, double *,
		    double *, double *, double *);
int     CheckSurfOpen(ParSurf *, int *, int *);
void	Circumcenter(double *, double *, double *, double, double *);
double	Circumradius(double *, double *, double *, double);
ParCurv	*clip_bezier_curve(ParCurv *, double, double);
int	close_enough(double);
void	CloseLog(void);
int	com_knot_vector(ParSurf *, ParSurf *, double []);
int	common_interval(double, double, double, double, double *, double *);
void	compare_curves(ParCurv **, int, int *, double *, int *, int *);
void	compare_patches(ParSurf ***, int, int, int *, int *, double *,
		double *, int *, int *, int *, int *);
int	ComparAdj(const void *, const void *);
void	CompareBoundingBox(float *, float *);
void    ComparePatches(ParSurf ***, int, int, int, double *, double *,
		       double *, double *, double *, double *, double *,
		       double *, double *, double *, double *, double *);
double	**compute_transfm(double *);
ParCurv	*construct_bspl(ParCurv **, int, int, int);
ParSurf	*construct_jbspl(ParSurf ***, int, int, int, int, int);
int	ContourFacets(double **, int, int *, int *, double,
		ParCurv ***, int ***, int ***, int **, int *);
void	ConvertToCos(FILE *, int, ParCurv **);
double	ConvexHullSurf(ParSurf *, ParSurf *, int, int);
double	ConvexHullTest(ParCurv *, ParCurv *);
void	Copy(short, float *, float *);
void    copyAdj(int , int *, int *, int, int **, int **);
void	copy_contpts(vector ***, vector ***, short, short);
void	copy_knots(double *, double *, short);
ParCurv *copyegeom_leo(ParCurv *, ParCurv *);
FacetSurf *CopyFacetSurf(FacetSurf *, FacetSurf *);
GridSurf *CopyGridSurf(GridSurf *, GridSurf *);
ListCurv *CopyListCurv(ListCurv *, ListCurv *);
ParUv  *CopyParUv(ParUv *, ParUv *);
TrimSurf *CopyTrimSurf(TrimSurf *, TrimSurf *);
double *CopyUniqueKnots(double *, short, short *, short);
ParSurf	*corner_offset(ParSurf **, int, double);
void	Cross(float *, float *, float *);
int	CrossEdge(double **, int, int, double, double *);
ParCurv *cross_link(double, ParSurf **, ParCurv **, short, short, double,
		double);
int	CrossTriangle(double **, double, int *, int *, int *);
vector	*curv_cont(vector **, vector *, vector **, vector *);
void	curv_contpts0(ParCurv *, vector *);
void	curv_contpts1(ParCurv *, vector *, short);
void	curve_oslo_per(short, int, short, double *, double *, double **,
		double **);
void	curve_oslo2(short, int, short, double *, double *, double **,
		double **);
void	curve_osloper_lt(short, int *, short, double *, double *,
		double **, double **, double **);
void	curve_osloper_rt(short, int *, short, double *, double *,
		double **, double **, double **);
void	CurveConvexHull(ParCurv *, double *, double *, double *, double *);
ParCurv *CurvDeviation(ParCurv *, ParCurv *, int, int, int, int);
int	CurveOrientation(ParCurv *, double, int);
int	cvx_hull_axis_intrsct(int, double *, double *, int, double *,
		double *, double *, double, double, double *, double *);
ParSurf	*cyl_canal_to_rat(ParCurv *, short, short, double, double *,
		ParCurv **, int *, short, double, double);
vector	**cyl_frenet_tr(ParCurv *, double, double);
ParCurv	*cyl_generatrix(ParCurv *, double, double, double);
ParCurv	*cyl_generatrix_lk(double, double);
short	cyl_sample_gencyl(ParCurv *, short, short, double, ParSurf *,
		ParCurv **, double *, int *, double, double);
void	decasteljau_curve(ParCurv *, double, ParCurv **, ParCurv **);
void	decasteljau_curve_contpts(vector **, double, int, vector **,
		vector **, double, double);
short	Delaunay(int, double **, int *, int *, int **, int);
void    DelaunayToEdgeU(int, double **, int *, double ***, int **, int *, int **);
void    DelaunayToEdgeV(int, double **, int *, double ***, int **, int *, int **);
int     InBdryEdge(int, int, int **, int);
int	DelaunayCompare(const void *, const void *);
int	deriv_sign(double *, int, double);
double	Det(double, double, double, double);
void	dhsv2rgb(double, double, double, float*, float*, float*);
void	disc_bsp(short, short, short, double *, double *, double **);
float	Dot(float *, float *);
vector	*eval_LE2D(double, int);
vector	*eval_LE3D(double, int);
vector	*EvaluateCos(ParCurv *, float, short, ParSurf *, vector **,
		vector **, vector **, double *, double *, double *,
		double *);
vector	*EvaluateCurv(ParCurv *, float, short, vector **, vector **,
		double *);
vector	*EvaluateSurf(ParSurf *, float, float, vector **, vector **,
		vector **, double *, double *, double *, double *);
vector	*eval_curve_bounded(ParCurv *, double, int);
vector	*eval_vert(vector *, vector *, double *, vector *);
ParCurv	*extract_cv(ParCurv *, int, int *);
void	extract_knots(ParSurf *, int, int, int *, int *, ParSurf *);
void	extract_knots_cv(ParCurv *, int, ParCurv *, int *);
ParSurf	*extract_patch(ParSurf *, int, int, int *, int *);
ParCurv *ExtractLeadingEdge2D(ParSurf *, double *, short, short);
ParCurv *ExtractLeadingEdge3D(ParSurf *, double *, short, short);
void	fair_per_knot(ParCurv *, short);
ParSurf ***fgeom_array2(int, int);
int	FilterNegative(double, ParCurv *, double *);
void	find_deviation(ParCurv *, ParCurv *, double *);
int	FindBot(double **, int, double ***, int);
int	FindLft(double **, int, double ***, int);
int	FindRht(double **, int, double ***, int);
int	FindThirdPoint(int *, int *, int, int);
int	FindTop(double **, int, double ***, int);
double	*FindUniqueKnots(double *, short, short *);
ParCurv *FitFunction(int, double *, double *, double *, double *, int);
float	***flt_array3(unsigned, unsigned, unsigned);
void	free_earray1(ParCurv **, unsigned);
void	free_earray2(ParCurv ***, short, short);
void	free_egeom_array1(ParCurv **, int);
void	free_fgeom_array2(ParSurf ***, int, int);
void	free_farray3(float ***);
void	FreeBezierContpts(vector *****, short, short, short, short);
void    FreeFacetSurf(FacetSurf *);
void    free_trim(TrimSurf *);
vector	**frenet_tr(ParSurf *, double, double, int);
void	functcv(double *, double *, double *);
double	functcv_sqr(double);
double	*fundamentals(vector *, vector **);
ParCurv	*generatrix(ParSurf *, double, double, int, double);
ParCurv	*generatrix_lk(double, double);
void	get_cv_removal(ParCurv *, int, ParCurv *, int *);
double	get_error(ParSurf ***, int, int, ParSurf *, int, int);
double	get_error_cv(ParCurv **, int, ParCurv *, int, int);
void	get_iso_removal(int, int, ParSurf *, ParSurf *, int, int *);
int	GetAdj(int *, int *, int, int);
void	GetCurvBoundBox(ParCurv *, float *);
vector	*GetCurvTang(ParCurv *, float, int);
char	*GetFileName(char *);
void	GetGridBoundBox(GridSurf *, float *);
int	GetIntersect(double *, double *, double *, double *, double *,
		double *);
void	GetListBoundBox(ListCurv *, float *);
int	GetLogDebug(void);
short	GetLogKey(void);
float	GetMaxRange(float *);
char	*GetPathName(char *);
void	GetSurfBoundBox(ParSurf *, float *);
void	GetUvBoundBox(ParUv *ugeom, float *);
void	GetVectBoundBox(double **, float *, int);
void	implicitize_ray(double *, double *, double *);
void	init_tmat(double *, double, double);
void	InsertTriangle(double **, AdjType *, int, int, int);
ParCurv *int_curve(double **, short, short, short);
int	int_number(ParCurv *, double * , double *, double);
int	int_sample_gencyl(ParCurv *, short, short, double, ParSurf *,
		ParCurv **, double *, double *, int *, int *, double,
		double);
void	intcyl_calc_bounds_offs(ParSurf *, double *, double *, int *, int *);
int	integral_build_offset(ParSurf *, ParSurf *, double);
void	integral_calc_bounds_offs(ParSurf *, double *, double *,
			       int *, int *);
short	integral_offset(ParSurf *, ParSurf **, double, int *, double, int *);
ParSurf *integral_offset_surf(ParSurf *, double, double **, int *,
		double,	int *);
int	integral_sample_offset(ParSurf *, ParSurf *, short, short,
		double, double, double *);
void	interp_points(ParCurv *, double **, double *, double *, double *,
		double *, double *, double **, double *, int);
void	interp_points_per(ParCurv *, double **, double *, double **,
		double **, double **, double **,  double **, double *,
		int);
double	InterpolateZ(double **, double, double, int, int, int);
int	intersect_curve_to_axis(ParCurv *, int, int *, double *);
TrimSurf *InvertTrimSurf(TrimSurf *, TrimSurf *);
int	IsAbove(ParCurv *, double **, int, int *, int *, double);
int     IsCollinear2(float, float, float, float, float, float);
double	iso_curv0(ParCurv *);
double	iso_curv1(ParCurv *);
void	IsoLoopIntersect(int, double w, int, TrimLoop *, int,
		double *, int *, int);
double	*knot_factor0(short, double *);
double	*knot_factor1(short, short, double *);
short	knot_normalize(double *, short);
double	knot_rem_and_pert(ParSurf ***, int, int, ParSurf *, double, int);
void	knot_rem_cv(ParCurv **, int, ParCurv *);
/*** obsolete ***
void	lefunct1(int *, double *, double *);
*****************/
/***** for NAG Marks 19,20 **********/
void	lefunct1(int *, double *, double *, int *, double *);
/************************************/
void	linear_interpolate_vector(double, vector *, double, vector *,
		double, vector *, double, double, double);
void	linear_interpolate_w(double, double *, double, double, double, double);
vector	*LinearCosNormal(ParCurv *, double, double);
ParCurv *link_alloc(short, short, short *);
int	LocateTri(double **, int *, int *, double, double,
		int, int, int, int *, int *, int *);
int	LocateTri2(double **, int *, int *, double, double,
		   int, int, int, int *, int *, int *, int *);
ParSurf *loft_integral(ParCurv **, double *, double *, short, short, short);
int	loop_rayintersect(ParCurv **, int, double *, double *);
void    MakeFileName(const char *, const char *, char *);
int	MakeTriangle(double **, int *, int *, int *, double, double,
		double, int *, AdjType *, int *);
double	max_curve_curvature(ParCurv *, short);
double	max_pert(ParSurf ***, int, int, ParSurf *, double, double);
double	max_pert_cv(ParCurv **, int, ParCurv *, double, double);
short	merge_can_off(ParSurf **, ParSurf **);
short	merge_canals(ParSurf *, ParSurf **, ParCurv ***, double **,
		double **, int *, int *, double);
void	merge_tol_edges(ParSurf *, ParSurf **, ParSurf **, double);
void	MergeKnots(ParCurv *, ParCurv *, short, ParCurv **, int *, int *);
vector	*midpoint(vector *, vector *, vector *);
int	min_bez_curve_to_axis_intrsct(double, double, double *, int *);
void	nagparam_(int *);	
double	next_u(ParCurv *, double);
int	NextAdjB(int *, int *, int, int);
short	NextTriangle(double **, int, int *, int *, int *, int *,
		double, int, int *, int *, int);
int	NoIntersects(double **, int *, int *, int, int, int, int);
int	NonConvexEdge(double **, int *);
double	normal_curvature(double *, double *);
void	Normalize(short, float *);
void	normalize_egeom(ParCurv *);
void	OpenLog(char *, char *);
void	opt_param(ParCurv *, double **, double *, short);
void	opt_param_per(ParCurv *, double **, double *, short);
int	Overlap(double *, double *, double *, double *);
SurfRev *ParSurf_to_SurfRev(ParSurf *, SurfRev *);
TabCyl	*ParSurf_to_TabCyl(ParSurf *, TabCyl *);
double	**partial(short, double *);
double	*partan_dir(vector *, vector *, vector *);
double	partl_u_functcv_sqr(double);
ParSurf	*patch_with_com_knots(ParSurf *, double [], int);
void	perturb_knots(ParSurf *, int, int, int *, int *, double);
void	perturb_knots_cv(ParCurv **, int, ParCurv *, double);
double	point_bspl(ParCurv *, double);
double	position_err(ParCurv *, ParSurf *, double, double);
void	post4x4(vector *, vector *);
void	principal_curv(double *, double *);
void	PrintLog(char *);
void	prodknot(ParCurv *, double *);
void	PruneFacets(int, double **, int **, int **, int, int, int **);
void	rainbow(double, double, double , float *, float *, float *);
void    RaiseAndMerge(ParSurf *, ParSurf *, ParSurf **, ParSurf **);
FacetSurf *ReadFacetSurf(FILE *, FacetSurf *);
GridSurf *ReadGridSurf(FILE *, GridSurf *);
ListCurv *ReadListCurv(FILE *, ListCurv *);
ParUv	*ReadParUv(FILE *, ParUv *);
TrimSurf *ReadTrimSurf(FILE *, TrimSurf *);
ParCurv	*recover_cv(ParCurv *, int, int *, int);
ParSurf	*recover_patch(ParSurf *, int, int, int *, int *, int);
int	RecursiveInt(ParCurv *, double, double, int);
ParCurv *reflectegeom(ParCurv *, ParCurv *, int);
ParSurf *reflectfgeom(ParSurf *, ParSurf *, int);
GridSurf *ReflectGridSurf(GridSurf *, GridSurf *, int);
ListCurv *ReflectListCurv(ListCurv *, ListCurv *, int);
double	**ReflectPts(double **, double **, int, int);
TrimSurf *ReflectTrimSurf(TrimSurf *, TrimSurf *, int);
void	rem_1st_knot(ParSurf *, int, int, int *, int *, ParSurf *);
void	renew_knots(ParSurf ***, int, int);
void	RevCurvParam(ParCurv *);
void	RevSurfUParam(ParSurf *);
void	RevSurfVParam(ParSurf *);
int	RightOf(double, double, double, double, double, double);
void	Rotate(float *, float *, float, float *);
void	rotate_section(vector **, double);
ParSurf	*RulSurf_to_ParSurf(RulSurf *, ParSurf *);
short	sample_blend(ParCurv *(*)(double, ParSurf **,ParCurv **,
		short, short, double, double), ParSurf *, double *,
		double *, short *, short, short, short, ParSurf **,
		ParCurv **, double, double);
short	sample_gencyl(ParSurf *, double, int, int, int, double, ParSurf *,
		ParCurv **, double *, int *, double);
short	sample_mod_surf(ParSurf *, ParSurf **, ParSurf **, int *,
		double, double);
void    SampleContinuity(ParSurf *, ParSurf *, int, double ***, int *,
			 double **);
void	Scale(short, float *, float, float *);
int     ScramblePoints(int, float *, float *, int, int *);
void	SearchBot(FILE *fp, double **, int, double, double, double,
		double,	int);
void	SearchLft(FILE *fp, double **, int, double, double, double,
		double,	int);
void	SearchRht(FILE *fp, double **, int, double, double, double,
		double,	int);
void	SearchTop(FILE *fp, double **, int, double, double, double,
		double,	int);
void	SegmentMinMax(double *, double *, double *, double *,
		double *, double *);
void	SetLogDebug(int);
void	SetLogKey(short);
int	signchange_ex(double *, int, double);
double	solve_int(ParCurv *, double *);
double	solve_u(ParCurv *, double, double);
void	spline_par(short, short *, short *);
void	split_bspl(ParCurv *, ParCurv *, ParCurv *);
void	split_it(double, ParCurv *, ParCurv *, ParCurv *);
void	SplitSurf(ParSurf *, int, double, ParSurf **, ParSurf **);
double	split_val(double);
void    stitch_facet (int, double **, int **, int *, int);
void	Sub(short, float *, float *, float *);
void	subdiv_bez(ParSurf *, ParSurf *, int, int);
void	subdiv_cv(ParCurv *, ParCurv *, int);
void	SubdivSurfNxM(ParSurf *, int, int, double *, double *, ParSurf ***);
vector	*surf_data(ParCurv *, ParSurf *, double, vector **);
vector	*surf_normal(vector *, vector *, vector *);
void    SurfDeviation(ParSurf *, ParSurf *, int, int, int, int *, double ***,
		      int **, int *, int **);
int	surfrev_compare(double *, double *);
ParSurf *SurfRev_to_ParSurf(SurfRev *, ParSurf *);
ParSurf	*swapsurf(ParSurf *);
void    SwitchFacetSurf(FacetSurf *);
void    SwitchFacetSurf2(int, double **, int *, int *);
void	SwitchSurfParam(ParSurf **);
ParSurf *TabCyl_to_ParSurf(TabCyl *, ParSurf *);
double	test_curvature0(ParCurv *, ParCurv *);
double	test_curvature1(ParCurv *, ParCurv *);
double	test_tangent(ParSurf *, vector *, double, double);
void    toFacet(int, int *, int, int *, int *, int *, int **, int **);
int	ToLeft(double *, double *, double *);
void	tol_calc_bandwidth(double *, double *, int, int, int *, int *);
ParSurf	*tol_loft_rational(ParCurv **, double *, short, int);
short	tol_sample_offsets(ParSurf *, ParSurf **, int, int, double, double);
short	tolerance_region(ParSurf *, ParSurf **, ParSurf **, ParSurf **,
		double, int *, int *, double, double, int *, int *, int *);
int	TraceContour(FILE *, double **, int, int *, int *, int,
		int *, double, int *);
void	TraceToBoundary(FILE *, double **, int, int *, int *,
		int *, int *, double, int *, int *, int *);
int	TraceToNext(FILE *, double **, int, int *, int *, int *,
		int *, double, int, int *, int *, int);
void	TraceToStart(FILE *, double **, int, int *, int *, int *,
		int *, double, int, int, int, int);
void	trans_per(short, short, short, double *, double *, double **);
ParCurv *transformegeom(ParCurv *, ParCurv *, double *);
ParSurf *transformfgeom(ParSurf *, ParSurf *, double *);
GridSurf *TransformGridSurf(GridSurf *, GridSurf *, double *);
ListCurv *TransformListCurv(ListCurv *, ListCurv *, double *);
double	**TransformPts(double **, double **, double *, int);
TrimSurf *TransformTrimSurf(TrimSurf *, TrimSurf *, double *);
int	TriFunction(int, int *, int *, int *, int, int, int);
int	TrimCompare(const void *, const void *);
int	TrimCompareXY(const void *, const void *);
void    TrimPoints(TrimSurf *, int, int, int, char *, int *, int *, int);
void    trmesh_(int *, float *, float *, int *, int *, int *, int *, int *);
void    trprnt (int, int *, int, float *, float *, int *, int *, int *);
void	unify_knots(ParSurf ***, int, int, ParSurf *);
void	unify_knots_cv(ParCurv **, int, ParCurv *);
void    updateAdj (int, int, int, int, int, int, int, int *, int *);
double	vec_angle(vector *, vector *);
void	vec_arr(int, vector **, double **);
int	vec_neq(vector *, vector *);
vector	*vec_rot(double, double, vector **);
vector	*vec_rot_fix(double, double, vector *, vector *);
int	vect_colinear(vector *, vector *);
double	vectang(vector *, vector *);
vector	*vector_dircos(vector *, vector *);
vector	*vector_rotate(vector *, double, double, double, vector *);
void	WriteFacet(FILE *, int, double **, int *, int, int *, short);
void	WriteFacetSurf(FILE *, FacetSurf *, short);
void	WriteGridSurf(FILE *, GridSurf *, short);
void    WriteInsp(FILE *, int, double **, double **, short);
void	WriteListCurv(FILE *, ListCurv *, short);
void	WriteParCurv2(FILE *, ParCurv *, short);
void	WriteParCurv2_Per(FILE *, ParCurv *, short);
void	WriteParSurf2(FILE *, ParSurf *, short);
void	WriteParUv(FILE *, ParUv *, short);
void    WriteTrimSurf(FILE *, TrimSurf *, short);

#ifdef __cplusplus
}
#endif

#endif

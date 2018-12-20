/* Copyright (C) 1992 Massachusetts Institute of Technology, Cambridge, MA */
/* All rights reserved */

/* praxiteles.h */

#ifndef PRAX_H
#define PRAX_H

#include <stdio.h>
#include "deslab.h"
#include "gen.h"
#include "editor.h"
#include "iges.h"

#ifdef __cplusplus
extern "C" {
#endif

/*********** 11.1 begin **********/
#define MAX_STRING 512
/*********** 11.1 end ************/

/* add prototypes not in strict ANSI library headers */

char *strdup(const char *);
int strcasecmp(const char *, const char *);
int strncasecmp(const char *, const char *, size_t);

#ifdef DEFINE_H

#define EXTERN
#define INIT(a)			= a
#define INIT2(a,b)		= {a,b}
#define INIT5(a,b,c,d,e)	= {a,b,c,d,e}
#define INIT6(a,b,c,d,e,f)	= {a,b,c,d,e,f}
#define INIT7(a,b,c,d,e,f,g)	= {a,b,c,d,e,f,g}
#define INIT8(a,b,c,d,e,f,g,h)	= {a,b,c,d,e,f,g,h}

#else

#define EXTERN extern
#define INIT(a)
#define INIT2(a,b)
#define INIT5(a,b,c,d,e)
#define INIT6(a,b,c,d,e,f)
#define INIT7(a,b,c,d,e,f,g)
#define INIT8(a,b,c,d,e,f,g,h)

#endif

/* message types */

#define PRAX_INFO       DESLAB_INFO
#define PRAX_WARNING	DESLAB_WARNING
#define PRAX_ERROR	DESLAB_ERROR

/* entity types */

#define PRAX_NONE	DESLAB_NONE
#define PRAX_LIST	DESLAB_LIST
#define PRAX_GRID	DESLAB_GRID
#define PRAX_CURV	DESLAB_CURV
#define PRAX_SURF	DESLAB_SURF
#define PRAX_COS	DESLAB_COS
#define PRAX_IGES	DESLAB_IGES
#define PRAX_MULT	DESLAB_MULT
#define PRAX_LOGO       DESLAB_LOGO
#define PRAX_DIST	DESLAB_DIST
#define PRAX_VECT	DESLAB_VECT
#define PRAX_FUNC	DESLAB_FUNC
#define PRAX_UV		DESLAB_UV
#define PRAX_LOCAL	DESLAB_LOCAL
#define PRAX_FACET	DESLAB_FACET
#define PRAX_TRIM	DESLAB_TRIM
#define PRAX_GROUP      DESLAB_GROUP
#define PRAX_HULL	DESLAB_HULL
#define PRAX_DIST_INSP  DESLAB_DIST_INSP
#define PRAX_LOCAL_INSP DESLAB_LOCAL_INSP
#define PRAX_SCRIPT     DESLAB_SCRIPT
#define PRAX_LOG        DESLAB_LOG
#define PRAX_PS         DESLAB_PS
#define PRAX_INSP       DESLAB_INSP

/* primitive curves */

#define PRAX_LINE       0
#define PRAX_CIRCLE     1
#define PRAX_ARC        2
#define PRAX_CONIC      3

/* primitive surfaces */

#define PRAX_PLANE	0
#define PRAX_CYLINDER	1
#define PRAX_CONE	2
#define PRAX_SPHERE	3
#define PRAX_TORUS	4
	
/* view pad labels */

#define PRAX_VIEW_ROTX	0
#define PRAX_VIEW_ROTY	1
#define PRAX_VIEW_ROTZ	2
#define PRAX_VIEW_PX	3
#define PRAX_VIEW_PY	4
#define PRAX_VIEW_PZ	5
#define PRAX_VIEW_DIST	6
#define PRAX_VIEW_FOV	7
#define PRAX_VIEW_NEAR	8
#define PRAX_VIEW_FAR	9

/* view tools */

#define PRAX_TOOL_NONE	-1
#define PRAX_TOOL_ROTX	0
#define PRAX_TOOL_ROTY	1
#define PRAX_TOOL_ROTZ	2
#define PRAX_TOOL_TRXY	3
#define PRAX_TOOL_TRZ	4
#define PRAX_TOOL_DIST	5
#define PRAX_TOOL_FOV	6
#define PRAX_TOOL_NEAR	7
#define PRAX_TOOL_FAR	8
#define PRAX_TOOL_RESET	9
#define PRAX_TOOL_WORKING 9

/* constants */

#define PRAX_DEG_TO_RAD 0.017453292519943295
#define PRAX_RAD_TO_DEG 57.295779513082323
#define PRAX_EMAX	1.0e30
#define PRAX_EMIN	-1.0e30
#define PRAX_EPSILON	0.001
#define PRAX_MINDIST    0.01
#define PRAX_ZERO6      1.0e-6
#define PRAX_ZERO10	1.0e-10

/* floating point format */

#define PRAX_F_FORMAT	0
#define PRAX_G_FORMAT	1

/* point font */

#define PRAX_DEFAULT_FONT	0
#define PRAX_POINT_FONT		1
#define PRAX_POINT		"\001"

/* line styles */

#define PRAX_SOLID	0
#define PRAX_DASHED	1

#define PRAX_YES	1
#define PRAX_NO		0

/* flags */

#define PRAX_SELECTED	0x0001
#define PRAX_IS_SAVED	0x0002
#define PRAX_IS_IGES	0x0004
#define PRAX_IS_COPY    0x0008
#define PRAX_IS_OPEN	0x0010
#define PRAX_IS_POWER	0x0020
#define PRAX_IS_RULED	0x0040
#define PRAX_IS_TAB_CYL	0x0080
#define PRAX_SURF_REV	0x0100
#define PRAX_IS_OFFSET	0x0200
#define PRAX_INTEGRAL	0x0400
#define PRAX_IS_LOFTED  0x0800
#define PRAX_ISOPHOTES	0x1000
#define PRAX_REFLECTION	0x2000
#define PRAX_GEODESICS	0x4000
#define PRAX_CURV_LINES	0x8000
#define PRAX_IS_KMIN	0x10000
#define PRAX_FUNCTION	0x20000
#define PRAX_IS_CLOSED  0x40000
#define PRAX_IS_PLANAR  0x80000
#define PRAX_UV_EDIT    0x100000
#define PRAX_CLOSED_U   0x200000
#define PRAX_CLOSED_V   0x400000

/* style flags */

#define PRAX_SELF	0x00000001
#define PRAX_POLYGON	0x00000002
#define PRAX_PORCUPINE	0x00000004
#define PRAX_FOCAL	0x00000008
#define PRAX_RADIAL	0x00000010
#define PRAX_EVALUATE	0x00000020
#define PRAX_OFFSET	0x00000040
#define PRAX_SHADED	0x00000080
#define PRAX_CURVATURE	0x00000100
#define PRAX_MIN_DIST	0x00000200
#define PRAX_MIN_PROJ	0x00000400
#define PRAX_ORTHOTOMIC	0x00000800
#define PRAX_TORSION	0x00001000
#define PRAX_HEIGHT	0x00002000
#define PRAX_PIVOT	0x00004000
#define PRAX_MAP_TO_SURF 0x0008000
#define PRAX_SPLIT	0x00010000
#define PRAX_TWO_SIDED	0x00020000
#define PRAX_OFFSET_1	0x00040000
#define PRAX_OFFSET_2	0x00080000
#define PRAX_ORDINAL    0x00100000
#define PRAX_POLYMESH   0x00200000
#define PRAX_SLOPE      0x00400000
#define PRAX_NON_CONVEX 0x00800000
#define PRAX_KNOTS      0x01000000
#define PRAX_CONTINUITY 0x02000000

#define PRAX_COS_Z      1
#define PRAX_COS_RANGE  2

/* colors */

#define PRAX_BLACK	0x00000000
#define PRAX_RED	0x000000ff
#define PRAX_GREEN	0x0000ff00
#define PRAX_YELLOW	0x0000ffff
#define PRAX_BLUE	0x00ff0000
#define PRAX_MAGENTA	0x00ff00ff
#define PRAX_CYAN	0x00ffff00
#define PRAX_WHITE	0x00ffffff
#define PRAX_GRAY	0x007f7f7f

/* lighting model materials */

#define PRAX_LIGHT1	1
#define PRAX_MAT_BLACK	1
#define PRAX_MAT_RED	2
#define PRAX_MAT_GREEN	3
#define PRAX_MAT_YELLOW	4
#define PRAX_MAT_BLUE	5
#define PRAX_MAT_MAGENTA 6
#define PRAX_MAT_CYAN	7
#define PRAX_MAT_WHITE	8

/* curvatures */

#define PRAX_GAUSSIAN	0
#define PRAX_MEAN	1
#define PRAX_ABSOLUTE	2
#define PRAX_NORMAL_U	3
#define PRAX_NORMAL_V	4
#define PRAX_MAX_PRINCIPAL 5
#define PRAX_MIN_PRINCIPAL 6
#define PRAX_MIN_RADII	7
#define PRAX_MAX_RADII	8
#define PRAX_RMS        9

/* contour ramping */

#define PRAX_LINEAR	0
#define PRAX_QUADRATIC	1
#define PRAX_CUBIC	2
#define PRAX_QUARTIC	3
#define PRAX_QUINTIC    4

/* curve/surface fitting */

#define PRAX_APPROXIMATE 0
#define PRAX_INTERPOLATE 1
#define PRAX_ISOPARAMETER 2
#define PRAX_PROJECTION  3
#define PRAX_LEADING_EDGE 4
#define PRAX_INTERP_LOFT 5

#define PRAX_CHORD	0
#define PRAX_UNIFORM	1
#define PRAX_FOLEY	2
#define PRAX_HARTLEY    3

#define PRAX_U 0
#define PRAX_V 1

#define PRAX_FROM_KEYBOARD	1
#define PRAX_FROM_FILE		2

#define PRAX_X	0
#define PRAX_Y	1
#define PRAX_Z	2

#define PRAX_ADD	0
#define PRAX_INSERT	1
#define PRAX_DELETE	2
#define PRAX_MOVE	3
#define PRAX_RESET	4

#define PRAX_DEFAULT_STORE	50

#define PRAX_IGES_100   1
#define PRAX_IGES_102   2
#define PRAX_IGES_104   3
#define PRAX_IGES_106	4
#define PRAX_IGES_108   5
#define PRAX_IGES_110   6
#define PRAX_IGES_112	7
#define PRAX_IGES_114	8
#define PRAX_IGES_116   9
#define PRAX_IGES_118  10
#define PRAX_IGES_120  11
#define PRAX_IGES_122  12
#define PRAX_IGES_126  13
#define PRAX_IGES_128  14
#define PRAX_IGES_132  15
#define PRAX_IGES_142  16
#define PRAX_IGES_144  17

/* system directory definition */

#define PRAX_COMPILED_PATH     1
#define PRAX_ENVIRONMENT_PATH  2
#define PRAX_COMMAND_LINE_PATH 3

/* zeros */

#define PRAX_ZERO_CURVATURE 1
#define PRAX_ZERO_DISTANCE  2
#define PRAX_ZERO_KNOT      3

/* errors */

#define PRAX_ERR_MAIN_LOOP	-1
#define PRAX_ERR_COMMAND_LINE	-2

/* system options */

#define PRAX_NUM_OPTIONS 1   /* number of system options */
#define PRAX_ORIENTED    1   /* allow oriented distance in localization */

/* tool bar button structure */

typedef struct {
  char		*label;		/* button label */
  short		iTool;		/* identifier */
} ButtonStruct;

/* external data */

/* function prototypes */

/* praxiteles program */

void	Background(void);
void    CamberLineDeviations(double, double, double, double, double, double,
			     double, double);
void    CamberLineResults(char *, vector *, vector *, vector *, double,
			  double, double, double, double, double, double,
			  double, double, double, double, double, int, int,
			  int, int, double);
void	CheckCosContinuity(int, int);
void	CheckCurvContinuity(int, int);
void	CheckSurfContinuity(int, int);
void    ClearScriptFile(void);
void	ClearView(void);
short	ColorToMaterial(long);
void	ConlocResults(int, int, int, int, double, int, double,
		      int, double, int, double, int, double *, int);
long	ContourRGB(float, float, float, short, short);
int     ConvertNonOpen(void);
void 	CountIters(int);
void	CycleColor(long *);
void	DeleteCos(short);
void	DeleteContourDialogs(void);
void	DeleteContourObjs(void);
void	DeleteCurv(short);
void	DeleteEvalObjs(void);
void	DeleteFacet(short);
void	DeleteGrid(short);
void	DeleteGroup(short);
void	DeleteHull(short);
void    DeleteInsp(short);
void	DeleteList(short);
void	DeleteStore(void);
void	DeleteSurf(short);
void	DeleteTrim(short);
void	DeleteUv(short);
void	DeleteVect(short);
void	DisplayBox(char *, char *, int, char **, int, char ***, int *);
void	DrawCos(short, float **);
void	DrawCosEvaluate(ParCurv *, short, float, ParSurf *, float);
void	DrawCosFocal(short, float, float **, float **, double *);
void	DrawCosKnots(ParCurv *, ParSurf *, int);
void	DrawCosPolygon(ParCurv *, ParSurf *);
void	DrawCosPorcupine(short, float, float **, float **, double *);
void	DrawCurv(short, float **);
void	DrawCurvEvaluate(ParCurv *, int, int,  float, float);
void	DrawCurvFocal(short, float, float **, float **, double *);
void	DrawCurvFunction(ParCurv *, short, short, short, char *);
void	DrawCurvKnots(ParCurv *, int);
void    DrawCurvOffset(ParCurv *, short, float, float **,float **);
void	DrawCurvOrthotomic(short, float, float *, float **, float **);
void	DrawCurvPolygon(ParCurv *);
void	DrawCurvPorcupine(short, float, float **, float **, double *);
void	DrawCurvRadial(short, float, float *, float **, double *);
void	DrawCurvTorsion(ParCurv *, short, short, float, float **, float **);
void	DrawFacetEvaluate(double, double, double, double);
void	DrawFacetEvaluateMap(double, double, double, double, ParSurf *);
void	DrawFacetHeight(int, double **, int *, int, int *, int,
		int, int, double, double, int, double *, double *,
		double, double, int);
void	DrawFacetHeightMap(int, double **, int *, int, int *,
		 int, int, int, double, double, int, double *,
		 double *, ParSurf *, double, double, int);
void	DrawFacetShaded(int, double **, int *, int, int *);
void	DrawFacetShadedMap(int, double **, int *, int, int *,
		ParSurf *);
void	DrawFacetSlope(int, double **, int *, int, int *, int,
		int, int, double, double, int, double *, double *,
		double, double, int);
void	DrawFacetSlopeMap(int, double **, int *, int, int *,
		 int, int, int, double, double, int, double *,
		 double *, ParSurf *, double, double, int);
void	DrawFacetWireframe(int, double **, int *, int, int *);
void	DrawFacetWireframeMap(int, double **, int *, int, int *,
		 ParSurf *);
void	DrawGrid(GridSurf *);
void	DrawGridHeight(int, int, double ***, short, short, short,
		double, double, short, double *, double *, double,
		double, short);
void	DrawGridPoints(double ***, int, int);
void	DrawGridPolygon(double ***, int, int);
void	DrawHullShaded(double **, int, int **);
void	DrawHullWireframe(double **, int, int **, int *, int);
void    DrawInsp(int, double **, double **, double);
void    DrawInspPolygon(int, double **);
void	DrawList(ListCurv *);
void	DrawListHeight(int, double **, short, short, short, double,
		double, short, double *, double *, double, double, short);
void	DrawListOrdinal(int, double **, short, short, short, double,
		double, short, double *, double *, double, double, short);
void	DrawListPoints(double **, int);
void	DrawListPolygon(double **, int);
void	DrawLocalAxis(float *, float, float, short);
void	DrawSurfCurv(short, short, short, short, float ***, double **,
		double **, double **, double **, short, double, double,
		short, double *, double *, double, double, short);
void	DrawSurfEvaluate(ParSurf *, float, float, float);
void	DrawSurfFocal(ParSurf *, short, short, short, short, float,
		float ***, float ***, double **, short);
void	DrawSurfHeight(short, short, short, float ***, short, short,
		double, double, short, double *, double *, double,
		double, short);
void    DrawSurfKnots(ParSurf *, int, int, int, int);
void	DrawSurfOffset(ParSurf *, short, short, short, short, float,
		float ***, float ***);
void	DrawSurfOffset1(short, short, float, float ***, float ***);
void	DrawSurfOffset2(short, short, float, float ***, float ***,
		long, float *, float, float, float, int);
void	DrawSurfOrthotomic(ParSurf *, short, short, short, short nsubv,
		float scal, float *, float ***, float ***);
void	DrawSurfPivot(ParSurf *, float, float, float, float);
void	DrawSurfPolygon(ParSurf *);
void	DrawSurfShaded(short, short, float ***, float ***);
void	DrawSurfSplits(ParSurf *, int, int, double *, double *, int, int);
void	DrawSurfTwoSided(short, short, float ***, float ***, long, float *,
		float, float, float, int);
void	DrawSurfWireframe(ParSurf *, short, short, short, short, float ***);
void	DrawTrimCurv(int, int *, int *, int, int, float **,
		double *, double *, double *, double *, int, double,
		double, int, double *, double *, double, double, int);
void	DrawTrimHeight(int, int *, int *, int, float **,
		int, int, double, double, int, double *, double *,
		double, double, int);
void	DrawTrimShaded(int, int *, int *, float **, float **);
void	DrawTrimWireframe(TrimSurf *, int, int, int);
void	DrawUv(ParUv *, ParSurf *);
void	DrawUvPolygon(ParUv *, ParSurf *);
void	DrawVect(int, double **, char *);
void	DrawViewingAxis(void);
void	DrawWorldAxis(void);
void	DrawLogo(void);
void    EditContSurfFile(int);
void    EditContSurfOk(void);
void    EditContSurfScript(char *);
void    EditContSurfUpdate(double);
void    EvalCosScript(char *);
void    EvalCurvScript(char *);
void    EvalFacetScript(char *);
void    EvalSurfScript(char *);
void    EvalTrimScript(char *);
void	Facet(int, int, int, char **);
void	FacetToEdge(int, int, int, char **);
short	FacetTriangles(int, double **, int *, int, int *);
char   *FileBox(char *, int, int *);
double	FloatBox(char *, char *, double, char *, double, char *, double,
		int *);
void    FormBox(char *, int, char **);
float	GetAspectRatio(long *, long *);
long	GetBackground(void);
int     GetBatchMode(void);
char   *GetColorName(long);
short	GetContourScale(void);
vector *GetCosValues(ParCurv *, float, short, ParSurf *, vector **,
		double *, double *);
short	GetCurrentTool(void);
/*
vector *GetCurvValues(ParCurv *, float, int, int, vector **, vector **,
		double *, double *);
*/
vector *GetCurvValues(ParCurv *, double, int, int, vector **, vector **,
		double *, double *);
char   *GetCwd(void);
int	GetDefaultIgesEntity(int);
short	GetDisplayType(void);
void	GetEntitiesScale(short *, float *, float *, float *, float *,
		float *, float *, short *, float *);
char   *GetFilter(short);
int     GetFirstExpose(void);
char   *GetFitName(short);
long	GetForeground(void);
short	GetFormat(void);
char   *GetHostName(void);
short	GetLocalAxis(void);
char   *GetLogFile(void);
short	GetLogoIters(void);
short	GetMaxCoss(void);
double	GetMaxContour(void);
short	GetMaxCurvs(void);
double	GetMaxDistance(void);
short	GetMaxFacets(void);
short	GetMaxGrids(void);
short	GetMaxGroups(void);
short	GetMaxHulls(void);
double	GetMaxHeight(void);
short   GetMaxInsps(void);
short	GetMaxLists(void);
double	GetMaxOrdinal(void);
double	GetMaxSlope(void);
short	GetMaxSurfs(void);
short	GetMaxTrims(void);
short	GetMaxUvs(void);
short	GetMaxVects(void);
double	GetMinContour(void);
double	GetMinDistance(void);
double	GetMinHeight(void);
void	GetMinMaxContour(double *, double *);
void	GetMinMaxDistance(double *, double *);
void	GetMinMaxHeight(double *, double *);
void	GetMinMaxOrdinal(double *, double *);
void	GetMinMaxSlope(double *, double *);
double	GetMinOrdinal(void);
double	GetMinSlope(void);
short	GetNumCoss(void);
short	GetNumCurvs(void);
short	GetNumFacets(void);
short	GetNumGrids(void);
short	GetNumGroups(void);
short	GetNumHulls(void);
short   GetNumInsps(void);
short	GetNumLists(void);
short	GetNumSurfs(void);
short	GetNumTrims(void);
short	GetNumUvs(void);
short	GetNumVects(void);
char   *GetParamName(short);
short   GetPostScriptImage(void);
char   *GetPostScriptName(void);
char   *GetPraxDirectory(void);
char   *GetPrimCurvName(short);
char   *GetPrimSurfName(short);
double	GetRadiusTolerance(void);
char   *GetRealName(void);
short	GetRevContour(void);
void	GetSampleSize(int *, int *);
float	GetScreenToWorld(void);
void    GetScriptDialog(char *, int);
char   *GetScriptFile(void);
FILE   *GetScriptFp(void);
int     GetScriptInt(void);
double  GetScriptReal(void);
void    GetScriptString(char *);
short	GetSelectedCoss(void);
short	GetSelectedCurvs(void);
short	GetSelectedFacets(void);
short	GetSelectedGrids(void);
short	GetSelectedGroups(void);
short	GetSelectedHulls(void);
short   GetSelectedInsps(void);
short	GetSelectedLists(void);
short	GetSelectedSurfs(void);
short	GetSelectedTrims(void);
short	GetSelectedUvs(void);
short	GetSelectedVects(void);
short	GetShowLogo(void);
short	GetShowMotd(void);
vector *GetSurfValues(ParSurf *, float, float, vector **, vector **,
		vector **, double *, double *, double *, double *);
int     GetSystemOptions(int);
char   *GetToolName(short);
int	GetTransformValues(char *, double *);
int     GetTrimSurfTry(void);
char   *GetUserName(void);
int     GetVectorScale(void);
short	GetViewAxis(void);
float	GetViewDist(void);
float	GetViewFar(void);
short	GetViewFov(void);
float	GetViewMaxRange(void);
float	GetViewNear(void);
float	GetViewPadInc(short);
char   *GetViewPadName(short);
float	GetViewRotX(void);
float	GetViewRotY(void);
void	GetViewRP(float *);
short	GetViewTwist(void);
void	GetViewUP(float *);
void	GetViewVP(float *);
void	GetViewVR(float *);
float	GetViewX(void);
float	GetViewY(void);
float	GetViewZ(void);
short	GetWorldAxis(void);
float	GetWorldToScreen(void);
double  GetZero(int);
void    GlxInput(int, double);
void	InitializeDefault(void);
void	InitializeFilters(void);
void	InitializeLighting(void);
void	InitializeLogo(void);
void	InitializePraxiteles(void);
void	InitializeStore(void);
void	InitializeWindows(unsigned, char **);
void    InitSystemOptions(void);
long	IntBox(char *, char *, long, char *, long, char *, long, int *);
void    InterrCosScript(char *);
void    InterrCurvScript(char *);
void    InterrFacetScript(char *);
void    InterrGridScript(char *);
void    InterrGroupScript(char *);
void    InterrHullScript(char *);
void    InterrInspScript(char *);
void    InterrListScript(char *);
void    InterrSurfScript(char *);
void    InterrTrimScript(char *);
void    InterrUvScript(char *);
short	IsViewInitialized(void);
void	JoinCoss(int, int);
void	JoinCurvs(int, int);
void	JoinSurfs(int, int);
void	LogCommandLine(int, int, int);
void	LogPraxDirectory(int, char *);
/* int	main(unsigned, char **); */
void	MainEventLoop(void);
void	MeasuredPoints(int, int);
void	MessageBox(short, char *, char *);
void	MessageOfTheDay(void);
void	MinDistance(int, int, int);
void	MinDistToList(int, int, int, char **);
void    MinDistTrim(int, int, int);
void	MultiFloatBox(char *, char *, int, char *[], double *, int, char **,
		      int *, int *);
void    MultiStringBox(char *, char *, int, char **, char **, int *);
void	OrthoDistance(int, int);
int	ParseCommandLine(int, char **);
void	PostScriptCos(short, float **);
void	PostScriptCosEvaluate(ParCurv *, short, float, ParSurf *, float);
void	PostScriptCosFocal(short, float, float **, float **, double *);
void	PostScriptCosKnots(ParCurv *, ParSurf *, int);
void	PostScriptCosPolygon(ParCurv *, ParSurf *);
void	PostScriptCosPorcupine(short, float, float **, float **, double *);
void	PostScriptCurv(short, float **);
void	PostScriptCurvEvaluate(ParCurv *, int, int, float, float);
void	PostScriptCurvFocal(short, float, float **, float **, double *);
void	PostScriptCurvFunction(ParCurv *, short, short, short, char *);
void	PostScriptCurvKnots(ParCurv *, int);
void    PostScriptCurvOffset(ParCurv *, short, float, float **, float **);
void	PostScriptCurvOrthotomic(short, float, float *, float **, float **);
void	PostScriptCurvPolygon(ParCurv *);
void	PostScriptCurvPorcupine(short, float, float **, float **, double *);
void	PostScriptCurvRadial(short, float, float *, float **, double *);
void	PostScriptCurvTorsion(ParCurv *, short, short, float, float **,
		float **);
void	PostScriptEntities(void);
void	PostScriptFacet(int, double **, int *, int, int *);
void	PostScriptFacetEvaluate(double, double, double, double);
void	PostScriptFacetEvaluateMap(double, double, double, double, ParSurf *);
void	PostScriptFacetMap(int, double **, int *, int, int *,
		ParSurf *);
void	PostScriptFrame(short, short, short, short);
void	PostScriptGrid(GridSurf *);
void	PostScriptGridHeight(int, int, double ***, short, short, short,
		double, double, short, double *, double *, double,
		double, short);
void	PostScriptGridPoints(double ***, int, int);
void	PostScriptGridPolygon(double ***, int, int);
void    PostScriptImage(char *);
void    PostScriptInsp(int, double **, double **, double);
void    PostScriptInspPolygon(int, double **);
void	PostScriptList(ListCurv *);
void	PostScriptListHeight(int, double **, short, short, short, double,
		double, short, double *, double *, double, double, short);
void	PostScriptListOrdinal(int, double **, short, short, short, double,
		double, short, double *, double *, double, double, short);
void	PostScriptListPoints(double **, int);
void	PostScriptListPolygon(double **, int);
void	PostScriptSurfEvaluate(ParSurf *, float, float, float);
void	PostScriptSurfFocal(ParSurf *, short, short, short, short, float,
		float ***, float ***, double **, short);
void	PostScriptSurfKnots(ParSurf *, int, int, int, int);
void	PostScriptSurfOffset(ParSurf *, short, short, short, short, float,
		float ***, float ***);
void	PostScriptSurfOrthotomic(ParSurf *, short, short, short, short,
		float, float *, float ***, float ***);
void	PostScriptSurfPolygon(ParSurf *);
void	PostScriptSurfWireframe(ParSurf *, short, short, short, short,
		float ***);
void	PostScriptTrimShaded(int, int *, int *, float **);
void	PostScriptTrimWireframe(TrimSurf *, int, int, int);
void	PostScriptUv(ParUv *, ParSurf *);
void	PostScriptUvPolygon(ParUv *, ParSurf *);
void	PostScriptVect(short, double **, char *);
void	PostScriptView(char *);
void	PostScriptWorldAxis(void);
char	*PromptBox(char *, char *, char *, int *);
short	RadioBox(char *, char *, int, char **, int, int *);
short	ReadDeslabCos(char *, char *, short, char *);
short	ReadDeslabCurv(char *, char *, short);
short	ReadDeslabFacet(char *, char *, char *);
short	ReadDeslabGrid(char *, char *);
short	ReadDeslabHull(char *, char *);
short   ReadDeslabInsp(char *, char *);
short	ReadDeslabList(char *, char *);
short	ReadDeslabLocal(char *, char *, int);
short	ReadDeslabMinDist(char *, char *, int);
void	ReadDeslabMult(char *);
short	ReadDeslabSurf(char *, char *);
short	ReadDeslabTrim(char *, char *);
short	ReadDeslabUv(char *, char *, char *);
short	ReadDeslabVect(char *, char *, short, short);
int     ReadIgesArc(char *, FILE *, char *, int *, char, char, char *,
		    Type124 *);
int     ReadIgesComposite(char *, FILE *, char *, int *, char, char, char *,
			  Type124 *);
int     ReadIgesConic(char *, FILE *, char *, int *, char, char, char *,
		      Type124 *);
int     ReadIgesConnectPoint(char *, FILE *, char *, int *, char, char, char *,
		      Type124 *);
int	ReadIgesCos(char *, FILE *, char *, int *, char, char, char *,
		    Type124 *);
int	ReadIgesCurv(char *, FILE *, char *, int *, char, char, char *,
		     Type124 *);
int	ReadIgesCos(char *, FILE *, char *, int *, char, char, char *,
		    Type124 *);
int	ReadIgesLine(char *, FILE *, char *, int *, char, char, char *,
		     Type124 *);
int	ReadIgesList(char *, FILE *, char *, int *, char, char, char *,
		     Type124 *);
int     ReadIgesPoint(char *, FILE *, char *, int *, char, char, char *,
		      Type124 *);
int	ReadIgesPowerCurv(char *, FILE *, char *, int *, char, char, char *,
			  Type124 *);
int	ReadIgesPowerSurf(char *, FILE *, char *, int *, char, char, char *,
			  Type124 *);
int	ReadIgesRuledSurf(char *, FILE *, char *, int *, char, char, char *,
			  Type124 *);
int	ReadIgesSurf(char *, FILE *, char *, int *, char, char, char *,
		     Type124 *);
int	ReadIgesSurfRev(char *, FILE *, char *, int *, char, char, char *,
			Type124 *);
int	ReadIgesTabCyl(char *, FILE *, char *, int *, char, char, char *,
		       Type124 *);
int	ReadIgesTrim(char *, FILE *, char *, int *, char, char, char *,
		     Type124 *);
int	ReadIgesUv(char *, FILE *, char *, int *, char, char, char *,
		   Type124 *);
void	ReadInitFile(char *, short);
int	ReadPraxDirectory(char *);
void	ReportSurfContinuity(int, int);
void	SaveEntity(int, int, char *);
int     SaveIgesMult(char *, int);
char   *ScrollBox(char *, char *, int);
void    SelectBox(char *, char *, int, char **, int *, int *);
void	SelectGroupEntities(int);
void	SetBackground(long);
void    SetBatchMode(int);
void	SetContourScale(short);
void	SetDisplayType(short);
void	SetFilter(short, char *);
void	SetForeground(long);
void	SetFormat(short);
void	SetCurrentTool(short);
void	SetCwd(char *);
void	SetDefaultIgesEntity(int, int);
void	SetLocalAxis(short);
void    SetLogFile(char *);
void	SetLogoIters(short);
void	SetMaxCoss(short);
void	SetMaxContour(double);
void	SetMaxCurvs(short);
void	SetMaxDistance(double);
void	SetMaxFacets(short);
void	SetMaxGrids(short);
void	SetMaxGroups(short);
void	SetMaxHulls(short);
void	SetMaxHeight(double);
void    SetMaxInsps(short);
void	SetMaxLists(short);
void	SetMaxOrdinal(double);
void	SetMaxSlope(double);
void	SetMaxSurfs(short);
void	SetMaxTrims(short);
void	SetMaxUvs(short);
void	SetMaxVects(short);
void	SetMinContour(double);
void	SetMinDistance(double);
void	SetMinHeight(double);
void	SetMinOrdinal(double);
void	SetMinSlope(double);
void	SetNumCoss(short);
void	SetNumCurvs(short);
void	SetNumFacets(short);
void	SetNumGrids(short);
void	SetNumGroups(short);
void	SetNumHulls(short);
void    SetNumInsps(short);
void	SetNumLists(short);
void	SetNumSurfs(short);
void	SetNumTrims(short);
void	SetNumUvs(short);
void	SetNumVects(short);
void	SetPerspective(void);
void    SetPostScriptImage(short);
void    SetPostScriptName(char *);
void	SetPostScriptTransform(short, float, float, float *, float,
		float, float, short);
void	SetPraxDirectory(char *);
void	SetRadiusTolerance(double);
void	SetRevContour(short);
void    SetScriptFile(char *);
void    SetScriptFp(FILE *);
void	SetSelectedCoss(short);
void	SetSelectedCurvs(short);
void	SetSelectedFacets(short);
void	SetSelectedGrids(short);
void	SetSelectedGroups(short);
void	SetSelectedHulls(short);
void    SetSelectedInsps(short);
void	SetSelectedLists(short);
void	SetSelectedSurfs(short);
void	SetSelectedTrims(short);
void	SetSelectedUvs(short);
void	SetSelectedVects(short);
void	SetShowLogo(short);
void	SetShowMotd(short);
void    SetSystemOptions(int);
void	SetViewAxis(short);
void	SetViewDist(float);
void	SetViewFar(float);
void	SetViewFov(short);
void	SetViewMaxRange(float);
void	SetViewNear(float);
void	SetViewPadInc(short, float);
void	SetViewPadText(short, float);
void	SetViewRotX(float);
void	SetViewRotY(float);
void	SetViewingAxis(void);
void	SetViewingTransform(short, float, float, float *, float,
		float, float, short);
void	SetViewTwist(short);
void	SetViewX(float);
void	SetViewY(float);
void	SetViewZ(float);
void	SetWorldAxis(short);
void    SplitScript(char *);
void	SplitSurface(ParSurf *, int, int, double *, double *);
void	StoreCosValues(ParCurv *, int, int, ParSurf *, float ***,
		float ***, double **, double **, short *);
void	StoreCurvValues(ParCurv *, int, int, int, float ***, float ***,
		double **, short *);
void	StoreSurfValues(ParSurf *, short, short, float ****, float ****,
		double ***, double ***, double ***, double ***, short *,
		short *, int);
void	SubDistance(int, int, int);
void    ToolBar(int, int);
void	TriangulateStructured(int, int, int, char **);
void	TriangulateStructuredToEdge(int, int, int, char **);
void	TriangulateUnstructured(int, int, int, char **);
void	TriangulateUnstructuredToEdge(int, int, int, char **);
void	TrimUnconloc(int, int, int);
void	Unconloc(int, int, int);
void	UnconlocBlade(int, int, int, int);
void	UnconlocResults(int, int, int, int, int, double, double, double *,
			int, char *);
void	Unconloc2dResults(int, int, int, int, double, double,
		double *, int);
void	UnselectAllCoss(void);
void	UnselectAllCurvs(void);
void	UnselectAllFacets(void);
void	UnselectAllGrids(void);
void	UnselectAllGroups(void);
void	UnselectAllHulls(void);
void    UnselectAllInsps(void);
void	UnselectAllLists(void);
void	UnselectAllSurfs(void);
void	UnselectAllTrims(void);
void	UnselectAllUvs(void);
void	UnselectAllVects(void);
void	UpdateEntities(void);
void	UpdateScale(double, double, short, short, double, double, short);
void	UpdateView(void);
void	UpdateViewingAxis(void);
void    UvScript(char *);
vector *VectorBox(char *, char *, int *);
void	ViewEntities(void);
void    ViewPadText(int, float);
void	WinsetGlxWindow(void);
void	WinsetViewingAxis(void);
int	WriteIgesStart(FILE *, char *);
short	YesNoBox(char *, char *, int, int *);

#ifdef __cplusplus
}
#endif

#endif

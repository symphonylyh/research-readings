/* Copyright (C) 1997 Massachusetts Institute of Technology, Cambridge, MA
   All rights reserved
*/

/* iges.h */

#ifndef IGES_H
#define IGES_H

/* #define POWER_DETAIL; */
#define POWER_DETAIL

#include "gen.h"

#ifdef __cplusplus
extern "C" {
#endif

#define PRAX_IGES_BINARY	'B'
#define PRAX_IGES_COMPRESSED	'C'
#define PRAX_IGES_START		'S'
#define PRAX_IGES_GLOBAL	'G'
#define PRAX_IGES_DIRECTORY	'D'
#define PRAX_IGES_PARAMETER	'P'
#define PRAX_IGES_TERMINATE	'T'

#define PRAX_ERR_IGES_BINARY		-1
#define PRAX_ERR_IGES_COMPRESSED	-2
#define PRAX_ERR_NOT_IGES		-3
#define PRAX_ERR_EOF			-4

#define PRAX_IGES_INTEGER	0
#define PRAX_IGES_REAL		1
#define PRAX_IGES_STRING	2

#define PRAX_IGES_ZERO	1.0e-10

/********** 11.1 begin **********/
#define CANCEL_SAVE -1000
/********** 11.1 end ************/

typedef struct {
  char   parmDelim;         /* parameter delimeter */
  char   recDelim;          /* record delimeter */
  char   *sendID;           /* name of the model */
  char   *fileName;         /* file name */
  char   *systemID;         /* system id */
  char   *preprocVers;      /* version id of system */
  int    bits;              /* number of bits per integer */
  int    singleMagnitude;   /* max magnitude for single precision */
  int    singlePrecision;   /* significant digits for single precision */
  int    doubleMagnitude;   /* max magnitude for double precision */
  int    doublePrecision;   /* significant digits for double precision */
  char   *receiveID;        /* id of receiving system */
  double scale;             /* ratio of model space units to cartesian space */
  int    unit;              /* flag specifying system unit */
                            /* flag   Measuring_system
                             *  ----   -------------------------------------
                             *   1    inches
                             *   2    millimeters
                             *   3    units named by units character string
                             *   4    feet
                             *   5    miles
                             *   6    meters
                             *   7    kilometers
                             *   8    mils
                             *   9    microns
                             *  10    centimeters
                             *  11    microinches */
  char   *units;            /* name of the units */
                            /* units  Units
                             * -----  ----------
                             * 'IN'   inches
                             * 'INCH' inches
                             * 'MM'   millimeters
                             * 'FT'   feet
                             * 'MI'   miles
                             * 'M'    meters
                             * 'KM'   kilometers
                             * 'MIL'  mils
                             * 'UM'   microns
                             * 'CM'   centimeters
                             * 'UIN'  microinches */
  int    gradations;        /* number of distinct line thicknesses */
  double weight;            /* actual width of thickest possible line */
  char   fileDate[14];      /* date of file creation */
  double resolution;        /* smallest discenible distance in model space */
  double coordinate;        /* absolute value of maximum coordinate value */
  char   *author;           /* author's name */
  char   *org;              /* author's organization */
  int    version;           /* IGES version */
                            /* version  Version
                             * -------  -------------------------
                             *    1     1.0
                             *    2     ANSI Y14.26M - 1981
                             *    3     2.0
                             *    4     3.0
                             *    5     ASME/ANSI Y14.26M - 1987
                             *    6     4.0
                             *    7     ASME Y14.26M - 1989
                             *    8     5.0 */
  int    standard;          /* drafting standard */
             /* standard Standard
              *     0     None
              *     1     ISO   International Organization for Standardization
              *     2     AFNOR French Association for Standardization
              *     3     ANSI  American National Standards Institute
              *     4     BSI   British Standards Institute
              *     5     CSA   Canadian Standards Association
              *     6     DIN   German Institute for Standardization
              *     7     JIS   Japanese Institute for Standardization */
  char   modelDate[14];     /* model creation date */
} GlobalStruct;

typedef struct {
  int   type;             /* entity type number */
  int   parameter;        /* sequence number of first line of parameter data */
  int   structure;        /* unsupported, defaults to 0 */
  int   pattern;          /* unsupported, defaults to 0 */
  int   level;            /* unsupported, defaults to 0 */
  int   view;             /* unsupported, defaults to 0 */
  int   matrix;           /* unsupported, defaults to 0 */
  int   associativity;    /* unsupported, defaults to 0 */
  int blank;              /* unsupported, defaults to 0 */
  int subordinate;        /* unsupported, defaults to 0 */
  int use;                /* unsupported, defaults to 0 */
  int hierarchy;          /* unsupported, defaults to 0 */
  int   weight;           /* unsupported, defaults to 0 */
  int   color;            /* unsupported, defaults to 0 */
  int   count;            /* number of lines in parameter section */
  int   form;             /* format number for entity-specific options */
  char  label[9];         /* eight character descriptive label */
  int   subscript;        /* unsupported, defaults to 0 */
  int DE;                 /* entity sequence number */
} DirStruct;

typedef struct entry_struct {
  DirStruct *dir;            /* directory entry structure */
  void *igeom;               /* the geometric entity of this entry */
  char *tmpfile;             /* reserved */
  struct entry_struct *next; /* next entry in a linked list */
                             /* if NULL, then this is the last in the list */
} EntryStruct;

typedef struct {
  int	start;            /* line count of Start section */
  int	global;           /* line count of Global section */
  int	directory;        /* line count of Directory section */
  int parameter;          /* line count of Parameter section */
} TermStruct;

typedef struct {
  double zt;              /* common z value */
  double *x;              /* array of x values */
  double *y;              /* array of y values */
} CopiousPair;

typedef struct {
  double *x, *y, *z;      /* arrays of x,y,z values */
} CopiousTriple;

typedef struct {
  double *x, *y, *z;      /* arrays of x,y,z values */
  double *i, *j, *k;      /* arrays of i,j,k values */
} CopiousSextuple;

typedef union {
  CopiousPair pair;
  CopiousTriple triple;
  CopiousSextuple sextuple;
} CopiousUnion;

typedef struct {     /* circular arc */
  double zt;         /* displacement from x,y plane */
  double x1, y1;     /* center point */
  double x2, y2;     /* start point */
  double x3, y3;     /* end point */
} Type100;

typedef struct {     /* composite curve */
  int n;             /* number of entities */
  int *de;           /* array of DEs to entities */
} Type102;

typedef struct {     /* conic arc */
  double a,b,c,d,e,f;/* conic coefficients */
  double zt;         /* z coordinate of plane of definition */
  double x1, y1;     /* start point */
  double x2, y2;     /* end point */
} Type104;

typedef struct {     /* copious data */
  int ip;            /* interpretation flag */
                     /* ip   Interpretation
                      * --   ---------------------------------
                      *  1   Coordinate pairs (x,y)
                      *  2   Coordinate triples (x,y,z)
                      *  3   Coordinate sextuples (x,y,z,i,j,k) */
  int n;             /* number of coordinate points */
  CopiousUnion pts;  /* the data points (either pairs, triples, or sextuples */
} Type106;

typedef struct {     /* plane */
  double a, b, c, d; /* coefficients of plane, Ax + By + Cz = D */
  int ptr;           /* pointer to DE of bounding curve */
  double x, y, z;    /* location of display symbol */
  double size;       /* size of display symbol */
} Type108;

typedef struct {     /* line */
  double x1;	     /* start point */
  double y1;
  double z1;
  double x2;	     /* terminate point */
  double y2;
  double z2;
} Type110;

typedef struct {     /* point */
  double x, y, z;    /* coordinates */
  int ptr;           /* pointer to DE of of symbol */
} Type116;

typedef struct {     /* ruled surface */
  int de1;	     /* pointer to first boundary curve */
  int de2;	 /* pointer to second boundary curve */
  int dirflg;	 /* direction flag */
  int devflg;
} Type118;

typedef struct { /* surface of revolution */
  int de1;	 /* pointer to axis of revolution */
  int de2;	 /* pointer to generatrix */
  double start; 	/* start angle (radians) */
  double term;	 /* terminate angle (radians) */
} Type120;

typedef struct { /* tabulated cylinder */
  int de;	 /* pointer to directrix curve */
  double lx;	 /* terminate point coordinates */
  double ly;
  double lz;
} Type122;

typedef struct {  /* transformation matrix */
  double r[3][3]; /* rotation matrix */
  double t[3];    /* translation vector */
} Type124;

typedef struct {  /* connect point */
  double x, y, z; /* coordinates of connection point */
  int ptr;        /* pointer to DE of display symbol */
  int tf;         /* type flag */
  int ff;         /* function flag */
  char *cid;      /* connection point function identifier */
  int pttcid;     /* pointer to DE of text display template entity for CID */
  char *cfn;      /* connection point function name */
  int pttcfn;     /* pointer to DE of text display template entity for CFN */
  int cpid;       /* unique connect point identifier */
  int fc;         /* connect point function code */
  int sf;         /* swap flag */
  int psfi;       /* pointer to DE of network subfigure instance entity */
} Type132;

typedef struct { /* curve on parametric surface */
  int crtn;	 /* creation flag */
  int sptr;	 /* surface pointer */
  int bptr;	 /* surface curve */
  int cptr;	 /* space curve */
  int pref;	 /* representational preference flag*/
                 /* pref  Preference
                  * ---   --------------------------
                  *   0   Unspecified
                  *   1   Parametric representation
                  *   2   Space curve representation
                  *   3   Either */
} Type142;

typedef struct { /* trimmed surface */
  int pts;	 /* surface pointer */
  int n1; 	 /* boundary flag */
                 /* if n1 = 0, then the outer boundary of the untrimmed
		  *            surface is a trimming loop of the trimmed
		  *            surface, else
		  * if n1 = 1, then not */
  int n2;	 /* number of trimming loops, not including the outer  */
                 /* boundary if n1 = 0 */
  int pt0;	 /* outer curve pointer */
  int *pti;	 /* inner curve pointers */
} Type144;

/* function prototypes */

void	AddNextInteger(FILE *, int, char, char *, int, char,
		int *, int);
void	AddNextReal(FILE *, double, char, char *, int, char,
		 int *, int);
void	AddNextString(FILE *, char *, char, char *, int, char,
		int *, int);
vector  *addv(vector *, vector *); 
Type106 *alloc_copious(int, int);
DirStruct *AllocIgesDirectory(int, int, int, int, int, int, int,
		int, int, int, int, int, int, int, int,
		int, char *, int, int DE);
GlobalStruct *AllocIgesGlobal(char, char, char *, char *, char *,
		char *, int, int, int, int, int, char *,
		double, int, char *, int, double, char *, double,
		double, char *, char *, int, int, char *);
Type124	*alloc_tmatrix(void);
double  berntomono1(int, int, double *);
double  berntomono2(int, int, int, int, double **);
void	calc_points_curv(double *, PowCurv *, double **, int);
void	calc_points_surf(double *, double *, PowSurf *, double **,
		double *, double *, int);
int	CheckCurv2d(ParCurv *);
int	CheckCurvClosed(ParCurv *);
int	CheckCurvIntegral(ParCurv *);
int	CheckCurvPlanar(ParCurv *, double *, double *, double *);
int	CheckIgesFile(FILE *, char *);
int	CheckSurfClosedU(ParSurf *);
int	CheckSurfClosedV(ParSurf *);
int	CheckSurfIntegral(ParSurf *);
double  comb(int, int);
void	DecodeIgesDirectory(char *, int *, DirStruct **);
void    DecodeIgesGlobal(FILE *, char *, GlobalStruct **);
TermStruct *DecodeIgesTerminate(char *);
void    FileCat(FILE *, FILE *);
void	free_copious(Type106 *);
void	free_pgeom(PowCurv *);
void	free_sgeom(PowSurf *);
void	free_tmatrix(Type124 *);
void	FreeIgesDirList(EntryStruct *);
void    FreeIgesDirectory(DirStruct *);
void	FreeIgesGlobal(GlobalStruct *);
EntryStruct *GetIgesDirEntry(EntryStruct *, int);
EntryStruct *GetIgesDirList(FILE *, char *);
char	*GetIgesEntity(int);
int	GetNextInteger(FILE *, char *, int *, char, char, int, int *);
int	GetNextReal(FILE *, char *, int *, char, char, int, double *);
int	GetNextString(FILE *, char *, int *, char, char, int, char **);
void	make_knots(double *, double *, int, int);
double  monotobern1(int, int, double *);
double  monotobern2(int, int, int, int, double **);
vector  *multv(double, vector *); 
char    *ParseIgesRecord(FILE *, char *, int *, char, char, int, int);
PowCurv *pgeomalloc(int);
PowCurv *ParCurv_to_PowCurv(ParCurv *, PowCurv *, double *);
PowCurv *ParCurv_to_PowCurv1(ParCurv *, PowCurv *, double *);
PowSurf *ParSurf_to_PowSurf(ParSurf *, PowSurf *, double *);
PowSurf *ParSurf_to_PowSurf2(ParSurf *, PowSurf *, double *);
void    PowCurv_error(PowCurv *, double *);
vector	*PowCurv_eval(PowCurv *, double, int);
ParCurv *PowCurv_to_ParCurv (PowCurv *, ParCurv *, double *, int);
ParCurv *PowCurv_to_ParCurv1 (PowCurv *, ParCurv *, double *, int);
void	PowParCurv_compare(ParCurv *, PowCurv *, double *, int);
void	PowParCurv_compare1(ParCurv *, PowCurv *, double *, int);
void	PowParSurf_compare(ParSurf *, PowSurf *, double *, int);
void	PowParSurf_compare2(ParSurf *, PowSurf *, double *, int);
void	PowSurf_error(PowSurf *, double *);
vector	*PowSurf_eval(PowSurf *, double u, double, int, int);
PowCurv *PowSurf_iso(PowSurf *, int, double, PowCurv *);
ParSurf *PowSurf_to_ParSurf(PowSurf *, ParSurf *, double *, int, int);
ParSurf *PowSurf_to_ParSurf2(PowSurf *, ParSurf *, double *, int);
void	PowSurf_to_ParSurf_int(PowSurf *, ParSurf *, double *, int);
void	PowSurf_to_ParSurf_loft(PowSurf *, ParSurf *,double *, int);
void	PrintIgesDirectory(DirStruct *, int);
void	PrintIgesGlobal(GlobalStruct *);
void    PrintIgesParameter(DirStruct *, FILE *, char *, int *, char, char);
void	PrintIgesTerminate(TermStruct *);
Type100 *ReadIges100(FILE *, char *, int *, char, char);
Type102 *ReadIges102(FILE *, char *, int *, char, char);
ParCurv *ReadIges102Curv(FILE *, char *, int *, char, char, Type102 *);
Type104 *ReadIges104(FILE *, char *, int *, char, char);
Type106 *ReadIges106(FILE *, char *, int *, char, char);
Type108 *ReadIges108(FILE *, char *, int *, char, char);
Type110 *ReadIges110(FILE *, char *, int *, char, char);
PowCurv *ReadIges112(FILE *, char *, int *, char, char);
PowSurf *ReadIges114(FILE *, char *, int *, char, char);
Type116 *ReadIges116(FILE *, char *, int *, char, char);
Type118 *ReadIges118(FILE *, char *, int *, char, char);
Type120 *ReadIges120(FILE *, char *, int *, char, char);
Type122 *ReadIges122(FILE *, char *, int *, char, char);
Type124 *ReadIges124(FILE *, char *, int *, char, char);
ParCurv *ReadIges126(FILE *, char *, int *, char, char);
ParSurf *ReadIges128(FILE *, char *, int *, char, char);
Type132 *ReadIges132(FILE *, char *, int *, char, char);
Type142 *ReadIges142(FILE *, char *, int *, char, char);
Type144 *ReadIges144(FILE *, char *, int *, char, char);
char	ReadNextIgesRecord(FILE *, char *, int);
int     RemoveCurveKnot(int, int, double *, vector **, double, int, int, int);
int     RemoveSurfKnotU(ParSurf *);
int     RemoveSurfKnotV(ParSurf *);
PowSurf *sgeomalloc(int, int);
void    SkipToNext(FILE *, char *, int *, int, char);
int	SkipToThis(FILE *, char *, int *, int, char);
vector  *subv(vector *, vector *); 
void	TransformIgesArray(Type124 *, int, int, vector ***);
void	TransformIgesList(Type124 *, int, double **);
void	TransformIgesUV(Type124 *, int, double **);
void	TransformIgesVect(Type124 *, int, vector **);
int	WriteIges106(FILE *, Type106 *, char, char, int, int *);
int	WriteIges110(FILE *, Type110 *, char, char, int, int *);
int	WriteIges112(FILE *, PowCurv *, char, char, int, int *);
int	WriteIges114(FILE *, PowSurf *, char, char, int, int *);
int	WriteIges118(FILE *, Type118 *, char, char, int, int *);
int	WriteIges120(FILE *, Type120 *, char, char, int, int *);
int	WriteIges122(FILE *, Type122 *, char, char, int, int *);
int	WriteIges126(FILE *, ParCurv *, char, char, int, int *);
int	WriteIges128(FILE *, ParSurf *, char, char, int, int *);
int	WriteIges142(FILE *, Type142 *, char, char, int, int *);
int	WriteIges144(FILE *, Type144 *, char, char, int, int *);
void	WriteIgesDirectory(FILE *, DirStruct *, int);
int	WriteIgesGlobal(FILE *, GlobalStruct *);
void	WriteIgesRecord(FILE *, char *, int, char, int);
void	WriteIgesTerminate(FILE *, int, int, int, int);

#ifdef __cplusplus
extern "C" {
#endif

#endif

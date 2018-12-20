/************************************************************************
 *									*
			 Copyright (C) 1989, 1990
     Massachusetts Institute of Technology, Cambridge, Massachusetts.
			   All rights reserved.

     This portion of source code was prepared in the Ocean Engineering
     Design Laboratory by Bradley A. Moran.
 *									*
 ************************************************************************/

#ifndef GEN_H			/* this is the header file gen.h */
#define GEN_H

#include <stdio.h>		/* For the FILE data type */

#ifdef __cplusplus
extern "C" {
#endif

typedef double double2;		/* for portability to old IRIS 3030 */

#if 0
typedef float real;		/* This is going unless somebody bitches */
#endif

#ifndef FALSE			/* generic conditionals */
#define FALSE 0
/***** this is inconsistent with definition in "compiler.h" *****
#define TRUE (!FALSE)
****** this is consistent with "compiler.h" (8/9/95) ***********/
#define TRUE 1
#endif

#define ENOGEOM (-1)		/* error code from dispatch() */

#define ZERO 1.0e-06		/* standard tolerance levels */

#if mips
#define MACHPREC 1.1102230246251567E-13	  /* changed on 25 JUN 90 */
#endif
#if vax
#define MACHPREC 1.3877787807814E-14
#endif

#define ACurve 192613		/* Implicit Polynomial Curve */
#define ASurface 192629		/* Implicit Polynomial Surface */
#define PCurveOpen 126		/* IGES NURBS curve */
#define PCurvePer 261326	/* Parametric Curve Periodic */
#define PSurfaceOpen 128	/* IGES NURBS surface */
#define PSurfacePerU 262931	/* Parametric Surface periodic U */
#define PSurfacePerV 262932	/* Parametric Surface periodic V */
#define PSurfacePerUV 263132	/* Parametric surface periodic UV */
#define PCurvePow 112		/* IGES entity number */
#define PSurfacePow 114		/* IGES entity number */

/* useful macros */
#define MAX(x,y)  (((x) >= (y)) ? (x): (y))  /* maximum of x and y */
#define MIN(x,y)  (((x) <= (y)) ? (x): (y))  /* minimum of x and y */


#define FABS(x) fabs ((double2) (x)) /* for portability */
#define SQRT(x) sqrt ((double2) (x))
#define UP_I(x,y) pow ((double2) (x),(double2) (y))
#define LPOW(x,y) pow ((double2) (x),(double2) (y))

typedef struct vector {		/* homogeneous coordinates vector */
     double2 x, y, z, w;
} vector;

typedef struct egeom {		/* parametric B-Spline curve segment */
     int type;
     int   order,         /* order of curve */
           ncontpts;      /* number of control points */
     int   kmem,          /* size of memory allocated for knot vector */
           pmem;          /* size of memory allocated for control points */
     double2 *knots;      /* knot vector */
     vector **contpts;    /* control points */
} ParCurv;
	
typedef struct fgeom {		/* parametric B-Spline surface patch */
     int type;
     int   uorder,        /* u order */
           vorder;        /* v order */
     int   ucontpts,      /* number of u control points */
           vcontpts;      /* number of v control points */
     int   ukmem,         /* size of memory allocated for u knot vector */
           vkmem,         /* size of memory allocated for v knot vector */
           upmem,         /* size of memory allocated for u control points */
           vpmem;         /* size of memory allocated for v control points */
     double2 *uknots,     /* u knot vector */
             *vknots;     /* v knot vector */
     vector ***contpts;   /* control points */
} ParSurf;

struct pgeom {			/* parametric power basis curve segment */
     int type;
     int   order,         /* order of curve */
           nsegmts;       /* number of polynomial segments */
     int   kmem,          /* size of memory allocated for break points */
           pmem;          /* size of memory allocated for control points */
     double *knots;       /* break points */
     vector ***contpts;   /* control points */
};

struct sgeom {			/* parametric power basis surface segment */
     int type;
     int   uorder,        /* u order */
           vorder,        /* v order */
           usegmts,       /* number of u polynomial segments */
           vsegmts;       /* number of v polynomial segments */
     int   ukmem,         /* size of memory allocated for u break points */
           vkmem,         /* size of memory allocated for v break points */
           upmem,         /* size of memory allocated for u control points */
           vpmem;         /* size of memory allocated for u control points */
     double *uknots,      /* u break points */
            *vknots;      /* v break points */
     vector ***contpts[16];	/* hardwired to be bicubic patch */
};

				/* define POWER_DETAIL before include */
				/* otherwise the details will be suppressed */
#ifdef POWER_DETAIL
typedef struct pgeom PowCurv;
typedef struct sgeom PowSurf;
#else
typedef struct { char _x[sizeof(struct pgeom)];} PowCurv;
typedef struct { char _x[sizeof(struct sgeom)];} PowSurf;
#endif				/* POWER_DETAIL */

#define BITS_PER_CHAR (unsigned char)8

/* Function prototypes follow */

/* alloc.c */

vector *vectalloc(void);
void vectfree(vector *);

/* avec02b.c */

vector *avec02b (void);

/* bit_array1.c */

unsigned char *bit_array1(unsigned int);
void free_barray1(unsigned char *);
void set_bit(unsigned char *, unsigned int, int);
int test_bit(unsigned char *, unsigned int);

/* copy.c */

ParCurv *copyegeom (ParCurv *, ParCurv *);
ParSurf *copyfgeom (ParSurf *, ParSurf *);
vector *copyvector (vector *, vector *);

/* dbl_array.c */

double2 *dbl_array1 (unsigned);
double2 **dbl_array2 (unsigned, unsigned);
double2 ***dbl_array3 (unsigned, unsigned, unsigned);
void free_darray1 (double2 *);
void free_darray2 (double2 **);
void free_darray3 (double ***);

/* dispatch.c */

int dispatch (char *);

/* distance.c */

double linearDistance(ParCurv *, vector *, double *);
double pointToLine(vector *, vector*, vector *, double *);

/* dynamic_memory.c */

#define GEN_B_PER_BLOCK 512.0
#define GEN_KB_PER_B (1.0/1024.0)
#define GEN_MB_PER_B (1.0/1048576.0)

double	DynamicMemoryUsage(double *, double *, int);
double	LogicalSwapSpace(double *, double *, double *);
char   *MemoryStatus(void);

/* egeom.c */

ParCurv *egeomalloc1 (int, int);
ParCurv *egeomalloc2 (int, int);
void free_egeom (ParCurv *);

/* error.c */

int errormsg (int, char *);

/* fgeom.c */

ParSurf *fgeomalloc1 (int, int, int, int);
ParSurf *fgeomalloc2 (int, int, int, int);
void free_fgeom (ParSurf *);

/* fileopen.c */

FILE *fileopen (char *, int, char *);

/* flt_array.c */

float *flt_array1 (unsigned);
float **flt_array2 (unsigned, unsigned);
void free_farray1(float *);
void free_farray2(float **);

/* gen_array.c */

char *gen_array1 (unsigned, unsigned);
char **gen_array2 (unsigned, unsigned, unsigned);
void free_garray1(char *);
void free_garray2(char **);

/* int_array.c */

int *int_array1 (unsigned);
int **int_array2 (unsigned, unsigned);
void free_iarray1 (int *);
void free_iarray2(int **);

/* matrix.c */

void mult4x4 (double **, vector *, vector *);
void matrix_mult(float[],float[],float[],
		 int,int,int);
void matrix_transpose(int,int,double2[],int,
		      double2[],int);
void Print_matrix(int,int,double2[],int);
void Read_matrix(FILE *,double2 *,int,int);
void matrixmu(int,int,int, int,int,int, 
	      double2[], double2[], double2[]);

/* nag_funs.c */

void set_fun1_ptr (void (*userfun1)(int *, double *, double *));
void set_fun2_ptr (void (*userfun2)(int *, double *, double *, double *));
void set_hes2_ptr (void (*userhes2)(int *, double *, double *, int *,
				    double *));
void funct1_ (int *, double *, double *);
void funct2_ (int *, double *, double *, double *);
void hess2_ (int *, double *, double *, int *, double *);

void save_fun1_ptr(void);
void save_fun2_ptr(void);
void restore_fun1_ptr(void);
void restore_fun2_ptr(void);

/* ptr_array.c */

char **ptr_array1 (unsigned);
char ***ptr_array2 (unsigned, unsigned);
void free_parray1 (char **);
void free_parray2 (char ***);

/* read.c */

ParCurv *ReadParCurv (FILE *, ParCurv *);
ParSurf *ReadParSurf (FILE *, ParSurf *);
int GetNextToken(FILE *, char *);

/* readper.c */

ParCurv *ReadParCurv_Per (FILE *, ParCurv *);
ParSurf *ReadParSurf_Peru (FILE *, ParSurf *);

/* sht_array.c */

short *sht_array1 (unsigned);
short **sht_array2 (unsigned, unsigned);
void free_sarray1 (short *);
void free_sarray2 (short **);

/* string.c */

void stripwhite (FILE *);
char *fgetstring (char *, int, FILE *);

/* vec_array.c */

vector **vec_array1 (unsigned);
vector ***vec_array2 (unsigned, unsigned);
void free_varray1 (vector **, unsigned);
void free_varray2 (vector ***, unsigned, unsigned);

/* vect_trans.c */

void rotx1(double2, vector *, vector *);
void roty1(double2, vector *, vector *);
void rotz1(double2, vector *, vector *);
void translate1(double2, double2, double2, vector *, vector *);

/* vectarith.c */

vector *add_vect (vector *, vector *);
void add_vect1 (vector *, vector *, vector *);
vector *cross (vector *, vector *);
void cross1 (vector *, vector *, vector *);
double2	distance (vector *, vector *);
double	distance1 (vector *, double, double, double);
double2 dot (vector *, vector *);
vector *glue_vector (double2, double2, double2);
void glue_vector1 (double2, double2, double2, vector *);
double2	mag (vector *);
void scale4 (double2, vector *);
vector *scale_vect (double2, vector *);
void scale_vect1 (double2, vector *, vector *);
void sethomogeq (vector *, double2);
vector *sub_vect (vector *, vector *);
void sub_vect1 (vector *, vector *, vector *);
vector *unitvector (vector *);
void unitvector1 (vector *, vector *);

/* vectarith2.c */

vector *midpoint(vector *, vector *, vector *);

/* write.c */

void WriteParSurf(FILE *, ParSurf *);
void WriteParCurv(FILE *, ParCurv *);
void Writevector(FILE *, vector *);

/* writeper.c */

void WriteParCurv_Per(FILE *, ParCurv *);
void WriteParSurf_Peru (FILE *, ParSurf *);

#ifdef __cplusplus
}
#endif

#endif

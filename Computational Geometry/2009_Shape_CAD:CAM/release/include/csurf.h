/* Copyright (C) 1994 Massachusetts Institute of Technology, Cambridge, MA
   All rights reserved
*/

/* csurf.h */

#ifndef CSURF_H
#define CSURF_H

#ifdef __cplusplus
extern "C" {
#endif

#define MACHPREC 1.1102230246251567E-13
#define MACHPR2 MACHPREC*MACHPREC
#define M_PI	3.14159265358979323846
#define PI05 (M_PI/2.0)
#define PI2  (M_PI*2.0)
#define NDFPTS  50    /* number of intermediate points on the camber line */
#define NCMPTS  300   /* maximum number of arc-length parameter values    */
#define MAXCLL  10000 /* max number of calls to routine providing the rhs */
#define CMBFC   1.0e5 /* factor multiplied with the tight tolerance       */
#define MAXSCT  301   /* maximum number of control polygon vertices       */
#define MAXCMB  150   /* Maximum number of control polygon vertices       */
#define CHRFC 0.505050505 /* dist of the midchord from the leading edge   */

/* prototypes not in standard headers under -ansi option */

char *strdup(const char *);
char *tempnam(const char *, const char *);

#include "editor.h"

/* prototypes for NAG functions */

void c05azf_(double *, double *,double *, double *, int *, double *w, int *,
	     int *);
void c05nbf_(void (int *, double *, double *, int *), int *, double *,
	     double *, double *, double *, int *, int *);
void d02chf_(double *, double *, int *, double *, double *, int *, double *,
	     void (double *, double *, double *), double (double *, double *),
	     double *, int *);
void d02cjf_(double *, double *, int *, double *, void (double *, double *,
							double *), double *,
	     char *, void (double *, double *), double (double *, double *),
	     double *, int *);
void e04jaf_(int *, int *, double *, double *, double *, double *, int *,
	     int *, double *, int *, int *);
void e04kaf_(int *, int *, double *, double *, double *, double *, double *,
	     int *, int *, double *, int *, int *);

/* prototypes for csurf */

void	bone_intersect3(ParCurv *, double, double *, double *, vector *,
		vector *);
int	camber_line_3d(ParCurv *, double *, double, double *, int, int,
		ParCurv **, double *);
vector	*camber_surf_3d(vector *);
void	camber_surface(ParSurf *, int, double *, double *, double *,
		int *, int *, ParSurf **, ParCurv **, double *,
		double **, short *, short *, int, int, int);
vector	*Camber2d_3d(vector *);
int	CamberBisect3d(ParCurv *, ParCurv **, double *, double);
void	critical_bone3(ParCurv *, double, double *, int *);
ParCurv	*cyl_section_param3(ParSurf *, double, double, int);
vector	*eval_cylsect_param3(double, int);
void	EvalMaxCamber3d(int *, double *, double *);
vector	*EvalSurfAndCyl3d(double, int);
void	fcn_crit3(int *, double *, double *, int *);
vector	*fn_camber3(double, int);
vector	*fn_thckns3(double, int);
void	fun_camber3(double *, double *, double *);
double	FuncSurfAndCyl3d(ParSurf *, double, double, double, double **,
		short *, short *, short *, int);
double	g_camber3(double *, double *);
ParCurv	*IntersectSurfAndCyl3d(ParSurf *, double, double, int, int);
ParCurv	*MakeCamberLine3d(char *, int);
ParCurv	*MakeThickness3d(char *, int);
double	max_camber3(void);
void	max_cmbr3(int *, double *, double *);
double	MaxCamber3d(ParCurv *);
void	output_ld3(double *, double *);
void	output_trl3(double *, double *);
int	sample_csurf(ParSurf *, ParCurv **, ParCurv **, double **,
		ParSurf *, double *, double *, double *, double *, int *,
		int *, int *, int, int);
void	set_camber_surface_offset(double);
void	set_constants3(ParSurf *);
void	SetConstants3d(ParSurf *, double **, short *, short *, short *, int);
void	SetConstantsR(double);
int	trim_section3(ParSurf *, double);
void	trim_start_par3(ParCurv *, double *, double *);
int	TrimSection3d(ParSurf *);

#ifdef __cplusplus
}
#endif

#endif

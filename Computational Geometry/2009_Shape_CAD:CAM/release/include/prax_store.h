/* Copyright (C) 1997 Massachusetts Institute of Technology, Cambridge, MA */
/* All rights reserved */

/* prax_store.h */

#ifndef PRAX_STORE_H
#define PRAX_STORE_H

#include <stdio.h>
#include <Xm/Xm.h>
#include "gen.h"
#include "editor.h"

#ifdef __cplusplus
extern "C" {
#endif

/* necessary to avoid conflict between X and GL */

typedef long GlObject;

/* view structure */

typedef struct {
  short		fov;		/* field-of-view, in tenths of degrees */
  short		twist;		/* rotation about z,in tenths of degrees */
  float		dist;		/* distance, in world space */
  float		far;		/* distance to far clip plane */
  float		max_range;	/* maximum x or y range, in world space */
  float		near;		/* distance to near clip plane */
  float		rotx;		/* rotation about x, in degrees */
  float		roty;		/* rotation about y, in degrees */
  float		rp[3];		/* reference point, in world space */
} ViewStruct;

/* full uv points structure */

typedef struct {
  char		*name;		/* identifier */
  char		*file;		/* file name */
  char		*copy_of;	/* identifier of copied curve */
  char		*surface;	/* surface identifier */
  ParUv		*ugeom;		/* uv points */
  float		box[6];		/* bounding box: xmin,xmax,...,zmin,zmax */
  long		flags;		/* flags */
  long		style;		/* style flags */
  long		color;		/* color */
  ViewStruct	view;		/* ViewStruct structure */
  GlObject	obj;		/* display list object */
  GlObject	polygonObj;	/* control polygon object */
  Widget	properties;	/* properties widget */
  Widget	interrogate;	/* interrogate widget */
  Widget	colorButton;	/* color button */
  Widget	uv[2];		/* uv coordinates */
  Widget	xyz[3];		/* xyz coordinates */
  Widget	uvWindow;	/* uv window */
  Widget	uvSpace;	/* uv space */
  Widget	tolerance;	/* pick tolerance */
  short		op;		/* point operation */
  float		tol;		/* pick tolerance */
  short		insert_point;	/* insert point */
} FullUv;

/* minimum distance structure */

typedef struct {
  char		*file;		/* file name */
  short		nSurfs;		/* number of surfaces */
  char		**surfs;	/* list of surfaces */
  short		*isIges;	/* iges flag */
  short		toCurve;	/* curve/surface/trimmed surf flag */
  char		**files;	/* list of files */
  int		nPoints;	/* number of points */
  double	*minDists;	/* minimum distances */
  double	**minUVs;	/* surface point UVs */
  double	**minProjs;	/* surface points */
  int		*minSurfs;	/* minimum surfaces */
  double	minDist;	/* minimum magnitude */
  double	maxDist;	/* maximum magnitude */
  double	minSign;	/* minimum signed distance */
  double	maxSign;	/* maximum signed distance */
  double	tr[7];		/* transformation */
  short		ntr;		/* number of tr */
  double	rmsDist;	/* rms distance */
} MinDistStruct;

/* full list structure */

typedef struct {
  char		*name;		/* identifier */
  char		*file;		/* file name */
  char		*copy_of;	/* identifier of copied curve */
  ListCurv	*list;		/* ListCurv structure */
  float		box[6];		/* bounding box: xmin,xmax,...,zmin,zmax */
  long		flags;		/* flags */
  long		style;		/* style flags */
  long		color;		/* color */
  ViewStruct	view;		/* ViewStruct structure */
  GlObject	obj;		/* display list object */
  GlObject	minDistObj;	/* minimum distance object */
  GlObject	minProjObj;	/* minimum projection object */
  GlObject	heightObj;	/* height function object */
  GlObject	polygonObj;	/* control polygon object */
  GlObject	ordinalObj;	/* ordinal object */
  Widget	properties;	/* properties widget */
  Widget	interrogate;	/* interrogate widget */
  Widget	colorButton;	/* color button */
  Widget	contourValue;	/* contour value */
  Widget	contourScale;	/* contour scale */
  Widget	minContValue;	/* min curvature value */
  Widget	maxContValue;	/* max curvature value */
  MinDistStruct minDist;	/* minimum distance structure */
  double	minToler;	/* minimum distance tolerance */
  short		minFlag;	/* minimum distance flag */
  short		heightcoord;	/* height function coordinate */
  short		fullrange;	/* contour scale range */
  short		contourRamp;	/* contour ramp */
  short	   	scalePosition;	/* scale thumb position */
  double	minContour;	/* minimum contour value */
  double	maxContour;	/* maximum contour value */
  double	atContour;	/* current value */
  double	contourMn;	/* minimum contour scale value */
  double	contourMx;	/* maximum contour scale value */
  short		discrete;	/* discrete vs continuous contours */
} FullList;

/* full grid structure */

typedef struct {
  char		*name;		/* identifier */
  char		*file;		/* file name */
  char		*copy_of;	/* identifier of copied grid */
  GridSurf	*grid;		/* GridSurf structure */
  float		box[6];		/* bounding box */
  long		flags;		/* flags */
  long		style;		/* style flags */
  long		color;		/* color */
  ViewStruct	view;		/* ViewStruct structure */
  GlObject	obj;		/* display list object */
  GlObject	heightObj;	/* height function object */
  GlObject	polygonObj;	/* control polygon object */
  Widget	properties;	/* properties widget */
  Widget	interrogate;	/* interrogate widget */
  Widget	colorButton;	/* color button */
  Widget	contourValue;	/* contour value */
  Widget	contourScale;	/* contour scale */
  Widget	minContValue;	/* min curvature value */
  Widget	maxContValue;	/* max curvature value */
  short		heightcoord;	/* height function coordinate */
  short		fullrange;	/* contour scale range */
  short		contourRamp;	/* contour ramp */
  short	   	scalePosition;	/* scale thumb position */
  double	minContour;	/* minimum contour value */
  double	maxContour;	/* maximum contour value */
  double	atContour;	/* current value */
  double	contourMn;	/* minimum contour scale value */
  double	contourMx;	/* maximum contour scale value */
  short		discrete;	/* discrete vs continuous contours */
} FullGrid;

/* full curve structure */

typedef struct {
  char		*name;		/* identifier */
  char		*file;		/* file name */
  char		*copy_of;	/* identifier of copied curve */
  ParCurv	*egeom;		/* ParCurv structure */
  float		box[6];		/* bounding box */
  long		flags;		/* flags */
  long		style;		/* style flags */
  long		color;		/* color */
  short		funcoord;	/* function coordinate */
  ViewStruct	view;		/* ViewStruct structure */
  GlObject	obj;		/* geometry object */
  GlObject	polygonObj;	/* control polygon object */
  GlObject	porcupineObj;	/* curvature porcupine object */
  GlObject	focalObj;	/* focal curve object */
  GlObject	radialObj;	/* radial curve object */
  GlObject	evaluateObj;	/* evaulate point object */
  GlObject	orthoObj;	/* orthotomic curve object */
  GlObject	torsionObj;	/* torsion map object */
  GlObject	knotsObj;	/* knots object */
  GlObject      offsetObj;      /* offset object */
  double	umin;		/* minimum u value */
  double	umax;		/* maximum u value */
  short		nsegs;		/* number of segments used to draw curve */
  short		nsegspan;	/* number of segments per knot span */
  Widget	properties;	/* properties widget */
  Widget	interrogate;	/* interrogate widget */
  Widget	evaluate;	/* evaluate widget */
  Widget	colorButton;	/* color button */
  Widget        evaluateLabel;  /* evaluate t */
  Widget        evaluateScale;  /* evaluate scale */
  Widget        evaluatePoint[3]; /* evaluate x,y,z */
  Widget	evaluateTangent[3]; /* evaluate tangent */
  Widget	evaluateNormal[3]; /* evaluate normal */
  Widget	curvatureK;	/* curvature k */
  Widget	torsion;	/* trosion */
  float		porcupineScale;	/* porcupine scale */
  float		focalScale;	/* focal curve scale */
  float		radialScale;	/* radial curve scale */
  float		radialPoint[3];	/* radial point x,y,z */
  float		evaluateT;	/* evaluate t */
  float		orthoPoint[3];	/* orthotomic point */
  float		orthoScale;	/* orthotomic scale */
  float		torsionScale;	/* torsion scale */
  float         offsetDist;	/* offset distance */
  float		**pts;		/* x,y,z points */
  float		**norms;	/* normals */
  double	*kcurv;		/* curvatures */
  char          *functionLabel[3]; /* function labels */
  FILE          *evalFp;        /* Evaluate file pointer */
  int           evalCount;      /* count of records in file */
} FullCurv;

/* full surf structure */

typedef struct {
  char		*name;		/* identifier */
  char		*file;		/* file name */
  char		*copy_of;	/* identifier of copied surface */
  char		*de1;		/* ruled boundary curve 1 */
  char		*de2;		/* ruled boundary curve 2 */
  ParSurf	*fgeom;		/* ParSurf structure */
  float		box[6];		/* bounding box */
  long		flags;		/* flags */
  long		style;		/* style flags */
  short		curvature;	/* curvature flags */
  long		color;		/* color */
  ViewStruct	view;		/* ViewStruct structure */
  GlObject	obj;		/* geometry object */
  GlObject	polygonObj;	/* control polygon object */
  GlObject	focalObj;	/* focal surface object */
  GlObject      offsetObj;	/* offset surface object */
  GlObject	evaluateObj;	/* evaulate point object */
  GlObject	shadedObj;	/* shaded object */
  GlObject	twoSidedObj;	/* two-sided shaded object */
  GlObject	curvatureObj;	/* curvature object */
  GlObject	orthoObj;	/* orthotomoc object */
  GlObject	heightObj;	/* height function object */
  GlObject	pivotObj;	/* pivot object */
  GlObject	*splitUobj;	/* split U object */
  GlObject	*splitVobj;	/* split V object */
  GlObject	knotsObj;	/* split V object */
  GlObject      continuityObj;  /* continuity object */
  double	umin;		/* minimum u value */
  double	umax;		/* maximum u value */
  double	vmin;		/* minimum v value */
  double	vmax;		/* maximum v value */
  short		nsegu;		/* number of segments used to draw u */
  short		nsegv;		/* number of segments used to draw v */
  short		nsubu;		/* number of subsegments used to draw u */
  short		nsubv;		/* number of subsegments used to draw v */
  short		nsegspanu;	/* number of segments per u knot span */
  short		nsegspanv;	/* number of segments per v knot span */
  short		nsegspanu2;	/* previous nsegspanu */
  short		nsegspanv2;	/* previous nsegspanv */
  short		nsubu2;		/* previous nsubu */
  short		nsubv2;		/* previous nsubv */
  Widget	properties;	/* properties widget */
  Widget	interrogate;	/* interrogate widget */
  Widget	evaluate;	/* evaluate widget */
  Widget	split;		/* split widget */
  Widget	colorButton;	/* color button */
  Widget        evaluatePoint[3]; /* evaluate x,y,z */
  Widget	evaluateTangentU[3]; /* evaluate tangent u */
  Widget	evaluateTangentV[3]; /* evaluate tangent v */
  Widget	evaluateNormal[3]; /* evaluate normal */
  Widget	evaluateLabelU;	/* evaluate u */
  Widget	evaluateLabelV;	/* evaluate v */
  Widget	evaluateScaleU; /* evaluate u scale */
  Widget	evaluateScaleV;	/* evaluate v scale */
  Widget	splitLabelU;	/* split u */
  Widget	splitLabelV;	/* split v */
  Widget	splitScaleU;	/* splti u scale */
  Widget	splitScaleV;	/* split v scale */
  Widget	styleButton;	/* shaded/wireframe button */
  Widget	contourScale;	/* contour scale */
  Widget	contourValue;	/* contour value */
  Widget	curvatureK[2];	/* min,max principal curvatures */
  Widget	curvatureNormal[2]; /* normal U,V curvatures */
  Widget	pivotScale;	/* pivot scale */
  Widget	pivotValue;	/* pivot value */
  Widget	pivotCurvature;	/* pivot curvature value */
  Widget	minContValue;	/* min curvature value */
  Widget	maxContValue;	/* max curvature value */
  Widget	offsetOne;	/* one-sided offset button */
  Widget	offsetTwo;	/* two-sided offset button */
  Widget	segScaleU;	/* segment scale for u */
  Widget	segScaleV;	/* segment scale for v */
  Widget	subScaleU;	/* sub-segment scale for u */
  Widget	subScaleV;	/* sub-segment scale for v */
  double   	atPivot;	/* current pivot value */
  short	   	pivotPosition;	/* scale thumb position */
  double	minContour;	/* min contour value */
  double	maxContour;	/* max contour value */
  double	contourMn;	/* min contour scale value */
  double	contourMx;	/* max contour scale value */
  double	atContour;	/* current value */
  short		contourRamp;	/* linear, quadratic, etc. */
  short	   	scalePosition;	/* scale thumb position */
  float		focalScale;	/* focal surface scale */
  short         focalRadii;     /* max/min radii of curvature */
  float         offsetDist;	/* offset distance */
  float		evaluateU;	/* evaluate u */
  float		evaluateV;	/* evaluate v */
  double	*splitsU;	/* split values u */
  double	*splitsV;	/* split values v */
  int		nSplitsU;	/* number of split values u */
  int		nSplitsV;	/* number of split values v */
  float		***pts;		/* x,y,z points */
  float		***norms;	/* normals */
  float		orthoPoint[3];	/* orthotomic point */
  float		orthoScale;	/* orthotomic scale */
  double	**kmin;		/* min principal curvature */
  double	**kmax;		/* max principal curvature */
  double	**normu;	/* normal u curvature */
  double	**normv;	/* normal v curvature */
  short		dirflg;		/* ruled direction flag */
  short		heightcoord;	/* height function coordinate */
  short		fullrange;	/* contour scale range */
  short		discrete;	/* discrete vs continuous contours */
  float		source[3];	/* light source for two-sided shading */
  double        continuityU;    /* continuity u */
  double        continuityV;    /* continuity v */
  double        posCont[3];     /* position */
  double        tangCont[3];    /* normal vector */
  double        curvCont;       /* normal curvature */
  FILE          *evalFp;        /* Evaluate file pointer */
  FILE          *pivotFp;       /* Pivot file pointer */
  int           evalCount;      /* count of records in file */
  int           pivotCount;     /* count of records in file */
} FullSurf;

/* full COS structure */

typedef struct {
  char		*name;		/* identifier */
  char		*file;		/* file name */
  char		*copy_of;	/* identifier of copied curve */
  char          *surface;       /* COS surface */
  ParCurv	*egeom;		/* ParCurv structure */
  float		box[6];		/* bounding box */
  long		flags;		/* flags */
  long		style;		/* style flags */
  long		color;		/* color */
  ViewStruct	view;		/* ViewStruct structure */
  GlObject	obj;		/* geometry object */
  GlObject	polygonObj;	/* control polygon object */
  GlObject	porcupineObj;	/* curvature porcupine object */
  GlObject	focalObj;	/* focal curve object */
  GlObject	evaluateObj;	/* evaulate point object */
  GlObject	knotsObj;	/* knots object */
  double	umin;		/* minimum u value */
  double	umax;		/* maximum u value */
  short		nsegs;		/* number of segments used to draw curve */
  short		nsegspan;	/* number of segments per knot span */
  Widget	properties;	/* properties widget */
  Widget	interrogate;	/* interrogate widget */
  Widget	colorButton;	/* color button */
  Widget        evaluateLabel;  /* evualate t */
  Widget        evaluateScale;  /* evualte scale */
  Widget        evaluateUV[2];  /* evaluate u,v */
  Widget        evaluatePoint[3]; /* evaluate x,y,z */
  Widget	evaluateTangentU[3]; /* evaluate tangent U */
  Widget	evaluateTangentV[3]; /* evaluate tangent V */
  Widget	evaluateNormal[3]; /* evaluate normal */
  Widget	curvatureK[2];	/* min,max principal curvature */
  Widget	curvatureNormal[2]; /* normal U,V curvature */
  Widget	evaluate;	/* evaluate widget */
  float		porcupineScale;	/* porcupine scale */
  float		focalScale;	/* focal curve scale */
  float		evaluateT;	/* evaluate t */
  float		**pts;		/* x,y,z points */
  float		**tangu;	/* tangent u */
  float		**tangv;	/* tangent v */
  float		**norms;	/* normals */
  double	*kmin;		/* min curvature */
  double	*kmax;		/* max curvature */
  double	*normu;		/* normal u curvature */
  double	*normv;		/* normal v curvature */
  short		focalRadii;	/* max radii? */
  short		porcupineRadii;	/* max radii? */
  short		orient;		/* counter-clockwise orientation */
  FILE          *evalFp;        /* evaluate file pointer */
  int           evalCount;      /* count of records in file */
} FullCos;

/* full vector structure */

typedef struct {
  char		*name;		/* identifier */
  char		*file;		/* file name */
  char		*copy_of;	/* identifier of copied curve */
  float		box[6];		/* bounding box */
  int		nPts;		/* number of points */
  double	**pts;		/* coordinates */
  char		*flg;		/* move/draw flag */
  long		flags;		/* flags */
  long		style;		/* style flags */
  long		color;		/* color */
  ViewStruct	view;		/* ViewStruct structure */
  GlObject	obj;		/* geometry object */
  Widget	properties;	/* properties widget */
} FullVect;

/* full facet structure */

typedef struct {
  char		*name;		/* identifier */
  char		*file;		/* file name */
  char		*copy_of;	/* identifier of copied curve */
  char          *surface;       /* facet surface */
  float		box[6];		/* bounding box */
  int		nPts;		/* number of points */
  double	**pts;		/* coordinates */
  int		*iEnd;		/* index matrix */
  int		nAdj;		/* number of adjacencies */
  int		*iAdj;		/* adjacency matrix */
  int		nTriang;	/* number of triangles */
  long		flags;		/* flags */
  long		style;		/* style flags */
  long		color;		/* color */
  ViewStruct	view;		/* ViewStruct structure */
  GlObject	obj;		/* geometry object */
  GlObject	shadedObj;	/* shaded object */
  GlObject	heightObj;	/* height object */
  GlObject	slopeObj;	/* slope object */
  GlObject	evaluateObj;	/* evaluate point object */
  Widget	colorButton;	/* color button */
  Widget	properties;	/* properties widget */
  Widget	interrogate;	/* interrogate widget */
  Widget	contourScale;	/* contour scale */
  Widget	contourValue;	/* contour value */
  Widget	minContValue;	/* min curvature value */
  Widget	maxContValue;	/* max curvature value */
  Widget	evaluate;	/* evaluate widget */
  Widget	evaluateLabelX;	/* evaluate x */
  Widget	evaluateLabelY;	/* evaluate y */
  Widget	evaluateLabelZ;	/* evaluate z */
  Widget	evaluateScaleX; /* evaluate y scale */
  Widget	evaluateScaleY;	/* evaluate y scale */
  Widget        evaluateXY[2];  /* evaluate x,y */
  double	minContour;	/* min contour value */
  double	maxContour;	/* max contour value */
  double	atContour;	/* current value */
  double	contourMn;	/* min contour scale value */
  double	contourMx;	/* max contour scale value */
  double	evaluateX;	/* evaluate x value */
  double	evaluateY;	/* evaluate y value */
  double	evaluateZ;	/* evaluate z value */
  short		contourRamp;	/* linear, quadratic, etc. */
  short	   	scalePosition;	/* scale thumb position */
  short		heightcoord;	/* height function coordinate */
  short		slopecoord;	/* slope coordinate */
  short		fullrange;	/* contour scale range */
  short		mapToSurface;	/* map to surface */
  short		lastA;		/* last node a visited by evaluation */
  short		lastB;		/* last node b visited by evaluation */
  short		lastC;		/* last node c visited by evaluation */
  short		discrete;	/* discrete vs continuous contours */
  FILE          *evalFp;        /* evaluate file pointer */
  int           evalCount;      /* count of records in file */
} FullFacet;

/* full trim structure */

typedef struct {
  char		*name;		/* identifier */
  char		*file;		/* file name */
  char		*copy_of;	/* identifier of copied trim */
  TrimSurf	*trim;		/* trimmed surface */
  float		box[6];		/* bounding box */
  long		flags;		/* flags */
  long		style;		/* style flags */
  short		curvature;	/* curvature flags */
  long		color;		/* color */
  ViewStruct	view;		/* ViewStruct structure */
  GlObject	obj;		/* wireframe object */
  GlObject	shadedObj;	/* shaded object */
  GlObject	polygonObj;	/* control polygon object */
  GlObject	focalObj;	/* focal surface object */
  GlObject	offsetObj;	/* focal surface object */
  GlObject	evaluateObj;	/* evaulate point object */
  GlObject	curvatureObj;	/* focal surface object */
  GlObject	orthoObj;	/* focal surface object */
  GlObject	heightObj;	/* focal surface object */
  double	umin;		/* minimum u value */
  double	umax;		/* maximum u value */
  double	vmin;		/* minimum v value */
  double	vmax;		/* maximum v value */
  short		nsegu;		/* number of segments used to draw surf */
  short		nsegv;		/* number of segments used to draw surf */
  short		nsegs;		/* number of segments used to draw curve */
  Widget	interrogate;	/* interrogate widget */
  Widget	properties;	/* properties widget */
  Widget	colorButton;	/* color button */
  Widget        evaluateUV[2];  /* evaluate u,v */
  Widget        evaluatePoint[3]; /* evaluate x,y,z */
  Widget	evaluateTangentU[3]; /* evaluate tangent U */
  Widget	evaluateTangentV[3]; /* evaluate tangent V */
  Widget	evaluateNormal[3]; /* evaluate normal */
  Widget	evaluateLabelU;	/* evaluate u */
  Widget	evaluateLabelV;	/* evaluate v */
  Widget	evaluateScaleU; /* evaluate u scale */
  Widget	evaluateScaleV;	/* evaluate v scale */
  Widget	styleButton;	/* shaded/wireframe button */
  Widget	contourScale;	/* contour scale */
  Widget	contourValue;	/* contour value */
  Widget	curvatureK[2];	/* min,max principal curvature */
  Widget	curvatureNormal[2]; /* normal U,V curvature */
  Widget	evaluate;	/* evaluate widget */
  Widget	minContValue;	/* min curvature value */
  Widget	maxContValue;	/* max curvature value */
  double	minContour;	/* min contour value */
  double	maxContour;	/* max contour value */
  double	atContour;	/* current value */
  double	contourMn;	/* min contour value */
  double	contourMx;	/* max contour value */
  short		contourRamp;	/* linear, quadratic, etc. */
  short	   	scalePosition;	/* scale thumb position */
  float		focalScale;	/* focal curve scale */
  short         focalRadii;     /* max/min radii of curvature */
  float         offsetDist;	/* offset distance */
  float		evaluateU;	/* evaluate u */
  float		evaluateV;	/* evaluate v */
  float		**pts;		/* x,y,z points */
  float		**norms;	/* normals */
  float		orthoPoint[3];	/* orthotomic point */
  float		orthoScale;	/* orthotomic scale */
  double	*kmin;		/* min curvature */
  double	*kmax;		/* max curvature */
  double	*normu;		/* normal u curvature */
  double	*normv;		/* normal v curvature */
  short		heightcoord;	/* height function coordinate */
  short		fullrange;	/* contour scale range */
  int		nPts;		/* number of trimmed points */
  int		nTriang;	/* number of triangles */
  int		*iEnd;		/* index matrix */
  int		*iAdj;		/* adjacency matrix */
  short		discrete;	/* discrete vs continuous contours */
  FILE          *evalFp;        /* evaluate file pointer */
  int           evalCount;      /* count of records in file */
} FullTrim;

/* full hull structure */

typedef struct {
  char		*name;		/* identifier */
  char		*file;		/* file name */
  char		*copy_of;	/* identifier of copied hull */
  float		box[6];		/* bounding box */
  long		flags;		/* flags */
  long		style;		/* style flags */
  long		color;		/* color */
  int		nPts;		/* number of points */
  int		nTris;		/* number of faces */
  int		nAdjs;		/* number of adjacencies */
  int		nEdges;		/* number of edges */
  int		nonConvex;	/* number of non-convex edges */
  double	**pts;		/* array of points */
  int		**tris;		/* array of points */
  int		**adjs;		/* array of points */
  int		*nc;		/* array of points */
  ViewStruct	view;		/* ViewStruct structure */
  GlObject	obj;		/* geometry object */
  GlObject	shadedObj;	/* shaded object */
  Widget	colorButton;	/* color button */
  Widget	properties;	/* properties widget */
  Widget	interrogate;	/* interrogate widget */
} FullHull;

/* full insp structure */

typedef struct {
  char		*name;		/* identifier */
  char		*file;		/* file name */
  char		*copy_of;	/* identifier of copied insp */
  float		box[6];		/* bounding box */
  long		flags;		/* flags */
  long		style;		/* style flags */
  long		color;		/* color */
  int		nPts;		/* number of points */
  double	**pts;		/* array of points */
  double	**norms;        /* array of normals */
  double        scale;          /* normal scale */
  ViewStruct	view;		/* ViewStruct structure */
  GlObject	obj;		/* geometry object */
  GlObject      polygonObj;     /* polygon object */
  Widget	colorButton;	/* color button */
  Widget	properties;	/* properties widget */
  Widget	interrogate;	/* interrogate widget */
  Widget	normalScale;	/* normal scale */
} FullInsp;

/* full group structure */

typedef struct {
  char		*name;		/* identifier */
  long		flags;		/* flags */
  long		style;		/* style flags */
  short		curvature;	/* curvature flags */
  long		color;		/* color */
  short		nLists;		/* # lists in group */
  short		nGrids;		/* # grids in group */
  short		nCurvs;		/* etc. */
  short		nSurfs;
  short		nCoss;
  short		nUVs;
  short		nVects;
  short		nFacets;
  short		nTrims;
  short		nHulls;
  short         nInsps;
  FullList	**lists;	/* lists in group */
  FullGrid	**grids;	/* grids in group */
  FullCurv	**curvs;	/* etc. */
  FullSurf	**surfs;
  FullCos	**coss;
  FullUv       	**uvs;
  FullVect	**vects;
  FullFacet	**facets;
  FullTrim	**trims;
  FullHull	**hulls;
  FullInsp      **insps;
  Widget	interrogate;	/* dialog box widget */
  Widget	properties;	/* properties widget */
  Widget	colorButton;	/* color button */
  Widget	contourValue;	/* contour value */
  Widget	contourScale;	/* contour scale */
  Widget	minContValue;	/* min curvature value */
  Widget	maxContValue;	/* max curvature value */
  double	minToler;	/* minimum distance tolerance */
  short		minFlag;	/* minimum distance flag */
  short		heightcoord;	/* height function coordinate */
  short		slopecoord;	/* slope coordinate */
  short		fullrange;	/* contour scale range */
  short		contourRamp;	/* contour ramp */
  short	   	scalePosition;	/* scale thumb position */
  double	minContour;	/* minimum contour value */
  double	maxContour;	/* maximum contour value */
  double	atContour;	/* current value */
  double	contourMn;	/* minimum contour scale value */
  double	contourMx;	/* maximum contour scale value */
  double	minDist;	/* minimum magnitude */
  double	maxDist;	/* maximum magnitude */
  double	minSign;	/* minimum signed distance */
  double	maxSign;	/* maximum signed distance */
  short		mapToSurface;	/* map to surface */
  short		discrete;	/* discrete vs continuous contours */
} FullGroup;

/* function prototypes */

short	CheckCos(ParCurv *);
short	CopyCos(short, FullCos *);
short	CopyCurv(short, FullCurv *);
short	CopyFacet(short, FullFacet *);
short	CopyGrid(short, FullGrid *);
short	CopyList(short, FullList *);
void	CopyMinDist(MinDistStruct *, MinDistStruct *);
short	CopySurf(short, FullSurf *);
short	CopyTrim(short, FullTrim *);
short	CopyUv(short, FullUv *);
short	CopyVect(short, FullVect *);
void    CosEvaluateFile(FullCos *, int);
void	CosEvaluateUpdate(FullCos *);
void	CreateUvWindow(FullUv *);
void    CurvEvaluateFile(FullCurv *, int);
void	CurvEvaluateUpdate(FullCurv *);
void	DeleteCosObjs(FullCos *);
void	DeleteCosPts(FullCos *);
void	DeleteCurvObjs(FullCurv *);
void	DeleteCurvPts(FullCurv *);
void    DeleteFacetObjs(FullFacet *);
void	DeleteSurfObjs(FullSurf *);
void	DeleteSurfPts(FullSurf *);
void	DeleteTrimObjs(FullTrim *);
void	DeleteTrimPts(FullTrim *);
void	DrawListMinDist(ListCurv *, MinDistStruct, short, double, short,
		short, double, double, short, double *, double *,
		double, double, short);
void	DrawListMinProj(ListCurv *, MinDistStruct);
void	DrawUvSpace(Widget, FullUv *);
void	EditKnotsCos(FullCos *);
void	EditKnotsCurv(FullCurv *);
void	EditKnotsSurf(FullSurf *);
Widget	EvalCos(FullCos *);
void    EvalCosOk(FullCos *);
Widget	EvalCurv(FullCurv *);
void    EvalCurvOk(FullCurv *);
Widget	EvalFacet(FullFacet *);
void    EvalFacetOk(FullFacet *);
Widget	EvalSurf(FullSurf *);
void    EvalSurfOk(FullSurf *);
Widget	EvalTrim(FullTrim *);
void    EvalTrimOk(FullTrim *);
void    FacetEvaluateFile(FullFacet *, int);
void	FacetEvaluateUpdate(FullFacet *);
void	FairCos(FullCos *);
void    FairCosScript(FullCos *);
void	FairCurv(FullCurv *);
void    FairCurvScript(FullCurv *);
void	FairSurf(FullSurf *);
void	FairSurfScript(FullSurf *, int);
FullCurv *FindCurv(char *);
FullSurf *FindSurf(char *);
FullTrim *FindTrim(char *);
FullCos	*GetCos(short);
FullCurv *GetCurv(short);
void	GetEntityView(ViewStruct *, short *, float *, float *, float *,
		float *, float *, float *, short *, float *);
FullFacet *GetFacet(short);
FullGrid *GetGrid(short);
FullGroup *GetGroup(short);
FullHull *GetHull(short);
FullInsp *GetInsp(short);
FullList *GetList(short);
FullSurf *GetSurf(short);
FullTrim *GetTrim(short);
FullUv	*GetUv(short);
FullVect *GetVect(short);
void	InitFullCos(FullCos *);
void	InitFullCurv(FullCurv *);
void	InitFullFacet(FullFacet *);
void	InitFullGrid(FullGrid *);
void	InitFullHull(FullHull *);
void    InitFullInsp(FullInsp *);
void	InitFullList(FullList *);
void	InitFullSurf(FullSurf *);
void	InitFullTrim(FullTrim *);
void	InitFullUv(FullUv *);
void	InitFullVect(FullVect *);
void	InitMinDist(FullList *);
void	InitializeView(ViewStruct *, float *);
FullInsp *InspectCurves(int, ParSurf *, int, ParCurv **, FullFacet *, int,
			int);
void    InterrCosOk(FullCos *);
void    InterrCurvOk(FullCurv *);
void    InterrFacetOk(FullFacet *);
void    InterrGridOk(FullGrid *);
void    InterrGroupOk(FullGroup *);
void    InterrGroupUpdate(FullGroup *);
void    InterrHullOk(FullHull *);
void    InterrInspOk(FullInsp *);
void    InterrListOk(FullList *);
void    InterrSurfOk(FullSurf *);
void    InterrTrimOk(FullTrim *);
void    InterrUvOk(FullUv *);
void	InterrGroupSetValues(FullGroup *);
Widget	InterrogateCos(FullCos *);
Widget	InterrogateCurv(FullCurv *);
Widget	InterrogateFacet(FullFacet *);
Widget	InterrogateGrid(FullGrid *);
Widget	InterrogateGroup(FullGroup *);
Widget	InterrogateHull(FullHull *);
Widget  InterrogateInsp(FullInsp *);
Widget	InterrogateList(FullList *);
Widget	InterrogateSurf(FullSurf *);
Widget	InterrogateTrim(FullTrim *);
Widget	InterrogateUv(FullUv *);
void	PostScriptListMinDist(ListCurv *, MinDistStruct, short, double,
		short, short, double, double, short, double *, double *,
		double, double, short);
void	PostScriptListMinProj(ListCurv *, MinDistStruct);
int	ReflectCos(int, FullCos *, int);
int	ReflectCurv(int, FullCurv *, int);
int	ReflectFacet(int, FullFacet *, int);
int	ReflectGrid(int, FullGrid *, int);
int	ReflectList(int, FullList *, int);
int	ReflectSurf(int, FullSurf *, int);
int	ReflectTrim(int, FullTrim *, int);
int	ReflectVect(int, FullVect *, int);
short	SaveDeslabCos(FullCos *, char *, short);
short	SaveDeslabCurv(FullCurv *, char *, short);
short	SaveDeslabFacet(FullFacet *, char *, short);
short	SaveDeslabGrid(FullGrid *, char *, short);
short   SaveDeslabInsp(FullInsp *, char *, short);
short	SaveDeslabList(FullList *, char *, short);
short	SaveDeslabSurf(FullSurf *, char *, short);
short	SaveDeslabTrim(FullTrim *, char *, short);
short	SaveDeslabUv(FullUv *, char *, short);
int	SaveIgesCos(FullCos *, char *, int);
int	SaveIgesCurv(FullCurv *, char *, int);
int	SaveIgesGrid(FullGrid *, char *);
int     SaveIgesInsp(FullInsp *, char *);
int	SaveIgesList(FullList *, char *);
int	SaveIgesMinDist(FullList *, char *, int);
int	SaveIgesSurf(FullSurf *, char *, int);
int	SaveIgesTrim(FullTrim *, char *);
int	SaveIgesUv(FullUv *, char *);
FullCos	*SetCos(short, FullCos *);
FullCurv *SetCurv(short, FullCurv *);
FullFacet *SetFacet(short, FullFacet *);
void	SetFacetContourCallbacks(FullFacet *);
FullGrid *SetGrid(short, FullGrid *);
void	SetGridContourCallbacks(FullGrid *);
FullGroup *SetGroup(short, FullGroup *);
void	SetGroupContourCallbacks(FullGroup *);
FullHull *SetHull(short, FullHull *);
FullInsp *SetInsp(short, FullInsp *);
FullList *SetList(short, FullList *);
void	SetListContourCallbacks(FullList *);
FullSurf *SetSurf(short, FullSurf *);
void	SetSurfContourCallbacks(FullSurf *);
FullTrim *SetTrim(short, FullTrim *);
void	SetTrimContourCallbacks(FullTrim *);
FullUv	*SetUv(short, FullUv *);
void	SetUvCallbacks(FullUv *);
FullVect *SetVect(short, FullVect *);
void    SplitSurfCancel(FullSurf *);
Widget	SplitSurfDialog(FullSurf *);
void    SplitSurfOk(FullSurf *);
void	StoreTrimValues(TrimSurf *, int, int, float ***, float ***,
		double **, double **, double **, double **, int,
		int *, int *, int **, int **);
void	SubdivideCos(FullCos *);
void	SubdivideCurv(FullCurv *);
void	SubdivideSurf(FullSurf *);
void    SurfEvaluateFile(FullSurf *, int);
void	SurfEvaluateUpdate(FullSurf *);
void    SurfPivotFile(FullSurf *, int);
void	SurfPivotUpdate(FullSurf *);
short	TransformCos(short, FullCos *);
short	TransformCurv(short, FullCurv *);
short	TransformFacet(short, FullFacet *);
short	TransformGrid(short, FullGrid *);
short	TransformList(short, FullList *);
short	TransformSurf(short, FullSurf *);
short	TransformTrim(short, FullTrim *);
short	TransformVect(short, FullVect *);
void    TrimEvaluateFile(FullTrim *, int);
void	TrimEvaluateUpdate(FullTrim *);
void	UpdateFacetContourScale(FullFacet *);
void	UpdateGridContourScale(FullGrid *);
void	UpdateGroupContourScale(FullGroup *);
void	UpdateListContourScale(FullList *);
void	UpdateSurfContourScale(FullSurf *);
void	UpdateTrimContourScale(FullTrim *);
void    UvCancel(FullUv *);
void    UvOk(FullUv *);

#ifdef __cplusplus
}
#endif

#endif

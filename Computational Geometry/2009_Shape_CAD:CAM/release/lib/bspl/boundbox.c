/* ***************************************************************************
 Copyright (C) 1996 Massachusetts Institute of Technology all rights reserved 
	Programmer: George A. Kriezis
**************************************************************************** */
#include <math.h>
#include "bspl.h"

/*****************************************************************************
*                                  boundbox
******************************************************************************
* 
* 1   Purpose
*     This function determines the min/max bounding box for a NURBS surface.
* 
* 2   Specification
*     #include "bspl.h"
*     void boundbox(ParSurf *fgeom, double *hull)
* 
* 3   Description
*     This routine finds the min/max bounding box of a NURBS surface.
* 
* 5   Parameters
*       1.ParSurf * fgeom
*         On entry: NURBS surface data structure containing surface of which a 
*                	  bounding box is to be computed.
* 
*       2.double * hull
*         On entry: array of lenght at least 6.
*         On exit: bounding box of the NURBS surface:
*                  hull[0] = minimum in x direction
* 		 hull[1] = maximum in x direction
* 		 hull[2] = minimum in y direction
* 		 hull[3] = maximum in y direction
* 		 hull[4] = minimum in z direction
* 		 hull[5] = maximum in z direction
* 
* 6   Return Values, Error Indicators and Warnings
*     Not applicable.
* 
* 7   Accuracy
*     Not applicable.
* 
* 8   Further Comments
*     Not applicable.
* 
* 9   Functions referenced by boehm_surface() are:
* 
* 10  Functions that reference boehm_surface() are:
* 
*****************************************************************************/

void boundbox(ParSurf *fgeom, double hull[])
{
  double amin,amax;             /* small number and large number for
				 * comparison */
  double w;                     /* homogeneous coordinate */    
  int i,j,k,jj;                 /* i,k -- indices
				 * j,jj -- order in u,v direction,
				 *         respectively */

  /* set the small and the large number */
  amin = 1000000.;
  amax = -1000000.;
  /* get the orders */
  j = fgeom->ucontpts;  jj = fgeom->vcontpts;

  /* find xmin and xmax */
  for (i=0; i<j; i++){
      for (k=0; k<jj; k++){
          w = fgeom->contpts[i][k]->w;
          if (amin > fgeom->contpts[i][k]->x/w)
	     amin = fgeom->contpts[i][k]->x/w;
          if (amax < fgeom->contpts[i][k]->x/w)
	     amax = fgeom->contpts[i][k]->x/w;
          }
      }
  hull[0]=amin;
  hull[1]=amax;

  /* find ymin and ymax */
  amin = 10000000.;
  amax = -10000000.;
  for (i=0; i<j; i++){
      for (k=0; k<jj; k++){
          w = fgeom->contpts[i][k]->w;
          if (amin > fgeom->contpts[i][k]->y/w)
	     amin = fgeom->contpts[i][k]->y/w;
          if (amax < fgeom->contpts[i][k]->y/w)
             amax = fgeom->contpts[i][k]->y/w;
          }
      }
  hull[2]=amin;
  hull[3]=amax;

  /* find zmin and zmax */
  amin = 10000000.;
  amax = -10000000.;
  for (i=0; i<j; i++){
      for (k=0; k<jj; k++){
          w = fgeom->contpts[i][k]->w;
          if (amin > fgeom->contpts[i][k]->z/w)
             amin = fgeom->contpts[i][k]->z/w;
          if (amax < fgeom->contpts[i][k]->z/w)
             amax = fgeom->contpts[i][k]->z/w;
         }
     }
  hull[4]=amin;
  hull[5]=amax;
}

/****************************************************************************
*                             boundbox_c()
*****************************************************************************
* 
* 1   Purpose
*     This function determines the min/max bounding box for a NURBS curve.
* 
* 2   Specification
*     #include "bspl.h"
*     void boundbox(ParCurv *egeom, double *hull)
* 
* 3   Description
*     This routine finds the min/max bounding box of a NURBS curve.
* 
* 5   Parameters
*        1.ParCurv * egeom
*          On entry: NURBS curve data structure containing curve of which a 
* 	 bounding box is to be computed.
* 
*        2.double * hull
*          On entry: array of lenght at least 6.
*          On exit: bounding box of the NURBS curve:
*           hull[0] = minimum in x direction
*           hull[1] = maximum in x direction
*           hull[2] = minimum in y direction
*           hull[3] = maximum in y direction
*           hull[4] = minimum in z direction
*           hull[5] = maximum in z direction
* 
* 6   Return Values, Error Indicators and Warnings
*     Not applicable.
* 
* 7   Accuracy
*     Not applicable.
* 
* 8   Further Comments
*     Not applicable.
* 
* 9   Functions referenced by boehm_curve() are: None
* 
* 10  Functions that reference boehm_curve() are: None
* 
*****************************************************************************/

void boundbox_c(ParCurv *fgeom, double hull[])
{
  double amin,amax;   /* small and large number for comparison */
  double w;           /* homogeneous coordinate */
  int i,j;            /* i -- index
		       * j -- order */

  amin = 1000000.;
  amax = -1000000.;
  j = fgeom->ncontpts; 

  /* find xmin and xmax */
  for (i=0; i<j; i++){
    w = fgeom->contpts[i]->w;
    if (amin > fgeom->contpts[i]->x/w)
      amin = fgeom->contpts[i]->x/w;
    if (amax < fgeom->contpts[i]->x/w)
      amax = fgeom->contpts[i]->x/w;
  }
  hull[0]=amin;
  hull[1]=amax;

  /* find ymin and ymax */
  amin = 10000000.;
  amax = -10000000.;
  for (i=0; i<j; i++){
    w = fgeom->contpts[i]->w;
    if (amin > fgeom->contpts[i]->y/w)
      amin = fgeom->contpts[i]->y/w;
    if (amax < fgeom->contpts[i]->y/w)
      amax = fgeom->contpts[i]->y/w;
  }
  hull[2]=amin;
  hull[3]=amax;

  /* find zmin and zmax */
  amin = 10000000.;
  amax = -10000000.;
  for (i=0; i<j; i++){
    w = fgeom->contpts[i]->w;
    if (amin > fgeom->contpts[i]->z/w)
      amin = fgeom->contpts[i]->z/w;
    if (amax < fgeom->contpts[i]->z/w)
      amax = fgeom->contpts[i]->z/w;
  }
  hull[4]=amin;
  hull[5]=amax;
}

/**************************************************************************
*                                compare_boxes()
***************************************************************************
* 
* 1   Purpose
*     This function compares min/max bounding boxes of two geometric entities 
*     for possible intersection.
* 
* 2   Specification
*     #include "bspl.h"
*     int compare_boxes(double *hull1, double *hull2)
* 
* 3   Description
*     This routine compares for intersection the min/max bounding boxes of 
*     two NURBS curves or surfaces. The routine returns 1 if there is an 
*     intersection, else it returns 0.
* 
* 5   Parameters
*        1.double * hull1
*          On entry: array of length at least 6 containing the minimum and 
* 	           maxima extents of a min/max box for the two geometric 
* 		   objects.
*        2.double * hull2
*          On entry: array of length at least 6 containing the minimum and 
* 	           maxima extents of a min/max box for the two geometric 
* 		   objects.
* 
* 6   Return Values, Error Indicators and Warnings
*     Not applicable
* 
* 7   Accuracy
*     Not applicable
* 
* 8   Further Comments
*     Bounding boxes are computed using routines boundbox (bspl _ boundbox.c) 
*     or boundbox_c (bspl _ boundbox.c).
* 
* 9   Functions referenced by compare_boxes() are: None
* 
* 10  Functions that reference compare_boxes() are: None
* 
***************************************************************************/
                              
int compare_boxes(double hull1[],double hull2[])
{
  int i;         /* return value */
  double eps;    /* samll number for comparison */
  eps = -1.e-16;

  /* i=0 no overlap, i=1 overlap */
  if (hull1[1] < hull2[0]+eps || hull2[1] < hull1[0]+eps)
    i=0;
  else{
    if (hull1[3] < hull2[2]+eps || hull2[3] < hull1[2]+eps)
      i=0;
    else{
      if (hull1[5] < hull2[4]+eps || hull2[5] < hull1[4]+eps)
	i=0;
      else
	i=1;
    }
  }
  return(i);
}

/**************************************************************************
*                                absmaxarray1()
***************************************************************************
* 
* 1   Purpose
*     Find the absolute maximum value in an array.
* 
* 2   Specification
*     #include "bspl.h"
*     double absmaxarray1(double *array, int n);
* 
* 3   Description
*     This routine finds the maximum absolute value (maximum magnitude) in an 
*     array.
* 
* 4   References
*     Not applicable.
* 
* 5   Parameters
*       1.double * array
*         On entry: the address of the array.
*       2.int n
*         On entry: the number of elements in the array.
* 
* 6   Return Values, Error Indicators and Warnings
*     The value of the maximum absolute value (maximum magnitude) is returned.
* 
* 7   Accuracy
*     Not applicable.
* 
* 8   Further Comments
*     Not applicable.
* 
* 9   Functions referenced by absmaxarray() are: None
* 
* 10  Functions that reference absmaxarray() are: monotobern()
* 
***************************************************************************/

double absmaxarray1(double *array,int n)
{
  /* find absolute maximum in array */
  double maxim;
  int i;
  maxim = FABS(array[0]);
  for (i=1; i<n; i++)
    if (FABS(array[i]) > maxim)
      maxim = FABS(array[i]);
  return(maxim);
}

/*****************************************************************************
*                                   surface_scale_trans()
******************************************************************************
* 
* 1   Purpose
*     This routine translate and scale the control points of a B-spline
*     surface.
* 
* 2   Specification
*     #include "bspl.h"
*     void surface_scale_trans(ParSurf *fgeom, double size, double tx,
*                              double ty, double tz) 
* 
* 3   Description
*     This routine translate the control points of a B-spline surface by tx, 
*     ty, tz, and then scale down the control points by a factor. 
* 
* 4   References
*     Not applicable
* 
* 5   Parameters
*       1.ParSurf * fgeom
*         On entry: NURBS surface data structure containing a B-spline surface.
* 	On exit:  NURBS surface data structure of the surface after translation
* 	          and scaling.
*       2.double size
*         On entry: the scaling factor.
*       3.double tx,ty,tz
*         On entry: the x,y,z components of the translating vector.
* 
* 6   Return Values, Error Indicators and Warnings
*     Not applicable.
* 
* 7   Accuracy
*     Not applicable.
* 
* 8   Further Comments
*     Not applicable.
* 
* 9   Functions referenced by surface_scale_trans() are: None
* 
* 10  Functions that reference surface_scale_trans() are: None
* 
*****************************************************************************/

void surface_scale_trans(ParSurf *fgeom, double size, double tx, double ty,
		    double tz)
{
  int i,j;
  double w;
  for (i=0; i<fgeom->ucontpts; i++)
    for (j=0; j<fgeom->vcontpts; j++){
      w = fgeom->contpts[i][j]->w;
      fgeom->contpts[i][j]->x = (fgeom->contpts[i][j]->x+w*tx)/size;
      fgeom->contpts[i][j]->y = (fgeom->contpts[i][j]->y+w*ty)/size;
      fgeom->contpts[i][j]->z = (fgeom->contpts[i][j]->z+w*tz)/size;
    }
}

/*****************************************************************************
*                                 surface_translate()
******************************************************************************
* 
* 1   Purpose
*     This function determines the translation and the scaling needed for making
*     two surfaces lie in the -1 to 1 box.
* 
* 2   Specification
*     #include "bspl.h"
*     void surface_translate(double *hull1, double *hull2, double *tx,,
* 		             double *tydouble *tz, double *size)
* 
* 3   Description
*     This routine uses the min/max bounding boxes of the two surfaces to
*     determine the translating values and the scaling factor which must be
*     applied on the two surfaces to make these two surfaces lie in the -1
*     to 1 box.
* 
* 4   References
*     Nonapplicable
* 
* 5   Parameters
*       1.double *hull1
*         On entry: bounding box of the first surface.
* 	On exit:  bounding box of the first surface after translation and
*                 scaling.
*       2.double *hull2
*         On entry: bounding box of the second surface.
* 	On exit:  bounding box of the second surface after translation and
*                 scaling.
*       3.double *tx
*         On entry: a pointer to a double data.
* 	On exit:  address of the element which contain x component of the
* 	          translating vector.
*       4.double *ty
*         On entry: a pointer to a double data.
* 	On exit:  address of the element which contain y component of the
* 	          translating vector.
*       5.double *tz
*         On entry: a pointer to a double data.
* 	On exit:  address of the element which contain z component of the
* 	          translating vector.
*       6.double *size
*         On entry: a pointer to a double data.
* 	On exit:  address of the element which contain the scaling factor.
* 
* 6   Return Values, Error Indicators and Warning
*     Not applicable.
* 
* 7   Further Comments
*     Not applicable.
* 
* 9   Functions referenced by surface_translate() are: None
* 
* 10  Functions that reference surface_translate() are: None
* 
*****************************************************************************/

void surface_translate(double *hull1, double *hull2, double *tx, double *ty,
		       double *tz, double *size)
/* this routine determines the translation to add to each of the surfaces
 * and the scaling factor in order to make the two surfaces lie in the -1 to 1 
 * box. It uses the min,max bounding boxes to do this calculation */
{
  double maxi,mini,siz;  /* maxi -- the maximum in one direction 
			  * mini -- the minimum in one direction
			  * siz -- the half edge of the box in one direction */

  /* determine the minimum and the maximum in x direction */
  if (hull1[0] < hull2[0])
    mini = hull1[0];
  else
    mini = hull2[0];
  if (hull1[1] > hull2[1])
    maxi = hull1[1];
  else
    maxi = hull2[1];
  /* calculate the translation in x direction */
  *tx = -(maxi + mini)/2.;
  /* calculate the scaling factor */
  *size = (maxi-mini)/2.;
 
  /* determine the minimum and the maximum in y direction */
  if (hull1[2] < hull2[2])
    mini = hull1[2];
  else
    mini = hull2[2];
  if (hull1[3] > hull2[3])
    maxi = hull1[3];
  else
    maxi = hull2[3];
  /* calculate the translation in y direction */
  *ty = -(maxi + mini)/2.;
  /* update the scaling factor */
  siz = (maxi-mini)/2.;
  if (siz > *size)
    *size = siz;

  /* determine the minimum and the maximum in z direction */
  if (hull1[4] < hull2[4])
    mini = hull1[4];
  else
    mini = hull2[4];
  if (hull1[5] > hull2[5])
    maxi = hull1[5];
  else
    maxi = hull2[5];
  /* calculate the translation in z direction */
  *tz = -(maxi + mini)/2.;
  siz = (maxi-mini)/2.;
  /* update the scaling factor */
  if (siz > *size)
    *size = siz;

  /* update bounding boxes also */
  hull1[0] = (hull1[0]+ *tx)/ *size;
  hull1[1] = (hull1[1]+ *tx)/ *size;
  hull1[2] = (hull1[2]+ *ty)/ *size;
  hull1[3] = (hull1[3]+ *ty)/ *size;
  hull1[4] = (hull1[4]+ *tz)/ *size;
  hull1[5] = (hull1[5]+ *tz)/ *size;
  hull2[0] = (hull2[0]+ *tx)/ *size;
  hull2[1] = (hull2[1]+ *tx)/ *size;
  hull2[2] = (hull2[2]+ *ty)/ *size;
  hull2[3] = (hull2[3]+ *ty)/ *size;
  hull2[4] = (hull2[4]+ *tz)/ *size;
  hull2[5] = (hull2[5]+ *tz)/ *size;
}

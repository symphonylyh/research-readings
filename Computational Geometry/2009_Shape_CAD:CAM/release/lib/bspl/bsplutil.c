/* ***************************************************************************
 Copyright (C) 1996 Massachusetts Institute of Technology all rights reserved 
	Programmer: George A. Kriezis
**************************************************************************** */
#include <stdio.h>
#include <math.h>
#include "gen.h"
#include "bspl.h"

/*****************************************************************************
*                              writevect()
******************************************************************************
* 
* 1   Purpose
*     Write a vector structure to a file.
* 
* 2   Specification
*     #include "bspl.h"
*     void writevect(FILE *fp, vector *v)
* 
* 3   Description
*     This function writes the values of a vector structure to an open file.
* 
* 4   References
*     Not applicable.
* 
* 5   Parameters
*       1.FILE * fp
*         On entry: the file pointer of a file open for writing.
*       2.vector * v
*         On entry: the address of a vector structure whose contents will be 
* 	written to the file.
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
****************************************************************************/

void writevect(FILE *fp, vector *v)
{
  fprintf(fp,"%20.16lf %20.16lf %20.16lf %20.16lf\n",v->x,v->y,v->z,v->w) ;
}

/*****************************************************************************
*                             alloc_fgeompts()
******************************************************************************
* 
* 1   Purpose
*     Allocate NURBS surface control point array.
* 
* 2   Specification
*     #include "bspl.h"
*     void alloc_fgeompts(ParCurv *egeom)
* 
* 3   Description
*     This function allocates dynamic memory to contain the control point 
*     coordinates of a ParSurf NURBS surface structure.
* 
* 4   References
*     Not applicable.
* 
* 5   Parameters
*       1.ParSurf * fgeom
*         On entry: the address of a NURBS surface structure.
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
* 9   Functions referenced by alloc_fgeompts() are:
*     vectalloc()
*  
* 10  Functions that reference alloc_fgeompts() are:
*     integral_offset()
*     integral_offset_surf()
*     merge_can_off()
*     rational_offset()
* 
******************************************************************************/

void alloc_fgeompts(ParSurf *fgeom)
/* Allocate memory (if necessary) to the elements of the array containing
 * the coordinates of the control points of a rational B-spline surface 
 * (fgeom->contpts) */
{
  int i,j ;

  for (i=0;i<fgeom->ucontpts;i++)
      for (j=0;j<fgeom->vcontpts;j++)
          if (fgeom->contpts[i][j]==NULL)
             fgeom->contpts[i][j] = vectalloc() ;
}

/*****************************************************************************
*                            alloc_egeompts()
******************************************************************************
* 
* 1   Purpose
*     Allocate NURBS curve control point array.
* 
* 2   Specification
*     #include "bspl.h"
*     void alloc_egeompts(ParCurv *egeom)
* 
* 3   Description
*     This function allocates dynamic memory to contain the control point 
*     coordinates of a ParCurv NURBS curve structure.
* 
* 4   References
*     Not applicable.
* 
* 5   Parameters
*       1.ParCurv * egeom
*         On entry: the address of a NURBS curve structure.
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
* 9   Functions referenced by alloc_egeompts() are:
*     vectalloc()
* 
* 10  Functions that reference alloc_egeompts() are: None
* 
*****************************************************************************/

void alloc_egeompts(ParCurv *egeom)
/* allocate memory (if necessary) to the elements of the array containing
 * the coordinates of the control points of a rational B-spline 
 * (egeom->contpts) */
{
int i ;
for (i=0;i<egeom->ncontpts;i++)
  if (egeom->contpts[i]==NULL) egeom->contpts[i] = vectalloc() ;
}

/*****************************************************************************
*                                extract_edge()
******************************************************************************
* 
* 1   Purpose
*     This function extracts an edge isoparametric curve, defined as a NURBS 
*     curve, from a NURBS surface in either parametric direction.
* 
* 2   Specification
*     #include "bspl.h"
*     void extract_edge(ParCurve *egeom, ParSurf *fgeom, int bool, int k)
* 
* 3   Description
*     This routine extracts a NURBS curve from the edge of a NURBS surface for 
*     either parametric direction.
* 
* 4   References
*     [1] S. T. Tuohy. Sculptured shape creation, approximation, and 
*         interrogation. Engineer's Thesis, Massachusetts Institute of, 
* 	  Technology Department of Ocean Engineering, Cambridge,
*          Massachusetts, 1992.
* 
* 5   Parameters
*        1.ParCurv * egeom
*          On entry: a NURBS data structure of allocated for the correct order 
* 	           and number of control points coinciding with the desired
*                  egde.
*          On exit: geometry representing the edge of the NURBS surface
*        2.ParSurf * fgeom
*          On entry: the NURBS surface data structure containing the surface 
* 	           geometry from which an edge will be extracted.
*        3.int bool
*          On entry: specifies if an edge in the u or v direction is to be 
* 	           extracted according to:
* 		       bool = 1 if u is required
* 		       bool = 0 if v is required
*        4.int k
*          On entry: specifies which edge is to be extracted according to:
*                        k = 0 if bool = 0 and the edge required is u=uknot[0]
*                        k = 0 if bool = 1 and the edge required is v=vknot[0]
*                        k = ucontpts if bool = 0 and the edge required is 
* 		                                       u = uknot[ucontpts]
*                        k = vcontpts if bool = 1 and the edge required is 
* 		                                       v = vknot[vcontpts]
* 
* 6   Return Values, Error Indicators and Warnings
*     This routine does not check for the value of k.
* 
* 7   Accuracy
*     Not applicable.
* 
* 8   Further Comments
*     Not applicable.
* 
* 9   Functions referenced by extract_edge() are:
*     copyvector()
*     vectalloc()
* 
* 10  Functions that reference extract_edge() are:
*     calc_disc_surf()
*     corner_offset()
*     fair_cubic_surf()
*     merge_can_off()
*     merge_fgeom()
*     ParCurv_iso()
*     raise_surf()
* 
*****************************************************************************/

void extract_edge (ParCurv *eg, ParSurf *fgeom, int n0, int k)
{
  int i ;

  if (n0) {                          /* u direction is required */
     eg->order = fgeom->uorder ; 
     eg->ncontpts = fgeom->ucontpts ;
     for (i=0;i<fgeom->ucontpts+fgeom->uorder;i++)
       eg->knots[i] = fgeom->uknots[i];
     for (i=0;i<fgeom->ucontpts;i++){
         if (eg->contpts[i] == NULL) 
             eg->contpts[i] = vectalloc();
         copyvector(fgeom->contpts[i][k],eg->contpts[i]) ;
         }
     }
  else {                             /* v direction is required */
     eg->order = fgeom->vorder ; 
     eg->ncontpts = fgeom->vcontpts ;
     for (i=0;i<fgeom->vcontpts+fgeom->vorder;i++)
       eg->knots[i] = fgeom->vknots[i];
     for (i=0;i<fgeom->vcontpts;i++){
         if (eg->contpts[i] == NULL) 
            eg->contpts[i] = vectalloc();
         copyvector(fgeom->contpts[k][i],eg->contpts[i]) ;
         }
     }
}

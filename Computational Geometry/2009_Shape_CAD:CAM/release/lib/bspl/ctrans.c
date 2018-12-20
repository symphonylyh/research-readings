/* ***************************************************************************
 Copyright (C) 1996 Massachusetts Institute of Technology all rights reserved 
	Programmer: George A. Kriezis
**************************************************************************** */
#include <math.h>
#include <malloc.h>
#include "gen.h"

/* ROUTINES FROM PRAKASH TO TRANSFORM A SURFACE AND PLACE IT IN AN
 * ORIENTATION IN THREE SPACE  */

/****************************************************************************
*                            copymat4()
*****************************************************************************
* 
* 1   Purpose
*     Copy 4x4 matrix.
* 
* 2   Specification
*     #include "bspl.h"
*     void copymat4(double M1[][4], double M2[][4]);
* 
* 3   Description
*     This function copies the contents of a 4times4 matrix.
* 
* 4   References
*     Not applicable.
* 
* 5   Parameters
*       1.double M1[][4]
*         On entry: the matrix to be copied.
*       2.double M2[][4]
*         On exit: the copy of the original matirx.
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
* 9   Functions referenced by copymat4() are: None
* 
* 10  Functions that reference copymat4() are: None
* 
*****************************************************************************/

void copymat4(double M1[][4], double M2[][4])
{
  int i,j;
  
  for(i=0; i<4; i++)
    for(j=0; j<4; j++)
      M2[i][j] = M1[i][j];
}

/***************** conflict with gl functions ******************/
/****************************************************************
* void translate(vector *v, vector *r0)
* {
* double w2,w;
* 
* w = v->w;
* w2 = r0->w;
* 
* v->x = v->x + v->w*r0->x/w2;
* v->y = v->y + v->w*r0->y/w2;
* v->z = v->z + v->w*r0->z/w2;
* }
* 
* void rotate(vector *v, vector *u1, vector *u2, vector *u3)
* {
* double x,y,z,w;
* 
* x = v->x;
* y = v->y;
* z = v->z;
* 
* v->x = u1->x*x + u2->x*y + u3->x*z;
* v->y = u1->y*x + u2->y*y + u3->y*z;
* v->z = u1->z*x + u2->z*y + u3->z*z;
* }
***********************************************************************/

/**************************************************************************
*                                chgsys()
***************************************************************************
* 
* 1   Purpose
*     This routine realizes the transformation of the coordinates of a point.
* 
* 2   Specification
*     #include "bspl.h"
*     void chgsys(double M[][4], vector *r)
* 
* 3   Description
*     This function calculates the coordinates of a point after transformation,
*     by post-multiply the transformation matrix.
*     
* 4   References
*     Not applicable
* 
* 5   Parameters
*       1.double M[][4]
*         On entry: 4X4 transformation matrix.
*       2.vector *r
*         On entry: the homogeneous coordinates of a point.
*         On exit:  the coordinates of the point after transformation.
* 
* 6   Return Values, Error Indicators and Warnings
*     Not applicable
* 
* 7   Accuracy
*     Not applicable
* 
* 8   Further comments
*     Not applicable
* 
* 9   Functions referenced by chgsys() are: None
* 
* 10  Functions that reference chgsys() are: None
* 
******************************************************************************/

void chgsys(double M[][4], vector *r) 
{
  double x,y,z,w;

  x = r->x;
  y = r->y;
  z = r->z;
  w = r->w;

  r->x = M[0][0]*x + M[0][1]*y + M[0][2]*z + M[0][3]*w;
  r->y = M[1][0]*x + M[1][1]*y + M[1][2]*z + M[1][3]*w;
  r->z = M[2][0]*x + M[2][1]*y + M[2][2]*z + M[2][3]*w;
  r->w = M[3][0]*x + M[3][1]*y + M[3][2]*z + M[3][3]*w;
}

/**************************************************************************
*                                chgsyspre()
***************************************************************************
* 
* 1   Purpose
*     This routine realizes the transformation of the coordinates of a point.
* 
* 2   Specification
*     #include "bspl.h"
*     void chgsyspre(double M[][4], vector *r)
* 
* 3   Description
*     This function calculates the coordinates of a point after transformation,
*     by pre-multiply the transformation matrix.
*     
* 
* 4   References
*     Not applicable
* 
* 5   Parameters
*       1.double M[][4]
*         On entry: 4X4 transformation matrix.
*       2.vector *r
*         On entry: the homogeneous coordinates of a point.
*         On exit:  the coordinates of the point after transformation.
* 
* 6   Return Values, Error Indicators and Warnings
*     Not applicable
* 
* 7   Accuracy
*     Not applicable
* 
* 8   Further comments
*     Not applicable
* 
* 9   Functions referenced by chgsyspre() are: None
* 
* 10  Functions that reference chgsyspre() are: None
* 
******************************************************************************/

void chgsyspre(double M[][4], vector *r)
{
double x,y,z,w;

x = r->x;
y = r->y;
z = r->z;
w = r->w;

r->x = M[0][0]*x + M[1][0]*y + M[2][0]*z + M[3][0]*w;
r->y = M[0][1]*x + M[1][1]*y + M[2][1]*z + M[3][1]*w;
r->z = M[0][2]*x + M[1][2]*y + M[2][2]*z + M[3][2]*w;
r->w = M[0][3]*x + M[1][3]*y + M[2][3]*z + M[3][3]*w;
}

/******************************************************************************
*                                 transmat()
*******************************************************************************
* 
* 1   Purpose
*     This routine calculates transformation matrix.
* 
* 2   Specification
*     #include "bspl.h"
*     void transmat(vector *r0, vector *u1, vector *u2, vector *u3,
*                   double M[][4])
* 
* 3   Description
*     This function calculates the matrx of transformation which first
*     translates by a translating vector and then totates to coincide the
*     coordinate axes with three specified directions. 
* 
* 4   References
*     Not applicable
*  
* 5   Parameters
*     1.vector *r0
*       On entry: the translating vector, i.e., the origin will translate by
*                 this vector.
*     2.vector *u1
*       On entry: the direction that x-axis will rotate to.
*     3.vector *u2
*       On entry: the direction that y-axis will rotate to.
*     4.vector *u3
*       On entry: the direction that z-axis will rotate to.
*     5.double M[][4]
*       On exit:  the transformation matrix.
* 
* 6   Return Values, Error Indicators and Warnings
*     Not applicable
* 
* 7   Accuracy
*     Not applicable
* 
* 8   Further Comments
*     Not applicable
* 
* 9   Functions referenced by transmat() are:
*     cross()
*     dot()
*     errormsg()
*     glue_vector()
*     mag()
*     matrixmu()
*     unitvector()
* 
* 10  Functions that reference transmat() are: None
*     
*****************************************************************************/

void transmat( vector *r0, vector *u1, vector *u2, vector *u3, double M[][4])
{
  double Mr[4][4],Mt[4][4];       /* rotation matrix and translation matrix */
  int i,j;                        /* indices */
  struct vector *ux, *uy, *uz, *iu, *ju, *ku;
                                  /* unit vectors */
  double d1, d2, d3, one, zero;   /* d1,d2,d3 -- magnitudes of some vectors */

  /* define 1.0 and 0.0 for convenience. */
  one = 1.0;
  zero = 0.0;

  /* initialize the two matrices */
  for(i=0; i<4; i++)
      for(j=0; j<4; j++)
	  Mt[i][j] = Mr[i][j] = 0.0;

  /* set the global scaling factor in the translation matrix. */
  for(i=0; i<4; i++)
      Mt[i][i] = r0->w;
  /* set the elements which determine the translation. */
  Mt[0][3] = -r0->x;
  Mt[1][3] = -r0->y;
  Mt[2][3] = -r0->z;

  /* calculate the magnitudes of the three rotation vectors. */
  d1 = mag(u1);
  d2 = mag(u2);
  d3 = mag(u3);

  /* set the elements of the rotation matrix. */
  /* if the first one is 0 vector while the 2nd and the 3rd are not. */
  if (d1 == 0.0 && d2 != 0.0 && d3 != 0.0)
     {
      /* calculate the three orthogonal unit vectors. */
      uy = unitvector(u2);
      uz = unitvector(u3);
      /* the first one is cross product of the other two. */
      ux = unitvector(cross(uy, uz));
     }
  /* if the 2nd one is 0 vector while the other two are not. */
  else if(d1 != 0.0 && d2 == 0.0 && d3 != 0.0)
     {
      /* calculate the unit vectors */
      ux = unitvector(u1);
      uz = unitvector(u3);
      /* the 2nd one is the cross product of the other two. */  
      uy = unitvector(cross(uz, ux));
     }
  /* if the 3rd one is 0 vector while the other two are not. */
  else if(d1 != 0.0 && d2 != 0.0 && d3 == 0.0)
     {
      /* calculate the unit vectors. */
      ux = unitvector(u1);
      uy = unitvector(u2);
      /* the 3rd one is the cross product of the other two. */
      uz = unitvector(cross(ux, uy));
     }
  /* if the 1st and the 2nd are 0 vectors but the 3rd are not. */
  else if(d1 == 0.0 && d2 == 0.0 && d3 != 0.0)
     {
      /* calculate the unit vector */
      uz = unitvector(u3);
      /* create three othogonal unit vectors. */
      iu = glue_vector(one,zero,zero);
      ju = glue_vector(zero,one,zero);
      ku = glue_vector(zero,zero,one);
      /* construct the other two unit vectors. */ 
      if ( dot(uz,iu) != 1.0)
         ux = unitvector(cross(uz,iu));
      else if ( dot(uz,ju) != 1.0)
              ux = unitvector(cross(uz,ju));
           else if ( dot(uz,ku) != 1.0)
                   ux = unitvector(cross(uz,ku));
      uy = unitvector(cross(uz, ux));
     }
  /* none of the above cases. */
  else
     {
       errormsg(1,"ctrans() : Cases of directions cosines u1 and u2 being specified are not handled well");
     }

  /* set the rotation matrix */
  Mr[0][0] = ux->x;
  Mr[0][1] = ux->y;
  Mr[0][2] = ux->z;
  Mr[1][0] = uy->x;
  Mr[1][1] = uy->y;
  Mr[1][2] = uy->z;
  Mr[2][0] = uz->x;
  Mr[2][1] = uz->y; 
  Mr[2][2] = uz->z;
  Mr[3][3] = 1.0;

  /* obtain the transformation matrix */
  matrixmu(4,4,4,4,4,4,&Mr[0][0],&Mt[0][0],&M[0][0]);

  /* free memory */
  free((char *)ux);
  free((char *)uy);
  free((char *)uz);
}

/******************************************************************************
                                pretransmat()
*******************************************************************************
*
* 1   Purpose
*     This routine calculates transformation matrix.
* 
* 2   Specification
*     #include "bspl.h"
*     void pretransmat(vector *r0, vector *u1, vector *u2, vector *u3,
*                      double M[][4])
* 
* 3   Description
*     This function calculates the matrx of transformation which first rotates 
*     to coincide the coordinate axes with three specified directions and
*     then translates by a translating vector.
* 
* 4   References
*     Not applicable
*  
* 5   Parameters
*     1.vector *r0
*       On entry: the translating vector, i.e., the origin will translate by
*                 this vector.
*     2.vector *u1
*       On entry: the direction that x-axis will rotate to.
*     3.vector *u2
*       On entry: the direction that y-axis will rotate to.
*     4.vector *u3
*       On entry: the direction that z-axis will rotate to.
*     5.double M[][4]
*       On exit:  the transformation matrix.
* 
* 6   Return Values, Error Indicators and Warnings
*     Not applicable
* 
* 7   Accuracy
*     Not applicable
* 
* 8   Further Comments
*     Not applicable
* 
* 9   Functions referenced by pretransmat() are:
*     cross()
*     dot()
*     errormsg()
*     glue_vector()
*     mag()
*     matrixmu()
*     unitvector()
* 
* 10  Functions that reference pretransmat() are: None
*     
******************************************************************************/

void pretransmat(vector *r0, vector *u1, vector *u2, vector *u3,double M[][4])
{
  double Mr[4][4],Mt[4][4];     /* rotation and translation matrix,
				 * respectively. */
  int i,j;                      /* indices */
  struct vector *ux, *uy, *uz, *iu, *ju, *ku;    /* unit vectors */
  double d1, d2, d3, one, zero; /* d1,d2,d3 -- magnitudes of some vectors */

  /* define 1.0 and 0.0 for convenience */
  one = 1.0;
  zero = 0.0;

  /* initialize the two matrices as 0 */
  for (i=0; i<4; i++)
      for (j=0; j<4; j++)
          Mt[i][j] = Mr[i][j] = 0.0;

  /* set the global scaling factor in the translation matrix. */
  for (i=0; i<4; i++)
      Mt[i][i] = r0->w;
  /* set the elements which determine the translation. */
  Mt[3][0] = -r0->x;
  Mt[3][1] = -r0->y;
  Mt[3][2] = -r0->z;

  /* calculate the magnitudes of the three rotation vectors. */
  d1 = mag(u1);
  d2 = mag(u2);
  d3 = mag(u3);

  /* set the rotation matrix. */
  /* if only the first direction is undefined. */
  if (d1 == 0.0 && d2 != 0.0 && d3 != 0.0)
     {
      /* calculte the unit vectors */
      uy = unitvector(u2);
      uz = unitvector(u3);
      /* define the other direction as the cross product */
      ux = unitvector(cross(uy, uz));
     }
  /* if only the 2nd direction is undefined. */
  else if(d1 != 0.0 && d2 == 0.0 && d3 != 0.0)
     {
      /* calculate the unit vectors (direction cosines ) */
      ux = unitvector(u1);
      uz = unitvector(u3);
      /* define the 2nd direction */  
      uy = unitvector(cross(uz, ux));
     }
  /* if only the 3rd direction is undefined. */
  else if(d1 != 0.0 && d2 != 0.0 && d3 == 0.0)
     {
      /* calculate the unit vectors (direction cosines ) */
      ux = unitvector(u1);
      uy = unitvector(u2);
      /* define the 3rd direction */
      uz = unitvector(cross(ux, uy));
     }
  /* if only the 3rd direction is defined. */
  else if(d1 == 0.0 && d2 == 0.0 && d3 != 0.0)
     {
      /* calculate the unit vector (direction cosine ) of the 3rd direction. */
      uz = unitvector(u3);
      /* create three unit vectors. */
      iu = glue_vector(one,zero,zero);
      ju = glue_vector(zero,one,zero);
      ku = glue_vector(zero,zero,one);
      /* construct the other two direction cosines. */
      if ( FABS(dot(uz,iu)) != 1.0)
         ux = unitvector(cross(uz,iu));
      else if ( FABS(dot(uz,ju)) != 1.0)
              ux = unitvector(cross(uz,ju));
              else if ( FABS(dot(uz,ku)) != 1.0)
                      ux = unitvector(cross(uz,ku));
      uy = unitvector(cross(uz, ux));
     }
  /* if none of the above cases. */
  else
     {
       errormsg(1,"ctrans() : Cases of directions cosines u1 and u2 being specified are not handled well");
     }

  /* set the rotation matrix */
  Mr[0][0] = ux->x;
  Mr[1][0] = ux->y;
  Mr[2][0] = ux->z;
  Mr[0][1] = uy->x;
  Mr[1][1] = uy->y;
  Mr[2][1] = uy->z;
  Mr[0][2] = uz->x;
  Mr[1][2] = uz->y;
  Mr[2][2] = uz->z;
  Mr[3][3] = 1.0;

  /* finally construct the transformtion matrix. */
  matrixmu(4,4,4,4,4,4,&Mt[0][0],&Mr[0][0],&M[0][0]);

  /* free memory */
  free((char *)ux);
  free((char *)uy);
  free((char *)uz);
}

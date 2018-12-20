/************************************************************************
 *									*
			    Copyright (C) 1990
     Massachusetts Institute of Technology, Cambridge, Massachusetts.
			   All rights reserved.

     This portion of source code was prepared in the Ocean Engineering
     Design Laboratory by Bradley A. Moran.
 *									*
 ************************************************************************/

#include <math.h>
#include <stdio.h>
#include "gen.h"

/* Function: mult4x4()
 * Purpose: Multiple a vector by a 4x4 matrix
 * Arguments:
 *  m - the 4x4 matrix
 *  in - the input vector
 *  out - the output vector
 */
/* Functions referenced by mult4x4() are:
 *  copyvector()
 */
/* Functions that reference mult4x4() are:
 *  con_confun()
 *  cyl_frenet_tr()
 *  cyl_generatrix()
 *  generatrix()
 *  localize_diagnostic()
 *  localize_sumsq()
 *  localize_sumsq_opt()
 *  rotate_section()
 *  transformegeom()
 *  transformfgeom()
 *  TransformGridSurf()
 *  TransformIgesArray()
 *  TransformIgesVect()
 *  TransformListCurv()
 *  TransformPts()
 */

void mult4x4 (double **m, vector *in, vector *out)
{
     vector tmp;
     tmp.x = m[0][0]*in->x + m[0][1]*in->y + m[0][2]*in->z + m[0][3]*in->w;
     tmp.y = m[1][0]*in->x + m[1][1]*in->y + m[1][2]*in->z + m[1][3]*in->w;
     tmp.z = m[2][0]*in->x + m[2][1]*in->y + m[2][2]*in->z + m[2][3]*in->w;
     tmp.w = m[3][0]*in->x + m[3][1]*in->y + m[3][2]*in->z + m[3][3]*in->w;
     copyvector (&tmp, out);
}

/* Function: matrix_mult:
 * Purpose: Multiply two matrices
 * Arguments:
 *  input1 - address of first matrix
 *  input2 - address of second matrix
 *  output - address of new matrix
 *  d1 - number of columns of first matrix
 *  d2 - number of rows of second matrix
 *  d3 - number of columns of new matrix
 */

void matrix_mult(float input1[],float input2[],float output[],
		 int d1,int d2,int d3)
{ 
  int i,j,k;
  for(i = 0; i < d1; i++)
  for(j = 0; j < d3; j++) {
    output[i*d3+j] = 0.0;
    for(k = 0; k < d2; k++)
      output[i*d3+j] += *(input1+i*d2+k)**(input2+k*d3+j);
  }
}

/* Function: matrix_transpose()
 * Purpose: Transpose a matrix
 * Method: Transpose rows and columns, mat[i][j] becomes mat[j][i]
 * Arguments:
 *  m - number of rows
 *  n - number of columns
 *  mat - the address of the matrix
 *  d1 - number of columns of the original matrix
 *  mattransp - the address of the transposed matrix
 *  d2 - number of columns of the transpose matrix
 *       this is equal to the number of rows of the original matrix
 */
/* Functions that reference matrix_transpose() are:
 *  berntomono()
 *  monotobern()
 */

void matrix_transpose(int m,int n,double2 mat[],int d1,
		      double2 mattransp[],int d2)

{   
  int i,j;
  for(i = 0; i < m; i++)
  for(j = 0; j < n; j++)
    mattransp[j*d2+i] = mat[i*d1+j];
}

/* Function: Print_matrix()
 * Purpose: Print a matrix row by row
 * Method: The contents of the m x n matrix are printed to the standard
 *         output unit row by row
 * Arguments:
 *  m - number of rows
 *  n - number of columns
 *  mat - address of the matrix
 *  d1 - the number of columns
 */

void Print_matrix(int m,int n,double2 mat[],int d1)
{  
  int i,j;
  for(i = 0; i < m; i++) {
    for(j = 0; j < n; j++)
      printf(" %lf ",*(mat+i*d1+j));
    printf("\n");
  }
}

/* Function: Read_matrix()
 * Purpose: Read matrix from a file
 * Arguments:
 *  file - file pointer of open file containing the matrix
 *  mat - the address of the matrix
 *  row - number of rows
 *  col - number of columns
 */

void Read_matrix(FILE *file,double2 *mat,int row,int col)
{
  int i,j;
  for(i = 0; i < row; i++)
  for(j = 0; j < col; j++)
    fscanf(file, "%lf", mat + i*col + j);
}

/* Function: matrixmu()
 */
/* Functions that reference matrixmu() are:
 *  berntomono()
 *  curve_oslo1()
 *  curve_oslo2()
 *  monotobern()
 *  pretransmat()
 *  transmat()
 */

void matrixmu(int m,int n,int p, int d1,int d2,int d3, 
	      double2 in1[], double2 in2[], double2 out[])
/* Here are the old declarations for the routine; I don't pretend to 
   understand them */
/* int m,n,p;      		 used members of the three matrices */
/* int d1,d2,d3;   		 column dimension of the three matrices */
/* double2 in1[],in2[],out[];    out = in1*in2  */
{
  int i,j,k;

  for(i = 0; i < m; i++) {
    for(j = 0; j < p; j++)
      out[d3*i+j] = 0.0;
    for(k = 0; k < n; k++)
      for(j = 0; j < p; j++)
	out[d3*i+j] += in1[d1*i+k]*in2[d2*k+j];
  }
}

/* ***************************************************************************
 Copyright (C) 1996 Massachusetts Institute of Technology all rights reserved 
	Programmer: George A. Kriezis
**************************************************************************** */
# include <math.h>
# include "gen.h"
# include "bspl.h"

/*****************************************************************************
*                                unionsrt()
******************************************************************************
* 
* 1   Purpose
*     This routine unions two sorted arrays into one sorted array.
* 
* 2   Specification
*     #include "bspl.h"
*     void unionsrt(double s1[], int ns1, double s2[], int ns2, double s3[],
* 	          int *ns3, double eps)
* 
* 3   Description
*     This routine takes two arrays sorted in non-descending order as input,
*     and combines them into one array in non-descending order. If two elements
*     from different arrays have same value, the value only occurs once in the
*     resulting array.
*     It assumes that the first and the last elements in two array are equal to
*     each other. 
* 
* 4   References
*     Not applicable
* 
* 5   Parameters
*       1.double s1[]
*         On entry: the first one dimensional array to be unioned.
*       2.int ns1
*         On entry: the length of the first array.
*       3.double s2[]
*         On entry: the second one dimensional array to be unioned.
*       4.int ns2
*         On entry: the length of the second array.
*       5.double s3[]
*         On exit:  the array which is a union of other two arrays.
*       6.double ns3
*         On exit:  the length of the union array.
*       7.double eps
*         On entry: the tolerance
* 
* 6   Return Values, Error Indicators and Warnings
*     The comparison whether two numbers are identical is based on a pre-
*     specified tolerance. If the difference of these two numbers is within,
*     the tolerance, they are considered as identical.
* 
* 7   Accuracy
* 
* 8   Further Comments
* 
* 9   Functions referenced by bernrs() are: None
* 
* 10  Functions that reference bernrs() are: None
* 
****************************************************************************/

void unionsrt(double s1[], int ns1, double s2[], int ns2, double s3[],
	      int *ns3, double eps)
{
  int i,j,k,flag;  /* indices and flag */

  #define MORE 1
  #define ENDIT 2

  flag = MORE;    /* flag indicates if the union will go on. */

  j = 0;
  k = 0;

  for (i=0; i<ns1; i++)
      { /* starts from s1[] */
       s3[k++] = s1[i];
       /* if the two elements are equal, skip to the next */
       if (flag == MORE){
          if ( FABS(s2[j]-s1[i]) < eps )
	     j++;
          }
       /* if s2[] reaches its end, then stop */
       if (j >= ns2)
	  flag = ENDIT;
       /* find the place for an element from s2[] */
       while (flag == MORE && (s1[i] < s2[j]) && (s2[j] < s1[i+1]))
             {
	      s3[k++] = s2[j];
	      j++;
	      /* if s2[] reaches its end, then stop */
	      if ( j >= ns2 )  
	 	 flag = ENDIT;
             }
      }

  *ns3 = k;
}

/*****************************************************************************
*                                 insertval()
******************************************************************************
* 
* 1   Purpose
*     This routine inserts a number into a sorted array.
* 
* 2   Specification
*     #include "bspl.h"
*     void insertval(double array[], int n, int ninsert, double val,
* 	           double newarray[])
* 
* 3   Description
*     This routine inserts a given value into a non-descending array such that 
*     the updated array remains non-descending.
* 
* 4   References
*     Not applicable
* 
* 5   Parameters
*       1.double array[]
*         On entry: one dimensional array sorted in non-descending order
*       2.int n
*         On entry: the length of the array
*       3.int ninsert
*         On entry: the number of times the value is to be inserted.
*       4.double val
*         On entry: the value to be inserted.
*       5.double newarray[]
*         On exit:  the updated array
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
* 9   Functions referenced by bernrs() are: None
* 
* 10  Functions that reference bernrs() are: None
* 
******************************************************************************/

void insertval(double array[], int n, int ninsert, double val,
	       double newarray[])
{
  int i=0,j;

  /* find the position to insert the value */
  for(i=0; i<n; i++)
     {
      newarray[i] = array[i];
      if( array[i] <= val && val < array[i+1] ) 
	break;
     }

  for(j=n-1; j>=i+1; j--)
	newarray[j+ninsert] = array[j];

  /* insert the value ninsert times */
  for(j=i+1; j<=(i+ninsert); j++)
      newarray[j] = val;

}

/****************************************************************************
*                                member()
*****************************************************************************
*
* 1   Purpose
*     This routine checks the multiplicity of a value in an array.
* 
* 2   Specification
*     #include "bspl.h"
*     int member(double array[], int n, double val, double eps)
* 
* 3   Description
*     If the number is an element of the array, the function returns its 
*     multiplicity in the array; if not, the multiplicity is 0.
* 
* 4   References
*     Not applicable
* 
* 5   Parameters
*       1.double array[]
*         On entry: one dimensional array
*       2.int n
*         On entry: the length of the array
*       3.double val
*         On entry: the number which is to be checked.
*       4.double eps
*         On entry: the tolerance
* 
* 6   Return Values, Error Indicators and Warnings
*     The function returns the number of occurence (multiplicity) of the value in 
*     the array. If the value is close to an element in the array within the 
*     tolerance, the multiplicity increases by 1.
* 
* 7   Accuracy
*     Not applicable
* 
* 8   Further Commments
*     Not applicable
* 
* 9   Functions referenced by bernrs() are: None
* 
* 10  Functions that reference bernrs() are: None
* 
******************************************************************************/

int member(double array[], int n, double val, double eps)
{
  int i,li,hi;   /* li,hi -- never used (glshen) */
  int multiplicity;

  multiplicity = 0;

  for(i=0; i<n; i++)
	{	
	if(FABS(array[i] - val) < eps)
		multiplicity++;	
	}

  return(multiplicity);
}

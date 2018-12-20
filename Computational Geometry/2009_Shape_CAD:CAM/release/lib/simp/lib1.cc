
#include <iostream>
#include "lib1.h"
#include "compiler.h"
#include "simpoly.h"
#include <cmath>
using namespace std;

ostream& operator<< (ostream& output, real_array& a)
{

     output << "\n";
     for (int i=0; i < a.n ; i++) {
	  output << " p[" ;
	  output << i ;
	  output << "] = " ;
	  output << a[i];
     }

     output << endl;

     return output;

}



int if_overlap(real_array& ra, real_array& rb)
{

     if (ra.n != rb.n) 
	  cout << "warning: wrong dimensional comparison in if_overlap()"
	       << endl;
     for (int i=0; i< ra.n ; i++) 
	  if ( if_overlap(ra[i], rb[i]) == FALSE) {
	 
	       return FALSE;
	  }

     return TRUE;
     
}

int if_overlap(double a, double b)
{
     if ( fabs(a - b) < 1.0e-10 )
	  return TRUE;
     else
	  return FALSE ;
}


int operator == (dbl_array &a, dbl_array &b) {
     if (a.n != b.n) return FALSE;
     int result = TRUE;
     for (int i= 0; i < a.n; i++) {
	  result = result && (a[i] == b[i]);
     }
     return result;
}


void dbl_array::init (int siz)           /* Initializer */
{
     if (siz == 0) p = NULL;
     else
	  p = alloc (siz);
     n = siz;
}

dbl_array::dbl_array (int siz)           /* Create a double array */
{ init (siz); }

dbl_array::dbl_array (const dbl_array &d)   /* Copy constructor */
{ init (d.size()); memcpy (p, (double *) d.p, n*sizeof(double)); }

double& dbl_array::operator[] (int i) const   /* Indexer */
    { return (((double *) p)[i]); }

void sht_array::init (int siz)           /* Initializer */
{
     if (siz == 0) p = NULL;
     else
	  p = alloc (siz);
     n = siz;
}

void int_array::init (int siz)           /* Initializer */
{     
     if (siz == 0) p = NULL;
     else
	  p = alloc (siz);
     n = siz;
}



void double_array2::init (int r, int c)      /* Initializer */
{ p = (double **) new double*[r];
  for (int i=0; i < r; i++)
    p[i] = new double[c];
  m = r; n = c;          /* Set m and n */
} 

void double_array2::init (int r, int c, int z) /* Initializer for array of 0's */
    { 

      p = (double **) new double*[r];
      for (int i=0; i < r; i++)
	p[i] = new double[c];
      int trash = z; trash++; /* not useful */
      m = r; n = c;            /* Set m and n */
    } 



double* & double_array2::operator[] (int i) const
{ return (((double**) p)[i]); }


#ifdef USE_INTERVAL

int operator == (Interval_array &a, Interval_array &b) {
     if (a.n != b.n) return FALSE;
     int result = TRUE;
     for (int i= 0; i < a.n; i++) {
	  result = result && (a[i] == b[i]);
     }
     return result;
}

Interval_array::Interval_array (int siz)           /* Create a Interval array */
{ init (siz); }

void Interval_array::init (int siz)           /* Initializer */
{ 
     if (siz == 0) p = NULL;
     else
	  p = alloc (siz);
     n = siz; 
}



/////////////////////////////////////////////////////////////////////



void Interval_array2::init (int r, int c)      /* Initializer */
{ 
  p = (Interval **) new Interval*[r];       /* Allocate row pointers */
  for (int i=0; i < r; i++)
    p[i] = new Interval[c];
  m = r; n = c;          /* Set m and n */
}

void Interval_array2::init (int r, int c, int z)  
{   
  int trash; trash = z;

  p = (Interval **) new Interval*[r];       /* Allocate row pointers */
  for (int i=0; i < r; i++)
    p[i] = new Interval[c];

  m = r; n = c;            /* Set m and n */
}


Interval* & Interval_array2::operator[] (int i) const
{ return (((Interval**) p)[i]); }

#endif

#ifdef USE_RAT
////////////////////////////////////////////////////////////////////////
void Interval_array2::init (int r, int c)      /* Initializer */
{ 
  p = (Rat **) new Rat*[r];       /* Allocate row pointers */
  for (int i=0; i < r; i++)
    p[i] = new Rat[c];
  m = r; n = c;          /* Set m and n */
}


#endif

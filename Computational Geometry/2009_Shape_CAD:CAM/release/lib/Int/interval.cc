// Copyright (C) Massachusetts Institute of Technology, 1997
// All rights reserved

/*****************************************************************************
   interval.cc                        last edit 10.3.97
 ****************************************************************************/

/* 10/3/97 - add new version of round(): correct and fast (SLA) */

extern "C" {
#include <stdlib.h>
#include <math.h>
}

#include "interval.h"
#define FLOATING  0
#define INTERVAL  1

int INPUT_STYLE = INTERVAL;
using namespace std;

typedef union du_sh {
  double d;
  unsigned short s[4];
} Double;

#define index_contain_exponent 3 // 0 for some UNIX platform
                                 // 3 for some PC
#define MSW index_contain_exponent

du_sh RU_DB;
du_sh ROUND;

/*** Takashi's original version; correct but slow **********************/
/*** use the latest version of round() instead *************************/
double ulp(double d)
{
 int expt;
 frexp(d, &expt);
 return ldexp(0.5, expt-52); // extract the ulp
}
/***********************************************************************/

double Interval:: max(double a, double b, double c, double d)
{
  double result = a;
  if(b > result) result = b;
  if(c > result) result = c;
  if(d > result) result = d;
  return result;
}

double Interval:: min(double a, double b, double c, double d)
{
  double result = a;
  if(b < result) result = b;
  if(c < result) result = c;
  if(d < result) result = d;
  return result;
}

Interval pow (Interval a, double n)
{
     Interval b = 1.0 ;
     if (n== 0.0) return b;

     
     for (int i=0; i< (int) n; i++)
     {
       
	  b = (b * a);
     }
     return b;
}

Interval pow (Interval a, int n)
{
     Interval b = 1.0;
     if (n== 0 ) return b;
     
     
     for (int i=0; i< n ; i++)
     {
	  b = b * a;
     }
     return b;
}

double pow (double a, int n)
{
     double b = 1.0;
     if (n== 0 ) return b;
     
     
     for (int i=0; i< n ; i++)
     {
	  b = b * a;
     }
     return b;
}

/*********************************************************************/
/*** return ulp of double precision number x (new version) ***********/
/*** correct and fast for normalized and denormalized ****************/

static unsigned short mask[16] = { /* bit masks for bits 0 - 15 */
  0x0001, 0x0002, 0x0004, 0x0008, 0x0010, 0x0020, 0x0040, 0x0080,
  0x0100, 0x0200, 0x0400, 0x0800, 0x1000, 0x2000, 0x4000, 0x8000};

double round(double x)
{
  Double U,   /* ulp of x */
    X;        /* working copy of x */
  int bit,    /* position of bit e-1 in 16-bit word */
    e1,       /* biased exponent - 1 */
    word;     /* index of 16-bit word containing bit e-1 */

  X.d = x;
  X.s[MSW] &= 0x7ff0;              /* isolate exponent in 16-bit word */

  /* X.s[0] now holds the exponent in bits 14-4 */

  U.d = 0.0;                       /* initialize exponent and mantissa to 0 */

  if (X.s[MSW] > 0x0340)           /* ulp is normalized number */
    U.s[MSW] = X.s[MSW]-0x0340;    /*  set exponent to e-52 */

  /* the value 0x0340 is 52 left-shifted 4 bits, i.e. 0x0340 = 832 = 52<<4 */

  else {                           /* ulp is denormalized number */
    e1 = (X.s[MSW]>>4) - 1;        /*  biased exponent - 1 */
    word = e1>>4;                  /*  find 16-bit word containing bit e-1 */
    if (MSW == 0) word = 3 - word; /* compensate for word ordering */
    bit  = e1%16;                  /*  find the bit position in this word */
    
    U.s[word] |= mask[bit];        /*  set the bit to 1 */
  }

  return U.d;                      /* return ulp */
}

/*************************************************************************/
/*** Hu's latest version; correct but slow for denormalized **************/
double round3(double d)
{
     RU_DB.d = d;

     RU_DB.s[index_contain_exponent] = 
     RU_DB.s[index_contain_exponent] & 0x7ff0;   
     //<--- extract the exponent part 

     if (RU_DB.s[index_contain_exponent] <= 0x340)
         return ulp(d);

     ROUND.d = 0.0;
	  
     RU_DB.s[index_contain_exponent] = RU_DB.s[index_contain_exponent] 
                                       - 0x0340;  // expt - 52 
     
     ROUND.s[index_contain_exponent] = RU_DB.s[index_contain_exponent];
     
     return ROUND.d ;
}

/***************************************************************************/
/*** Hu's original version; fast but incorrect for denormalized ************/
double round2( double  d)
{
     ROUND.d = 0.0;
     RU_DB.d = d;

     RU_DB.s[0] = RU_DB.s[0] & 0x7ff0;   // extract the exponent part 

      if ( !RU_DB.s[0]) {
	ROUND.s[0] = 0x3ca0;
	return ROUND.d;
      }
	  
     RU_DB.s[0] = RU_DB.s[0] - 0x0340;  // expt - 52 
     
     ROUND.s[0] = RU_DB.s[0];
     
     return ROUND.d ;
}

/***************************************************************************/
/*** Takashi's original version; correct but slow **************************/
double round1 (double d)
{
  double man, round;
  int expt;

  man = frexp(d, &expt);
  round = ldexp(0.5, expt-52);
  return round;
}

Interval merge(Interval a, Interval b)
  // Merge two Intervals
{
  double low;
  double upp;

  if (a.low < b.low)
    low = a.low;
  else
    low = b.low;

  if (a.upp > b.upp)
    upp = a.upp;
  else
    upp = b.upp;
  
  Interval inv; inv.low = low; inv.upp = upp; return(inv);
}

int if_overlap(Interval a, Interval b)
{
     if (a.low > b.upp) return FALSE;
     else if (a.upp < b.low) return FALSE;
     else return TRUE;

}


  
ostream& operator<<( ostream& os, Interval& inv)
{
   os << "[" << inv.low << ", " << inv.upp << "]";

//  os << "{" <<  ((inv.low + inv.upp) / 2.0) << ",}";
   return os;
}






istream& operator>>(istream& is, Interval& inv)
{
   double d;

   is >> d;

   inv.low = d;

   if (INPUT_STYLE == INTERVAL)
   is >> d;
   inv.upp = d;
 
   return is;
}


  // --------------------
  // Arithmetic Operators
  // --------------------

Interval operator + (Interval a, Interval b)
  {
    Interval c; add (a, b, c); return c;
  }

Interval operator - (Interval a, Interval b)
  {
    Interval c; sub (a, b, c); return c;
  }

Interval operator * (Interval a, Interval b)
  {
    Interval c; mul (a, b, c); return c;
  }

Interval operator / (Interval a, Interval b)
  {
    Interval c; div (a, b, c); return c;
  }



Interval operator + (Interval a, double d)
  {
    Interval b(d,d), c; add (a, b, c); return c;
  }

Interval operator - (Interval a, double d)
  {
    Interval b(d,d), c; sub (a, b, c); return c;
  }

  
Interval operator * (Interval a, double d)
  {
    Interval b(d,d), c; mul (a, b, c); return c;
  }

Interval operator / (Interval a, double d)
  {
    Interval b(d,d), c; div (a, b, c); return c;
  }

Interval operator + (double d, Interval b)
  {
    Interval a(d, d), c; add (a, b, c); return c;
  }

Interval operator - (double d, Interval b)
  {
    Interval a(d, d), c; sub (a, b, c); return c;
  }

  
Interval operator * (double d, Interval b)
  {
    Interval a(d, d), c; mul (a, b, c); return c;
  }

Interval operator / (double d, Interval b)
  {
    Interval a(d, d), c; div (a, b, c); return c;
  }

  
  //---------------------
  // Comparison Operators
  //---------------------

int operator > (Interval a, Interval b)
  {
    if (a.low > b.upp) return TRUE;
    else if(a.upp <= b.low) return FALSE;
    else return MAYBE;

//    else if (a.upp > b.upp) return TRUE;
//    else return FALSE;
  
  }


int comp_greater(Interval a, Interval b) // conservatively compare
  {
    if (a.get_low() > b.get_upp()) return TRUE;
    else if(a.get_upp() <= b.get_low()) return FALSE;
//    else return MAYBE;

    else if (a.get_upp() > b.get_upp()) return TRUE;
   else return FALSE;
  
  }

int operator < (Interval a, Interval b)
  {
    if (a.upp < b.low) return TRUE;
    else if(a.low >= b.upp) return FALSE;
    else return MAYBE;
//    else if (a.low < b.low) return TRUE;
//    else return FALSE;
  }


int comp_less(Interval a, Interval b) // conservatively compare
  {
    if (a.get_upp() < b.get_low()) return TRUE;
    else if(a.get_low() >= b.get_upp()) return FALSE;
//    else return MAYBE;
    else if (a.get_low() < b.get_low()) return TRUE;
    else return FALSE;
  }

int operator >= (Interval a, Interval b)
  {
    if (a.low >= b.upp) return TRUE;
    else if(a.upp <= b.low) return FALSE;
    else return MAYBE;
  }

int operator <= (Interval a, Interval b)
  {
    if (a.upp <= b.low) return TRUE;
    else if(a.low >= b.upp) return FALSE;
    else return MAYBE;
  }


int operator == (Interval a, Interval b)
  {
    return (a.low == b.low && a.upp == b.upp);
  }

int operator != (Interval a, Interval b)
  {
    return (a.low != b.low || a.upp != b.upp);
  }

int operator > (Interval a, double b)
  {
    if (a.low > b) return TRUE;
    else if(a.upp <= b) return FALSE;
    else return MAYBE;
  }

int operator < (Interval a, double b)
  {
    if (a.upp < b) return TRUE;
    else if(a.low >= b) return FALSE;
    else return MAYBE;
  }

int operator >= (Interval a, double b)
  {
    if (a.low >= b) return TRUE;
    else if(a.upp <= b) return FALSE;
    else return MAYBE;
  }

int operator <= (Interval a, double b)
  {
    if (a.upp <= b) return TRUE;
    else if(a.low >= b) return FALSE;
    else return MAYBE;
  }


int operator == (Interval a, double b)
  {
    return (a.low == b && a.upp == b);
  }

    int operator != (Interval a, double b)
  {
    return (a.low != b || a.upp != b);
  }

    int operator > (double a, Interval b)
  {
    if (a > b.upp) return TRUE;
    else if(a <= b.low) return FALSE;
    else return MAYBE;
  }

    int operator < (double a, Interval b)
  {
       cout << " in < double, Interval " << endl;
    if (a < b.low) {cout << "out < d , I ,true" << endl;return TRUE;}
    else if(a >= b.upp) {cout << "out < d , I ,false" << endl;return FALSE;}
    else { cout << "out < d , I ,maybe" << endl;return MAYBE;}
  }

    int operator >= (double a, Interval b)
  {
    if (a >= b.upp) return TRUE;
    else if(a <= b.low) return FALSE;
    else return MAYBE;
  }

    int operator <= (double a, Interval b)
  {
    if (a <= b.low) return TRUE;
    else if(a >= b.upp) return FALSE;
    else return MAYBE;
  }


    int operator == (double a, Interval b)
  {
    return (a == b.low && a == b.upp);
  }

    int operator != (double a, Interval b)
  {
    return (a != b.low || a != b.upp);
  }


 Interval Interval::lower()
  {
    Interval inv; inv.low = low; inv.upp = low; return(inv);
  }

 Interval Interval::upper()
  {
    Interval inv; inv.low = upp; inv.upp = upp; return(inv);
  }

   Interval Interval::operator += (Interval inv)
  {
    add (inv, *this, *this);  return *this;
  }

   Interval Interval::operator -= (Interval inv)
  {
    sub (*this, inv, *this);  return *this;
  }

   Interval Interval::operator *= (Interval inv)
  {
    mul (*this, inv, *this);  return *this;
  }

   Interval Interval::operator /= (Interval inv)
  {
    div (*this, inv, *this);  return *this;
  }

   Interval Interval::operator += (double d)
  {
    Interval inv(d, d); add (inv, *this, *this);  return *this;
  }

   Interval Interval::operator -= (double d)
  {
    Interval inv(d, d); sub (*this, inv, *this);  return *this;
  }

   Interval Interval::operator *= (double d)
  {
    Interval inv(d, d); mul (*this, inv, *this);  return *this;
  }

   Interval Interval::operator /= (double d)
  {
    Interval inv(d, d); div (*this, inv, *this);  return *this;
  }


Interval::Interval (const Interval& inv)
  {
    low = inv.low;   upp = inv.upp;
  }



Interval::~Interval (){}

Interval::Interval()
{ 
     low = 0.; upp=0.;
}

//Interval Interval::operator = (Interval inv)
//{
//     low = inv.low; upp = inv.upp; 
//     return *this;
//}

// assignment operator overloading
Interval& Interval::operator = (const Interval& inv)
{
     low = inv.low; upp = inv.upp; 
     return *this;
}



// defined in interval.h
//double Interval::center()
//{ 
//     return (low + upp)/2.;
//}

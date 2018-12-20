// Copyright (C) 1995 Massachusetts Institute of Technology, Cambridge, MA
// All rights reserved
/****************************************************************************
  interval.h                last edit 10.8.97
 ***************************************************************************/

#ifndef INTERVAL_H
#define INTERVAL_H

#include <iostream>
#include <cmath>
#include <fpu_control.h>

#ifndef FALSE
#define TRUE 1
#define FALSE 0
#endif

#ifndef MAYBE
#define MAYBE -1
#endif



class Interval{

 private:
     double low;
     double upp;
 public:
     Interval();
     Interval (double);
     Interval (double, double);
     Interval (const Interval&);
     ~Interval ();
     
     double max(double, double, double, double);
     double min(double, double, double, double);
     double center();
     double range();
     double magnitude();
     Interval lower();
     Interval upper();
     double get_low() {return low;}
     double get_upp() { return upp;}

  // Assignment Operators
 
//     Interval operator = (Interval);
     Interval operator += (Interval);
     Interval operator -= (Interval);
     Interval operator *= (Interval);
     Interval operator /= (Interval);   

     Interval operator = (double );
     Interval operator += (double);
     Interval operator -= (double);
     Interval operator *= (double);
     Interval operator /= (double);   

  // assignment operator overloading
     Interval& operator=( const Interval& );


 // Arithmetic Operators
	  
     friend Interval operator + (Interval, Interval);
     friend Interval operator - (Interval, Interval);
     friend Interval operator * (Interval, Interval);
     friend Interval operator / (Interval, Interval);
     
     friend Interval operator + (Interval, double);
     friend Interval operator - (Interval, double);
     friend Interval operator * (Interval, double);
     friend Interval operator / (Interval, double); 

     friend Interval operator + (double, Interval);
     friend Interval operator - (double, Interval);
     friend Interval operator * (double, Interval);
     friend Interval operator / (double, Interval);

  // Comparison Operators

     friend int operator > (Interval, Interval);
     friend int operator < (Interval, Interval);
     friend int operator >= (Interval, Interval);
     friend int operator <= (Interval, Interval);
     friend int operator == (Interval, Interval);
     friend int operator != (Interval, Interval);

     friend int operator > (Interval, double);
     friend int operator < (Interval, double);
     friend int operator >= (Interval, double);
     friend int operator <= (Interval, double);
     friend int operator == (Interval, double);
     friend int operator != (Interval, double);

     friend int operator > (double, Interval);
     friend int operator < (double, Interval);
     friend int operator >= (double, Interval);
     friend int operator <= (double, Interval);
     friend int operator == (double, Interval);
     friend int operator != (double, Interval);

   // Other useful functions  

     friend void add (Interval, Interval, Interval&); // Addition
     friend void sub (Interval, Interval, Interval&); // Subtraction
     friend void mul (Interval, Interval, Interval&); // Multiplication
     friend void div (Interval, Interval, Interval&); // Division
     friend Interval pow (Interval, double);           // Raise to a power
     friend Interval pow (Interval, int);           // Raise to a power
     friend Interval sqrt(Interval);                   // square root
     friend double round (double);                      // Round off error
     friend double round1 (double);                     // Takashi
     friend double round2 (double);                     // Hu (incorrect)
     friend double round3 (double);                     // Hu
     friend Interval merge(Interval, Interval);
     friend int if_overlap(Interval, Interval);
     friend std::istream &operator >> (std::istream&, Interval&);
     friend std::ostream &operator << (std::ostream&, Interval&);
     friend int comp_greater(Interval, Interval); // conservatively compare >
     friend int comp_less(Interval, Interval); // conservatively compare <
};





  inline Interval::Interval (double a)
  //  Constructor from a double
  {
    low = a;     upp = a;
  }

  inline Interval::Interval(double a, double b)
  // Constructor from two doubles
  {
    low = a; upp = b;
  }

  inline double Interval::center()
  { 
    return (low + upp)/2.;
  }

  inline double Interval::range()
  {
    return(fabs(upp - low));
  }

  inline double Interval::magnitude()
  {
    return(fabs(upp) > fabs(low) ? fabs(upp) : fabs(low));
  }

  //--------------------
  // Assigment Operators
  //--------------------

  inline Interval Interval::operator = (double d)
  {
    low = d; upp = d; 
    return *this;
  }

inline  void add (Interval a, Interval b, Interval& c)
{

  short save = 0x037f;
  short  cwlow = 0x077f;
  short  cwhigh = 0x0b7f;

  _FPU_SETCW(cwlow);
  c.low = a.low + b.low;

  _FPU_SETCW(cwhigh);
  c.upp = a.upp + b.upp;

  _FPU_SETCW(save);
}

inline void sub (Interval a, Interval b, Interval& c)
{
  
  short save = 0x037f;
  short  cwlow = 0x077f;
  short  cwhigh = 0x0b7f;

  _FPU_SETCW(cwlow);
  c.low = a.low - b.upp;
  _FPU_SETCW(cwhigh);
  c.upp = a.upp - b.low;
  _FPU_SETCW(save);
}

inline void mul (Interval a, Interval b, Interval& c)
{

  short save = 0x037f;
  short  cwlow = 0x077f;
  short  cwhigh = 0x0b7f;

  _FPU_SETCW(cwlow);
  c.low = c.min(a.low*b.low, a.low*b.upp, a.upp*b.low, a.upp*b.upp);


  _FPU_SETCW(cwhigh);
  c.upp = c.max(a.low*b.low, a.low*b.upp, a.upp*b.low, a.upp*b.upp);
  _FPU_SETCW(save);
}

inline void div (Interval a, Interval b, Interval& c)
{
  short save = 0x037f;
  short  cwlow = 0x077f;
  short  cwhigh = 0x0b7f;

  if(b.low==0 || b.upp==0){
	std::cerr << "error: attempt to divide by zero";
  }
  else {
    _FPU_SETCW(cwlow);
    c.low = c.min(a.low/b.low, a.low/b.upp, a.upp/b.low, a.upp/b.upp);
    _FPU_SETCW(cwhigh);
    c.upp = c.max(a.low/b.low, a.low/b.upp, a.upp/b.low, a.upp/b.upp);
    _FPU_SETCW(save);
  }
}


inline Interval sqrt(Interval a)
{

  short save = 0x037f;
  short  cwlow = 0x077f;
  short  cwhigh = 0x0b7f;
  
  Interval b;
  if(a.low < 0.)
     std::cerr << "error: negative value inside sqrt"; 
  _FPU_SETCW(cwlow);
  b.low = sqrt(a.low);
  _FPU_SETCW(cwhigh);
  b.upp = sqrt(a.upp);
  _FPU_SETCW(save);

  return b;
}  


#endif

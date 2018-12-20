#ifndef _RAT_H_
#define _RAT_H_

#include "Int.h"   /* Depends on arbitrary length ints */

/* Header file rat.h -- supporting rational numbers. */

/* This code uses some of the algorithms originally found in Knuth's
   The Art of Computer Programming, vol. 2. */

class Rat {
  Int num;			/* Numerator, denominator */
  Int den;
public:
  Rat ();			/* Constructor functions */
  Rat (long, long);
  Rat (Int&, Int&);
  Rat (const Rat&);
  Rat (double);
  Rat (int);   // added
  ~Rat ();
  operator double();
  operator Int();

  char sign ();			/* For changing and getting signs */
  char sign (char);
  void simplify ();		/* Make rel. prime */
  //  friend void copy_Rat (Rat&,	/* Copy rationals */
  //			Rat&);
  friend void make_Rat (const Rat&,	/* Like copy_Rat, but for the constructor */
			Rat&);
  friend int comp (Rat&, Rat&);	/* Comparison function */
  friend void add (Rat&,
		   Rat&,
		   Rat&);	/* Addition */
  friend void sub (Rat&,
		   Rat&,
		   Rat&); 	/* Subtraction */
  friend void mul (Rat&,
		   Rat&,
		   Rat&); 	/* Multiplication */
  friend void div (Rat&,
		   Rat&,
		   Rat&); 	/* Division */
  friend void pow (Rat&,	/* Raise to a power */
		   long,
		   Rat&);
  friend void copy_Rat (Rat &,	/* Copy rationals */
			Rat &);
  friend int sign(Rat&);	/* Sign of Rat: -1, 0, 1 */
  friend std::istream &operator >> (std::istream &, Rat &);
  friend std::ostream &operator << (std::ostream &, Rat &);
  friend Rat fabs (Rat&);      // absolute value   : added

/* Operators */

/* Arithmetic */

  friend Rat operator + (Rat&, Rat&);
  friend Rat operator - (Rat&, Rat&);
  friend Rat operator * (Rat&, Rat&);
  friend Rat operator / (Rat&, Rat&);
  friend Rat operator % (Rat&, Rat&);
  friend Rat operator ^ (Rat&, long);	/* Exponentiation, not xor */
  friend Rat operator - (Rat&, int);  // added
  friend Rat operator * (Rat&, int);  // added
  friend Rat operator / (Rat&, int);  // added

/* Comparison */

  friend int operator > (Rat&, Rat&);
  friend int operator < (Rat&, Rat&);
  friend int operator >= (Rat&, Rat&);
  friend int operator <= (Rat&, Rat&);
  friend int operator == (Rat&, Rat&);
  friend int operator != (Rat&, Rat&);

/* Assignment */

  Rat operator = (Rat&);
  Rat operator += (Rat&);
  Rat operator -= (Rat&);
  Rat operator *= (Rat&);
  Rat operator /= (Rat&);
  Rat operator %= (Rat&);
  Rat operator ^= (long);

  Rat operator = (double);  //added
  
};

inline char Rat::sign()

{
  return (num.sign());
}

inline char Rat::sign(char s)

{
  return (num.sign(s));
}

inline Rat::Rat (): num(), den()

/* Constructor with no args */

{}

inline Rat::Rat (long n, long d): num(n), den(d)

/* Constructor with n/d */

{
  simplify();
}

inline Rat::Rat (Int &I1, Int &I2): num(), den()

/* Constructor with two Ints. */

{
  copy_Int (I1, num);
  copy_Int (I2, den);
  simplify();
}

inline Rat::Rat (const Rat &R)

/* Constructor from a Rat */

{
  make_Rat (R, *this);
}


inline Rat::Rat (int i): num((long) i), den(1)
{}


inline Rat::~Rat ()

/* Destructor */

{}

/* Operators */

inline Rat operator + (Rat &a, Rat &b)

{
  Rat c; add (a,b,c); return c;
}

inline Rat operator - (Rat &a, Rat &b)

{
  Rat c; sub (a,b,c); return c;
}

inline Rat operator * (Rat &a, Rat &b)

{
  Rat c; mul (a,b,c); return c;
}

inline Rat operator / (Rat &a, Rat &b)

{
  Rat c; div (a,b,c); return c;
}

inline Rat operator - (Rat &a, int l)  //added

{
  Rat c; Rat b = Rat(l); sub (a,b,c); return c; 
}

inline Rat operator * (Rat &a, int l)  //added

{
  Rat c; Rat b = Rat(l); mul (a,b,c); return c;
}

inline Rat operator / (Rat &a, int l)  //added

{
  Rat c; Rat b = Rat(l); div (a,b,c); return c;
}


inline Rat operator ^ (Rat &a, long b)

{
  Rat c; pow (a,b,c); return c;
}

inline int operator > (Rat &a, Rat &b)

{
  return (comp (a,b) > 0);
}

inline int operator < (Rat &a, Rat &b)

{
  return (comp (a,b) < 0);
}

inline int operator >= (Rat &a, Rat &b)

{
  return (comp (a,b) >= 0);
}

inline int operator <= (Rat &a, Rat &b)

{
  return (comp (a,b) <= 0);
}

inline int operator == (Rat &a, Rat &b)

{
  return (comp (a,b) == 0);
}

inline int operator != (Rat &a, Rat &b)

{
  return (comp (a,b) != 0);
}

inline Rat Rat::operator = (Rat & R)

{
  copy_Rat (R, *this); return *this;
}

inline Rat Rat::operator += (Rat &R)

{
  add (R, *this, *this); return *this;
}

inline Rat Rat::operator -= (Rat &R)

{
  sub (*this, R, *this); return *this;
}

inline Rat Rat::operator *= (Rat &R)

{
  mul (*this, R, *this); return *this;
}

inline Rat Rat::operator /= (Rat &R)

{
  div (*this, R, *this); return *this;
}

inline Rat Rat::operator ^= (long l)

{
  pow (*this, l, *this); return *this;
}

inline Rat Rat::operator = (double d)      //added

{
  Rat R(d); copy_Rat (R, *this); return *this;
}


#endif


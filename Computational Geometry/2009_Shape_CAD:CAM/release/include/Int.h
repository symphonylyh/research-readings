#ifndef _INT_H_
#define _INT_H_

#include <iostream>

/* Header file int.h -- supporting arbitrary length integers. */

/* This code uses some of the algorithms originally found in Knuth's
   The Art of Computer Programming, vol. 2. */

typedef unsigned short ushort;

#define LONGBYTES 4
#define SHORTBYTES 2
#define LONGBITS (8 * LONGBYTES)
#define SHORTBITS (8 * SHORTBYTES)
#define SHORTPERLONG 2
#define BASE 10			/* For output, base 10 */
#define MAXSHORT ((ushort)65535)   /* Largest ushort */
#define HIGHBIT ((ushort)32768)	  /* Binary 1000000000000000 */
#define HIGHLBIT 2147483648UL
#define RADIX 65536UL		/* Radix of our computations */
#define MAXLONG 4294967295UL

#define MAX(x,y) ((x>y) ? (x) : (y))
#define MIN(x,y) ((x<y) ? (x) : (y))

struct intstruct {
  ushort maxlength;		/* Maximum length it can grow */
  ushort length;		/* Current length */
  char sign;			/* 1 means positive, 0 means negative */
  ushort *rep;			/* Representation of the integer */
};

class Int {
  intstruct *internal;
  ushort maxlength();		/* Data abstractions */
  ushort maxlength (ushort);
  ushort length (ushort);
  ushort *rep (ushort *);

  void rehash (ushort, int);	/* Re-allocate space */
  void resize (ushort);		/* For resizing the Int */
  void reinit (ushort);		/* Re-initialize int to this length */
  void killzeros();		/* Kill leading zeros */
  void zero();			/* Make the int into zero */

public:

  char sign();			/* For getting and setting signs */
  char sign (char);		/* Users should use the friend sign func. */
  ushort length();		/* These must be public, but no changing */
  ushort *rep();		/* is allowed unless func is a friend */
  int iszero();			/* Is number = 0? */
  int isone();			/* Is number = 1? */
  Int ();			/* Constructor functions */
  Int (long);
  Int (const Int&);
  ~Int ();			/* Destructor */
  operator double();		/* Convert Int to double */
  operator long();		/* Convert Int to long */

  friend void uadd (Int&, Int&, Int&);	 /* Unsigned add */
  friend void usub (Int&, Int&, Int&);	 /* Unsigned subtract */
  friend void gcd_recur (Int&, Int&, Int&);   /* Recursive gcd */
  friend void copy_Int (Int&,	/* Copy integers */
			Int&);
  friend void make_Int (const Int&,	/* Like copy, but for a constructor */
			Int&);
  friend char *Itoa (Int&);	/* Integer to string */
  friend void atoI (char *, Int&);   /* String to int */
  friend int ucomp (Int&, Int&);   /* Unsigned comparison */
  friend int comp (Int&, Int&);	  /* Signed comparison */

  friend void add (Int&,	/* Add function */
		   Int&,
		   Int&);
  friend void sub (Int&,	/* Subtraction function */
		   Int&,
		   Int&);
  friend void mul (Int&,	/* Multiply */
		   Int&,
		   Int&);
  friend void divmod (Int&,	/* Divide and mod */
		      Int&,
		      Int&,
		      Int&);
  friend void div (Int&,	/* Divide */
		   Int&,
		   Int&);
  friend void mod (Int&,	/* Modulo */
		   Int&,
		   Int&);
  friend void gcd (Int&,	/* Greatest common divisor */
		   Int&,
		   Int&);
  friend void pow (Int&,	/* Raising to a power */
		   unsigned long,
		   Int&);

/* Operators */

/* Arithmetic */

  friend Int operator + (Int&, Int&);
  friend Int operator - (Int&, Int&);
  friend Int operator * (Int&, Int&);
  friend Int operator / (Int&, Int&);
  friend Int operator % (Int&, Int&);
  friend Int operator ^ (Int&, unsigned long);	 /* Exponentiation, not xor */

/* Comparison */

  friend int operator > (Int&, Int&);
  friend int operator < (Int&, Int&);
  friend int operator >= (Int&, Int&);
  friend int operator <= (Int&, Int&);
  friend int operator == (Int&, Int&);
  friend int operator != (Int&, Int&);

/* Assignment */

  Int operator = (Int&);
  Int operator += (Int&);
  Int operator -= (Int&);
  Int operator *= (Int&);
  Int operator /= (Int&);
  Int operator %= (Int&);
  Int operator ^= (unsigned long);
};

void uadd (ushort *, ushort *, ushort *, ushort, ushort, ushort &);
void usub (ushort *, ushort *, ushort *, ushort, ushort, ushort &);
void mul (ushort *, ushort *, ushort *, ushort, ushort);
std::istream &operator >> (std::istream &, Int &);
std::ostream &operator << (std::ostream &, Int &);
ushort basereduce (ushort *, ushort, ushort &);
int ucomp (ushort *, ushort *, ushort);

inline Int::Int ()		/* Constructor with no arguments */

{
  internal = 0;
}

inline Int::Int (const Int &I)

{
  make_Int (I, *this);
}

inline ushort Int::maxlength ()

{
  return internal->maxlength;
}

inline ushort Int::maxlength (ushort v)

{
  internal->maxlength = v;
  return internal->maxlength;
}

inline ushort Int::length ()

{
  return internal->length;
}

inline ushort Int::length (ushort v)

{
  internal->length = v;
  return internal->length;
}

inline char Int::sign ()

{
  return internal->sign;
}

inline char Int::sign (char v)

{
  internal->sign = v;
  return internal->sign;
}

inline ushort *Int::rep ()

{
  return internal->rep;
}

inline ushort *Int::rep (ushort *v)

{
  internal->rep = v;
  return internal->rep;
}

inline int Int::iszero()

/* Is the integer zero? */

{
  if ((length() == 1) &&
      *(rep()) == 0)
    return 1;
  else
    return 0;
}

inline int Int::isone()

/* Is the integer one? */

{
  if ((length() == 1) &&
      *(rep()) == 1)
    return 1;
  else
    return 0;
}

inline void div (Int &I1, Int &I2, Int &I3)

/* Divide I1 by I2. Result goes in I3. */

{
  Int I4;			/* Temporary int for modulo */
  divmod (I1, I2, I3, I4);	/* Get div and mod */
}

inline void mod (Int &I1, Int &I2, Int &I3)

/* Take I1 % I2. Result in I3. */

{
  Int I4;			/* Temporary int for quotient */
  divmod (I1, I2, I4, I3);	/* Get div and mod */
}

/* Operators */

inline Int operator + (Int &a, Int &b)

{
  Int c; add (a,b,c); return c;
}

inline Int operator - (Int &a, Int &b)

{
  Int c; sub (a,b,c); return c;
}

inline Int operator * (Int &a, Int &b)

{
  Int c; mul (a,b,c); return c;
}

inline Int operator / (Int &a, Int &b)

{
  Int c; div (a,b,c); return c;
}

inline Int operator % (Int &a, Int &b)

{
  Int c; mod (a,b,c); return c;
}

inline Int operator ^ (Int &a, unsigned long b)

{
  Int c; pow (a,b,c); return c;
}

inline int operator > (Int &a, Int &b)

{
  return (comp (a,b) > 0);
}

inline int operator < (Int &a, Int &b)

{
  return (comp (a,b) < 0);
}

inline int operator >= (Int &a, Int &b)

{
  return (comp (a,b) >= 0);
}

inline int operator <= (Int &a, Int &b)

{
  return (comp (a,b) <= 0);
}

inline int operator == (Int &a, Int &b)

{
  return (comp (a,b) == 0);
}

inline int operator != (Int &a, Int &b)

{
  return (comp (a,b) != 0);
}

inline Int Int::operator = (Int & I)

{
  copy_Int (I, *this); return *this;
}

inline Int Int::operator += (Int &I)

{
  add (I, *this, *this); return *this;
}

inline Int Int::operator -= (Int &I)

{
  sub (*this, I, *this); return *this;
}

inline Int Int::operator *= (Int &I)

{
  mul (*this, I, *this); return *this;
}

inline Int Int::operator /= (Int &I)

{
  div (*this, I, *this); return *this;
}

inline Int Int::operator %= (Int &I)

{
  mod (*this, I, *this); return *this;
}

inline Int Int::operator ^= (unsigned long l)

{
  pow (*this, l, *this); return *this;
}

#endif

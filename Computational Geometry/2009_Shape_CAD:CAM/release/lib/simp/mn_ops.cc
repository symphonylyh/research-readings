#include "lib1.h"
#include "multinom.h"
using namespace std;

real rabs (real) ; /* finding the absolute value of real */

#ifndef NO
#define NO 0
#endif 

#ifndef YES
#define YES 1
#endif

#if 0
multinomial *multinomial::copy ()

/* Allocates a copy of the multinomial and returns it. */

{
  multinomial *mn2 =            /* The new multinomial */
    new multinomial (ndim, ordlist);
  for (int i=0; i < cosize; i++)
    mn2->bp[i] = bp[i];         /* Copy all the coefficients */
  return (mn2);
}
#endif

multinomial::multinomial (const sht_array &ol) : bp (ol),
cosize(bp.size()),ordlist (ol) 
{   
} 
/* Construct from ndim, ordlist */


multinomial::~multinomial () {}            /* Destructor */

multinomial::multinomial (multinomial &mn) : bp (mn.bp),
cosize(mn.cosize), ordlist (mn.ordlist) {} /* Copy constructor */


istream &operator >> (istream &stream, multinomial *&mn)

/* Read a multinomial */

{
  short nd;
  stream >> nd;                 /* Read ndim */
  sht_array olist(nd);          /* Prepare to read ordlist */
  int i;                        /* Loop counter */
  for (i=0; i < nd; i++)        /* Read in orders */
    stream >> olist[i];
  mn = new multinomial (olist); /* Allocate the multinomial */
  real *bp = mn->bp;            /* For speed */
  for (i=0; i < mn->cosize; i++, bp++)
    stream >> (*bp);            /* Call real inserter */
  return (stream);
}


int 
multinomial::every_points_not_far_away_from_axis_by_eps(real eps)
{
#ifndef USE_RAT
  double multi = 100.0;
#else
  real   multi = Rat(100,1);
#endif
     for (int i=0 ; i < cosize; i++) { 
       //       cout << " bp[" << i<< "] = " << bp[i] << endl;
       real temp = rabs(bp[i]);
       //       cout <<  "rabs(" << bp[i] <<") = " << temp << endl;
       real temp1 = eps * multi ;

       //       cout << " eps * " <<  multi <<" = " << temp1 << endl;
       //       cout << "rabs (bp[i]) > eps * " << multi <<"  = ";
       //      cout<< (rabs(bp[i]) > temp1)<<endl;
	  if (rabs(bp[i]) > eps * multi) return NO;  
     }
     return YES;
}


int multinomial::convex_hull_cross_every_axis()
{

  const int plus = 1;
  const int minus = -1;
//  static int counter = 1;

//  counter++;
 
  int sign1, sign2;
#ifndef USE_RAT
  if ( bp[0] <= 0.0 ) sign1 = minus;
#else
  if ( bp[0] <= Rat(0,1) ) sign1 = minus;
#endif
  else sign1 = plus;

//  if (counter ==1 ) cout << "-----------------------------------"<<endl;

     for (int i=1 ; i < cosize; i++) {
#ifndef USE_RAT
       if (bp[i] >= 0.0) sign2 = plus;
#else
       if (bp[i] >= Rat(0,1)) sign2 = plus;
#endif
       else sign2 = minus;

/*       if (counter == 1) {
	 cout << bp[0] << endl;
	 cout << bp[i] << endl; 
	 cout << " sign1 = " << sign1 <<endl;
	 cout << " sign2 = " << sign2 <<endl;
       }
*/
       
       if ( (sign1 * sign2) == -1) {
//	 if (counter == 1 ) cout << " return YES" <<endl;
	 return YES;
       }
       
       
       
     }
  
//  if (counter == 1) cout << " return NO " << endl;
     return NO;

}

real rabs (real r) 
{
  real temp;
  real result ;
#ifndef USE_RAT
  if ( r > 0) temp = r;
#else
  if ( r > Rat(0,1)) temp = r;
#endif
  else
#ifndef USE_RAT
  temp = ((-1.0) * r);
#else
  temp = (Rat(-1,1) * r);
#endif
#ifdef USE_INTERVAL  
  result = temp.center();
#else
  result = temp;
#endif
  return result;
}

ostream &operator << (ostream &stream, multinomial &mn)

/* Write a multinomial */

{
  int i;                        /* Loop counter */
  int j=0;
  
  stream << " in multinomial output " << endl;
//  stream << " mn.ordlist.size is " << mn.ordlist.size() << endl;

//  for (int ii = 0 ; ii < mn.ordlist.size(); ii++)
//    stream << ii + 1 << "-th degree is " << mn.ordlist[ii] << endl;

  int factor = mn.ordlist[j]+1; /* Tell us when to skip */
  
//  stream << " facotr = " << factor << endl;
  
  stream << mn.bp[0] << ' ';    /* Print out first element */
  for (i=2; i <= mn.cosize; i++) {
 //   stream << " in the for loop ,i = " << i << endl; 
    stream << mn.bp[i-1] << ' ';
    if (i % factor == 0) {      /* Do we need to skip lines? */
      while ((i % factor == 0)  && (j < mn.ordlist.size() - 1) ) {
	/* I changed the above line from Evan's routine to avoid coredump */
	
//	stream << " in while loop j = " << j << " i = " << i << endl;
  
	factor *=               /* Multiply factor */
	  mn.ordlist [++j]+1;
	stream << endl;         /* Skip a line */
      }
      j = 0;                    /* Restore j and factor */
      factor = mn.ordlist[j]+1;
    }
  }
  stream << endl;
//  stream << " out  multinomial output " << endl;
  return (stream);
}

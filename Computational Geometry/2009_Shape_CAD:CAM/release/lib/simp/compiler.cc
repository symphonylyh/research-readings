#include "compiler.h"
#include <iostream>
using namespace std;


#ifdef USE_INTERVAL
real_array::real_array():Interval_array(0){
}

real_array::real_array (int siz) : Interval_array(siz) {

}

real_array::real_array (const real_array& d)
{ Interval_array::init (d.size());
 int i;
 for(i=0;i<n;i++) 
   p[i] = d.p[i];
}
//  memcpy (p, (real *) d.p, n*sizeof(real)); }

real_array real_array::operator=(real_array ra)  
{
  int i;
  if(n == ra.n) {
    for(i=0;i<n;i++)
      p[i] = ra.p[i];
  } else {
    resize(ra.n);
    for(i=0;i<n;i++)
      p[i] = ra.p[i];
  }
  /*
  if (n == ra.n) memcpy (p, (real *) ra.p, n*sizeof(real));
  else {
    resize(ra.n);
    memcpy (p, (real *) ra.p, n*sizeof(real));
    }*/
  return *this;
}
#else
#ifdef USE_RAT
real_array::real_array():Rat_array(0){}

real_array::real_array (int siz) : Rat_array(siz) {}

real_array::real_array (const real_array& d)
{ Rat_array::init (d.size());
  memcpy (p, (real *) d.p, n*sizeof(real)); }
#else

real_array::real_array (const real_array& d)
{ dbl_array::init (d.size());
  memcpy (p, (real *) d.p, n*sizeof(real)); }

real_array::real_array():dbl_array(0){

}
#endif
#endif


real_array::~real_array(){}

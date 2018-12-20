/***************************************************************************
 si_phull.cc           last edit 11.2.97
 **************************************************************************/

#include <cmath>
#include "simpoly.h"
#define SCALE 1.e-2
#define EPS_PHULL 1.e-6
using namespace std;

extern int szplug;

/* Use the projected-polyhedron algorithm to generate a new box.
   si_phull can be called in exactly the same way as si_lp, so
   that the two drivers may be interchanged if desired. */

void print_points (real_array x, real_array y)

{
  int n = x.size();
  cout << "\n";
  for (int i=0; i < n; i++)
    cout << x[i] << " " << y[i] << "\n";
  cout << "\n";
}

static int nonleft (real x0, real y0, real x1, real y1,
			   real x2, real y2)

/* Checks to see if p0, p1, p2 constitute a non-left turn; 1 if yes */

{
  real dx1 = x1-x0;
  real dx2 = x2-x0;
  real dy1 = y1-y0;
  real dy2 = y2-y0;

#ifndef USE_RAT
  if ((dx1*dy2 - dx2*dy1 < 0.) != 0) // explicitly ignore MAYBE (SLA, 12/1/97)
#else
  if (dx1*dy2 - dx2*dy1 < Rat(0,1))
#endif
    return 1;
  else
    return 0;
}

 static real rootof (real x1, real y1, real x2, real y2)

/* Find intersection of the line segment with y=0 */

{
  return (y2*x1-y1*x2) / (y2-y1);
}

static void si_convexify (real_array &x, real_array &y)

/* Traverses the (x,y) points and throws away those which are not in
   the convex hull. The arrays are resized to reflect the new size */

{
  int i;                                /* Loop counter */
  int n=x.size();                       /* Number of points */
  int top=1;                            /* Length of processed points */
  for (i=2; i < n; i++) {
    while ((top>0) && (nonleft (x[top-1],          /* If a non-left turn, throw away */
		 y[top-1], x[top], y[top],
				x[i], y[i])))
      top--;
    top++;
    x[top] = x[i];
    y[top] = y[i];
  }
  x.resize(top+1);
  y.resize(top+1);
}

static void si_traverse (real_array x, real_array y, real &a, real &b)
  
/* Traverse the (x,y) points (already in counterclockwise order)
   to find the interval subtended by the convex polygon */

{                  
  int i;
  int n = y.size();
  int zeroflag = 0;                     /* Whether we're on 0 or not */
#ifndef USE_RAT
  a = -1.;                               /* Signifies a is unassigned */
  b = -1.;                               /* and b too */
  if ((y[0]*y[n-1] <= 0.) != 0) {
    if (y[0]==0.)
      zeroflag = 1;
    a = 0.;                             /* Check if 0 case is true */
#else
  a = Rat(-1,1);                         /* Signifies a is unassigned */
  b = Rat(-1,1);                         /* and b too */
  if ((y[0]*y[n-1] <= Rat(0,1)) != 0) {
    if (y[0]==Rat(0,1))
      zeroflag = 1;
    a = Rat(0,1);                       /* Check if 0 case is true */
#endif
  }
  for (i=0; i < n-1; i++) {
    if (x[i] == x[i+1] && y[i] == y[i+1]) // if 2 adjacent points are the same
      i = i+1;                            // ignore the first one (SLA, 12/1/97
#ifndef USE_RAT
    if ((y[i]*y[i+1] <= 0.) != 0)
#else
    if ((y[i]*y[i+1] <= Rat(0,1)) != 0)
#endif
      if (!zeroflag) {
#ifndef USE_RAT
	if (y[i]*y[i+1] == 0.) {
#else
	if (y[i]*y[i+1] == Rat(0,1)) {
#endif
	  if (y[i] == y[i+1]) { 
	                                    /* Check for a nasty tangency */
	    a = x[i];
	    b = x[i+1];
	    break;
	  } else {                      /* Just passing through */
#ifndef USE_RAT
	    if (a == -1.)
#else
	    if (a == Rat(-1,1))
#endif
	      a = x[i+1];
	    else {
	      b = x[i+1];
	      break;
	    }
	  }
	  zeroflag = 1;
	} else {
	  zeroflag = 0;
#ifndef USE_RAT
	  if (a == -1.)
#else
	  if (a == Rat(-1,1))
#endif
	    a = rootof (x[i], y[i], x[i+1], y[i+1]);
	  else {
	    b = rootof (x[i], y[i], x[i+1], y[i+1]);
	    break;
	  }
	}
      } else zeroflag = 0;
      }
#ifndef USE_RAT
  if (b == -1.) {                        /* VERY rare case -- tangency */
#else
  if (b == Rat(-1,1)) {                  /* VERY rare case -- tangency */
#endif
    cerr << "Warning! Tangency in convex hull traversal.\n";
    b = a;
  }
  if (b < a) {                          /* Swap needed */
    //cout << " swap needed " << " a= " << a << " b = " << b << endl;
    real tmp = a;
    a = b;
    b = tmp;
  }

#ifdef USE_INTERVAL  
 
  a = a.lower();
  
  b = b.upper();
   
#endif

}

static int si_interv (multinomial &mn, int j, real &a, real &b)

/* Determine the intersection of the jth projected convex hull with
   the x-axis. */

{
//  cout << " in si_interv " << endl;
//  cout << " j = " << j << endl;
  int ii, io;                           /* Outer and inner loops */
  int i;                                /* Other loop counter */
  int lastii=1, lastio=1;               /* Their limits */
  int ndim = mn.ndim();
  real_array x(2*mn.ordlist[j] + 2);   /* The arrays of points */
  real_array y(2*mn.ordlist[j] + 2);
/*  {
    cout << "ndim = " << ndim << endl;
    for (i=0; i < ndim; i++)
      cout << i+1 <<"-th degree is " << mn.ordlist[i] << endl;
    
  }
*/
  for (i=0; i < j; i++)           /* Set up inner and outer limits */
    lastio *= mn.ordlist[i] + 1;
  for (i=j+1; i < ndim; i++)
    lastii *= mn.ordlist[i] + 1;

/*
  { static int counter_si_inter = 0;
    cout << " counter_si_inter = " << counter_si_inter << endl;
    counter_si_inter++;
  }
*/
  //  cout << setprecision(5) << mn << endl;

/* The main loop; for each x, get a min and max y */
  int pinc = lastii*(mn.ordlist[j]+1);
  int negflg = 0;                       /* Are there y's < 0? */
  int posflg = 0;                       /* "   "     "   > "? */
  for (i=0; i <= mn.ordlist[j]; i++) {
       real iii = i;
       real orderlist = mn.ordlist[j];
       x[i] = x[x.size()-i-1] =            /* Set up x */
	    iii / orderlist;
       //        cout << " i= " << i << " x=" << x[i] << endl;
#ifndef USE_RAT
       real max ;max = -1.e30;
       real min ; min = 1.e30;
#else
       real max ;max = (Rat)-1.e30;
       real min ; min = (Rat)1.e30;
#endif
       real *p = &mn.bp[i*lastii];
       for (io=0; io < lastio; io++, p += pinc) {
	    real *q = p;
	    //   cout << " *q = " << *q<< endl;
	    
	    for (ii=0; ii < lastii; ii++, q++) {
	      //	  cout << " *q = " << *q<< endl;
	      //	  cout << " min = " << min ;
	      //   cout << " max = " << max << endl;
	      //	  int compare;
	      //	  compare = (*q > max);
	      //	  cout << " compare for" << *q ;
	      //   cout << " > " << max << " is ";
	      //	  cout << compare << endl;
#ifdef USE_INTERVAL
		 if (comp_greater(*q , max))
#else
                 if ( *q > max)
#endif     
		      max = *q;
		 //	  cout << " after compare, max is " << max << endl;
		 
		 
		 //	 compare = (*q < max);
		 //  cout << " compare for" << *q << " < " << min << " is ";
		 // cout << compare << endl;
#ifdef USE_INTERVAL
		 if (comp_less(*q,  min))
#else
		 if (*q < min)
#endif
		      min = *q;
		 //	  cout << " after compare, min is " << min << endl;
	    }
	    
       }
       /*	 real scale = min.magnitude();
	 if(scale < max.magnitude())
	 scale = max.magnitude();
	 min -= EPS_PHULL*scale;
	 max -= EPS_PHULL*scale;
	 */
 //         cout << "before " << " min= " << min << " max= " << max << endl;
#ifdef USE_INTERVAL       
       min = min.lower(); max = max.upper();
#endif
 //         cout << " min= " << min << " max= " << max << endl;
#ifndef USE_RAT
       if (min <= 0.) negflg = 1;	        /* Set flags */
       if (max >= 0.) posflg = 1;
#else
       if (min <= Rat(0,1)) negflg = 1;	        /* Set flags */
       if (max >= Rat(0,1)) posflg = 1;
#endif
       y[i] = min;                     /* Put max and min in array */
       y[y.size()-i-1] = max;      
     }
  
#ifndef USE_RAT
  x[0] = x[x.size()-1] = 0.;
  x[mn.ordlist[j]] = x[x.size()-mn.ordlist[j]-1] = 1.0;
#else
  x[0] = x[x.size()-1] = Rat(0,1);
  x[mn.ordlist[j]] = x[x.size()-mn.ordlist[j]-1] = Rat(1,1);
#endif

  // cout << "before convex hull " << endl;  
  // print_points (x,y);
                                   
  if (!(negflg && posflg)){   /* All above/below axis? */
    //      cout << " all above/below axis " << endl;
    //   cout <<"si_interv j= " << j << endl;
    //   print_points(x,y);
    return (1);
  }


  if(ndim == 1){
       real_array xl(mn.ordlist[j]+1);
       real_array yl(mn.ordlist[j]+1);
       real_array xu(mn.ordlist[j]+1);
       real_array yu(mn.ordlist[j]+1);

       for(i=0; i<mn.ordlist[j]+1; i++){
	 xl[i] = x[i];
	 yl[i] = y[i];
	 xu[i] = x[i+mn.ordlist[j]+1];
	 yu[i] = y[i+mn.ordlist[j]+1];
	 //cout << " xl = " << xl[i] << " yl = " << yl[i] << " xu = " << xu[i] << " yu = " << yu[i] << endl;
       }
       
       si_convexify(xl, yl);
       si_convexify(xu, yu);
       //for(i=0;i<xu.size();i++)
	 //cout << " xu = " << xu[i] << " yu = " << yu[i] << endl;
	 

       int nl = xl.size();
       int nu = xu.size();
       int nt = nl + nu;

       x.resize(nt);
       y.resize(nt);

       for(i=0; i<nl; i++){
	 x[i] = xl[i];
	 y[i] = yl[i];
       }

       for(i=nl; i<nt; i++){
	 x[i] = xu[i-nl];
	 y[i] = yu[i-nl];
       }
  }
  else  
       si_convexify (x,y);         /* Discard unimportant points */
  //  cout << "after convex hull " << endl;
  //  print_points (x,y);
  si_traverse (x,y,a,b);                /* Traverse to get a and b */

  //  cout << " a= " << a << " b= " << b << endl;
  
  return (0);
}

void si_phull (mn_array &mn_list, Bbox &B, int_array &check)

/* Find a new box */

{
//  cout << " in si_phull " << endl;
 

  int ndim = mn_list[0].ndim(); /* # of dimension and # of unknown */
  int n_eq = mn_list.size();        // number of eqtaions
  //  cout << "mn_list[0] is " << endl;
  //  cout << mn_list[0] << endl;
  //  if (n_eq > 1) {  cout << "mn_list[1] is " << endl;
  //		   cout << mn_list[1] << endl;}
  int i,j;

  //  cout << "ndim = " << ndim << endl;
  // cout << "n_eq = " << n_eq << endl;

  for (i=0; i < ndim; i++) {            /* For each dimension */
    if (!check[i]) continue;            /* Skip over it if small enough */
#ifndef USE_RAT
    B.a[i] = 0.;                         /* Init to [0,1] */
    B.b[i] = 1.;
#else
    B.a[i] = Rat(0,1);                         /* Init to [0,1] */
    B.b[i] = Rat(1,1);
#endif
    for (j=0; j < n_eq; j++) {  /* for each eqation */
      if (!mn_list[j].ordlist[i])       /* Only if ith variable appears */
	continue;

      real aij, bij;                    /* Will hold interval */

      //    cout << " before si_interv, i = " << i << " j= " << j<< endl;

      
      int flag = si_interv              /* Call interval driver */
	(mn_list[j], i, aij, bij);

      //   cout << " box is " << endl;
      // cout << B << flush;
      // cout << "aij = " << aij << ", bij = " << bij << endl;
      // cout << "\n si_phull flag= " << flag << endl;

/*    real delta = bij-aij;
      aij -= delta*SCALE;
      bij += delta*SCALE;
*/


#ifdef USE_INTERVAL
      aij = aij.lower();	/* Tolerancing */
      bij = bij.upper();
#endif
      if (flag) {
	B.stat = 1;
	//	cout << " stat become 1 for flag is 1 " << endl;
	return;
      }
      if ((aij > B.b[i]) ||             /* Make sure intervals intersect */
	  (bij < B.a[i])) {
	//  cout << " i = " << i << endl;
	//  cout << " aij = " << aij << " B.b[" << i << "]="<<B.b[i]<<endl;
	//  cout << " bij = " << bij << " B.a[" << i << "]=" <<B.a[i]<<endl;
	  B.stat = 1;
	  
	  //	cout << " stat become 1 for interva not intersect " << endl;
	return;
      }
#ifdef USE_INTERVAL
      if (aij > B.a[i])
	B.a[i] = aij.lower();
      if (bij < B.b[i])
	B.b[i] = bij.upper();
#else
      if (aij > B.a[i])
	B.a[i] = aij;
      if (bij < B.b[i])
	B.b[i] = bij;
#endif
    }  // for each equation
  } // for each dimension
  B.stat = 0;
}


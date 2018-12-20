/*****************************************************************************
  si_inter.cc            last edit 10.22.92
 ****************************************************************************/


#include "simpoly.h"
#include "compiler.h"
/* The heart of the intersection algorithm */

#ifndef NO 
#define NO 0
#endif

#ifndef YES
#define YES 1
#endif
using namespace std;

int PRINT_BOX = NO;

extern int si_reset_flag;       /* Do we reset everything? */
int szplug;

void si_multipush (Rstack &stack, int_array &l, Bbox &B)

/* Push unresolved problems onto the stack */

{
//  cout << " in si_multipush " << endl;
  int i,j;                      /* Loop counters */
  int n = B.n;                  /* Number of dimensions */
  int m = 1;                    /* Number of recursive subproblems */
  for (i=0; i < n; i++)
    if (l[i]) m *= 2;           /* Get new value of m */
  Bbox **B_list =               /* List of boxes */
    new Bbox*[m];
  for (i=0; i < m; i++)         /* Allocate B_list elements */
    B_list[i] =
      new Bbox(n);

//  cout << "original box is " << B << endl;
  for (i=0; i < m; i++) {       /* Set up boxes */
    int k=1;                    /* k = 2^(number of 1's we have passed) */
    for (j=0; j < n; j++)
      if (l[j]) {
//		cout << j << "-th dir has  been subdivided " << endl;
	if (i & k) {            /* Position equal to 1? */
#ifdef USE_INTERVAL
	  B_list[i]->a[j] =     /* Second half of the interval */
	    (0.5 * (B.b[j]+B.a[j])).lower();
	  B_list[i]->b[j] =
	    B.b[j].upper();
#else
#ifdef USE_RAT
	  B_list[i]->a[j] =     /* Second half of the interval */
	    (Rat(1,2) * (B.b[j]+B.a[j]));
	  B_list[i]->b[j] =
	    B.b[j];
#else
	  B_list[i]->a[j] =     /* Second half of the interval */
	    (0.5 * (B.b[j]+B.a[j]));
	  B_list[i]->b[j] =
	    B.b[j];
#endif
#endif
	  
	} else {
#ifdef USE_INTERVAL
	  B_list[i]->a[j] =     /* First half of interval */
	    B.a[j].lower();
	  B_list[i]->b[j] =
	    (0.5 * (B.b[j]+B.a[j])).upper();
#else
#ifdef USE_RAT
	  B_list[i]->a[j] =     /* First half of interval */
	    B.a[j];
	  B_list[i]->b[j] =
	    (Rat(1,2) * (B.b[j]+B.a[j]));
#else
	  B_list[i]->a[j] =     /* First half of interval */
	    B.a[j];
	  B_list[i]->b[j] =
	    (0.5 * (B.b[j]+B.a[j]));
#endif
#endif
	}
	k *= 2;                 /* Increment k's exponent; we passed a 1 */
      } else {                  /* No binary subdivision in this index */
//	cout << j << "-th dir has not been subdivided " << endl;
	B_list[i]->a[j] =
	  B.a[j];
	B_list[i]->b[j] =
	  B.b[j];
      }
  }
//  cout<<"one multinomial has been subdivided into " << m << "multis"<< endl;
  
  for (i=0; i < m; i++) {
    stack.push (B_list[i]);     /* Push onto the stack */
//    cout << *(B_list[i]) << endl;
  }
  delete B_list;                /* Free array of pointers (but not the boxes!) */

//  cout << " out si_multipush " << endl;
}

void si_sub (mn_array &mns, mn_array &mn_list, Bbox &B)

/* Subdivide mns, making mn_list correspond to the subdivision of mns
   according to the box B. Calls mnr_sub as a driver. */

{
  int i,j;
  int n = mns.size();   // n is number of equations
  int n_dim;  /* number of unknows */
  n_dim = mns[0].ndim();



//  cout << "inside si_sub " << endl;
//  cout << " B = " << B << endl;
//  cout << " no of eqs = " << n << endl;
//  cout << " no of unknowns = " << n_dim << endl;

/*  {
    for (int ii = 0; ii < n ; ii++) {
      cout << "mn_list[" << ii << "] is " << endl;
      cout << mn_list[ii] << endl;
    }
  }
*/
  for (i=0; i < n; i++) {  // number of equations
       
//       cout << " i = " << i << endl;
//       cout << B.a[0].lower() << endl;
//       cout << B.b[0].upper() << endl;
//       cout << mn_list[i] << endl;

#ifdef USE_INTERVAL
    mns[i].sub                 /* Subdivide in first direction */
      (0, B.a[0].lower(), B.b[0].upper(), &mn_list[i]); 
#else
    mns[i].sub                 /* Subdivide in first direction */
      (0, B.a[0], B.b[0], &mn_list[i]); 
#endif

    
    for (j=1; j < n_dim; j++) {  /* Subdivide in other directions */

//	 cout << " B.a[" << j << "] = " << B.a[j] ;
//	 cout << " B.b[" << j << "] = " << B.b[j] << endl;
#ifdef USE_INTERVAL
      mn_list[i].sub
            (j, B.a[j].lower(), B.b[j].upper(), &mn_list[i]);
#else
      mn_list[i].sub
	(j, B.a[j], B.b[j], &mn_list[i]);
#endif
     }
   }

/*  cout << " after sub " << endl;
  {
    for (int ii = 0; ii < n ; ii++) {
      cout << "mn_list[" << ii << "] is " << endl;
      cout << mn_list[ii] << endl;
    }
  }
    cout << " multi.sub " << endl; 
*/
}

rootlist *si_pinter (mn_array &mns, real eps)

/* Intersect the multinomials iteratively using the
   Projected-Polyhedron algorithm. Start with the box [0,1] x ... x [0,1]
   bounding all possible roots, and then subdivide the multinomials,
   eliminating areas where no roots can occur.  This elimination results
   from checking intersections of the convex hulls of the multinomials
   with one another. If at any point the size of the new box is more than
   80% of the current interval size, split the problem into two
   subproblems of approximately equal size by splitting the interval in
   half. Use an event stack to hold all the pending subproblems.

   When all intervals are small enough, we stop splitting and splice
   the root onto the list of roots which already exists. Return 0
   if something very wrong occurs.

   We assume without checking that the size of mn_list is equal to the
   number of dimensions of any multinomial in the list. */

{
  rootlist *points;		/* The list we return */
  int  i;               /*  loop counters */
  int ndim = mns[0].ndim();  /* number of unknowns */
  
  mn_array mn_list(mns);     /* Major list for subdivision -- make a copy */
  Bbox *B = new Bbox (ndim);    /* The current bounding box */
  Rstack stack;                 /* The event stack */
  int nbinary = 0 ;           /* the number of binary subdivision */
  int counter = 0 ;           /* the number of Iterations */
  int lastcheck=0;            /* Added 1/3/94; this flag checks to see 
    if we are doing a final verification check
      on a sufficiently small box. */
  
  //     cout << " just in si_pinter " << endl;
  //     cout << " ndim = " << ndim << endl;
  
  //     cout << " mn_list[0] = \n " << setprecision(15) << mn_list[0] << endl;
  
  //     if (ndim >1)
  //     cout << " mn_list[1] = \n " << setprecision(15) << mn_list[1] << endl;


  /* Initialize points and mn_list */
  
  points = new rootlist;
  
  stack.push (B);             /* Push governing data onto the event stack */
  
  /* Now begin the intersection process */
  
  while (stack.pop (&B)) {    /* While there are things to process */
    
//    cout << "*******************************************" <<endl;
//    cout << " counter = " << counter << endl;//    counter++;
    
    counter++;
    /* Reduce the box in question to the necessary size.
      Use a bounding box approach to generate a new box. */
    if (PRINT_BOX) {
      //      cout << "try box " << setprecision(15) << *B << flush;
      //      cout << endl;
    }

    int_array check(ndim);
    Bbox *B2 = new Bbox(ndim);  /* Allocate B2 */
    for (i=0; i < ndim; i++)
      if (B->b[i] - B->a[i] < eps) {
	
	//		     cout << " i= " << i << " a= " << B->a[i] ;
	//		     cout << " b= " << B->b[i] <<endl;
	real b_a = B->b[i] - B->a[i];
	//		     cout << " b - a= " << b_a << endl;
	real mid =              
#ifndef USE_RAT
	  (B->b[i] + B->a[i]) / 2.;
	
	B->a[i] = mid-eps/2.01;
	/* This makes some artificial correction */
	B->b[i] = mid+eps/2.01;        
#else
	  (B->b[i] + B->a[i]) / Rat(2,1);
	
	B->a[i] = mid-eps/Rat(201,100);
	/* This makes some artificial correction */
	B->b[i] = mid+eps/Rat(201,100);        
#endif
	
	check[i] = 0;
#ifndef USE_RAT
	B2->a[i] = 0.0;
	B2->b[i] = 1.0;
#else
	B2->a[i] = Rat(0,1);
	B2->b[i] = Rat(1,1);
#endif
//	cout << " check[" << i <<"] = " << check[i] << endl;
      } else {
	check[i] = 1;
//	cout << " check[" << i <<"] = " << check[i] << endl;
      }
//    cout << "inside si_pinter just before si_sub" << endl;
//    cout << " B = " << setprecision(8) <<  *B << endl;
    
    //	  cout << " mn_list[0] = \n " << mn_list[0] << endl;
    //         if (ndim >1)
    //	  cout << " mn_list[1] = \n " << mn_list[1] << endl;
    
    
   si_sub (mns, mn_list, *B);  /* Subdivision or Bezier clipping */
    si_phull(mn_list, *B2, check); /* Generate it */
    // 	  cout << " after su_sub" << endl;
    //	  cout << " last check = " << lastcheck << endl;
    //  	  cout << " mn_list[0] = \n " << mn_list[0] << endl;
    //        if (ndim >1) cout << " mn_list[1] = \n " << mn_list[1] << endl;
    
    if (B2->flag()) {           /* Nothing in the box? */
    //  cout << " No root" << endl;
      delete B2;		/* Added */
      delete B;		/* Added */
      continue;                 /* Break away */
    }
    B2->scale (*B);      /* Scale the second box with the first */
    
//    cout << " after scale " << endl;
//    cout << " B2 = " << setprecision(8) <<  *B << endl;
    
    /* Now, check the new box. We may need to do some binary */
    /* subdivision. If so, perform the subdivision and push the other */
    /* subproblems onto the event stack. */
	  if (lastcheck) {
	    real_array root(ndim);
	    for (i=0; i < ndim; i++) {
#ifdef USE_INTERVAL		
	      Interval low(0.0);
	      Interval high(0.0);
	      Interval epsilon_correction(0.0); 
	      // due to artificial correction (above)
	      epsilon_correction = eps * ( 1.0 / 2.0 - 1.0 / 2.01);
	      //		    cout << " ep_cor = "<< epsilon_correction << endl;
	      //		    cout << " B2->a[i] = " << B2->a[i] << endl;
	      //		    cout << " B2->b[i] = " << B2->b[i] << endl;
	      
	      low = B2->a[i]  - 1.0 * epsilon_correction;
	      high = B2->b[i] + 1.0 * epsilon_correction;
	      
	      //		    cout << " low = " << low << endl;
	      //		    cout << " high = " << high << endl;
	      
	      root[i] = merge(low, high);
	      
	      //		    root[i] = merge(B2->a[i], B2->b[i]);  
		    /* Put points in */ 
#else
#ifdef USE_RAT
	      root[i] = (Rat(1,2) * (B2->a[i] + B2->b[i]));
#else
	      root[i] = (0.5 * (B2->a[i] + B2->b[i]));
#endif
#endif
      //   cout << root[i] << " , ";
	    }
//cout << endl;
//	    cout << " before final check " << endl;


	    int yes_or_no;
	    yes_or_no = mn_list.convex_hull_cross_every_axis();

	    if( yes_or_no == YES)  /* Added by Hu 2/4/94 */
	    {
	      //cout << "Add to points!" << endl;
	      points->insert(root);
		 }
	    delete B2;
	    lastcheck = 0;
	  }
	  else if (B2->small(eps)) {
	    /* Small enough? */ 
	    //	  cout << "B2 small enough" << setprecision(15) << *B2 << endl;
	    //	         cout << " eps = " << eps << endl;
	    //	       for (i=0; i<B2->n; i++)
	    //	         cout <<"b[i] - a[i] = "<< B2->b[i] - B2->a[i] << endl;
	    //	         cout << " root" << endl;
	    lastcheck = 1;
	    stack.push(B2);
	  } else {
	    int_array flags (ndim);
	    B->subcheck               /* Subdivison flags */
	      (*B2, check, flags);
	    if (B->flag()) {          /* Binary subdivision needed? */
	      si_multipush            /* Push new boxes onto the stack */
		(stack, flags, *B2);
	      delete B2;		/* Added */
	      nbinary++;
//	      	 	cout << " Binary sub." << endl;
	    } else {                  /* No binary subdivision needed */
	      stack.push (B2);        /* Push it on for reconsideration */
	      //	cout << endl;
	    }
	  }
    delete B;                   /* Trash B */
  }                             /* End of while loop */
     //cout << " number of roots = " << points->size() << endl;
     //cout << " number of iterations = " << counter << endl;
     //cout << " number of binary subdivision = " << nbinary << endl;
	  
  points->head();		/* Move back to head of list */
  return (points);
}                               /* End of si_pinter */

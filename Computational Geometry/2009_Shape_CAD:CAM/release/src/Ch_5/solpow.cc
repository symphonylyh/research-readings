// Copyright (C) Massachusetts Institute of Technology, 2008
// All rights reserved

/* ------------------------------------------------------------------------
                            solpow.cc
  
   Program to illustrate the Projected Polyhedron nonlinear system solver
   (w/ input as power basis polynomial equations) 

   How to make:
   prompt> make
   How to run:
       Example 5.6.1 in page 127,
       prompt> ex.5.6 ex.5.6.1.in
       Example 5.6.2 in page 133,
       prompt> ex.5.6 ex.5.6.2.in

       Note for "Example 5.6.2":

       In order to use the IPP solver, we need a re-parametrization such that:
       x = 4t - 2, y = 2s -1 where 0 <= s,t <= 1
       and re-formulate the given equations f and g to get the corresponding 
       coefficients data (as shown in ex.5.6.2.in)

       We also note that a lot of root boxes are generated during solution process
       due to the tangential intersection at (s,t)=(0.5,1) i.e. (x,y)=(2,0).
       Such root boxes are merged and consolidated as one root through the 
       root consolidation process.

       Also note you will need to substitute the resulting roots in (s,t) into:
       x = 4t - 2 and y = 2s -1 to have the roots in (x,y) as mentioned above.
------------------------------------------------------------------------- */

#include <iostream>
#include <fstream>
using namespace std;


#include "lexPP.h"
#include "monoPP.h"
#include "pp_solver.h"
#include "consolidate.h"
#ifdef USE_INTERVAL
#define ARITH "(RIA)"
#else
#define ARITH "(FPA)"
#endif

int main(unsigned argc, char *argv[])
{
  if (argc < 2) {
    cerr << "usage: " << argv[0] << " <input file>" << endl;
    exit(-1);
  }

  fstream fin;                                    // open input data file
  fin.open(argv[1], ios::in); 
  if(!fin) {
	cerr << "Can't open an input file: " << argv[1] << endl;
      exit(-1);
	}
   
  real **bp;
  short **ordlists;
  int numeq, numvar;
  short i, j, size;                                

  fin >> numeq >> numvar;

  bp = (real**) new real*[numeq];
  ordlists = (short**) new short*[numeq];
  for(i=0;i<numeq;i++)
    ordlists[i] = new short[numvar];

  // ifstream inFile("six-dimension.inp");
  int t2,tsum;
  double ttin;
  for(i=0;i<numeq;i++) {
    // fin >> temp;
    tsum = 1;
    for(j=0;j<numvar;j++) {
      fin >> t2;
      /* ordlists[i][j] = t2+1; this is for solbern version which needs "order" */  
      ordlists[i][j] = t2; /* solpow needs "degree" */
      tsum *= (t2+1);
    }

    bp[i] = new real[tsum];

    for(j=0;j<tsum;j++) {
      fin >> ttin;
      bp[i][j] = ttin;
    }
  }

  fin.close();

  cout << "Parse equations and convert to Bernstein basis " << ARITH << endl;
  cout << "Number of equations = " << numeq  << endl;
  cout << "Number of variables = " << numvar << endl;

  monoToBern(bp,ordlists,numeq,numvar);    // convert power basis to Bernstein basis for solver

  for (i=0; i<numeq; i++) {                       // coefficients
    cout << "Equation " << i+1 << ':';
    size = 1;
    for (j=0; j<numvar; j++) {
      size *= ordlists[i][j];
      cout << ' ' << ordlists[i][j];
    }
    cout << endl;
    
    for (j=0; j<size; j++)
      cout << ' ' << bp[i][j];
    cout << endl;
  }

  real **roots;                                   // roots of system
  real eps = 1.0e-8;                              // solver tolerance
  int numroots;                                   // number of roots
  sim_solve(bp, ordlists, numeq, numvar, eps, &roots, &numroots);

  cout << "Solution for system of " << numeq << " equations with " << numvar;
  cout << " variables " << ARITH << endl;
  cout << "Number of roots = " << numroots << endl;
  cout.precision(15);              // set printout precision to 15 digits
  for (i=0; i<numroots; i++) {
    cout << "Root " << i+1 << ':';
    for (j=0; j<numvar; j++)
      cout << ' ' << roots[i][j];
    cout << endl;
  }

  cout << endl << "Consolidating the roots... " << endl;
  real epsCon = 0.1;                           // consolidate tolerance
  consolidate(epsCon, eps, bp, ordlists, numeq, numvar, roots, numroots);

  cout << "Root consolidation result: " << endl;
  cout << "Number of roots = " << numroots << endl;
  cout.precision(15);              // set printout precision to 15 digits
  for (i=0; i<numroots; i++) {
    cout << "Root " << i+1 << ':';
    for (j=0; j<numvar; j++)
      cout << ' ' << roots[i][j];
    cout << endl;
  }

#ifdef USE_FLOAT
cout << endl;
cout << "Note: " << ARITH << " solver may miss the roots and the consolidation procedure may miss the roots further so, we recommend to use the Interval Projected Polyhedron algorithm i.e. (RIA) version, instead."
     << endl;
#endif


  for (i=0; i<numeq; i++)                         // deallocate memory
    delete [] ordlists[i];
  delete [] ordlists;

  for (i=0; i<numeq; i++)
    delete [] bp[i];
  delete [] bp;

  for (i=0; i<numroots; i++)
    delete [] roots[i];
  delete[] roots;
  return 0;
}

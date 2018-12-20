// Copyright (C) Massachusetts Institute of Technology, 2008
// All rights reserved

/* ------------------------------------------------------------------------
                            solbern.cc

   Program to illustrate the Projected Polyhedron nonlinear system solver
   (w/ input as Bernstein basis polynomial equations) 

   How to make:
   prompt> make
   How to run:
   (1) for floating point arithmetic:
       prompt> solbern-fpa <input file name>
   (2) for rounded iterval arithmetic:
       prompt> solbern-ria <input file name>

   Examples:
 
      (1) intersection between a circle x^2 + y^2 = 1 and 
         an ellipse (x^2)/4 + 4(y^2) = 1 in the first quadrant: 
       prompt> solbern-fpa solbern.in
       prompt> solbern-ria solbern.in

      (2) degree 20 Wilkinson polynomial in Section 4.9 of the hyperbook
       prompt> solbern-fpa wil.in
       prompt> solbern-ria wil.in
       Note solbern-fpa misses some roots due to robustness issues described in 
       section 4.6 of the hyperbook

      (3) another example: finding the extrema of the squared distance between 
          two spheres, which results in 6 equations w/ 6 unknowns (xi, i=1,2,..,6):
          (x1-0.2)^2 + (y1-0.2)^2 + (z1-0.2)^2 - 0.04 = 0
          (x2-0.2)^2 + (y2-0.2)^2 + (z2-0.8)^2 - 0.04 = 0 
          (x1-x2)(z1-0.2) - (z1-z2)(x1-0.2) = 0
          (y1-y2)(z1-0.2) - (z1-z2)(y1-0.2) = 0
          (x1-x2)(z2-0.8) - (z1-z2)(x2-0.2) = 0
          (y1-y2)(z2-0.8) - (z1-z2)(y1-0.2) = 0
       prompt> solbern-fpa solbern6.in
       prompt> solbern-ria solbern6.in  
------------------------------------------------------------------------- */

#include <iostream>
#include <fstream>
using namespace std;


#include "lexPP.h"
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
      ordlists[i][j] = t2+1;
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
  real epsCon = 1.0e-1;                           // consolidate tolerance
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
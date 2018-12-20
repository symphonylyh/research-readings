// Copyright (C) Massachusetts Institute of Technology, 2008
// All rights reserved

/* ------------------------------------------------------------------------
                            sollex.cc
  
   Program to illustrate the Projected Polyhedron nonlinear system solver
   (w/ input as lexical representation of power basis polynomial equations) 

   How to make:
   prompt> make
   How to run:
   (1) for floating point arithmetic:
       prompt> sollex-fpa <input file name>
   (2) for rounded iterval arithmetic:
       prompt> sollex-ria <input file name>

   Example: 3 eqns, 3 unknowns:
            10x - 20y + 30z = 14
            20x + 10y - 40z = -2
           -30x + 40y - 10z = -2
       prompt> sollex-fpa sollex.in
       prompt> sollex-ria sollex.in   
------------------------------------------------------------------------- */


#include <iostream>
#include <fstream>

#include "lexPP.h"
#include "pp_solver.h"
#include "consolidate.h"
#ifdef USE_INTERVAL
#define ARITH "(RIA)"
#else
#define ARITH "(FPA)"
#endif

using namespace std;

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
  fin.close();

  real **bp;
  short **ordlists;
  int numeq, numvar;
  lexConv(argv[1], bp, ordlists, numeq, numvar);  // parse and convert file

  cout << "Parse equations and convert to Bernstein basis " << ARITH << endl;
  cout << "Number of equations = " << numeq  << endl;
  cout << "Number of variables = " << numvar << endl;

  short i, j, size;                               // print Bernstein
  for (i=0; i<numeq; i++) {                       // coefficients
    cout << "Equation " << i+1 << ':';
    size = 1;
    for (j=0; j<numvar; j++)
      size *= ordlists[i][j];
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

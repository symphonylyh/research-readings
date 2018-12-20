// Copyright (C) Massachusetts Institute of Technology
// All rights reserved

#include "interval.h"
#include "conv.h"
#include "monoPP.h"
#include "pp_solver.h"

using namespace std;

#ifdef USE_INTERVAL
typedef Interval real;
#else
typedef double real;
#endif

void consolidate(real epsCon, real eps, real **bp, short **ordlists,
		 int numeq, int numvar, real **&roots, int &numroots)
{
  int i, j;
  double low, upp;

  if( numroots == 0 ) return;

  else if( numroots == 1 ) {
     for( j = 0; j < numvar; j++ ) {
#ifdef USE_INTERVAL
          low = roots[0][j].get_low();
          upp = roots[0][j].get_upp();
#else
          low = roots[0][j];
          if (low > 1.0)
	      low = 1.0;
	    else if (low < 0.0)
	      low = 0.0;
          else ;
#endif
          
#ifdef USE_INTERVAL
          real temp(low, upp);
#else
          real temp = low;
#endif
          roots[0][j] = temp;

    }//------> "for( j = ..."

  return;
 }//---------> "else if( numroots == 1 )"



  int numRoots = numroots; 
  int row;
#ifdef USE_INTERVAL
  real dummy = eps;
  real **root = roots;
#else
  Interval** root = new Interval*[numroots];
  for (row=0; row<numroots; row++)
    *(root + row) = new Interval [numvar];

  for (i=0; i<numroots; i++)
    for (j=0; j<numvar; j++) {
      Interval temp( roots[i][j] - eps/2.0, roots[i][j] + eps/2.0); 
      root[i][j] = temp;
    }
#endif

  int preNumRoots = 0;
  int k;
  while (preNumRoots != numroots) {
    preNumRoots = numroots;

    for (i=0; i<numroots-1; i++)
      for (j=i+1; j<numroots; j++) {
	short decision = 1;

	for (k=0; k<numvar && decision; k++)
	  decision *= ((root[j][k].get_low() >= root[i][k].get_low() &&
			root[j][k].get_low() <= root[i][k].get_upp()) ||
		       (root[i][k].get_low() >= root[j][k].get_low() &&
			root[i][k].get_low() <= root[j][k].get_upp()));

	if (decision ) {
	  numroots--;

	  for (k=0; k<numvar; k++) {
            if (root[i][k].get_low() <= root[j][k].get_low())
	      low = root[i][k].get_low();
            else
	      low = root[j][k].get_low();

            if (root[i][k].get_upp() >= root[j][k].get_upp())
	      upp = root[i][k].get_upp();
            else
	      upp = root[j][k].get_upp();

	    Interval temp(low, upp); 
	    root[i][k] = temp;

            for (int l=j; l<numroots; l++)
	      root[l][k] = root[l+1][k];
	  }
          j--;
	}//---------> "if (decision)"
      }//----------> double "for"
  }//------------> "while (preNumRoots == numroots)"


  cout << "Result from merging intersecting root-boxes: " << endl;
  cout << "Number of roots = " << numroots << endl;
  cout.precision(15);              // set printout precision to 15 digits
  for (i=0; i<numroots; i++) {
    cout << "Root " << i+1 << ':';
    for (j=0; j<numvar; j++)
      cout << ' ' << root[i][j];
    cout << endl;
  }


  //*
  if (numroots != 1) {
    cout << endl
         << "Consolidating further by checking if a root exists in each merged root-box..." << endl;
    
    short **ordLists = new short*[numeq];
    for (i=0; i<numeq; i++) {
      ordLists[i] = new short [numvar];
      for (j=0; j<numvar; j++)
	ordLists[i][j] = ordlists[i][j] - 1;
    }

    bernToMono(bp, ordLists, numeq, numvar);

    //__________________________NEW_________________
    real **tmpBp = new real*[numeq];

    for (i=0; i<numeq; i++) {
      int cosize = 1;
      for (j=0; j<numvar; j++)
	cosize *= ordLists[i][j]; 

      tmpBp[i] = new real [cosize];
    }
    //__________________________NEW_________________

    real **intv = new real*[numvar];
    for (i=0; i<numvar; i++) 
      intv[i] = new real [2];

    real **dumRoots;
    int dumNumRoots;

    for (i=0; i<numroots; i++) {
      for (j=0; j<numvar; j++) {
	intv[numvar-1-j][0] = root[i][j].get_low();
	intv[numvar-1-j][1] = root[i][j].get_upp();
      }

      //__________________________NEW_________________
      for (int ii=0; ii<numeq; ii++) {
	int cosize = 1;
	for (j=0; j<numvar; j++)
	  cosize *= ordLists[ii][j]; 

	for (j=0; j<cosize; j++)
	  tmpBp[ii][j] = bp[ii][j];
      }
      //__________________________NEW_________________

      convPolyBox(tmpBp, ordLists, intv, numvar, numeq);

      for (j=0; j<numeq; j++)
	for (int k=0; k<numvar; k++)
	  ordLists[j][k] = ordLists[j][k] - 1;

      monoToBern(tmpBp, ordLists, numeq, numvar);

      sim_solve(tmpBp, ordLists, numeq, numvar, epsCon, &dumRoots,
		&dumNumRoots);

      if (dumNumRoots == 0) {
	numroots--;

	for (j=i; j<numroots; j++)
	  for (int k = 0; k < numvar; k++)
	    root[j][k] = root[j+1][k];
	i--;
      }
    } //----------> "for (i=0; i<numroots; i++)"

    for (i=0; i<dumNumRoots; i++)
      delete [] dumRoots[i];
    delete [] dumRoots;
      
    for (i=0; i<numvar; i++)
      delete [] intv[i];
    delete [] intv;

    for (i=0; i<numeq; i++)
      delete [] ordLists[i];
    delete [] ordLists;
      
    for (i=0; i<numeq; i++)
      delete [] tmpBp[i];
    delete [] tmpBp;
  } //----------> "if (numroots != 1)"

  //*/

  real** tmpRoots = new real*[numroots];
  for (row=0; row<numroots; row++)
    *(tmpRoots + row) = new real [numvar];

  for (i=0; i<numroots; i++)
    for (j=0; j<numvar; j++) {
#ifdef USE_INTERVAL
      low = roots[i][j].get_low();
      upp = roots[i][j].get_upp();
#else
      low = root[i][j].center();
      if (low > 1.0)
	    low = 1.0;
      else if (low < 0.0)
	    low = 0.0;
	else ;
#endif
      
#ifdef USE_INTERVAL
      real temp(low, upp);
#else
      real temp = low;
#endif
      tmpRoots[i][j] = temp;
    }

  for (i=0; i<numRoots; i++)
    delete [] roots[i];
  delete [] roots;

  roots = tmpRoots;

  return;
}











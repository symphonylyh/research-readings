// Copyright (C) Massachusetts Institute of Technology, 1995
// All rights reserved

#include "simpoly.h"
using namespace std;

typedef real *realptr;

void sim_solve (real **bp, short **ordlists, int numeq, real eps,
		real **&roots, int &numroots)

/* A call to si_pinter, masked */

{
  int i,j,k;
  mn_array mna (numeq);
  for (i=0; i < numeq; i++) {
    sht_array ol(numeq);
    for (j=0; j < numeq; j++)
      ol[j] = ordlists[i][j]-1;
    multinomial *mn = new multinomial (ol);

    for(k=0;k<mn->cosize;k++)
      mn->bp[k] = bp[i][k];
    //    memcpy ((real *) mn->bp, bp[i], mn->cosize*sizeof(real));
    ((multinomial **) mna)[i] = mn;
    //    cout << mna[i];
  }
  rootlist *l = si_pinter (mna, eps);
  numroots = 0;
  l->head();
  numroots = l->size();
  roots = new realptr[numroots];
  for (i=0; i < numroots; i++) {
    roots[i] = new real[numeq];
    for (j=0; j < numeq; j++) {
      roots[i][j] = l->root()[j];
      //      cout << roots[i][j] << ' ';
    }
    //    cout << endl;
    l->next();
  }
}


void sim_solve (real **bp, short **ordlists, int numeq, int numvar, real eps,
		real ***roots, int *numroots)

/* A call to si_pinter, masked */

{
  int i,j,k;
  mn_array mna (numeq);
  for (i=0; i < numeq; i++) {
    sht_array ol(numvar);
    for (j=0; j < numvar; j++)
      ol[j] = ordlists[i][j]-1;
    multinomial *mn = new multinomial (ol);
    for(k=0;k<mn->cosize;k++)
      mn->bp[k] = bp[i][k];
    //    memcpy ((real *) mn->bp, bp[i], mn->cosize*sizeof(real));
    ((multinomial **) mna)[i] = mn;
    //   cout << mna[i];
  }
  //  cout << "End of mna...\n";cout.flush();
  rootlist *l = si_pinter (mna, eps);
  //  cout << "End of si_pinter...\n";cout.flush();
  *numroots = 0;
  l->head();
  (*numroots) = l->size();
  (*roots) = new realptr[*numroots];
  for (i=0; i < *numroots; i++) {
    (*roots)[i] = new real[numeq];
    for (j=0; j < numeq; j++) {
      (*roots)[i][j] = l->root()[j];
   //   cout << (*roots)[i][j] << ' ';
    }
//    cout << endl;
    l->next();
  }
}

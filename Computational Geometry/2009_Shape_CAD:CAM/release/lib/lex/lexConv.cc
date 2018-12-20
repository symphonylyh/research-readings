//#include <stdio.h>
#include <iostream>
#include <iomanip>
#include "lexPP.h"
using namespace std;


FILE *fn;
#define _cplusplus

int lexConv(char *filename, 
	    real **&bp, short **&ordlists, int &numeq, int &numvar)
{
  int i,j; 
  mn_array *mn_list;


  fn = fopen(filename, "r");
  lexCode(mn_list);

  //  numeq = ((*mn_list)[0]).ndim();
    numeq = (*mn_list)[0].ndim();
  for(i=0; i<numeq; i++){
    ((*mn_list)[i]).monotobern(&((*mn_list)[i]));
  }

  bp =  new real*[numeq];
  for(i=0; i<numeq; i++){
    int cosize =( (*mn_list)[i] ).cosize;
    bp[i] = new real [cosize];
    for(j=0; j<cosize; j++){
      bp[i][j] = ((*mn_list)[i]).bp[j];
    }
  }

  numvar = numeq;
  ordlists = new short*[numeq];
  for(i=0; i<numeq; i++){
    ordlists[i] = new short [numvar];
    for(j=0; j<numvar; j++){
      ordlists[i][j] =  ((*mn_list)[i]).ordlist[j]+1;
    }
  }


  return 0;
}


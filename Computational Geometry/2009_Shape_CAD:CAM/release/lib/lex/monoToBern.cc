//#include <stdio.h>
#include <iostream>
#include <iomanip>
#include "monoPP.h"

using namespace std;

void monoToBern(real **bp, short **ordlists, int numeq, int numvar)
{
  int i,j, print=0; 
  multinomial *mn;
  sht_array ol(numeq);


  for(i=0; i<numeq; i++){
    for(j=0; j<numvar; j++){
      ol[j] = ordlists[i][j];
    }

    mn = new multinomial (ol);
    int cosize=mn->cosize;
    for(j=0; j<cosize; j++)     
      mn->bp[j] = bp[i][j];

    if(print==1)
      cout<<*mn<<endl;

    mn->monotobern(mn);

    if(print==1)
      cout<<*mn<<endl;

    for(j=0; j<numvar; j++)
      ordlists[i][j] = mn->ordlist[j]+1;
    

    for(j=0; j<cosize; j++)
      bp[i][j] = mn->bp[j];
    
    delete mn;
  }

}


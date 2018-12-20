
//extern "C"{
//  #include <stdlib.h>
//}


#include "simpoly.h"
using namespace std;


#define GREATER 1
#define LESS -1 


int WHICH_CORDNT_TO_SORT;
int compar_rlist_element(const void *, const void *);

ostream& operator<<(ostream& output, rlist_element*& a)
{
     output << a->u ;
     return output;
}

ostream& operator<<(ostream& output, rlist_element& a)
{
     output << a.u ;
     return output;
}

ostream& operator<<(ostream& output, rootlist* a)
{
     a->clp = a->hd;
     int num = 0 ;

     while (a->clp) {
	  output << a->clp;
	  a->next();
	 num++;
     }

  
     output << " root number is " << a->size() << endl;
     return output;

}


ostream& operator<<(ostream& output, rootlist& a)
{
     a.clp = a.hd;
//     int num = 0 ;

/*     while (a.clp) {
	  output << a.clp;
	  a.next();
	 num++;
     }
  */
     output << " root number is " << a.size() << endl;
     return output;

}


ostream& operator<<(ostream& output, bagroot_list*& a)
{
     a->clp = a->hd;
     int num = 0 ;

     while (a->clp) {
	  bagrootlist_element* brle;

	  brle =   a->curent_pointer();
	  output << brle;
	  a->next();
	 num++;
     }  
     output << " bagroot number is " << num << endl;
     return output;

}


void integerlist::insert (int item) /* Insert item after clp */
{ integer_element *el = new integer_element(item);
  if (!clp) {
       hd = clp = el;
       clp->next = 0;
  } else {
       el->next = clp->next;	/* Set next pointer of el */
       clp->next = el;		/* Set next pointer of clp to el */
       clp = el;		/* Clp now points to just inserted el */
  }
  n++; 
}


ostream& operator<<(ostream& output , integerlist* &il) {
     il->head();
     int num = 0;
     while(il->clp) {
	  output <<(*il->clp);
	  il->next();
	  num++; 
     }
     
     output << " region number is " << num<< endl;
     return output;
}


rlist_element::rlist_element(int siz):u(siz){ next = NULL;}
rlist_element::rlist_element (const real_array &item) : u(item) {next = NULL;} 

real& rlist_element::operator[](int i){
     if (i >= u.size()) {
	  cout << " error: in rlist_element[]" << endl;
	  exit(1);
     }
     return u[i];
}

real rlist_element::operator()(int i){
     if (i >= u.size()) {
	  cout << " error: in rlist_element()" << endl;
	  exit(1);
     }
     return u[i];
}

rlist_element& rlist_element::operator=(rlist_element& rl){

  u.n = rl.u.size();
  delete u.p;
  u.p = new real[u.n];
  //  u.p = realloc(u.p, u.n * sizeof(real));


  for(int i=0;i<u.n;i++) {
    u.p[i] = rl.u[i];
  }
  //  memcpy(u.p, (real *) (rl.u), u.n * sizeof(real));
  next = rl.next;

  return *this;
}


root_array::root_array( rootlist &rl) 
{
  int i;
     n = rl.size();
     rl.head();
     p = alloc(n * sizeof(rlist_element));
     for (i=0; i< n; i++) {  // initialize each rlist_element
	  (((rlist_element *) p)[i]).u.init_size();
	  (((rlist_element *) p)[i]).u.init_pointer();
	  
     }
     for ( i=0; i < n; i++) {
	  ((rlist_element *) p)[i] = (*(rl.curent_pointer()));
	  rl.next();	       
     }
     
}

rlist_element &root_array::operator[](int i) 
{
     if (i >= n) {
	  cout << " error: in rlist_element[]" << endl;
	  exit(1);
     }
     return (((rlist_element *)p)[i]);
}

ostream& operator<<(ostream& output, root_array& a)
{
  int num;
     output << " in cout << root_array " << endl;
     for (num=0; num < a.size() ; num++) {
	  
	  output << a[num];
     }

     output << " root number is " << num << endl;
     return output;

}




void root_array::init(int siz) {
     p = alloc(siz * sizeof(rlist_element)); 
     n=siz;
     for (int i=0; i< n; i++) {  // initialize each rlist_element
	  (((rlist_element *) p)[i]).u.init_size();
	  (((rlist_element *) p)[i]).u.init_pointer();
	  
	  }
}


root_array* rootlist::sort_wrt_cordnt(int cordnt)
{
  WHICH_CORDNT_TO_SORT = cordnt;

  root_array*  ra = new root_array(*this);

  qsort( ((rlist_element *) (*ra)), n, sizeof(rlist_element), 
	compar_rlist_element); 
  /* p is the pointer to the array of points */
  /* n is the number of elements in p */

  
  return ra;

}



int compar_rlist_element(const void *pa, const void *pb)
{

     rlist_element* p1;
     rlist_element* p2;
     
     p1 = (rlist_element *) pa;
     p2 = (rlist_element *) pb;

  if ( WHICH_CORDNT_TO_SORT >= p1->size() ||
      WHICH_CORDNT_TO_SORT >= p2->size() ) {
    cout << " running error: cordnt to compare > size ";
    cout << "in compare_rlist_element "<< endl;
    exit(1);
  }  
  
     real p1_int, p2_int;

     p1_int = (*p1)[WHICH_CORDNT_TO_SORT];
     p2_int = (*p2)[WHICH_CORDNT_TO_SORT];
     
#ifdef USE_INTERVAL
  if (p1_int.get_low() <
      p2_int.get_low())

      return LESS;

  else if (p1_int.get_upp() >
	   p2_int.get_upp())
  
    return GREATER;

  return LESS;

#else

  if (p1_int < p2_int)

      return LESS;

  else if (p1_int > p2_int)
  
    return GREATER;

  return LESS;
#endif
}


void rootlist_array::init(int siz) {
     p = alloc(siz * sizeof(rootlist)); 
     n=siz;
     for (int i=0; i< n; i++) {  // initialize each rootlist
	  (((rootlist *) p)[i]).hd = NULL;
	  (((rootlist *) p)[i]).clp = NULL;
	  (((rootlist *) p)[i]).n = 0;
	  	  
	  }
}

rootlist_array::~rootlist_array() 
{
     for (int i=0; i < n; i++) {
	  rootlist *rl = &(((rootlist *) p)[i]);
	  rl->clp = rl->hd;	/* Start at beginning */
	  while (rl->clp) {	  /* While elements left */
	       delete (rl->clp);
	       rl->clp = rl->clp->next;
	  }
     }
     dealloc();
}


ostream& operator<<(ostream& output, rootlist_array& a)
{
  int num;
     output << " in cout << rootlist_array " << endl;
     for (num=0; num < a.size() ; num++)
        output << a[num];


     output << " root bag number is " << num << endl;
     return output;

}


rootlist::~rootlist()			/* Free list; frees EVERYTHING */
{ // cout << " Hi I am in rootlist deletor " << endl;
     clp = hd;			/* Start at beginning */
     while (clp) {		/* While elements left */
	  delete clp;
	  clp = clp->next;
     }
     
}

/* return the i-th value of */
/* the current element.     */
real rootlist::operator()(int i )
{
 return (*clp)(i);
 } 


void rootlist::insert (real_array &item) /* Insert item after clp */
{
     rlist_element *el =
	  new rlist_element (item);
     if (!clp) {
	  hd = clp = el;
	  clp->next = 0;
     } else {
	  el->next = clp->next; /* Set next pointer of el */
	  clp->next = el;	/* Set next pointer of clp to el */
	  clp = el;		/* Clp now points to just inserted el */
     }
     n++; 
}


bagrootlist_element::bagrootlist_element(int siz):rlist_element(siz){}

bagrootlist_element::~bagrootlist_element(){ }

bagroot_list::~bagroot_list()
{
//     cout << " in bagroot_list deletor " << endl;
     clp = hd;			/* Start at beginning */
     while (clp) {		/* While elements left */
	  delete (bagrootlist_element *) clp;
	  clp = clp->next;
     }
}

root_array::operator rlist_element* () const { /* Conversion operator */
     return ((rlist_element *) p);
}

ostream& operator<<(ostream& output, bagrootlist_element* b)
{
     output << " the covered region is " << endl;
     output << b->u << endl;
//     output << "root patches are" << endl;
//     output << b->bucket << endl;

     return output;
}


void bagroot_list::insert(bagrootlist_element* &el) {
     /* before call this routine 
	  make sure 'el' is created by new */
     if (!clp) {
	  hd = clp = el;
	  clp->next = 0;
     } else {
	  el->next = clp->next;	/* Set next pointer of el */
	  clp->next = el;		/* Set next pointer of clp to el */
	  clp = el;		/* Clp now points to just inserted el */
     }
     n++; 
}  

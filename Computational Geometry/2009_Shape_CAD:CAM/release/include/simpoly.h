// Copyright (C) 1995 Massachusetts Institute of Technology, Cambridge, MA
// All rights reserved

#ifndef SIMPOLY_H		/* Header file for simultaneous poly. */
#define SIMPOLY_H

#include <iostream>
#include "multinom.h"
#include "bbox.h"

struct stackel {
  void *item;
};

class rootlist;    /* Forward declaration */
class root_array;  /* Forward declaration */

class Rstack {

  stackel *base;                /* Base of the stack */
  int maxsiz;                   /* Maximum size it can attain */
  int n;                        /* Number of elements */

public:

  Rstack();                     /* Constructor */
  ~Rstack();                    /* Destructor */
  void push(void *);            /* Generic push and pop */
  int pop(void **);
  void push(Bbox *);
  int pop(Bbox **);
};

class rlist_element {
 protected:
  real_array u;			/* N-dimensional coordinate */
public:
  rlist_element() {}
  rlist_element *next;		/* Pointer to next element */
  rlist_element(int);
  rlist_element (const real_array &);
  ~rlist_element () {/* cout <<"In rlist_element deletor" << endl; */}

  real_array root () 	/* Abstraction; returns copy of u */
    { return u; }
  real& operator[](int );
  real operator()(int );
  
  rlist_element& operator=(rlist_element&);
  friend int operator == (rlist_element& a, rlist_element& b) {
#ifndef USE_RAT
       return ( a.u == b.u);
#else
       return ((Rat_array)a.u == (Rat_array)b.u);
#endif
  }

  void modify_root( real_array &ra) { u = ra;}


  friend int if_overlap(rlist_element* a, rlist_element* b)
  { return if_overlap( a->u, b->u);}

  int size() {return u.size();}
  friend std::ostream& operator<<(std::ostream&, rlist_element* &);
  friend std::ostream& operator<<(std::ostream&, rlist_element &);

  friend class rootlist;
  int if_rlist_elem_in_rootlist(rootlist* );
  friend class root_array;
};

class rootlist {
 protected:
     rlist_element *hd;		/* Head */
     rlist_element *clp;		/* Current list pointer */
     int n;			/* Number of elements */
 public:
     rootlist() {hd = NULL; clp= NULL; n=0;}
     /* Make a new null list */
     ~rootlist();
     void head()  { clp = hd; /* Move clp to head */ }

     void next()      { clp = clp->next; }
     void insert (real_array &); /* Insert item after clp */
     real operator()(int i ); 
     /* return the i-th value of */
     /* the current element.     */
     
     real_array root () 	/* Return root at clp */
     { return (clp->root()); }
     int size () const
     { return n; }
     rlist_element* curent_pointer() const {return clp;}
     root_array* sort_wrt_cordnt(int);
     
     friend class rootlist_array;
  friend std::ostream& operator<<(std::ostream& , rootlist* );
  friend std::ostream& operator<<(std::ostream& , rootlist &);
};

class bagrootlist_element: public rlist_element {
 protected:
     rootlist bucket;
 public:
     bagrootlist_element(int);
     ~bagrootlist_element();

     rootlist& bag(){return bucket;}

   friend std::ostream& operator<<(std::ostream&, bagrootlist_element* );
};

class bagroot_list: public rootlist {
     
 public:
     bagroot_list(): rootlist(){}
     ~bagroot_list();
     
     void insert(bagrootlist_element* &); 
     /* before call this routine 
	make sure 'el' is created
	by new */

     bagrootlist_element* curent_pointer() const
     {
	  return ((bagrootlist_element *) clp);
     }
     
     friend std::ostream& operator<<(std::ostream& output, bagroot_list*&);
};

/* array of rlist_element */
class root_array: public generic_array<rlist_element> {
 protected:
     root_array(){ n=0; p=NULL;}
     void init(int);
 public:
     root_array(int siz) { init(siz);}
     root_array(const root_array &ra) {
	  int siz = ra.size();
	  p = alloc(siz);
	  n = siz;
	  for(int i=0;i<n;i++)
	    p[i] = ra[i];
	  //	  memcpy(p, (rlist_element *)ra, n*sizeof(rlist_element));
     }
     root_array( rootlist &);

     ~root_array() {dealloc();}
     operator rlist_element* () const;  /* Conversion operator */

     
     rlist_element &operator[](int );
     void resize (int siz) {p=ralloc(siz);n=siz;}

  friend std::ostream& operator<<(std::ostream& , root_array &);
};

class integer_element {
  int u; 
public:
  integer_element *next;		/* Pointer to next element */
  integer_element (int &item) : u(item), next(0) {} /* Constructor */
  ~integer_element () {};		/* Destructor */
  int value () const	/* Abstraction; returns copy of u */
    { return u; }
  int& operator()(){
       return u;
  }
  
  integer_element& operator=(integer_element& ie) {
       next = ie.next;
       u = ie.u;
       return *this;
  }

  friend std::ostream& operator<<(std::ostream& output, integer_element &ie){
       output << " u = " << ie.u << std::endl;
       return output;
  }
};

class integerlist {
  integer_element *hd;		/* Head */
  integer_element *clp;		/* Current list pointer */
  int n;			/* Number of elements */
public:
  integerlist() :hd(0),clp(0),n(0) {} /* Make a new null list */
  ~integerlist()		      /* Free list; frees EVERYTHING */
    { // cout << " Hi I am in integerlist deletor " << endl;
	 clp = hd;			/* Start at beginning */
      while (clp) {		/* While elements left */
	delete clp;
	clp = clp->next;
      }
    }
  void head()			/* Move clp to head */
    { clp = hd; }
  void next()
    { clp = clp->next; }
  void insert (int ); /* Insert item after clp */

  int value () const	/* Return value at clp */
    { return (clp->value()); }

  int size () const
    { return n; }
  integer_element* curent_pointer() const {return clp;}
  friend std::ostream& operator<<(std::ostream& , integerlist* &);
};

class rootlist_array:  public generic_array<rootlist> {
 protected:
     rootlist_array(){ n=0; p=NULL;}
     void init(int);
 public:
     rootlist_array(int siz) { init(siz);}
     ~rootlist_array();
     operator rootlist * () const { /* Conversion operator */
	  return ((rootlist *) p);
     }
     
     rootlist &operator[](int i) {
	  if (i >= n) {
	       std::cout << " error: in rootlist_array[]" << std::endl;
	       exit(1);
	  }
	  return (((rootlist *)p)[i]);
     }
     void resize (int siz) {p=ralloc(siz);n=siz;}

  friend std::ostream& operator<<(std::ostream& , rootlist_array &);
};

extern "C" {
  void e04kcf_(int&, int&, real *, real *, real *, real&, real *,
	       int *, int&, real *, int&, int&);
/*  void funct2_ (int *, real *, real *, real *); */
  void e04mbf_(int&, int&, int&, int&, int&, int&, real *,
	       real *, real *, real *, int&, real *, int *,
	       real&, real *, int *, int&, real *, int&, int&);
  void simplx (real **, int, int, int, int, int, int *, int *, int *);
}

void sim_solve(real**, short**, int, real, real***, int*);
void sim_solve(real**, short**, int, real, real**&, int&);
rootlist *si_pinter (mn_array &, real);
rootlist *si_pinter_more_unknown_than_eq(mn_array &, real, Bbox*);
rootlist *si_pinter_more_eq_than_unknown(mn_array &, real, Bbox*);
rootlist *si_linter (mn_array &mn_list, real);
rootlist *si_hybinter (mn_array &mn_list, real, real);

void si_lpinit (mn_array &);
void si_lp (mn_array &, Bbox &, int_array &);
void si_phull (mn_array &, Bbox &, int_array &);
void si_multipush (Rstack &, int_array &, Bbox &);
void si_sub (mn_array &, mn_array &, Bbox &);

void si_ginit (int);
void si_show_graph (real, real, real, real, int);
void si_show_graph_2 (real, real);
void si_reset_graph (int);
int if_overlap(double , double );
bagroot_list* separate_roots(rootlist *, real);
bagroot_list* separate_roots(rootlist *);
rootlist_array *bucket_each_coord(root_array*, int, real);
rootlist_array *bucket_each_coord(root_array*, int);

extern int si_graph;

inline Rstack::Rstack ()

/* Constructor for stacks */

{
  base = (stackel *)            /* Allocate the base */
    malloc (sizeof (stackel));
  maxsiz = 1;                   /* Initial max size */
  n = 0;                        /* Empty */
}

inline Rstack::~Rstack ()

/* Destructor */

{
  free ((void *) base);
}

inline void Rstack::push (void *p)

/* Generic push of one item */

{
  if (n >= maxsiz) {
    maxsiz *= 2;                /* Resize it */
    base = (stackel *)
      realloc (base,
	       maxsiz * sizeof (stackel));
  }
  base[n].item = p;                     /* Store the argument */
  n++;
}

inline int Rstack::pop (void **p)

/* Generic pop of one item */

{
  if (n==0)
    return (0);
  else {
    n--;
    *p = base[n].item;
    return (1);
  }
}

inline void Rstack::push (Bbox *B)

/* Push a Bbox */

{
  push ((void *) B);
}

inline int Rstack::pop (Bbox **B)

/* Pop a Bbox */

{
  return (pop ((void **) B));
}

#define LP                              /* Defined to use LP boxes */
#endif

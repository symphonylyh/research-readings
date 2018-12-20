// Copyright (C) 1995 Massachusetts Institute of Technology, Cambridge, MA
// All rights reserved

#ifndef LIB1_H
#define LIB1_H
#ifdef USE_INTERVAL
#include "interval.h"
#else
#ifdef USE_RAT
#include "rat.h"
#endif
#endif
//#include <stdlib.h>             /* For malloc, free */
#include <string>             /* For memcpy */
#include <iostream>

template<class base_type>class generic_array {           /* This is the class of generic arrays */
protected:
  base_type *p;    
  int n;           
  base_type *alloc (int siz) { return new base_type[siz]; }
  void dealloc () { if(p==NULL) return; else delete[] p; }

  base_type *ralloc (int siz) { 
    register int i;
    base_type *temp;
    siz = siz > n ? siz : n;
    temp = new base_type[siz];
    for(i = 0;i < n;i++)
      temp[i] = p[i];
    delete[] p;
    return temp;
  }
 
  base_type *zalloc (int siz)        /* Do a calloc */
    {return new base_type[siz]; }
public:
  int size () const { return n; }    /* Size of array */
  void init_size() {n= 0;}
  void init_pointer(){p = NULL;}
};

class dbl_array : public generic_array<double> {   /* A double array */
protected:
  dbl_array () {}               /* Default constructor for inheritance */
  void init (int);           /* Initializer */

public:
  dbl_array (int );           /* Create a double array */
    
  dbl_array (const dbl_array &);   /* Copy constructor */
    
  ~dbl_array ()
    { dealloc (); }
  double &operator[] (int ) const ;  /* Indexer */
   

  void resize (int siz)
    { p = ralloc (siz); n = siz; }

  dbl_array operator=(dbl_array &ia)
  {
       if (n == ia.n) memcpy (p, (double *) ia.p, n*sizeof(double));
       else {
	    resize(ia.n);
	    memcpy (p, (double *) ia.p, n*sizeof(double));
       }
       return *this;
  }

  friend int operator == (dbl_array &, dbl_array &);
};

class sht_array : public generic_array<short> {   /* A short array */
protected:
  sht_array () {}               /* Default constructor for inheritance */
  void init (int);           /* Initializer */
public:
  sht_array (int siz)           /* Create a double array */
    { init (siz); }
  sht_array (const sht_array &d)   /* Copy constructor */
    { init (d.size()); 
    for(int i=0;i<n;i++)
      p[i]=d.p[i];
    //memcpy (p, (short *) d.p, n*sizeof(short)); }
    }
  ~sht_array ()
    { dealloc (); }
  short &operator[] (int i) const   /* Indexer */
    { return (((short *) p)[i]); }
  void resize (int siz)
    { p = ralloc (siz); n = siz; }
};

class int_array : public generic_array<int> {   /* An integer array */
protected:
  int_array () {}               /* Default constructor for inheritance */
  void init (int);           /* Initializer */
public:
  int_array (int siz)           /* Create a double array */
    { init (siz); }
  int_array (const int_array &d)   /* Copy constructor */
    { init (d.size()); 
    for(int i=0;i<n;i++)
      p[i] = d.p[i];
    //memcpy (p, (int *) d.p, n*sizeof(int)); }
    }
  ~int_array ()
    { dealloc (); }
  int &operator[] (int i) const   /* Indexer */
    { return (((int *) p)[i]); }

  void resize (int siz)
    { p = ralloc (siz); n = siz; }
};

#ifdef USE_RAT
class Rat_array : public generic_array<Rat> {   /* A rational array */
protected:
  Rat_array () {}               /* Default constructor for inheritance */
  void init (int siz)           /* Initializer */
    { p = zalloc (siz); n = siz; }
public:
  Rat_array (int siz)           /* Create a double array */
    { init (siz); }
  Rat_array (const Rat_array &d)   /* Copy constructor */
    { init (d.size()); for (int i=0; i < n; i++) p[i] = d.p[i]; }
  ~Rat_array ()
    { 
      //            for (int i=0; i < n; i++)
      //      	(*this)[i].Rat::~Rat(); 
    dealloc (); }
  Rat &operator[] (int i) const   /* Indexer */
    { return (((Rat *) p)[i]); }
  operator Rat* () const     /* Conversion operator */
    { return ((Rat *) p); }
  void resize (int siz) 
    { 
      p = ralloc(siz);n = siz;
    }
 Rat_array operator=(Rat_array &ia)
  {
    int i;
    if (n == ia.n) {
      for(i=0;i<n;i++)
	p[i] = ia.p[i];
    } else {
      resize(ia.n);
      for(i=0;i<n;i++)
	p[i] = ia.p[i];
    }


    return *this;
  }
};
#endif

#ifdef USE_INTERVAL
class Interval_array : public generic_array<Interval> {   /* A Interval array */
protected:
  Interval_array () {}              /* Default constructor for inheritance */
  void init (int);           /* Initializer */
public:
  Interval_array (int);
  Interval_array (const Interval_array &d)   /* Copy constructor */
    { init (d.size()); 
    for(int i=0;i<n;i++)
      p[i] = d.p[i];
    //memcpy (p, (Interval *) d.p, n*sizeof(Interval)); }
    }
  ~Interval_array ()
    { dealloc (); }
   Interval &operator[] (int i) const  /* Indexer */
   {
     return p[i]; }

    operator Interval* () const     /* Conversion operator */
     { return (Interval *) p; }
  void resize (int siz)
    { p = ralloc (siz); n = siz; }
  Interval_array operator=(Interval_array &ia)
  {
    int i;
    if (n == ia.n) {
      for(i=0;i<n;i++)
	p[i] = ia.p[i];
    } else {
      resize(ia.n);
      for(i=0;i<n;i++)
	p[i] = ia.p[i];
    }


      /*
       if (n == ia.n) memcpy (p, (Interval *) ia.p, n*sizeof(Interval));
       else {
	    resize(ia.n);
	    memcpy (p, (Interval *) ia.p, n*sizeof(Interval));
       }
      */

    return *this;
  }

  friend int operator == (Interval_array &, Interval_array &);
  
};
#endif


template<class base_type>class generic_array2 {          /* This is the class of 2D generic arrays */
protected:
  base_type **p;                     /* Pointer to data */
  int m,n;                      /* Number of elements in each direction */
  base_type *alloc (int siz)         /* Allocate n bytes */
    {return new base_type*[siz];}
  void dealloc ()               /* Deallocate p and p[0] */
    { int i,j;
    for(i=0;i<m;i++)
      delete[] p[i];
    delete[] p;
    }
  //    { free ((char *) p[0]); free ((char *) p); }
  //  void *zalloc (int siz)        /* Do a calloc */
  //    { return (calloc (siz, 1)); }
public:
  int rows () const             /* # rows */
    { return m; }
  int cols () const             /* # cols */
    { return n; }
};

class double_array2: public generic_array2<double> {       /* 2-D double arrays */
protected:
  void init (int , int );      /* Initializer */
  void init (int , int ,int );   /* Initializer for array of 0's */  

  double_array2 () {}              /* Default constructor -- inheritance only */
public:
  double_array2 (int r,int c)      /* Create a 2-D double array */
    { init (r, c); }
  double_array2 (int r,int c,int z)/* Create a 2-D double array of 0's */
    { init (r, c, z); }
  ~double_array2 ()                /* Destructor */
    { dealloc (); }
  double* &operator[] (int) const;
  operator double** () const     /* Conversion operator */
    { return ((double **) p); }
};


#ifdef USE_INTERVAL
class Interval_array2: public generic_array2<Interval> {       /* 2-D Interval arrays */
protected:
  void init (int , int );      /* Initializer */
  void init (int , int , int );   /* Initializer for array of 0's */  
  Interval_array2 () {}         /* Default constructor -- inheritance only */
public:
  Interval_array2 (int r,int c)      /* Create a 2-D Interval array */
    { init (r, c); }
  Interval_array2 (int r,int c,int z)/* Create a 2-D Interval array of 0's */
    { init (r, c, z); }
  Interval_array2 (const Interval_array2 &d)
    { init (d.rows(), d.cols());
    int i,j;
    for(i=0;i<d.rows();i++) {
      for(j=0;j<d.cols();j++) {
	p[i][j] = d[i][j];
      }
    }
    }
  ~Interval_array2 ()                /* Destructor */
    { dealloc (); }
  Interval* &operator[] (int) const;
  operator Interval** () const     /* Conversion operator */
    { return ((Interval **) p); }
};
#endif



#ifdef USE_RAT

class Rat_array2: public generic_array2<Rat> {     
protected:
  void init (int r, int c);     
  /*    { p = (void **) alloc      
	(r * sizeof (Rat *));
      p[0] = (Rat *) zalloc
	(r * c * sizeof (Rat));
      for (int i=1; i < r; i++)
	p[i] = (void *)    
	  ((Rat *) p[i-1] + c);
	  m = r; n = c; }    */
  Rat_array2 () {}         
public:
  Rat_array2 (int r,int c) 
    { init (r, c); }
  Rat_array2 (const Rat_array2 &d)
    { init (d.rows(), d.cols());
    /*      Rat *r1 = (Rat *) p[0];
      Rat *r2 = (Rat *) &d[0][0];
      for (int i=0; i < m*n; i++, r1++, r2++)
      *r1 = *r2; */
    register int i,j;
    for(i=0;i<d.rows();i++)
      for(j=0;j<d.cols();j++)
	p[i][j] = d.p[i][j];
    }
  ~Rat_array2 ()           
    { 
      //Rat *q = (Rat *) p[0];
      //      for (int i=0; i < m*n; i++, q++)
      //	(*q).Rat::~Rat();
      dealloc (); }
  Rat* &operator[] (int i) const
    { return (((Rat**) p)[i]); }
  operator Rat** () const  
    { return ((Rat **) p); }
    };
#endif




#if 0
class generic_narray {          /* A generic array in n dimensions */
protected:
  void *p;                      /* The data */
  int n;                        /* Number of dimensions */
  unsigned short *m;            /* Array of dimension sizes */
  void *alloc (int siz)         /* Allocate n bytes */
    { return (malloc (siz)); }
  void dealloc (void *q)        /* Deallocate p */
    { free (q); }
  void *ralloc (int siz)        /* Reallocate p */
    { return (realloc (p, siz)); }
  void *zalloc (int siz)        /* Do a calloc */
    { return (calloc (siz, 1)); }
};

class dbl_narray: public generic_narray {       /* Double N-dim arrays */
  int size;                     /* This is internally useful */
public:
  dbl_narray (int nd, unsigned short *md);
  dbl_narray (int nd, ...);
  ~dbl_narray ();
  double &operator () (int i, ...);
};
#endif
#endif


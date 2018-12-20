// Copyright (C) Massachusetts Institute of Technology, 1995
// All rights reserved

// geom.h

#ifndef GEOM_H
#define GEOM_H

#include <iostream>
#include <iomanip>
#include "gen.h"

class geom {
 public:
  geom(int);                       // class constructor
  geom(const int, int *);
  geom(int, int *, double *, vector **);
  geom(geom &);
  ~geom();                                // class destructor

  friend void     add           (geom &, geom &, geom &);
  friend geom    *deri_geom     (geom *, int);
  friend geom   **dist_obj      (geom *); 
  friend vector **value_of_geom (geom *);
  friend void     make_geom     (geom &, geom &);
  friend void     sub           (geom &, geom &, geom &);

  friend geom *operator + (geom &, geom &);
  friend geom *operator * (geom &, geom &);

  geom  operator -= (geom &);
  geom  operator += (geom &);
  geom  operator =  (geom &);
  geom *operator =  (geom *);
  
  friend ostream& operator << (ostream & , const geom & );

  double *knots;           // array of knots
  int nd, nc, nk;          // num of dimensions
  int *ncontpts;           // list of num of contpts
  int *order;              // list of orders
  vector **contpts;        // array of points of contpts
};

void   add          (geom &, geom &, geom &);
void   add_base     (int, int *, int *);
double combine_cc   (int, int);
void   copy_array   (int, int *, int *);
geom  *deri_geom    (geom *, int);
double factorial_cc (int);
double get_coef0    (int, int *, int *);
double get_coef1    (int, int *, int *, int *, int *);
int    get_MIN_MAX  (int, int *, int *, int *, int *, int *);
int    get_total    (int, int *);
void   hash         (int, int, int *, int *);
void   rhash        (int, int, int *, int *);
int    rsort        (int, int *, int *);
int    sort         (int, int *, int *);
void   sub          (geom &, geom &, geom &);
void   sub_array    (int, int *, int *, int *);
void   sub_array1   (int, int *, int *, int *);

#endif

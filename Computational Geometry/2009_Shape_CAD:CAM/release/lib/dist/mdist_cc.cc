// Copyright (C) Massachusetts Institute of Technology, 1997
// All rights reserved

#include <iostream>
using namespace std;

#include <math.h>
#include "bspl.h"
#include "gen.h"
#include "geom.h"
#include "mdist_cc.h"
#include "multinom.h"
#include "simpoly.h"

// #define USE_INTERVAL

/////////////////////////////////////////////////////////////////////////
// C functions 
/////////////////////////////////////////////////////////////////////////

double pointToCurv(double x, double y, double z, ParCurv *egeom, double eps,
		   double *u0) 
{
  return robustMinDistance(x, y, z, egeom, eps, u0);
}

void pointToCurv_n(double x, double y, double z, ParCurv *egeom, double eps,
		   double *u0, double *d0) 
{
  *d0 = robustMinDistance(x, y, z, egeom, eps, u0);
}

double pointToSurf(double x, double y, double z, ParSurf *fgeom, double eps,
		   double *u0, double *v0)
{
  return robustMinDistance(x, y, z, fgeom, eps, u0, v0);
}

void pointToSurf_n(double x, double y, double z, ParSurf *fgeom, double eps,
		   double *u0, double *v0, double *d0)
{
  *d0 = robustMinDistance(x, y, z, fgeom, eps, u0, v0);
}

double curvToCurv(ParCurv *egeom1, ParCurv *egeom2, double eps, double *u1,
                  double *u2)
{
  return robustMinDistance(egeom1, egeom2, eps, u1, u2);
}

void curvToCurv_n(ParCurv *egeom1, ParCurv *egeom2, double eps, double *u1,
                  double *u2, double *d0)
{
  *d0 = robustMinDistance(egeom1, egeom2, eps, u1, u2);
}

double curvToSurf(ParCurv *egeom, ParSurf *fgeom, double eps, double *u0,
                  double *u1, double *v1)
{
  return robustMinDistance(egeom, fgeom, eps, u0, u1, v1);
}

void curvToSurf_n(ParCurv *egeom, ParSurf *fgeom, double eps, double *u0,
                  double *u1, double *v1, double *d0)
{
  *d0 = robustMinDistance(egeom, fgeom, eps, u0, u1, v1);
}

double surfToSurf(ParSurf *fgeom1, ParSurf *fgeom2, double eps, double *u1,
		  double *v1, double *u2, double *v2)
{
  return robustMinDistance(fgeom1, fgeom2, eps, u1, v1, u2, v2);
}

void surfToSurf_n(ParSurf *fgeom1, ParSurf *fgeom2, double eps, double *u1,
		  double *v1, double *u2, double *v2, double *d0)
{
  *d0 = robustMinDistance(fgeom1, fgeom2, eps, u1, v1, u2, v2);
}

/////////////////////////////////////////////////////////////////////////
// Point to curve
/////////////////////////////////////////////////////////////////////////

double robustMinDistance(double x, double y, double z, ParCurv *egeom,
			 double eps, double *u0)
{
  geom *a = dist_bez_pc(x, y, z, egeom);    // form equations
  rootlist *l = solveMinDistance(a, eps);   // solve equations
  delete a;                                 // delete geom

  vector *n, *p = vectalloc(), *q, *r;
  p->x = x;
  p->y = y;
  p->z = z;
  p->w = 1.0;

  *u0 = 0.0;
  double dist, d0 = 1.0e30, u;
  l->head();                                // starting at first root

  for (int i=0; i<l->size(); i++) {         // for each root in list
    real_array t = l->root();               // get current root

#ifndef USE_INTERVAL
    u = t[0];
#else
    u = t[0].center();
#endif

    if (u < egeom->knots[egeom->order-1])   // insure that root is inside 
      u = egeom->knots[egeom->order-1];     // the knot range
    if (u > egeom->knots[egeom->ncontpts])
      u = egeom->knots[egeom->ncontpts];

    r = evalbsp(egeom, u);                  // evaluate curve
    n = normalcurv(egeom, u, 1);            // curve normal

    dist = distance(p, r);                  // sqrt(dist^2), i.e. |dist|
    q = sub_vect(p, r);
    if (dot(q, n) < 0.0)                    // add sign to absolute value
      dist *= -1.0;

    if (fabs(dist) < fabs(d0)) {            // update minimum distance
      d0 = dist;
      *u0 = u;                              // save parameters of min dist
    }

    vectfree(r);
    vectfree(n);
    vectfree(q);

    l->next();                              // point to next root
  }

  l = 0;  
  delete l;                                 // delete root list

  double end[2];
  end[0] = egeom->knots[egeom->order-1]; 
  end[1] = egeom->knots[egeom->ncontpts];

  // Cho beg
  // for (i=0; i<2; i++) {                     // also check ends of curve
  for (int i=0; i<2; i++) {                     // also check ends of curve
  // Cho end
    u = end[i];

    r = evalbsp(egeom, u);
    n = normalcurv(egeom, u, 1);

    dist = distance(p, r);                  // sqrt(dist^2), i.e. |dist|
    q = sub_vect(p, r);
    if (dot(q, n) < 0.0)                    // add sign to absolute value
      dist *= -1.0;

    if (fabs(dist) < fabs(d0)) {            // update minimum distance
      d0 = dist;
      *u0 = u;                              // save parameters of min dist
    }

    vectfree(r);
    vectfree(n);
    vectfree(q);
  }
  vectfree(p);

  return d0;
}


//////////////////////////////////////////////////////////////////////////
// Point to surface
//////////////////////////////////////////////////////////////////////////

double robustMinDistance(double x, double y, double z, ParSurf *fgeom,
			 double eps, double *u0, double *v0)
{
  geom *a = dist_bez_ps(x, y, z, fgeom);    // form equations
  rootlist *l = solveMinDistance(a, eps);   // solve equations
  delete a;                                 // delete geom

  vector *n, *p = vectalloc(), *q, *r;
  p->x = x;
  p->y = y;
  p->z = z;
  p->w = 1.0;

  *u0 = *v0 = 0.0;
  double dist, d0 = 1.0e30, u, v;
  l->head();                                // starting at first root
  for (int i=0; i<l->size(); i++) {         // for each root in list
    real_array t = l->root();               // get current root
#ifndef USE_INTERVAL
    u = t[0];
    v = t[1];
#else
    u = t[0].center();
    v = t[1].center();
#endif
    if (u < fgeom->uknots[fgeom->uorder-1]) // insure that root is inside 
      u = fgeom->uknots[fgeom->uorder-1];   // the knot range
    if (u > fgeom->uknots[fgeom->ucontpts])
      u = fgeom->uknots[fgeom->ucontpts];

    if (v < fgeom->vknots[fgeom->vorder-1]) // insure that root is inside 
      v = fgeom->vknots[fgeom->vorder-1];   // the knot range
    if (v > fgeom->vknots[fgeom->vcontpts])
      v = fgeom->vknots[fgeom->vcontpts];

    r = evalsurf(fgeom, u, v);              // evaluate surface
    n = normalsurf(fgeom, u, v);            // surface normal

    dist = distance(p, r);                  // sqrt(dist^2), i.e. |dist|
    q = sub_vect(p, r);
    if (dot(q, n) < 0.0)                    // add sign to absolute value
      dist *= -1.0;

    if (fabs(dist) < fabs(d0)) {            // update minimum distance
      d0 = dist;
      *u0 = u;                              // save parameters of min dist
      *v0 = v;
    }

    vectfree(r);
    vectfree(n);
    vectfree(q);

    l->next();                              // point to next root
  }

  l = 0; 
  delete l;                                 // delete root list
  vectfree(p);

  //also check the boundary curves
  int k;                                    // for boundary extraction
  ParCurv *bcurve = egeomalloc1(fgeom->vorder, fgeom->vcontpts);
  //first, check the u=const. boundaries
  // Cho beg
  //  for ( i=0;i<2;i++ ) {
  for (int i=0;i<2;i++ ) {
  // Cho end
    u = (double)i;
    k = (i==0?0:fgeom->ucontpts-1);         // which boundary curve
    extract_edge(bcurve, fgeom, 0, k);      // extract the boundary curve

    dist = robustMinDistance(x, y, z, bcurve, eps, &v);
    if (fabs(dist) < fabs(d0)) {            // update minimum distance
      d0 = dist;
      *u0 = u;                              // save parameters of min dist
      *v0 = v;
    }
  }
  free_egeom(bcurve);
  //then, check the v=const. boundaries
  bcurve = egeomalloc1(fgeom->uorder, fgeom->ucontpts);
  // Cho beg
  //  for ( i=0;i<2;i++ ) {
  for (int i=0;i<2;i++ ) {
  // Cho end
    v = (double)i;
    k = (i==0?0:fgeom->vcontpts-1);
    extract_edge(bcurve, fgeom, 1, k);

    dist = robustMinDistance(x, y, z, bcurve, eps, &u);

    if (fabs(dist) < fabs(d0)) {            // update minimum distance
      d0 = dist;
      *u0 = u;                              // save parameters of min dist
      *v0 = v;
    }
  }

  free_egeom(bcurve);

  return d0;
}

/////////////////////////////////////////////////////////////////////////
// Curve to curve
/////////////////////////////////////////////////////////////////////////

double robustMinDistance(ParCurv *egeom1, ParCurv *egeom2, double eps,
			 double *u1, double *u2)
{
  geom *a = dist_bez_cc(egeom1, egeom2);    // form equations
  rootlist *l = solveMinDistance(a, eps);   // solve equations
  delete a;                                 // delete geom

  vector *r1, *r2;                          // evaluated points

  *u1 = 0.0; *u2 = 0.0;
  double dist, d0 = 1.0e30, t1, t2;
  l->head();                                // starting at first root
  for (int i=0; i<l->size(); i++) {         // for each root in list
    real_array t = l->root();               // get current root
#ifndef USE_INTERVAL
    t1 = t[0]; t2 = t[1];
#else
    t1 = t[0].center(); t2 = t[1].center(); 
#endif

    if (t1 < egeom1->knots[egeom1->order-1])   // insure that root is inside 
      t1 = egeom1->knots[egeom1->order-1];     // the knot range
    if (t1 > egeom1->knots[egeom1->ncontpts])
      t1 = egeom1->knots[egeom1->ncontpts];

    if (t2 < egeom2->knots[egeom2->order-1])   // insure that root is inside 
      t2 = egeom2->knots[egeom2->order-1];     // the knot range
    if (t2 > egeom2->knots[egeom2->ncontpts])
      t2 = egeom2->knots[egeom2->ncontpts];

    r1 = evalbsp(egeom1, t1);                // evaluate curve
    r2 = evalbsp(egeom2, t2);

    dist = distance(r1, r2);                 // sqrt(dist^2), i.e. |dist|

    if (fabs(dist) < fabs(d0)) {             // update minimum distance
      d0 = dist;
      *u1 = t1;                              // save parameters of min dist
      *u2 = t2;
    }

    vectfree(r1);
    vectfree(r2);

    l->next();                              // point to next root
  }

  l = 0; 
  delete l;                                 // delete root list

  //first, check the endpoints of the first curve to the 2nd one
  // Cho beg
  //  for (i=0; i<2; i++) {                     
  for (int i=0; i<2; i++) {                     
  // Cho end
    t1 = (double)i;
    r1 = evalbsp(egeom1, t1);               // evaluate the endpoints
    dist = robustMinDistance(r1->x*r1->w, r1->y*r1->w, r1->z*r1->w, 
			     egeom2, eps, &t2);
    if (fabs(dist) < fabs(d0)) {
      d0 = fabs(dist);                      // in case that dist<0
      *u1 = t1;
      *u2 = t2;
    }
  }
  //then, check the endpoints of the 2nd curve to the 1st one
  // Cho beg
  //  for (i=0; i<2; i++) {                     
  for (int i=0; i<2; i++) {                     
  // Cho end
    t2 = (double)i;
    r2 = evalbsp(egeom2, t2);
    dist = robustMinDistance(r2->x*r2->w, r2->y*r2->w, r2->z*r2->w, 
			     egeom1, eps, &t1);
    if (fabs(dist) < fabs(d0)) {
      d0 = fabs(dist);
      *u1 = t1;
      *u2 = t2;
    }
  }

  return d0;
}

////////////////////////////////////////////////////////////////////////
// Curve to surface
////////////////////////////////////////////////////////////////////////

double robustMinDistance(ParCurv *egeom, ParSurf *fgeom, double eps,
			 double *u0, double *u1, double *v1)
{
  geom *a = dist_bez_cs(egeom, fgeom);      // form equations
  rootlist *l = solveMinDistance(a, eps);   // solve equations
  delete a;                                 // delete geom

  vector *r1, *r2;                          //evaluated points
  *u0 = 0.0; 
  *u1 = 0.0; *v1 = 0.0;
  double dist, d0 = 1.0e30, t1, u, v;
  l->head();                                // starting at first root

  for (int i=0; i<l->size(); i++) {         // for each root in list
    real_array t = l->root();               // get current root

#ifndef USE_INTERVAL
    t1 = t[0]; 
    u = t[1]; v = t[2];
#else
    t1 = t[0].center(); 
    u = t[1].center(); v = t[2].center(); 
#endif

    if (t1 < egeom->knots[egeom->order-1])   // insure that root is inside 
      t1 = egeom->knots[egeom->order-1];     // the knot range
    if (t1 > egeom->knots[egeom->ncontpts])
      t1 = egeom->knots[egeom->ncontpts];

    if (u < fgeom->uknots[fgeom->uorder-1])   // insure that root is inside 
      u = fgeom->uknots[fgeom->uorder-1];     // the knot range
    if (u > fgeom->uknots[fgeom->ucontpts])
      u = fgeom->uknots[fgeom->ucontpts];

    if (v < fgeom->vknots[fgeom->vorder-1])   // insure that root is inside 
      v = fgeom->vknots[fgeom->vorder-1];     // the knot range
    if (v > fgeom->vknots[fgeom->vcontpts])
      v = fgeom->vknots[fgeom->vcontpts];

    r1 = evalbsp(egeom, t1);                // evaluate curve
    r2 = evalsurf(fgeom, u, v);             //evaluate surface

    dist = distance(r1, r2);                 // sqrt(dist^2), i.e. |dist|

    if (fabs(dist) < fabs(d0)) {             // update minimum distance
      d0 = dist;
      *u0 = t1;                              // save parameters of min dist
      *u1 = u; *v1 = v;
    }

    vectfree(r1);
    vectfree(r2);

    l->next();                              // point to next root
  }

  l = 0; 
  delete l;                                 // delete root list

  //first, check the endpoints of the curve to the surface
  // Cho beg
  //  for ( i=0;i<2;i++ ) {
  for (int i=0;i<2;i++ ) {
  // Cho end
    t1 = double(i);
    r1 = evalbsp(egeom, t1);                // evaluate the endpoints

    dist = robustMinDistance(r1->x*r1->w, r1->y*r1->w, r1->z*r1->w, 
			     fgeom, eps, &u, &v);

    if (fabs(dist) < fabs(d0)) {             // update minimum distance
      d0 = fabs(dist);                       // in case that dist<0
      *u0 = t1;                              // save parameters of min dist
      *u1 = u; *v1 = v;
    }
  }

  //then, check u = 0,1 boundary curves to the curve
  int k;                                     // for boundary extraction
  ParCurv* bcurve = egeomalloc1(fgeom->vorder, fgeom->vcontpts);

  // Cho beg
  //  for ( i=0;i<2;i++ ) {
  for (int i=0;i<2;i++ ) {
  // Cho end
    u = (double)i;
    k = (i==0?0:fgeom->ucontpts-1);          // which boundary curve

    extract_edge(bcurve, fgeom, 0, k);

    dist = robustMinDistance(egeom, bcurve, eps, &t1, &v);

    if (fabs(dist) < fabs(d0)) {
      d0 = dist;
      *u0 = t1;
      *u1 = u; *v1 = v;
    }
  }

  free_egeom(bcurve);

  //then, check v = 0,1 boundary curves to the curve
  bcurve = egeomalloc1(fgeom->uorder, fgeom->ucontpts);

  // Cho beg
  //  for ( i=0;i<2;i++ ) {
  for (int i=0;i<2;i++ ) {
  // Cho end
    v = (double)i;
    k = (i==0?0:fgeom->vcontpts-1);

    extract_edge(bcurve, fgeom, 1, k);

    dist = robustMinDistance(egeom, bcurve, eps, &t1, &u);

    if (fabs(dist) < fabs(d0)) {
      d0 = dist;
      *u0 = t1;
      *u1 = u; *v1 = v;
    }
  }

  free_egeom(bcurve);

  return d0;
}

//////////////////////////////////////////////////////////////////////////
// Surface to surface
//////////////////////////////////////////////////////////////////////////

double robustMinDistance(ParSurf *fgeom1, ParSurf *fgeom2, double eps,
			 double *u1, double *v1, double *u2, double *v2)
{
  geom *a = dist_bez_ss(fgeom1, fgeom2);    // form equations
  rootlist *l = solveMinDistance(a, eps);   // solve equations
  delete a;                                 // delete geom
  vector *r1, *r2;                          //evaluated points

  *u1 = 0.0; *v1 = 0.0;
  *u2 = 0.0; *v2 = 0.0;
  double dist, d0 = 1.0e30, s1, t1, s2, t2;
  l->head();                                // starting at first root
  for (int i=0; i<l->size(); i++) {         // for each root in list
    real_array t = l->root();               // get current root
#ifndef USE_INTERVAL
    s1 = t[0]; t1 = t[1]; 
    s2 = t[2]; t2 = t[3];
#else
    s1 = t[0].center(); t1 = t[1].center(); 
    s2 = t[2].center(); t2 = t[3].center(); 
#endif

    if (s1 < fgeom1->uknots[fgeom1->uorder-1])   // insure that root is inside 
      s1 = fgeom1->uknots[fgeom1->uorder-1];     // the knot range
    if (s1 > fgeom1->uknots[fgeom1->ucontpts])
      s1 = fgeom1->uknots[fgeom1->ucontpts];

    if (t1 < fgeom1->vknots[fgeom1->vorder-1])   // insure that root is inside 
      t1 = fgeom1->vknots[fgeom1->vorder-1];     // the knot range
    if (t1 > fgeom1->vknots[fgeom1->vcontpts])
      t1 = fgeom1->vknots[fgeom1->vcontpts];

    if (s2 < fgeom2->uknots[fgeom2->uorder-1])   // insure that root is inside 
      s2 = fgeom2->uknots[fgeom2->uorder-1];     // the knot range
    if (s2 > fgeom2->uknots[fgeom2->ucontpts])
      s2 = fgeom2->uknots[fgeom2->ucontpts];

    if (t2 < fgeom2->vknots[fgeom2->vorder-1])   // insure that root is inside 
      t2 = fgeom2->vknots[fgeom2->vorder-1];     // the knot range
    if (t2 > fgeom2->vknots[fgeom2->vcontpts])
      t2 = fgeom2->vknots[fgeom2->vcontpts];

    r1 = evalsurf(fgeom1, s1, t1);           // evaluate curve
    r2 = evalsurf(fgeom2, s2, t2);           //evaluate surface

    dist = distance(r1, r2);                 // sqrt(dist^2), i.e. |dist|

    if (fabs(dist) < fabs(d0)) {             // update minimum distance
      d0 = dist;
      *u1 = s1; *v1 = t1;                    // save parameters of min dist
      *u2 = s2; *v2 = t2;
    }

    vectfree(r1);
    vectfree(r2);

    l->next();                              // point to next root
  }

  l = 0; 
  delete l;                                 // delete root list

  //check the boundary curves
  //the boundary curves of the first surface to the second surface
  int k;
  ParCurv *bcurve;
  bcurve = egeomalloc1(fgeom1->vorder, fgeom1->vcontpts);
  // Cho beg
  //  for ( i=0;i<2;i++ ) {
  for (int i=0;i<2;i++ ) {
  // Cho end
    s1 = (double)i;
    k = (i==0?0:fgeom1->ucontpts-1);
    extract_edge(bcurve, fgeom1, 0, k);
    dist = robustMinDistance(bcurve, fgeom2, eps, &t1, &s2, &t2);
    if (fabs(dist) < fabs(d0)) {            
      d0 = dist;
      *u1 = s1; *v1 = t1;                   
      *u2 = s2; *v2 = t2;
    }
  }
  free_egeom(bcurve);
  bcurve = egeomalloc1(fgeom1->uorder, fgeom1->ucontpts);
  // Cho beg
  //  for ( i=0;i<2;i++ ) {
  for (int i=0;i<2;i++ ) {
  // Cho end
    t1 = (double)i;
    k = (i==0?0:fgeom1->vcontpts-1);
    extract_edge(bcurve, fgeom1, 1, k);
    dist = robustMinDistance(bcurve, fgeom2, eps, &s1, &s2, &t2);
    if (fabs(dist) < fabs(d0)) {            
      d0 = dist;
      *u1 = s1; *v1 = t1;                   
      *u2 = s2; *v2 = t2;
    }
  }
  free_egeom(bcurve);
  
  //the boundary curves of the second surface to the first surface
  bcurve = egeomalloc1(fgeom2->vorder, fgeom2->vcontpts);
  // Cho beg
  //  for ( i=0;i<2;i++ ) {
  for (int i=0;i<2;i++ ) {
  // Cho end
    s2 = (double)i;
    k = (i==0?0:fgeom2->ucontpts-1);
    extract_edge(bcurve, fgeom2, 0, k);
    dist = robustMinDistance(bcurve, fgeom1, eps, &t2, &s1, &t1);
    if (fabs(dist) < fabs(d0)) {            
      d0 = dist;
      *u1 = s1; *v1 = t1;                   
      *u2 = s2; *v2 = t2;
    }
  }
  free_egeom(bcurve);
  bcurve = egeomalloc1(fgeom2->uorder, fgeom2->ucontpts);
  // Cho beg
  //  for ( i=0;i<2;i++ ) {
  for (int i=0;i<2;i++ ) {
  // Cho end
    t2 = (double)i;
    k = (i==0?0:fgeom2->vcontpts-1);
    extract_edge(bcurve, fgeom2, 1, k);
    dist = robustMinDistance(bcurve, fgeom1, eps, &s2, &s1, &t1);
    if (fabs(dist) < fabs(d0)) {            
      d0 = dist;
      *u1 = s1; *v1 = t1;                   
      *u2 = s2; *v2 = t2;
    }
  }
  free_egeom(bcurve);

  return d0;
}

//////////////////////////////////////////////////////////////////////////

rootlist *solveMinDistance(geom *a, double eps)
{
  geom **geom_list = distObj(a);            // formulate equations
  mn_array *mn_list;
  mn_list = convertToMonomial(geom_list);   // convert to monmial basis
  for (int i=0; i<a->nd; i++) {
    delete geom_list[i];                    // delete each geom in array
  }
  delete [] geom_list;                      // delete geom array



  rootlist *l = si_pinter(*mn_list, eps);   // solve system of equations

  delete mn_list;                           // delete multinomial array

  return l;
}

// convertToMonomial - convert bernstein to monomial basis

mn_array *convertToMonomial(geom** geom_list)
{
  int i, j, n = geom_list[0]->nd;

  mn_array *mn_list = new mn_array(n);
  sht_array olist(n);

  for (i=0; i<n; i++) {
    for(j=0; j<n; j++)
      olist[j] = geom_list[i]->order[j] - 1;
    
    multinomial *mn = new multinomial (olist);
    
    for (j=0; j<geom_list[i]->nc; j++)
      mn->bp[j] = geom_list[i]->contpts[j]->x;
    ((multinomial **)*mn_list)[i] = mn;
  }

  return mn_list;
}

/************************************************************************
  This portion of routine creat the basic object function for calculating 
  dist between two rational bezier expressions
 ************************************************************************/

geom **geom_l;

geom **distObj(geom *a)
{
  geom w(*a);
  for (int i=0; i<a->nc; i++)
    w.contpts[i]->x = a->contpts[i]->w;

  geom *x_2 = (*a)*(*a);

  geom y_2(*x_2), z_2(*x_2);

  // Cho beg
  //  for(i=0; i<x_2->nc; i++) {
  for(int i=0; i<x_2->nc; i++) {
  // Cho end
    y_2.contpts[i]->x = y_2.contpts[i]->y;
    z_2.contpts[i]->x = z_2.contpts[i]->z;
  }

  geom *xy2 = *x_2 + y_2;
  geom *xyz = *xy2 + z_2;
  delete x_2;
  delete xy2;

  int nd = a->nd;
  geom_l = new geom*[nd];

  // Cho beg
  //  for (i=0; i<nd; i++)
  for (int i=0; i<nd; i++) {
    geom_l[i] = deri_geom(xyz, i);
  }
  // Cho end
  delete xyz;

  return geom_l;
}

/************************************************************************
  This portion of routine creat the basic function for calculating dist
  between two rational bezier curves

  The basic function is:

                           S(u) = P0 - P(u);

  contpts:    C_{i} = P_i - P0*w_i;			   
 ************************************************************************/

geom *dist_bez_pc(double x, double y, double z, ParCurv *egeom)
{
  int nd = MDIST_PC, *ncont = new int[nd];

  ncont[0] = egeom->ncontpts;
  
  int nc = 1;
  // Cho beg
  int i;
  //  for (int i=0; i<nd; i++)
  for (i=0; i<nd; i++)
    nc *= ncont[i];
  // Cho end

  int j = egeom->order + egeom->ncontpts;
    
  double *kn = new double[j];
  vector **cont = vec_array1(nc);

  int k;
  for (i=k=0; i<egeom->ncontpts; i++, k++) {
    cont[k]->x = egeom->contpts[i]->x - x*egeom->contpts[i]->w;
    cont[k]->y = egeom->contpts[i]->y - y*egeom->contpts[i]->w;
    cont[k]->z = egeom->contpts[i]->z - z*egeom->contpts[i]->w;
    cont[k]->w=egeom->contpts[i]->w;
  }

  for (i=j=0; i<egeom->order + egeom->ncontpts; i++, j++)
    kn[j] = egeom->knots[i];

  geom *c = new geom(nd, ncont, kn, cont);

  free_varray1(cont, nc);
  delete [] kn;
  delete [] ncont;

  return c;
}

/************************************************************************
  This portion of routine creat the basic function for calculating dist
  between ta point p0 and a Bezier surface:

  The basic function is:

                           S(u,v) = P0 - P(u,v);

  contpts:    C_{i} = P_i - P0*w_i;			   
 ************************************************************************/

geom *dist_bez_ps(double x, double y, double z, ParSurf *fgeom)
{
  int *ncont=new int[MDIST_PS];

  ncont[0] = fgeom->ucontpts;
  ncont[1] = fgeom->vcontpts;
  
  int nc = 1;
  // Cho beg
  int i;
  //  for (int i=0; i<MDIST_PS; i++)
  for (i=0; i<MDIST_PS; i++)
    nc *= ncont[i];
  // Cho end

  int j = fgeom->uorder + fgeom->ucontpts + fgeom->vorder + fgeom->vcontpts;
    
  double *kn = new double[j];
  vector **cont = vec_array1(nc);

  int k;
  for(i=k=0; i<fgeom->ucontpts; i++)
    for(j=0; j<fgeom->vcontpts; j++, k++) {
      cont[k]->x = fgeom->contpts[i][j]->x - x*fgeom->contpts[i][j]->w;
      cont[k]->y = fgeom->contpts[i][j]->y - y*fgeom->contpts[i][j]->w;
      cont[k]->z = fgeom->contpts[i][j]->z - z*fgeom->contpts[i][j]->w;
      cont[k]->w = fgeom->contpts[i][j]->w;
    }


  for (i=k=0; i<fgeom->uorder + fgeom->ucontpts; i++, k++)
    kn[k]=fgeom->uknots[i];
  for (i=0; i<fgeom->vorder + fgeom->vcontpts; i++, k++)
    kn[k]=fgeom->vknots[i];

  geom *c = new geom(MDIST_PS, ncont, kn, cont);

  free_varray1(cont, nc);
  delete [] kn;
  delete [] ncont;

  return c;
}

/************************************************************************
  This portion of routine creat the basic function for calculating dist
  between two rational bezier curves

  The basic function is:

                            P(u) - Q(v) = S(u,v)

    where P(u) and Q(v) are two rational bezier curves, S(u,V) is a 
    bezier surface with control points:

                  C_{i,j} =  w_j * P_i - w_i * Q_j   and weights :

		  w_i * w_j
 ************************************************************************/

geom *dist_bez_cc(ParCurv *egeom1, ParCurv *egeom2)
{
  int nd = MDIST_CC, *ncont = new int[nd];

  ncont[0] = egeom1->ncontpts;
  ncont[1] = egeom2->ncontpts;
  
  int nc = 1;
  // Cho beg
  int i;
  //  for (int i=0; i<nd; i++)
  for (i=0; i<nd; i++)
    nc *= ncont[i];
  // Cho end

  int j = egeom1->order + egeom1->ncontpts + egeom2->order + egeom2->ncontpts;
    
  double *kn=new double[j];
  vector **cont=vec_array1(nc);

  int k;
  for (i=k=0; i<egeom1->ncontpts; i++)
    for (j=0; j<egeom2->ncontpts; j++, k++) {
      cont[k]->x = egeom2->contpts[j]->w * egeom1->contpts[i]->x -
	           egeom1->contpts[i]->w * egeom2->contpts[j]->x;
      cont[k]->y = egeom2->contpts[j]->w * egeom1->contpts[i]->y -
	           egeom1->contpts[i]->w * egeom2->contpts[j]->y;
      cont[k]->z = egeom2->contpts[j]->w * egeom1->contpts[i]->z -
	           egeom1->contpts[i]->w * egeom2->contpts[j]->z;
      cont[k]->w = egeom1->contpts[i]->w * egeom2->contpts[j]->w;
    }

  for (i=j=0; i<egeom1->order + egeom1->ncontpts; i++, j++)
    kn[j] = egeom1->knots[i];
  for (i=0; i<egeom2->order + egeom2->ncontpts; i++, j++)
    kn[j] = egeom2->knots[i];

  geom *c = new geom(nd, ncont, kn, cont);

  free_varray1(cont, nc);
  delete [] kn;
  delete [] ncont;

  return c;
}

/************************************************************************
  This portion of routine creat the basic function for calculating dist
  between a rational bezier curve and a rational Bezier patch

  The basic function is:

                        P(s) - Q(u,v) = S(s, u,v)

    where P(s) and Q(u,v) are two rational bezier curves, S(u,V) is a 
    bezier surface with control points:

             C_{i,j} =  w_{j,k} * P_i - w_i * Q_{j,k}   and weights :

		  w_i * w_{j,k}
 ************************************************************************/

geom *dist_bez_cs(ParCurv *egeom, ParSurf *fgeom)
{
  int *ncont = new int[MDIST_CS];

  ncont[0] = egeom->ncontpts;
  ncont[1] = fgeom->ucontpts;
  ncont[2] = fgeom->vcontpts;
  
  int nc = 1;
  // Cho beg
  int i;
  //  for (int i=0; i<MDIST_CS; i++)
  for (i=0; i<MDIST_CS; i++)
    nc *= ncont[i];
  // Cho end

  int j = egeom->order + egeom->ncontpts+ fgeom->uorder +
          fgeom->ucontpts+fgeom->vorder + fgeom->vcontpts;
    
  double *kn = new double[j];
  vector **cont = vec_array1(nc);

  int k, l; 
  for(i=l=0; i<egeom->ncontpts; i++)
    for (j=0; j<fgeom->ucontpts; j++)
      for (k=0; k<fgeom->vcontpts; k++, l++) {
	cont[l]->x = fgeom->contpts[j][k]->w *egeom->contpts[i]->x -
	             egeom->contpts[i]->w * fgeom->contpts[j][k]->x;
	cont[l]->y = fgeom->contpts[j][k]->w * egeom->contpts[i]->y -
	             egeom->contpts[i]->w * fgeom->contpts[j][k]->y;
	cont[l]->z = fgeom->contpts[j][k]->w * egeom->contpts[i]->z -
	             egeom->contpts[i]->w * fgeom->contpts[j][k]->z;
	cont[l]->w = egeom->contpts[i]->w * fgeom->contpts[j][k]->w;
      }
      
  for (i=j=0; i<egeom->order + egeom->ncontpts; i++, j++)
    kn[j]=egeom->knots[i];
  for (i=0; i<fgeom->uorder + fgeom->ucontpts; i++, j++)
    kn[j]=fgeom->uknots[i];
  for (i=0; i<fgeom->vorder + fgeom->vcontpts; i++, j++)
    kn[j]=fgeom->vknots[i];

  geom *c = new geom(MDIST_CS, ncont, kn, cont);

  free_varray1(cont, nc);
  delete [] kn;
  delete [] ncont;

  return c;
}

/************************************************************************
  This portion of routine creat the basic function for calculating dist
  between two rational Bezier patches

  The basic function is:

                        P(s,t) - Q(u,v) = S(s,t, u,v)

    where P(s) and Q(u,v) are two rational bezier curves, S(u,V) is a 
    bezier surface with control points:

             C_{i,j,k,l} =  w_{k,l} * P_{i,j} - w_{i,j} * Q_{k,l}  
	     
    and weights :

		  w_{i,k} * w_{k,l}
 ************************************************************************/

geom *dist_bez_ss(ParSurf *fgeom1, ParSurf *fgeom2)
{
  int *ncont = new int[MDIST_SS];

  ncont[0] = fgeom1->ucontpts;
  ncont[1] = fgeom1->vcontpts;
  ncont[2] = fgeom2->ucontpts;
  ncont[3] = fgeom2->vcontpts;
  
  int nc=1;
  // Cho beg
  int i;
  //  for (int i=0; i<MDIST_SS; i++)
  for (i=0; i<MDIST_SS; i++)
    nc *= ncont[i];
  // Cho end

  int j = fgeom1->uorder + fgeom1->ucontpts + fgeom1->vorder +
          fgeom1->vcontpts + fgeom2->uorder + fgeom2->ucontpts +
	  fgeom2->vorder + fgeom2->vcontpts;
    
  double *kn = new double[j];
  vector **cont = vec_array1(nc);

  int k, l, m;
   for(i=m=0; i<fgeom1->ucontpts; i++)
    for(j=0; j<fgeom1->vcontpts; j++)
      for(k=0; k<fgeom2->ucontpts; k++)
	for(l=0; l<fgeom2->vcontpts; l++, m++){
	  cont[m]->x = fgeom2->contpts[k][l]->w * fgeom1->contpts[i][j]->x -
	               fgeom1->contpts[i][j]->w * fgeom2->contpts[k][l]->x;
	  cont[m]->y = fgeom2->contpts[k][l]->w * fgeom1->contpts[i][j]->y -
	               fgeom1->contpts[i][j]->w * fgeom2->contpts[k][l]->y;
	  cont[m]->z = fgeom2->contpts[k][l]->w * fgeom1->contpts[i][j]->z -
	               fgeom1->contpts[i][j]->w * fgeom2->contpts[k][l]->z;
	  cont[m]->w = fgeom1->contpts[i][j]->w * fgeom2->contpts[k][l]->w;
	}
      
  for (i=j=0; i<fgeom1->uorder + fgeom1->ucontpts; i++, j++)
    kn[j]=fgeom1->uknots[i];
  for (i=0; i<fgeom1->vorder + fgeom1->vcontpts; i++, j++)
    kn[j]=fgeom1->vknots[i];
  for (i=0; i<fgeom2->uorder + fgeom2->ucontpts; i++, j++)
    kn[j]=fgeom2->uknots[i];
  for (i=0; i<fgeom2->vorder + fgeom2->vcontpts; i++, j++)
    kn[j]=fgeom2->vknots[i];

  geom *c = new geom(MDIST_SS, ncont, kn, cont);

  free_varray1(cont, nc);
  delete [] kn;
  delete [] ncont;

  return c;
}

// geom.cc - geom class

geom::geom(int)
{}

geom::geom(const int n, int *ncont)
{
  nd = n;                        // num of dimensions 
  ncontpts = new int[nd];        // num of contpts in each dimension
  order = new int[nd];           // num of order in each dimension

  nc = 1;
  nk = 0;
  for (int i=0; i<nd; i++) {
    nk += ncont[i]+ncont[i];     // num of total knots
    nc *= ncont[i];              // num of total contpts
    ncontpts[i] = ncont[i];        
    order[i] = ncontpts[i];
  }

  knots = new double[nk];        // alloc memory for knots
  contpts = vec_array1(nc);      // alloc memory for contpts
}

geom::geom(int n, int *ncont, double *kn, vector **cont)  
{
  nd = n;                        // num of dimensions 
  ncontpts = new int[nd];        // num of contpts in each dimension
  order = new int[nd];           // num of order in each dimension

  nc = 1;
  nk = 0;
  // Cho beg
  int i;
  //  for (int i=0; i<nd; i++) {
  for (i=0; i<nd; i++) {
  // Cho end
    nk += ncont[i]+ncont[i];              // num of total knots
    nc *= ncont[i];                       // num of total contpts
    ncontpts[i] = ncont[i];        
    order[i] = ncontpts[i];
  }
  knots = new double[nk];                // alloc memory for knots
  for (i=0; i<nk; i++)
    knots[i] = kn[i];

  contpts = vec_array1(nc);
  for (i=0; i<nc; i++) {
    contpts[i]->x = cont[i]->x;
    contpts[i]->y = cont[i]->y;
    contpts[i]->z = cont[i]->z;
    contpts[i]->w = cont[i]->w;
  }
}

geom::geom(geom &a)
{
  int i=0;

  nd = a.nd;
  nk = a.nk;
  nc = a.nc;

  ncontpts = new int[nd];
  order = new int[nd];
  for (i=0; i<nd; i++) {
    ncontpts[i] = a.ncontpts[i];
    order[i] = a.order[i];
  }

  knots = new double[nk];
  for (i=0; i<nk; i++)
    knots[i] = a.knots[i];
  
  contpts = vec_array1(nc);
  for (i=0; i<nc; i++) {
    contpts[i]->x = a.contpts[i]->x;
    contpts[i]->y = a.contpts[i]->y;
    contpts[i]->z = a.contpts[i]->z;
    contpts[i]->w = a.contpts[i]->w;
  }
}

geom::~geom()                    // destructor of geom
{
  free_varray1(contpts, nc);
  
  delete[] knots;
  delete[] order;
  delete[] ncontpts;
}

geom geom::operator -= (geom &G)
{
  sub (*this, G, *this);
  return *this;
}

geom geom::operator += (geom &G)
{
  add (*this, G, *this);
  return *this;
}

geom* geom::operator=(geom *a)
{
  geom *c = new geom(a->nd, a->ncontpts);

  // Cho beg
  int i;
  //  for (int i=0; i<c->nk; i++)
  for (i=0; i<c->nk; i++)
    c->knots[i] = a->knots[i];
  // Cho end

  for (i=0; i<nc; i++) {
    c->contpts[i]->x = a->contpts[i]->x;
    c->contpts[i]->y = a->contpts[i]->y;
    c->contpts[i]->z = a->contpts[i]->z;
    c->contpts[i]->w = a->contpts[i]->w;
  }
  
  return c;
}

geom geom::operator = (geom &a)
{
  geom c(a);
 
  return c;
}

ostream &operator << (ostream &stream, const geom &g)
{
  int i, j=0;
  int factor = 5;

  stream << g.nd << ' ';
  stream << endl;
  
  for (i=0; i<g.nd; i++) {
    // Cho beg
    stream.precision (15);
    stream << g.order[i] << ' ';
  //    stream <<  setprecision (15) << g.order[i] << ' ';
    // Cho end
  }

  stream << endl;

  for (i = 0; i < g.nc; i++){
    // Cho beg
    stream.precision (15);
    stream << g.contpts[i]->x << " ";
    //    stream <<  setprecision (15) << g.contpts[i]->x << " ";
    // Cho end
    if ((i+1)/factor > j) {
      cout << endl;
      j++;
    }
  }
  stream << endl;

  for ( i=0;i<g.nk;i++ )
    cout<<g.knots[i]<<" ";
  cout<<endl;

  return (stream);
}

geom* operator * (geom &a, geom &b)
{
//This routine calculate multiplication of two bezier expressions
//Let the two original bezier expressions have order (or contpts):
// (m_1, m_2, ... , m_nd) and (n_1, n_2, ... , n_nd), then
// the new orders are (m_1+n_1-1, m_2+n_2-1, ..., m_nd+n_nd-1)
// the new contpts are of (m_1+n_1-1)*(m_2+n_2-1)*...*(m_nd+n_nd-1)
// the ith contpts is defined by:
//                1                          MIN(m_1-1, k_1)
//------------------------------------  *    ---------------  *
//  + m_1+n_1 +       + m_nd + n_nd +        \
// +           +...  +               +       /
// +           +     +               +       ---------------
//  +   k_1   +       +     k_nd    +     i_1=MAX(0,k_1-n_1+1)
//
//          MIN(m_nd-1, k_nd)
// ... *    -----------------      + m_1 +      + m_nd +
//          \                     +       +....+        + X(1)_{i1,..., i_nd}  
//          /                     +       +    +        +
//          -----------------      + i_1 +      + i_nd +
//       i_nd=MAX(0,k_nd-n_nd+1)
//
//   *X(2)_{i1,..., i_nd};
//

  int *ncont = new int[a.nd];
  int *a_ncont = new int[a.nd];
  int *b_ncont = new int[a.nd];
  int *ab_ncont = new int[a.nd];

  // Cho beg
  int i;
  //  for (int i=0; i<a.nd; i++) {
  for (i=0; i<a.nd; i++) {
  // Cho end
    if (a.ncontpts[i] != a.order[i]) {
      cerr << "a is not a Bezier expression in Bezier mul" << endl;
/*      exit(1);*/
    }
    
    ncont[i] = a.ncontpts[i]+b.ncontpts[i]-1;
    a_ncont[i] = a.ncontpts[i]-1;             // a's degree
    b_ncont[i] = b.ncontpts[i]-1;              // b's degree
    ab_ncont[i] = a_ncont[i]+b_ncont[i];       // a*b's degree
  }

  geom *c = new geom(a.nd, ncont);  

  int j = 0, *k, *MI, *MA, *MM;
  int count_i = 0, count_j = 0;

  for (i=0; i<c->nc; i++) {
    k = new int[a.nd];

    rhash(i, a.nd, k, c->ncontpts);

    MI = new int[a.nd];
    MA = new int[a.nd];
    MM = new int[a.nd];
   
    int ni = get_MIN_MAX(a.nd, k, a_ncont, b_ncont, MI, MA);

    sub_array1(a.nd, MA, MI, MM);   // MA[i]-MI[i] +1=MM[i]

    double tmp_x = 0.0, tmp_y = 0.0, tmp_z = 0.0, coef = 0.0;
    
    for (j=0; j<ni; j++) {
      int *i_a = new int[a.nd];
      int *i_b = new int[a.nd];

      rhash(j, a.nd, i_a, MM);  // the combination form for i1, i2, ..., in

      add_base(a.nd, i_a, MA); // where i1 is from MAX(0, k1-n1) to 
	                        // Min(m1, k1), etc.
      
      sub_array(a.nd, k, i_a, i_b);  // i_b[i]=k[i]-i_a[i]

      coef = get_coef1(a.nd, a_ncont, b_ncont, i_a, i_b);
     
      count_i = rsort(a.nd, i_a, a.ncontpts);
      count_j = rsort(a.nd, i_b, b.ncontpts);

     
      tmp_x += coef * a.contpts[count_i]->x * b.contpts[count_j]->x;
      tmp_y += coef * a.contpts[count_i]->y * b.contpts[count_j]->y;
      tmp_z += coef * a.contpts[count_i]->z * b.contpts[count_j]->z;

      //      cout<<coef<<" "<<tmp_x<<" "<<tmp_y<<" "<<tmp_z<<endl;

      delete[] i_a;
      delete[] i_b;
    }
    delete[] MI;
    delete[] MA;
    delete[] MM;
     
    coef = get_coef0(a.nd, ab_ncont, k);
    delete[] k;

    if (fabs(coef) < ZERO) 
      cout << " coef =0->0 in multiplying two Bezier Exp." << endl;

    //    cout<<"final = "<<coef<<" "<<tmp_x<<" "<<tmp_y<<" "<<tmp_z<<endl;
    
    c->contpts[i]->x = tmp_x/(double)coef;
    c->contpts[i]->y = tmp_y/(double)coef;
    c->contpts[i]->z = tmp_z/(double)coef;
    c->contpts[i]->w = 1.0;
  }

  delete[] ncont;
  delete[] a_ncont;
  delete[] b_ncont;
  delete[] ab_ncont;
 
  return c;
}

geom *operator + (geom &a, geom &b)            // addition of two integer geoms
{     
  geom *c = new geom(a.nd, a.ncontpts);   // result of two geoms addition

  for (int i=0; i<a.nk; i++)
    c->knots[i] = a.knots[i];
  
  // Cho beg
  //  for (i=0; i<a.nc; i++) {
  for (int i=0; i<a.nc; i++) {
  // Cho end
    c->contpts[i]->x=a.contpts[i]->x+b.contpts[i]->x;
    c->contpts[i]->y=a.contpts[i]->y+b.contpts[i]->y;
    c->contpts[i]->z=a.contpts[i]->z+b.contpts[i]->z;
    c->contpts[i]->w=1.0;
  }
  return c;
}

void add(geom &a, geom &b, geom &c)
{
  for (int i=0; i<a.nk; i++)
    c.knots[i] = a.knots[i];

  // Cho beg
  //  for (i=0; i<a.nc; i++) {
  for (int i=0; i<a.nc; i++) {
  // Cho end
    c.contpts[i]->x = a.contpts[i]->x+b.contpts[i]->x;
    c.contpts[i]->y = a.contpts[i]->y+b.contpts[i]->y;
    c.contpts[i]->z = a.contpts[i]->z+b.contpts[i]->z;
    c.contpts[i]->w = 1.0;
  }
}

void add_base(int nd, int* i_a, int* MI)
{
  for (int i=0; i<nd; i++)
    i_a[i] += MI[i];
}

double combine_cc(int a, int b)   // a must be larger than b
{
/*  if (a < b) {
    cerr<<" a is less than b in combine_cc "<<"\n";
    exit(1);
  }*/
  return factorial_cc(a)/(factorial_cc(b)*factorial_cc(a-b));
}

void copy_array(int nd, int *ps1, int *ps2)   // copy array ps1 to ps2
{
  for (int i=0; i<nd; i++)
    ps2[i] = ps1[i];
}

// friend function to calculate the first deri[Bvative of Bezier expression
//The new expression will have ncontpts[0]* .. (ncontpts[i]-1)*..*ncontpts[nd]
// contpts; and the j contpts will be:
//   (ncontpts[i]-1)(C_{i1, i2,...i_(i+1), ..., i_nd} -
//                   C_{i1, i2,...i_(i), ..., i_nd}) = V_{i1, i2, ..., i_nd}
//
//   In the routine, a is the original Bezier expression, i is the direction
// (i=0, 1, ..., nd). It can be for rational Bezier expression but 
// the returned geom is an integer Bezier expression.
//

geom *deri_geom(geom *a, int i)
{
  int *ncont=new int[a->nd];
  
  for (int j=0; j<a->nd; j++) {
    ncont[j] = a->ncontpts[j];
    if (j == i)
      ncont[j]--;
  }

  geom *b = new geom(a->nd, ncont);
  
  int j1 = 0, j2 = 0;
  // Cho beg
  //  for (j=0; j<b->nc; j++) {
  for (int j=0; j<b->nc; j++) {
  // Cho end
    int *ps1 = new int[b->nd];
    int *ps2 = new int[b->nd];
    rhash(j, a->nd, ps1, b->ncontpts);
    
    copy_array(a->nd, ps1, ps2);
    ps2[i] = ps2[i]+1;

    j1 = rsort(a->nd, ps1, a->ncontpts);
    j2 = rsort(a->nd, ps2, a->ncontpts);

    b->contpts[j]->x = b->order[i] * (a->contpts[j2]->x/a->contpts[j2]->w -
				      a->contpts[j1]->x/a->contpts[j1]->w);
    b->contpts[j]->y = b->order[i] * (a->contpts[j2]->y/a->contpts[j2]->w -
    				      a->contpts[j1]->y/a->contpts[j1]->w);
    b->contpts[j]->z = b->order[i] * (a->contpts[j2]->z/a->contpts[j2]->w -
    				      a->contpts[j1]->z/a->contpts[j1]->w);
    b->contpts[j]->w = 1.0;
    
    delete [] ps1;
    delete [] ps2;
  }

  // Cho beg
  //  j = 0;
  int j = 0;
  // Cho end

  for (j1=0; j1<b->nd; j1++) {
    for (j2=0; j2<b->order[j1]; j2++) { 
      b->knots[j] = 0.0;
      j++;
    }
    for (j2=b->order[j1]; j2<b->order[j1]+b->order[j1]; j2++) {
      b->knots[j] = 1.0;
      j++;
    }
  }

  delete [] ncont;
  return b;
}

double factorial_cc(int n)
{
  double prod = 1.0;
  for (int i=0; i<n; i++)
    prod = prod*(i+1.0);

  return prod;
}

double get_coef0(int n, int *a, int *b)
{
  if (n == 1)
    return combine_cc(a[0], b[0]);  
  else
    return combine_cc(a[n-1], b[n-1])*get_coef0(n-1, a, b);
}

double get_coef1(int n, int *a, int *b, int *i_a, int *i_b)
{
  return (double)get_coef0(n, a, i_a) * get_coef0(n, b, i_b);
}

  ///////////////////////////////////////////
  //       This routine calc the result:
  //     ++   ++   ++    ++        ++     ++
  //    ++ a1  ++ ++  a2  ++  ... ++  an   ++
  //    ++ b1  ++ ++  b2  ++  ... ++  bn   ++
  //     ++   ++   ++    ++        ++     ++
  ///////////////////////////////////////////

// get the min and max value for *MA, *MI and return total sum of 
// the number which is needed in 
//      ---   ---       ---
//      \     \    ...  \
//      /     /         /
//      ---   ---       ---

int get_MIN_MAX(int n, int *k, int *a, int *b, int *MI, int *MA)
{
  int j=1; int z=0;
  for (int i=0; i<n; i++) {
    MA[i] = MAX(z, k[i]-b[i]);
    MI[i] = MIN(a[i], k[i]);
    j *= MI[i] - MA[i] + 1;
  }
  
  return j;
}

// get the total number of n strings in the base
// It is: base[0]*base[1]* .. *base[n-1]

int get_total(int n, int *base)
{
  int j=1;

  for (int i=0; i<n; i++)
    j *= base[i];

  return j;
}

// hash is defined to determine the position of a point in the train
// of contpts with subscripts i1, i2, ... i_nd base

void hash(int i, int n, int *ps, int *base)
{
/*  if (n < 1) {
    cout<<" error in hash: n<1"<<endl;
    exit(1);
  }*/

  if (n == 1) {
    ps[0] = i;
    return;
  }
  else {
    int np = get_total(n-1, base);
    if (i >= np) {
      int rem  = i/np;

      ps[n-1] = rem;
      int next = i-rem*np;
      hash(next, n-1, ps, base);
    }
    else {
      ps[n-1]=0;
      hash(i, n-1, ps, base);
    }
  }
}

void rhash(int i, int n, int *ps, int *base)
{
  int *rbase=new int[n];
  for (int j=0; j<n; j++)
    rbase[n-j-1]=base[j];
  
  int *rps = new int[n];
  hash(i, n, rps, rbase);

  delete[] rbase;

  // Cho beg  
  //  for (j=0; j<n; j++)
  for (int j=0; j<n; j++)
    ps[n-j-1] = rps[j];
  // Cho end

  delete[] rps;
}

int rsort(int n, int *ps, int *base)
{
  int *rbase=new int[n];
  for (int i=0; i<n; i++)
    rbase[n-i-1]=base[i];

  int j = ps[n-1];

  // Cho beg
  //  for(i=1; i<n; i++)
  for(int i=1; i<n; i++)
    j += get_total(i, rbase) * ps[n-i-1];
  // Cho end

  delete [] rbase;
  return j;
}

// get the position in a string of contpts according to its *ps in base
// it is the reverse process of hase
// it is: 1+ps[0]+base[0]*ps[1]+base[1]*ps[2]+...+bsae[n-1]*ps[n]

int sort(int n, int *ps, int *base)
{
  int j = ps[0];
  for (int i=1; i<n; i++)
    j += get_total(i-1, base)*ps[i];

  return j;
}

void sub (geom &G1, geom &G2, geom &G3)
{
  for (int i=0; i<G1.nc; i++) {
    G3.contpts[i]->x = G1.contpts[i]->x - G2.contpts[i]->x;
    G3.contpts[i]->y = G1.contpts[i]->y - G2.contpts[i]->y;
    G3.contpts[i]->z = G1.contpts[i]->z - G2.contpts[i]->z;
    G3.contpts[i]->w = 1.0;
  }
}

void sub_array(int nd, int* MA, int* MI, int* MM)
{
  for (int i=0; i<nd; i++)
    MM[i] = MA[i] - MI[i];
}

void sub_array1(int nd, int* MA, int* MI, int* MM)
{
  for (int i=0; i<nd; i++)
    MM[i] = MI[i] - MA[i] +1 ;
}

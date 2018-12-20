
#include "simpoly.h"

using namespace std;

bagroot_list* bagrut_list_1d_from_rootlist_array(rootlist_array**);
bagroot_list* bagrut_list_2d_from_rootlist_array(rootlist_array**);
bagroot_list* bagrut_list_3d_from_rootlist_array_by_3rd_para(rootlist_array**);
bagroot_list* bagrut_list_3d_from_rootlist_array(rootlist_array**);
 
 
/*In this routine, if the two contiguous real number differ by only at 
most tolerance  then, treat them as one contiguous number */

bagroot_list*  separate_roots(rootlist *l, real tolerance)
{

  if (l->size() == 0) return NULL ;
  l->head();
  
  int dim = l->curent_pointer()->size();
  
  root_array** ra;
  rootlist_array** rla;

  ra = (root_array **) new root_array*[dim];
  /* In each dimension, ra[i] is an array containing all roots */
  
  rla = (rootlist_array **) new rootlist_array*[dim];
  /* In each dimension, rla[i] is an array containing several 
    rootlists, in which roots 
    in the same contiguous region are in the same rootlist */ 

  for (int cordnt=0; cordnt < dim; cordnt++) {
    ra[cordnt] = l->sort_wrt_cordnt(cordnt);
//     cout << "---------------------"<< cordnt<<"-th coordinate"<<endl;
//        cout<< *(ra[cordnt]) << endl;
    rla[cordnt] = bucket_each_coord(ra[cordnt], cordnt, tolerance);
//    for (int i=0; i < (rla[cordnt])->size(); i++) 
 //     cout << &((*(rla[cordnt]))[i]) ;
//	cout << "ra[ "<<cordnt<<" ] : " << *ra[cordnt] << endl;
//	cout << "rla[ "<<cordnt<<" ] : " << *rla[cordnt] << endl;

  }

  bagroot_list* brl;

  if (dim == 1) 
    brl = bagrut_list_1d_from_rootlist_array( rla);
  else if (dim == 2) 
    brl = bagrut_list_2d_from_rootlist_array( rla);
  else if (dim == 3) 
    brl = bagrut_list_3d_from_rootlist_array_by_3rd_para( rla);
  else
    { cout << " In separate_roots(): dim > 3 not supported! Use 'ria' instead of 'ria-sep'  " << endl;
      exit(1);
    }
  
//   for (int i=0; i < dim; i++) {
//     delete ra[i];
//     //delete rla[i];
//   }

  delete[] ra;
  delete[] rla;

  return brl;

}



bagroot_list*  separate_roots(rootlist *l)
{

  if (l->size() == 0) return NULL ;
  l->head();
  
  int dim = l->curent_pointer()->size();
  
  root_array** ra;
  rootlist_array** rla;

  ra = (root_array **) new root_array*[dim];
  /* In each dimension, ra[i] is an array containing all roots */
  
  rla = (rootlist_array **) new rootlist_array*[dim];
  /* In each dimension, rla[i] is an array containing several 
    rootlists, in which roots 
    in the same contiguous region are in the same rootlist */ 

  for (int cordnt=0; cordnt < dim; cordnt++) {
    ra[cordnt] = l->sort_wrt_cordnt(cordnt);
    // cout << "---------------------"<< cordnt<<"-th coordinate"<<endl;
    //    cout<< *(ra[cordnt]);
    rla[cordnt] = bucket_each_coord(ra[cordnt], cordnt);
//    for (int i=0; i < (rla[cordnt])->size(); i++) 
//      cout << &((*(rla[cordnt]))[i]) ;

  }


  bagroot_list* brl;

  if (dim == 1) 
    brl = bagrut_list_1d_from_rootlist_array( rla);
  else if (dim == 2) 
    brl = bagrut_list_2d_from_rootlist_array( rla);
  else if (dim == 3) 
    brl = bagrut_list_3d_from_rootlist_array_by_3rd_para( rla);
  else
    { cout << " In separate_roots(): dim > 3 not supported! Use 'ria' instead of 'ria-sep' " << endl;
      exit(1);
    }
  
//   for (int i=0; i < dim; i++) {
//     delete ra[i];
//     //delete rla[i];
//   }

  delete[] ra;
  delete[] rla;

  return brl;

}


int rlist_element::if_rlist_elem_in_rootlist(rootlist* rl) 
{
  rl->head();
  while(rl->curent_pointer()){
    if (rl->curent_pointer()->u == u) return TRUE;
    rl->next();
  }
  return FALSE;
}

/* if the two contiguous real number differ by only at most tolerance
then, treat them as one contiguous number */
rootlist_array *bucket_each_coord(root_array* ra, int cordnt, real tolerance)
{
  
  
  integerlist* il = new integerlist();
  
  if (ra->size() <= 0) return NULL;
  
  int counter_ra=1;
  il->insert(1);

  /* decide how many regions and  how many root in each region */
  while (counter_ra < ra->size() ) {
   real diff;
   
//   cout << "***************************************" << endl;
//   cout << " high is " << (*ra)[counter_ra][cordnt] << endl;
//   cout << " low is " <<  (*ra)[counter_ra -1][cordnt] << endl;    

   diff =  (*ra)[counter_ra][cordnt] - (*ra)[counter_ra -1][cordnt] ;    
   /* if the two contiguous real number differ by only at most tolerance
      then, treat them as one contiguous number */
//   cout << " differ is " << diff << endl;
//   cout << " tolerance is " << tolerance << endl;
//   cout << " compare diff < tolerance  is " << (diff < tolerance) << endl;

//   cout << "***************************************" << endl;
   if ( ((*ra)[counter_ra-1][cordnt] > (*ra)[counter_ra][cordnt] )
	 == MAYBE || diff < tolerance )
     {
       (*(il->curent_pointer()))() = il->curent_pointer()->value() + 1;
       
     }
   else {
     il->insert(1);
   }
   counter_ra++;      
 }
  
//  cout <<  il;

  /* put roots in the same region into the same rootlist */

  int counter = 0;
  rootlist_array*  rutare = new rootlist_array(il->size());
  il->head();
  
  for (int i=0; i < il->size(); i++) {
    for (int j=0; j < il->value(); j++) {
	 real_array tu(((*ra)[counter++]).root());
      
      (*rutare)[i].insert( tu);
    }
    il->next();
    if (!(il->curent_pointer())) break;
  }
  
//  for ( i=0; i < il->size(); i++) cout << &((*rutare)[i]) ;
  
  return rutare;
}



rootlist_array *bucket_each_coord(root_array* ra, int cordnt)
{
  
  
  integerlist* il = new integerlist();
  
  if (ra->size() <= 0) return NULL;
  
  int counter_ra=1;
  il->insert(1);

  /* decide how many regions and  how many root in each region */
  while (counter_ra < ra->size() ) {



   
   
   if ( ((*ra)[counter_ra-1][cordnt] > (*ra)[counter_ra][cordnt] )
	 == MAYBE)
     {
       (*(il->curent_pointer()))() = il->curent_pointer()->value() + 1;
       
     }
   else {
     il->insert(1);
   }
   counter_ra++;      
 }
  
//  cout <<  il;

  /* put roots in the same region into the same rootlist */

  int counter = 0;
  rootlist_array*  rutare = new rootlist_array(il->size());
  il->head();
  
  for (int i=0; i < il->size(); i++) {
    for (int j=0; j < il->value(); j++) {
	 real_array tu(((*ra)[counter++]).root());
      
      (*rutare)[i].insert( tu);
    }
    il->next();
    if (!(il->curent_pointer())) break;
  }
  
//  for ( i=0; i < il->size(); i++) cout << &((*rutare)[i]) ;
  
  return rutare;
}



bagroot_list* bagrut_list_1d_from_rootlist_array(rootlist_array** rla)
{

 bagroot_list* brl;
 rootlist_array* temp_rutlist;

 temp_rutlist = rla[0];

 int dim =  (((*temp_rutlist)[0]).root()).size();


  brl = NULL;
  brl = new bagroot_list();
  
  
  if (dim == 1) {
    rootlist_array* rla_ar = rla[0]; /* 0, since only one dimension */
    
    for (int i=0; i< rla_ar->size(); i++) { /* how many contiguous
					      regions  */
      rootlist* rl;
      rl = &((*rla_ar)[i]);
      rl->head();
      bagrootlist_element* be;
      be = new bagrootlist_element(dim);
      double low = 2.0;
      double high = -1.0;

      while(rl->curent_pointer()) { /* put roots in the same contiguous
				    region into the same bagroot */
#ifdef USE_INTERVAL
	if ((rl->root())[0].get_low() < low)
	  low = (rl->root())[0].get_low();
	if ((rl->root())[0].get_upp() > high)
	  high = (rl->root())[0].get_upp();
#else
	if ((rl->root())[0] < low)
	     low = (rl->root())[0];
	if ((rl->root())[0] > high)
	     high = (rl->root())[0];
#endif
	real_array tu(rl->root());
	be->bag().insert(tu);
	rl->next();
      }

      if (low > high) {
	cout << " running error in separate_roots()" << endl;
	exit(1);}

      /* assign to  bagrootlist_element (be) an as small as possible 
	region which cover all roots */
      Interval tol(low, high);
      (*be)[0] = tol;
      brl->insert(be); /* insert bagrootlist_elem into bagroot_list */ 
    }
  }

 return brl;

}



bagroot_list* bagrut_list_2d_from_rootlist_array(rootlist_array** rla)
{
  int i,j,i1;
 bagroot_list* brl;
 rootlist_array* temp_rutlist;

 temp_rutlist = rla[0];

 int dim =  (((*temp_rutlist)[0]).root()).size();

 brl = NULL;
 brl = new bagroot_list();

 
 
 if (dim == 2) {
   rootlist_array* rla_ar_x = rla[0];
   rootlist_array* rla_ar_y = rla[1];

   double *low, *high;
   low  = new double[dim];
   high = new double[dim];
   rootlist* rl_y;
   rootlist* rl_x; 
   
   for (i=0; i<rla_ar_x->size(); i++)
     for (j=0; j<rla_ar_y->size(); j++) { /* how many possible 
						contiguous regions  */
  
       rl_x = &((*rla_ar_x)[i]);
       rl_x->head();
       rl_y = &((*rla_ar_y)[j]);
       rl_y->head();						
       bagrootlist_element* be;
       be = new bagrootlist_element(dim);
       
       for (i1=0 ; i1 < dim ; i1++) {
	 low[i1] = 2.0;
	 high[i1] = -1.0;
       }
  
       rl_x->head();
       while(rl_x->curent_pointer()) { /* put roots in the same contiguous
				       region into the same bagroot */

	 if (rl_x->curent_pointer()->if_rlist_elem_in_rootlist(rl_y)) {
	   /* rl_x->curent_pointer is an common to rl_x and rl_y */
	   
	   /* take the covered region */
	   for (int ii = 0; ii < dim ; ii++) {
#ifdef USE_INTERVAL
	     if ((rl_x->root())[ii].get_low() < low[ii])
	       low[ii] = (rl_x->root())[ii].get_low();
	     if ((rl_x->root())[ii].get_upp() > high[ii])
	       high[ii] = (rl_x->root())[ii].get_upp();
#else 
	     if ((rl_x->root())[ii] < low[ii])
	       low[ii] = (rl_x->root())[ii];
	     if ((rl_x->root())[ii] > high[ii])
	       high[ii] = (rl_x->root())[ii];
#endif
	   }

	   for (i1=0; i1< dim ; i1++) 
	     if (low[i1] > high[i1]) {
	       cout << " running error ";
	       cout << "in bagrut_list_from_rootlist_array()" << endl;
	       exit(1);}
	   
	   /* assign to  bagrootlist_element (be) an as small as possible 
	     region which cover all roots */
	   for (i1=0 ; i1 < dim; i1++) {
	     real tol(low[i1], high[i1]);
	     (*be)[i1] = tol;
	   }
	   
	   real_array tu(rl_x->root());
	   be->bag().insert(tu);
	 }
              
	 rl_x->next();
       }// while (rl_x->)
      
       if (be->bag().size() > 0) 
	 brl->insert(be); /* insert bagrootlist_elem into bagroot_list */ 
     }
 delete low;
 delete high;
 }

 return brl;
}




bagroot_list* bagrut_list_3d_from_rootlist_array_by_3rd_para(
						  rootlist_array** rla)
{
  int i,j,k;
 bagroot_list* brl;
 rootlist_array* temp_rutlist;

 temp_rutlist = rla[0];

 int dim =  (((*temp_rutlist)[0]).root()).size();

 brl = NULL;
 brl = new bagroot_list();

 
 
 if (dim == 3) {
   rootlist_array* rla_ar_z = rla[2];
  


   double *low, *high;
   low  = new double[dim];
   high = new double[dim];
  
   rootlist* rl_z; 
   


       for ( k = 0; k < rla_ar_z->size(); k++) { /* how many possible 
						contiguous regions  */
  
	 			
	 rl_z = &((*rla_ar_z)[k]);
	 rl_z->head();						
	 bagrootlist_element* be;
	 be = new bagrootlist_element(dim);
	 
	 for (i=0 ; i < dim ; i++) {
	   low[i] = 2.0;
	   high[i] = -1.0;
	 }
	 
	 rl_z->head();
	 while(rl_z->curent_pointer()) { /* put roots in the same contiguous
					   region into the same bagroot */
	   
	   
	     
	     /* take the covered region */
	     for (int ii = 0; ii < dim ; ii++) {
#ifdef USE_INTERVAL
	       if ((rl_z->root())[ii].get_low() < low[ii])
		 low[ii] = (rl_z->root())[ii].get_low();
	       if ((rl_z->root())[ii].get_upp() > high[ii])
		 high[ii] = (rl_z->root())[ii].get_upp();
#else 
	       if ((rl_z->root())[ii] < low[ii])
		 low[ii] = (rl_z->root())[ii];
	       if ((rl_z->root())[ii] > high[ii])
		 high[ii] = (rl_z->root())[ii];
#endif
	     }
	     
	     for (i=0; i< dim ; i++) 
	       if (low[i] > high[i]) {
		 cout << " running error ";
		 cout << "in bagrut_list_from_rootlist_array()" << endl;
		 exit(1);}
	     
	     /* assign to  bagrootlist_element (be) an as small as possible 
	       region which cover all roots */
	     for (i=0 ; i < dim; i++) {
	       real tol(low[i], high[i]);
	       (*be)[i] = tol;
	     }
	     
	     real_array tu(rl_z->root());
	     be->bag().insert(tu);
	   
	   
	   rl_z->next();
	 }// while (rl_z->curent_pointer)
	 
	 if (be->bag().size() > 0) 
	   brl->insert(be); /* insert bagrootlist_elem into bagroot_list */ 
       }
 delete low;
 delete high;
 }
 
 return brl;
}






bagroot_list* bagrut_list_3d_from_rootlist_array( rootlist_array** rla)
{
  int i,j,k,i1;
 bagroot_list* brl;
 rootlist_array* temp_rutlist;

 temp_rutlist = rla[0];

 int dim =  (((*temp_rutlist)[0]).root()).size();

 brl = NULL;
 brl = new bagroot_list();

 
 
 if (dim == 3) {
   rootlist_array* rla_ar_x = rla[0];
   rootlist_array* rla_ar_y = rla[1];
   rootlist_array* rla_ar_z = rla[2];


   double *low, *high;
   low  = new double[dim];
   high = new double[dim];

   rootlist* rl_y;
   rootlist* rl_x; 
   rootlist* rl_z; 
   
   for (i=0; i<rla_ar_x->size(); i++)
     for ( j=0; j<rla_ar_y->size(); j++) 
       for ( k = 0; k < rla_ar_z->size(); k++) { /* how many possible 
						contiguous regions  */
  
	 rl_x = &((*rla_ar_x)[i]);
	 rl_x->head();
	 rl_y = &((*rla_ar_y)[j]);
	 rl_y->head();						
	 rl_z = &((*rla_ar_z)[k]);
	 rl_z->head();						
	 bagrootlist_element* be;
	 be = new bagrootlist_element(dim);
	 
	 for ( i1=0 ; i1 < dim ; i1++) {
	   low[i1] = 2.0;
	   high[i1] = -1.0;
	 }
	 
	 rl_x->head();
	 while(rl_x->curent_pointer()) { /* put roots in the same contiguous
					   region into the same bagroot */
	   
	   if (rl_x->curent_pointer()->if_rlist_elem_in_rootlist(rl_y) &&
	       rl_x->curent_pointer()->if_rlist_elem_in_rootlist(rl_z)) {
	     /* rl_x->curent_pointer is  common to rl_x, rl_y and rl_z */
	     
	     /* take the covered region */
	     for (int ii = 0; ii < dim ; ii++) {
#ifdef USE_INTERVAL
	       if ((rl_x->root())[ii].get_low() < low[ii])
		 low[ii] = (rl_x->root())[ii].get_low();
	       if ((rl_x->root())[ii].get_upp() > high[ii])
		 high[ii] = (rl_x->root())[ii].get_upp();
#else 
	       if ((rl_x->root())[ii] < low[ii])
		 low[ii] = (rl_x->root())[ii];
	       if ((rl_x->root())[ii] > high[ii])
		 high[ii] = (rl_x->root())[ii];
#endif
	     }
	     
	     for (i1=0; i1< dim ; i1++) 
	       if (low[i1] > high[i1]) {
		 cout << " running error ";
		 cout << "in bagrut_list_from_rootlist_array()" << endl;
		 exit(1);}
	     
	     /* assign to  bagrootlist_element (be) an as small as possible 
	       region which cover all roots */
	     for (i=0 ; i < dim; i++) {
	       real tol(low[i], high[i]);
	       (*be)[i] = tol;
	     }
	     
	     real_array tu(rl_x->root());
	     be->bag().insert(tu);
	   }
	   
	   rl_x->next();
	 }// while (rl_x->)
	 
	 if (be->bag().size() > 0) 
	   brl->insert(be); /* insert bagrootlist_elem into bagroot_list */ 
       }
	  delete low;
	  delete high;
 }
 
 return brl;
}

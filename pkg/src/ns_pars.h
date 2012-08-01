#ifndef NS_PARS_H
#define NS_PARS_H

#include <iostream>
#include <vector>
using namespace std;

class ns_pars
{
 public:
  vector<double> x;
  int n;
  
  double a;
  double b;
  vector<double> *pc, *pd;
  vector<size_t> *pwlevels;
  int l;
  double r;
  double rdelta;
  double rshape;
  double log_fwhh_prop_var;
  double log_fwhh_prior_mean;
  double log_fwhh_prior_var;
  double log_fwhh_re_prior_var;
  double log_fwhh_re_prop_var;
  double thresh;
  double steep;                  
  double freq;
  double shapestrleft;
  double shapestrright;

  int nlev;
  int p;
  int hs;
  vector<int> levsize;
  vector<double> h_vec;

 public:
  ns_pars& operator =( const ns_pars& A )
  {
      x.assign(A.x.begin(),A.x.end());
      n = A.n;
      a = A.a;
      b = A.b; 
      freq = A.freq;
      if(pc)
	  {
    	  delete pc;
    	  delete pd;
    	  delete pwlevels;
	  }
      if(A.pc)
	  {
    	  pc = new vector<double>(A.pc->size());
    	  *pc=*(A.pc);
    
    	  pd = new vector<double>(A.pd->size()); 
    	  *pd=*(A.pd);
    
    	  pwlevels = new vector<size_t>(A.pwlevels->size()); 
    	  *pwlevels=*(A.pwlevels);
	  }
      else
	  {
    	  pc=0;
    	  pd=0;
    	  pwlevels=0;
	  }
      l = A.l;
      r = A.r;
      rdelta = A.rdelta;
      rshape = A.rshape;
      log_fwhh_prior_mean=A.log_fwhh_prior_mean; 
      log_fwhh_prop_var=A.log_fwhh_prop_var;
      log_fwhh_prior_var=A.log_fwhh_prior_var;
      log_fwhh_re_prop_var=A.log_fwhh_re_prop_var;
      log_fwhh_re_prior_var=A.log_fwhh_re_prior_var;
      freq=A.freq;
      thresh = A.thresh;
      steep = A.steep;                  

      nlev=A.nlev;
      p=A.p;
      hs=A.hs;
      levsize.assign(A.levsize.begin(),A.levsize.end());
      h_vec.assign(A.h_vec.begin(),A.h_vec.end());;

      return * this;
  }   
  ns_pars( const ns_pars& A )
  {
      x.assign(A.x.begin(),A.x.end());
      n = A.n;
      a = A.a;
      b = A.b;
      pc = new vector<double>(A.pc->size());
      *pc=*(A.pc);
      pd = new vector<double>(A.pd->size()); 
      *pd=*(A.pd);
      pwlevels = new vector<size_t>(A.pwlevels->size()); 
      *pwlevels=*(A.pwlevels);
      l = A.l;
      r = A.r; 
      freq=A.freq;
      rdelta = A.rdelta;
      rshape = A.rshape;
      log_fwhh_prior_mean=A.log_fwhh_prior_mean; 
      log_fwhh_prop_var=A.log_fwhh_prop_var;
      log_fwhh_prior_var=A.log_fwhh_prior_var;
      log_fwhh_re_prop_var=A.log_fwhh_re_prop_var;
      log_fwhh_re_prior_var=A.log_fwhh_re_prior_var;
     
      thresh = A.thresh;
      steep = A.steep;                  

      nlev=A.nlev;
      p=A.p;
      hs=A.hs;
      levsize.assign(A.levsize.begin(),A.levsize.end());
      h_vec.assign(A.h_vec.begin(),A.h_vec.end());;    
  }  
  ns_pars()
  {
      pc=0;
      pd=0;
      pwlevels=0;
  }  
  ~ns_pars()
  {
      if(pc)
      {
    	  delete pc;
    	  delete pd;
    	  delete pwlevels;
	  }
  }
} ;
#endif

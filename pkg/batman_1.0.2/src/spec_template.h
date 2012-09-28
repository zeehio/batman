#ifndef SPEC_TEMPLATE_H
#define SPEC_TEMPLATE_H


#ifdef HPC
   #include <mkl_lapacke.h>
#include <mkl.h>

typedef MKL_INT MY_LAPACK_INT;
typedef double MY_LAPACK_DBL;
#else
   #include "my_lapack.h"
#endif

#include "ns_pars.h"
#include "metab_template.h"
#include "myheader.h"

class spec_template
{
 public:
  ns_pars pars;
  rngType *rng;
  bool Adapt;
  bool BurnIn;
  int tit;
 
  
  vector<double> datay;
  vector<double> sfitmean;
  int smean;
  int smeancount;
  
  vector<metab_template> FTems;
  matrix L;
  
  vector<double> dataWy;
  vector< vector<double> > B;
  vector< vector<double> > BBt_diags;
  vector< vector<double> > BBt_diags_tran;
  vector< vector<double> > vex;

  double log_fwhh_draw;
  vector<double> theta_draw;
  vector<double> theta_mean;
  
  double lambda_draw;
  vector<double> psi_draw;
  vector<double> phi_draw;
  vector<double> tau_draw;
  
  
  vector<double> log_fwhh_sam;
  matrix theta_sam;
  matrix sfit;
  matrix meta_sam;
  
  vector<double> lambda_sam;
  matrix phi_sam;
  matrix psi_sam;
  matrix tau_sam;
  vector<double> beta_draw;
  vector<double> beta_mean;
  matrix beta_sam;
  
  int it;
  vector<double> WtLbeta_sam;
  
  
  double shape_uplogsd;
  int shape_accepcount;
  double shape_adaptinterval;
  double noshadaptations;
  bool fixed;
  bool splitshape;
  bool largestore;
  double Temp;
  double repen;
  
 public:
  spec_template(vector<metab_template>* pFTems, ns_pars *ppars, vector<double> *py, rngType* rg)
  {    
      BurnIn=true;
      noshadaptations=0;
      datay = (*py);
      FTems=(*pFTems);
      pars = (*ppars);
      Adapt = false;
      largestore=true;
      splitshape=false;
      rng = rg;
      smean = 0;
      smeancount = 0;
      sfitmean.assign(pars.n,0.0);

      if(WAVELETS)
   	  { 
    	  const double ph[]={SYM6};
    	  pars.h_vec.assign(12,0);
    	  pars.hs=12;
    	  db_vit h_vit;
    	  unsigned int t=0;
    	 
    	  for(h_vit=pars.h_vec.begin();h_vit<pars.h_vec.end();h_vit++)
    	  {
    	      *h_vit=ph[t];
    	      t++;
    	  }
    	  
    	  pars.nlev=dwt_nlev(pars.n, pars.hs);
    	  pars.p=dwt_size(pars.n, pars.nlev, pars.hs);
    	  pars.pc=new vector<double>(pars.nlev+1);
    	  pars.pd=new vector <double>(pars.nlev+1);
    	  
    	  for(int t=0;t<pars.nlev+1; t++)
          {
    	      double mean=2*pow(1.05, (double)t);
    	      double  var=7*pow(1.05, (double)t);
    	      var=var*var/2;
    	      (*(pars.pd))[t]=mean*mean/var; 
    	      (*(pars.pc))[t]=mean/var;
    	  }
    	  pars.pwlevels=new vector <size_t>(pars.p);
    	      
    	  dataWy.assign(pars.p,0);
    	  pars.levsize.assign(pars.nlev+2,0);	 
    	 
    	  dwt_levvec(pars.n, pars.nlev,pars.hs, pars.levsize);
    	  
    	  int it1=0;
    	  
    	  for(unsigned int it2=0;it2<pars.levsize.size()-1;it2++)
  	      {
    	      for(int it3=0;it3<pars.levsize[it2];it3++)
    		  {
        		  (*(pars.pwlevels))[it1]=it2;
        		  it1++;
	          }
   	      }
    
    	  dwt(datay, pars.nlev, pars.h_vec, dataWy);
    	   
    	  theta_draw.assign(pars.p,0);
    	  theta_mean.assign(pars.p,0);
    	  
    	  theta_draw=dataWy;
    	  if(pars.p==pars.n)
    	  {
              printf("\nWavelet boundary error.\n");
              exit(1);
          }
    	  vector<double> unit(pars.p,0.0);
    	  vector<double> Bunit(pars.n,0.0);
    	  
    	  vector<double> bufvec(pars.n,0.0);
    	  
    	  vector<vector<double> > Mid; 
    	  vector<vector<double> > BtimesB(pars.p, unit);
    	      
    	  for(int lit=0;lit<pars.p;lit++)
    	  {
    	      if(lit>0)
    	    	unit[lit-1]=0.0;
    	      unit[lit]=1.0;
    	      
    	      Mid.push_back(unit);
    	      idwt(unit, pars.levsize, pars.h_vec, Bunit); 
    	      B.push_back(Bunit);
    	      for(int t=0;t<pars.n;t++)
    		    bufvec[t]=Bunit[t]*Bunit[t];
    	      BBt_diags.push_back(bufvec);
    	  }	  
    	  for(int lit=0;lit<pars.n;lit++)
    	  {
    	      BBt_diags_tran.push_back(unit); 
    	      for(int t=0;t<pars.p;t++)
    	   	     BBt_diags_tran[lit][t]=BBt_diags[t][lit];
    	  }
    	  Msq(&B, &BtimesB);
    	  MMdif(&Mid, &BtimesB, &Mid);
    	
          MY_LAPACK_DBL* pam = new MY_LAPACK_DBL[pars.p*pars.p]; 
          for(int it1=0;it1<pars.p;it1++)
    	    for(int it2=0;it2<pars.p;it2++)
  	        {
    		pam[it2+pars.p*it1]=Mid[it1][it2];
    	    }
    	  MY_LAPACK_INT* jpvt=new MY_LAPACK_INT[pars.p];
    	  MY_LAPACK_DBL* tau=new MY_LAPACK_DBL[pars.p];
    	  for(int it=0;it<pars.p;it++)
    	    jpvt[it]=0;
    	 
    #ifdef HPC
    	  LAPACKE_dgeqp3(LAPACK_COL_MAJOR, pars.p, pars.p, pam, pars.p, jpvt, tau); 
    	  LAPACKE_dorgqr(LAPACK_COL_MAJOR, pars.p, pars.p-pars.n, pars.p-pars.n, pam, pars.p, tau); 
    #else
    	  MY_LAPACK_INT lminus1=-1;
    	  MY_LAPACK_DBL size_work;
    	  MY_LAPACK_INT info=0;
    	  MY_LAPACK_INT lsize_work;
    	  MY_LAPACK_INT PARSP=pars.p;
    	  MY_LAPACK_INT DIFF=pars.p-pars.n;
    	  dgeqp3_(&PARSP, &PARSP, pam, &PARSP, jpvt, tau, &size_work, &lminus1, &info); 
    	  lsize_work=(int)size_work;

    	  MY_LAPACK_DBL* work = new MY_LAPACK_DBL[lsize_work];
    	  dgeqp3_(&PARSP, &PARSP, pam, &PARSP, jpvt, tau, work, &lsize_work, &info); 
    	  delete[] work;
    	  dorgqr_(&PARSP, &DIFF, &DIFF, pam, &PARSP, tau, &size_work, &lminus1, &info);   
    	  lsize_work=(int)size_work;

    	  work = new MY_LAPACK_DBL[lsize_work];
    	  dorgqr_(&PARSP, &DIFF, &DIFF, pam, &PARSP, tau, work, &lsize_work, &info);   
    	  delete[] work;
    #endif
    	  for(int it=0;it<pars.p-pars.n;it++)
    	    vex.push_back(vector<double> (pam+pars.p*it, pam + pars.p*(it+1))); 
    	  delete [] pam; 
   	  }
      else
	  {
    	  pars.p=pars.n;
    	  theta_draw.assign(pars.n,0);
    	  theta_mean.assign(pars.p,0);
    	  dataWy.assign(pars.p,0);
    	  dataWy=datay;
      }
      //initialize with mean

      psi_draw.assign((pars).p,1.0);
      phi_draw.assign((pars).p,1.0);
      lambda_draw=(pars.a/pars.b)/2.0;

      tau_draw.assign(pars.n,0);
      for (unsigned int i = 0; i<tau_draw.size(); i++)
	     tau_draw[i] = min(0.0,datay[i])-1.0/pars.steep;
     
      // start shape at mean of prior
      log_fwhh_draw=pars.log_fwhh_prior_mean;
      shape_uplogsd=-1;       
      
      beta_draw.assign((FTems).size(),0.0);
      beta_mean.assign((FTems).size(),0.0);

      Temp=1;
      fixed=false;

      vector<double> Ln (pars.p,0.0);
      L.resize(pars.l,Ln);   
      vector<double> dec(pars.p,0.0);
     
      
      if(pars.l>0)  
	  {
    	  for (int lit = 0; lit < pars.l; lit++)
          {      
    	      L[lit].assign(pars.p,0.0);
    	      FTems[lit].build_curve(pars.x, log_fwhh_draw, pars.freq); 
    	      FTems[lit].prop_accept_log_fwhh();
    	      if(WAVELETS)
    		  {  
        		  dwt(FTems[lit].curve, pars.nlev, pars.h_vec, dec); 
        		  L[lit] = dec; 
    		  }
    	      else
    		      L[lit] = FTems[lit].curve;
    	  }
	  }
    }

  void sample_metropolis_local(int dm, double re_fwhh);
  
  double calculate_metropolis_ratio_eta_local(matrix *prop);
  double gibbs_sample_lambda();

  void gibbs_sample_theta(vector<double> *theta_sample);
  void gibbs_sample_tau(vector<double> *tau_sample); 
  void gibbs_sample_tau_scale_free(vector<double> *tau_sample);
  double gibbs_sample_lambda_scale_free();
  
  void gibbs_sample_psi(vector<double> *psi_sample);

  void adapt_shape(double timepase);

  void sample_log_fwhh(); 
  void sample_log_fwhh_re(size_t t);
  void sample_pos(int t);
  void store_state();

  void sample_pos_theta_block(int t);
  void sample_theta_conditional_prop(matrix *Lprop, vector<double> *beta_prop, 
				     vector<double> *pthetaprop, double *plogthetaprprob);
  double sample_theta_conditional_revprop_logpdf();
  double calculate_metropolis_ratio_eta_theta(matrix *Lprop, vector<double> *thetaprop);
};
#endif

#ifndef CHAIN_TEMPLATE_H
#define CHAIN_TEMPLATE_H

#include "myheader.h"
#include <numeric>
#include <string>
#include <algorithm>
#include <iterator>

#include "parclass.h"
#include "spec_template.h"



class chain_template 
{
   public:
        rngType * rng;
        int it;
        parclass pars;
        matrix beta_sam;
        vector<double> beta_draw; 
        vector<double> beta_mean;
        double noadaptations;
        vector<spec_template> spectra;
        double Temp;
        double repen;
        bool BurnIn;
        double Adapt;
        ns_pars newspec_pars;
	    //rngType** p_rng;
    
    public:
	chain_template( vector<metab_template> * FTems,  matrix *data, rngType *rg, opt* options)
    {	    
        BurnIn = true;
        it = 0;
        rng = rg;

        pars.s = (*data).size()-1;      
        pars.l =(*FTems).size();
        pars.r = pow (10.0,(-8.0));
        pars.rdelta = (*options).rdelta;
        beta_draw.assign((*FTems).size(), pow(10.0,-6.0));
        beta_mean.assign((*FTems).size(), 0.0);

        Temp = 1;
        
        noadaptations = 0;
        
        newspec_pars.x = (*data)[0];
        newspec_pars.n = (*data)[0].size();

        newspec_pars.a=(*options).a;
        newspec_pars.b=(*options).b;

        newspec_pars.l = pars.l;
        newspec_pars.r=pars.r;
        newspec_pars.rdelta=(*options).rdelta;
        newspec_pars.rshape=(*options).rshape; 
		newspec_pars.log_fwhh_prior_mean=options->log_fwhh_prior_mean;
		newspec_pars.log_fwhh_prior_var=options->log_fwhh_prior_var;
		newspec_pars.log_fwhh_prop_var=options->log_fwhh_prop_var;
		newspec_pars.log_fwhh_re_prop_var=options->log_fwhh_re_prop_var;
		newspec_pars.log_fwhh_re_prior_var=options->log_fwhh_re_prior_var;
		newspec_pars.thresh=(*options).thresh;
        newspec_pars.steep=(*options).steep;
		newspec_pars.freq=options->spec_freq;

        spec_template st(FTems, &newspec_pars, &((*data)[1]), rng);
        spectra.resize(pars.s, st);

        for (int lit = 0; lit < pars.s; lit++)
        {
		  spec_template st(FTems, &newspec_pars, &((*data)[lit+1]), rng);                      
		  spectra[lit] = st; 
        }
        char strb[]="beta";
        up_date_children(strb);
    }
        
    void sample_metropolis(bool print, int dm, double re_fwhh);
    void gibbs_sample_beta(vector<double> * beta_sample);
    void sample_theta_beta_block(int t);

    double calculate_metropolis_ratio_eta_beta_theta(int s, matrix *Lprop, vector<double> *betaprop, matrix * thetaprop);
    double calculate_metropolis_ratio_beta_theta(vector<double> *betaprop, matrix *thetaprop);
    // this function sets the temperature parameter of the likelihood.
    void set_temp(double Temp);  
    void set_repen(double repen);    
    void set_adapt(double adapt); 
    void stop_burn();    
    void up_date_children(char par[]);
    void store_state();
    double sample_beta_conditional_revprop_logpdf(int t);
    void sample_eta_theta_beta_block(int s, int t);
    void sample_beta_noshift_conditional_prop(int t, double *beta_sample, double *logprob);
    void sample_beta_conditional_prop(int s, int t, matrix *WtLprop, double * beta_sample, 
    double *logprob);
};


#endif



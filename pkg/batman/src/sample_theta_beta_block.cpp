// written by Dr. Jie Hao, Dr William Astle
#include "chain_template.h"

void chain_template::sample_theta_beta_block(int t)
{      
    vector<double> betaprop(beta_draw);
    double logbetaprprob =  0 ;
    double logthetaprprob = 0;
    double sumlogthetaprprob = 0;
    double sumlogthetaprob = 0;
    
    sample_beta_noshift_conditional_prop(t, &(betaprop[t]), &logbetaprprob);
    //Wavelet
    vector<double> thetaprop (spectra[0].dataWy.size());
    matrix thetaprop_struct(pars.s, thetaprop);
    
    for (int sit = 0; sit < pars.s; sit++)
    {          
        spectra[sit].sample_theta_conditional_prop(&(spectra[sit].L), &betaprop, &thetaprop, &logthetaprprob);
        sumlogthetaprprob = sumlogthetaprprob+logthetaprprob;
        thetaprop_struct[sit] = thetaprop;
        sumlogthetaprob=sumlogthetaprob+spectra[sit].sample_theta_conditional_revprop_logpdf();
    }

    double logbetaprob = sample_beta_conditional_revprop_logpdf(t);
    double logratio = calculate_metropolis_ratio_beta_theta(&betaprop, &thetaprop_struct);

    //proposal
    double logproposal= logbetaprob- logbetaprprob ;
    logproposal=logproposal+(sumlogthetaprob-sumlogthetaprprob);

    logratio=logratio+logproposal;
     //logproposal
    double logacceptance=log(my_rand(rng));
    //to rig the position acceptance
    if (logratio>logacceptance)  
    {        
        //note storing here dont do this way get rid of
        beta_draw[t]=betaprop[t];

        for (int sit = 0; sit < pars.s; sit++)
        {
            for (unsigned int i = 0; i<spectra[sit].theta_draw.size(); i++)       
               spectra[sit].theta_draw=thetaprop_struct[sit];
        }
    }  
}

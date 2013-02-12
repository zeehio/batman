// written by Dr. Jie Hao, Dr William Astle
#include "spec_template.h"


void spec_template::sample_metropolis_local(int dm, double re_fwhh)
{
//main MCMC iterator for a spectrum
//gibbs sample the wavelet coefficients

    vector<double> theta_sample(theta_draw);
    gibbs_sample_theta(&theta_sample);
    theta_draw =  theta_sample;
  
    if(!fixed)
       sample_log_fwhh();
 
    if (re_fwhh < 0.0000001)
    { 
        for (int t = 0; t< pars.l; t++)
            sample_log_fwhh_re(t);
    } 
  
    if (dm == 0)
    {
        for (int t = 0; t< pars.l; t++)
           sample_pos_theta_block(t);
        
        for (int t = 0; t < pars.l; t++)
            sample_pos(t);
        /*if(BurnIn && Adapt)
        {  
            if(it % 50 == 0)
            {    
                for (int t = 0; t < pars.l; t++)
                    FTems[t].adapt_pos(50.0);
            }
        }*/
    }
 
    lambda_draw=gibbs_sample_lambda_scale_free();
    vector<double> psi_sample(psi_draw);
    gibbs_sample_psi(&psi_sample);
    psi_draw = psi_sample;
    vector<double> tau_sample(tau_draw);
    gibbs_sample_tau_scale_free(&tau_sample); 
    tau_draw=tau_sample;
}



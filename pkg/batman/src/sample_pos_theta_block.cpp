// written by Dr. Jie Hao, Dr William Astle
#include "spec_template.h"


void spec_template::sample_pos_theta_block(int t)
{
    double delta_prop = 0;
    double logacceptance = 0;
    double logthetaprprob = 0;
    double logproposal = 0;
    double logratio = 0;
    double logthetaprob = 0;

    vector<double> thetaprop(pars.p,0);

    vector<double> decon(pars.p,0.0);
    
    double frdelta = pars.rdelta;
    double delta_sigma = sqrt(Temp)*(frdelta)/5.0;

    for (unsigned int m = 0; m < FTems[t].no_mults; m++)
    {
        if (((FTems[t].multiplet_sites[m].fvar - 0)<0.0000001)&& ((FTems[t].multiplet_sites[m].fvar - 0)>-0.0000001))
        {
            continue;
        }
        else if (FTems[t].multiplet_sites[m].fvar > 0.0000001)
        {
            delta_sigma = sqrt(FTems[t].multiplet_sites[m].fvar);
    	}
        else if (FTems[t].multiplet_sites[m].frdelta > 0.0000001)
        {
            frdelta = FTems[t].multiplet_sites[m].frdelta;
            delta_sigma = sqrt(Temp)*(frdelta)/5.0;
        }
    	
        delta_prop=sample_truncated_t_twosided(FTems[t].multiplet_sites[m].delta_draw, 
        frdelta/2.0, 0.015, -frdelta, frdelta,rng);
        FTems[t].build_curve_prop_delta(pars.x,m,delta_prop, log_fwhh_draw, pars.freq);
        matrix Llocal(L);
       	if(WAVELETS)
    	{    
    	    dwt(FTems[t].curve_prop,pars.nlev, pars.h_vec, decon); 
    	    Llocal[t] = decon;  
    	}
    	else
        {
			Llocal[t]= FTems[t].curve_prop;
		}
    
        sample_theta_conditional_prop(&Llocal, &beta_draw, &thetaprop, &logthetaprprob);
        logthetaprob = sample_theta_conditional_revprop_logpdf();
        //ratio 
    
        logratio = calculate_metropolis_ratio_eta_theta(&Llocal, &thetaprop);
    
        //prior
        logratio=logratio+truncated_normal_twosided_logpdf(delta_prop, 0, 
        delta_sigma, -frdelta, frdelta);
        logratio=logratio-truncated_normal_twosided_logpdf(FTems[t].multiplet_sites[m].delta_draw, 
        0.0, delta_sigma, -frdelta,frdelta);
        //proposal
        logproposal=log(truncated_t_twosided_pdf(FTems[t].multiplet_sites[m].delta_draw, 
        delta_prop, frdelta/2, 0.015, -frdelta,frdelta));
        logproposal=logproposal-log(truncated_t_twosided_pdf(delta_prop, 
        FTems[t].multiplet_sites[m].delta_draw, frdelta/2, 0.015, 
        -frdelta,frdelta));
        logproposal=logproposal+(logthetaprob-logthetaprprob);
        logratio=logratio+logproposal;
    
        logacceptance = log(my_rand(rng));
        
        if (logratio>logacceptance)
        {
            // note storing here dont do this way get rid of
            FTems[t].prop_accept_delta(m,delta_prop, 0);
            L[t] = Llocal[t];
            theta_draw = thetaprop;
        }
    }

}

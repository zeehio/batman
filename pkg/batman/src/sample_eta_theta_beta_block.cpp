// written by Dr. Jie Hao, Dr William Astle
#include "chain_template.h"

void chain_template::sample_eta_theta_beta_block(int s, int t)
{
    int sit = 0;
    double delta_prop = 0;
    double logbetaprprob = 0;
    double logthetaprprob = 0;
    double sumlogthetaprprob = 0;
    double sumlogthetaprob = 0;
    double logbetaprob = 0;
    double logratio = 0;
    double logproposal = 0;
    double logacceptance = 0;
    vector<int> spec_indices(pars.s-1,0);
    double frdelta = pars.rdelta;
    double delta_sigma = sqrt(Temp)*(frdelta)/5.0;

    vector<double> decon(spectra[0].pars.p);
    for (unsigned int m = 0; m < spectra[s].FTems[t].multiplet_sites.size(); m++)
    { 
        if (spectra[s].FTems[t].multiplet_sites[m].fvar == 0)
        {
            continue;
        }
        else if (spectra[s].FTems[t].multiplet_sites[m].fvar > 0)
        {
            delta_sigma = sqrt(spectra[s].FTems[t].multiplet_sites[m].fvar);
        }
        else if (spectra[s].FTems[t].multiplet_sites[m].frdelta > 0)
        {
            frdelta = spectra[s].FTems[t].multiplet_sites[m].frdelta;
            delta_sigma = sqrt(Temp)*(frdelta)/5.0;
        }
        delta_prop=sample_truncated_t_twosided(spectra[s].FTems[t].multiplet_sites[m].delta_draw, 
        frdelta/2.0, 0.015, -frdelta, frdelta, rng);
        
        spectra[s].FTems[t].build_curve_prop_delta(spectra[s].pars.x, m, delta_prop, spectra[s].log_fwhh_draw, spectra[s].pars.freq);
        matrix Lprop(spectra[s].L);
    
        if(WAVELETS)  
    	{
    	  dwt(spectra[s].FTems[t].curve, spectra[s].pars.nlev, spectra[s].pars.h_vec, decon); 
    	  Lprop[t] = decon;
    	}
        else
        {
			Lprop[t] = spectra[s].FTems[t].curve_prop;
		}
    
        vector<double> betaprop(beta_draw);
        sample_beta_conditional_prop(s,t, &Lprop, &betaprop[t], &logbetaprprob);
        assign_ind(s, &spec_indices); 
        
        vector<double> thetaprop (spectra[0].dataWy.size(),0);
        matrix thetaprop_struct(pars.s, thetaprop);
    
        spectra[s].sample_theta_conditional_prop(&Lprop, &betaprop, &thetaprop, &logthetaprprob);
        sumlogthetaprprob = logthetaprprob;
        
        thetaprop_struct[s] = thetaprop;
        
        sumlogthetaprob = spectra[s].sample_theta_conditional_revprop_logpdf();
        for (unsigned int st = 0; st < spec_indices.size(); st++)
        {          
            sit = spec_indices[st];
            spectra[sit].sample_theta_conditional_prop(&(spectra[sit].L), &betaprop, &thetaprop, &logthetaprprob);
            sumlogthetaprprob = sumlogthetaprprob + logthetaprprob;
            thetaprop_struct[sit] = thetaprop;
            sumlogthetaprob = sumlogthetaprob + spectra[sit].sample_theta_conditional_revprop_logpdf();
        }
        logbetaprob = sample_beta_conditional_revprop_logpdf(t);
    
        logratio = calculate_metropolis_ratio_eta_beta_theta(s, &Lprop, &betaprop, &thetaprop_struct);
        
        //prior
        logratio=logratio+truncated_normal_twosided_logpdf(delta_prop, 0.0, 
        delta_sigma, -frdelta, frdelta);
        logratio=logratio-truncated_normal_twosided_logpdf(spectra[s].FTems[t].multiplet_sites[m].delta_draw,
        0.0, delta_sigma, -frdelta, frdelta);
      
        
        logproposal=log(truncated_t_twosided_pdf(spectra[s].FTems[t].multiplet_sites[m].delta_draw, 
        delta_prop, frdelta/2.0, 0.015,-frdelta,frdelta));
        logproposal=logproposal-log(truncated_t_twosided_pdf(delta_prop, spectra[s].FTems[t].multiplet_sites[m].delta_draw, 
        frdelta/2.0, 0.015,-frdelta,frdelta));
        logproposal=logproposal+(logbetaprob-logbetaprprob);
        logproposal=logproposal+(sumlogthetaprob-sumlogthetaprprob);
        logratio=logratio+logproposal;
        
        logacceptance=log(my_rand(rng));
        //logratio
        if (logratio>logacceptance)
        {   
           // note storing here dont do this way get rid of
            spectra[s].FTems[t].prop_accept_delta(m,delta_prop, 0);
            spectra[s].L[t] = Lprop[t];     
            beta_draw[t] = betaprop[t];
            for (sit = 0; sit< pars.s; sit++)
            {    
                 for (unsigned int i = 0; i<spectra[sit].theta_draw.size(); i++)
				 {
					 spectra[sit].theta_draw=thetaprop_struct[sit];
				 }
            }
            
        }
    }


}

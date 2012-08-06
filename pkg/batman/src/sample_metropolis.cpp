#include "chain_template.h"

void chain_template::sample_metropolis(bool print, int dm, double re_fwhh)
{
  it=it+1;
  char strt[]="it";
  up_date_children(strt);
  
  for (int lit = 0; lit < pars.s; lit++)
      spectra[lit].sample_metropolis_local(dm, re_fwhh);    

  vector<double> beta_sample(beta_draw);
  gibbs_sample_beta(&beta_sample);
  beta_draw =  beta_sample;
  char strb[]="beta";
  up_date_children(strb);
  
  for (int t = 0; t < pars.l; t++)
  {
      sample_theta_beta_block(t);
      up_date_children(strb);  
  }
  for (int s = 0; s < pars.s; s++)
  {
      for (int t = 0; t < pars.l; t++)
      {
    	  sample_eta_theta_beta_block(s,t);
    	  up_date_children(strb);
      }
  }
  if (!BurnIn || dm == 1)
  {        
    for (int s = 0; s < pars.s; s++)
    {
        vector<double> V(spectra[s].pars.n);
        matrix reconL(spectra[s].L.size());
      
    	VVsum (&(spectra[s].beta_mean), &(spectra[s].beta_draw), &(spectra[s].beta_mean));
    	VVsum (&(spectra[s].theta_mean), &(spectra[s].theta_draw), &(spectra[s].theta_mean));
    	spectra[s].smeancount ++;
    
        for (unsigned int d2 = 0; d2 < spectra[s].FTems.size();d2++)
        {
            for (unsigned int d3 = 0; d3 < spectra[s].FTems[d2].multiplet_sites.size(); d3++)
            {
            spectra[s].FTems[d2].multiplet_sites[d3].delta_draw_mean = 
            spectra[s].FTems[d2].multiplet_sites[d3].delta_draw_mean + spectra[s].FTems[d2].multiplet_sites[d3].delta_draw;
            }
        }
    }
    VVsum (&(beta_mean), &(beta_draw), &(beta_mean));
  }
}




#include "spec_template.h"


void spec_template::sample_pos(int t)
{
  double logratio = 0;
  double delta_prop = 0;
  double logproposal = 0;
  double logacceptance = 0;
  double delta_sigma2 = 0;
  double frdelta = pars.rdelta;
  double delta_sigma = sqrt(Temp)*(frdelta)/5.0;

  vector<double> decon(pars.p);
  for(unsigned int m=0;m<FTems[t].multiplet_sites.size();m++)
  {
      delta_sigma2 = exp(FTems[t].multiplet_sites[m].delta_uplogsd);
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
      
      delta_prop=sample_truncated_normal_gen(FTems[t].multiplet_sites[m].delta_draw, 
					     delta_sigma2, -frdelta, frdelta, rng);
      
      FTems[t].build_curve_prop_delta(pars.x,m,delta_prop, log_fwhh_draw, pars.freq);
      matrix Lprop(L);
      if(WAVELETS)
	  {
    	  dwt(FTems[t].curve_prop, pars.nlev, pars.h_vec, decon); 
    	  Lprop[t]= decon;
	  }
      else
	      Lprop[t]=FTems[t].curve_prop;
    
      logratio = calculate_metropolis_ratio_eta_local(&Lprop);

      // prior
      logratio=logratio+truncated_normal_twosided_logpdf(delta_prop, 0.0, delta_sigma, -frdelta, frdelta);
      logratio=logratio-truncated_normal_twosided_logpdf(FTems[t].multiplet_sites[m].delta_draw, 0.0, delta_sigma, -frdelta, frdelta);
      // logratio
      
      logproposal=truncated_normal_twosided_logpdf(FTems[t].multiplet_sites[m].delta_draw, 
						 delta_prop, delta_sigma2,-frdelta, frdelta);

      logproposal=logproposal-truncated_normal_twosided_logpdf(delta_prop, FTems[t].multiplet_sites[m].delta_draw, 
							     delta_sigma2,-frdelta, frdelta);
    
      logratio=logratio+logproposal;
    
      logacceptance=log(my_rand(rng));
    
      if (logratio>logacceptance)
      {   
	  FTems[t].prop_accept_delta(m,delta_prop, 1);
	  FTems[t].multiplet_sites[m].delta_accepcount=FTems[t].multiplet_sites[m].delta_accepcount+1;
	  L[t]= Lprop[t];
      }
   }
}

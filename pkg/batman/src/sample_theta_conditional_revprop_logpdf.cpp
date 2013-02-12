// written by Dr. Jie Hao, Dr William Astle
#include "spec_template.h"

double spec_template::sample_theta_conditional_revprop_logpdf()
{
    vector<double> Lbeta(L[0].size(),0);
    vector<double> meanvec(L[0].size(),0);

    MVprod(&L, &beta_draw, &Lbeta);
    double logprob=0;
    vector<double> var(psi_draw.size(),0);

    vector<double> meanbuf(pars.p);
    
    if(WAVELETS)
    { 
      for(int i = 0; i < pars.p; i++)
	  {
	    var[i]=(1/(1+psi_draw[i]));
	    meanbuf[i]=var[i]*(dataWy[i]-Lbeta[i]);
	  }
	  idwt(meanbuf, pars.levsize, pars.h_vec, meanvec);  
	  vector<double> y_draw(pars.n,0.0);
	
      idwt(theta_draw, pars.levsize, pars.h_vec, y_draw);  
	  for (int j = 0; j <pars.n; j++)
	  {  
	    logprob=logprob+truncated_normal_logpdf(y_draw[j], meanvec[j], sqrt(Temp/lambda_draw), tau_draw[j]);
	  } 
	
	  vector<double> addedpoints(pars.p-pars.n,0.0);
	  MVprod(&vex,&meanbuf,&addedpoints); 
	
	  vector<double> tmp2(pars.p-pars.n,0.0);
	  MVprod(&vex,&theta_draw,&tmp2);
	
	  VVdif(&tmp2,&addedpoints, &tmp2);
	  Vconstprod(&tmp2,sqrt(Temp/lambda_draw), &tmp2);
	
	  for(int j = 0; j <pars.p-pars.n; j++)
	    logprob+=normlogpdf(tmp2[j], 0, 1);
   }
   else
   {
      for(unsigned int i = 0; i<psi_draw.size(); i++)
  	  {
	    var[i] = (1/(1 + psi_draw[i]));
	    meanvec[i] = var[i]*(dataWy[i] - Lbeta[i]);
	  }      
  	  for (int j = 0; j < pars.p; j++) 
	    logprob=logprob+truncated_normal_logpdf(theta_draw[j], meanvec[j],sqrt(Temp*var[j]/lambda_draw), tau_draw[j]);
   }
   return logprob;
}



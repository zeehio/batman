#include "spec_template.h"

void spec_template::sample_theta_conditional_prop(matrix *Lprop, vector<double> *beta_prop,
 vector<double> *ltheta_draw, double *logprob)
{
  vector<double> Lbeta(pars.p,0.0);
  vector<double> var(pars.p,0.0);
  vector<double> meanvec(pars.n, 0.0);
  MVprod(Lprop,beta_prop, &Lbeta);
  vector<double> meanbuf(pars.p, 0.0);

  if(WAVELETS)
  {
    for(int i = 0; i < pars.p; i++)
	{
	  var[i]=(1/(1+psi_draw[i]));
	  meanbuf[i]=var[i]*(dataWy[i]-Lbeta[i]);
	}
    idwt(meanbuf, pars.levsize, pars.h_vec, meanvec);  
  }
  else
    for(int i = 0; i < pars.p; i++)
    {  
	var[i]=(1/(1+psi_draw[i]));
	meanvec[i]=var[i]*(dataWy[i]-Lbeta[i]);
    }
  *logprob=0;
  (*ltheta_draw).assign(pars.p,0);

  if(WAVELETS)
  { 
    vector<double> y_draw(pars.n,0.0);
    for (int j = 0; j <pars.n; j++)
	{
	  y_draw[j]=sample_truncated_normal_below(meanvec[j], sqrt(Temp/lambda_draw), tau_draw[j], rng);   
	  *logprob=*logprob+truncated_normal_logpdf(y_draw[j], meanvec[j], sqrt(Temp/lambda_draw), tau_draw[j]);
	}
    vector<double> addsam(pars.p-pars.n,0.0);
      
    for(int j = 0; j <pars.p-pars.n; j++)
	{
	  addsam[j]=my_normrnd(0, 1,rng); 
	  *logprob+=normlogpdf(addsam[j], 0, 1);
	}
    dwt(y_draw,pars.nlev,pars.h_vec, *ltheta_draw);
    vector<double> addedpoints(pars.p-pars.n,0.0);
  
    MVprod(&vex,&meanbuf,&addedpoints);
    vector<double> tmp2(pars.p-pars.n,0.0);
  
    MVprod(&vex,ltheta_draw,&tmp2);
    VVdif(&addedpoints, &tmp2, &tmp2);
    Vconstprod(&addsam,sqrt(Temp/lambda_draw), &addsam);
    VVsum(&addsam, &tmp2, &tmp2); 
    vector<double> tmp3(pars.p,0.0);
    MVprod(&vex,&tmp2,&tmp3);	  
      
    VVsum(&tmp3,ltheta_draw,ltheta_draw);   
    //      COME BACK HERE
  }
  else
  {  
    for (int j = 0; j <pars.p; j++)
	{  
	  (*ltheta_draw)[j]=sample_truncated_normal_below(meanvec[j], sqrt(Temp*var[j]/lambda_draw), tau_draw[j], rng); 
	  *logprob=*logprob+truncated_normal_logpdf((*ltheta_draw)[j], meanvec[j], sqrt(Temp*var[j]/lambda_draw), tau_draw[j]);
	}
 }   
}


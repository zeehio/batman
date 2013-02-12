// written by Dr. Jie Hao, Dr William Astle
#include "spec_template.h"

void spec_template::gibbs_sample_theta(vector<double> *theta_sample)
{
  double mean = 0;
  //gibbs sampling of wavelet coefficients
  double Inf = numeric_limits<double>::infinity();
  (*theta_sample) = theta_draw;
  
  int i =0;

  if(WAVELETS)
  {
      double toplim, botlim;
      vector<double> ratio(pars.n,0.0);
      
      bool anytoppers;    
      bool anybotters;
      int j;
      vector<double> diff(tau_draw);
      vector<double> wtrtmp(pars.n,0.0);
      idwt(theta_draw, pars.levsize, pars.h_vec, wtrtmp); 
      VVdif(&diff, &wtrtmp, &diff);
      
      vector<double> Lr(L.size(),0.0); 
      vector<double> prodbuf(pars.p,0.0);
      for(i = 0; i< pars.p; i++)
      {  
    	  anytoppers=false;
    	  anybotters=false;
    	  vector<bool> toppers(pars.n,false);
    	  vector<bool> botters(pars.n,false);	
        
    	  prodbuf=B[i];
    	  Vconstprod(&prodbuf,(*theta_sample)[i], &prodbuf);
    	  VVsum(&diff, &prodbuf,&diff);
        	  
    	  for(j=0;j<pars.n;j++)
          {
            if(B[i][j]<0)
    		{		  
                toppers[j]=true;
                if(!anytoppers)
                    anytoppers=true;
    		}
    	    if(B[i][j]>0)
    		{
                botters[j]=true; 
                if(!anybotters)
                    anybotters=true;
    		}
    	  }
    	  VVdivdot(&diff,&B[i],&ratio); 
    	  
    	  toplim=Inf;
    	  if(anytoppers)
    	    for(j=0;j<pars.n;j++)
    	      if((toppers[j])&&(toplim>ratio[j]))
        		toplim=ratio[j]; 
    	  botlim=-Inf;
    	  if(anybotters)
    	    for(j=0;j<pars.n;j++)
    	      if((botters[j])&&(botlim<ratio[j]))
        		botlim=ratio[j];
        		    
    	  if(toplim<=botlim)
    	  {	      
              double centre=0.5*(toplim+botlim);
              botlim=centre-numeric_limits<double>::min();
              toplim=centre+numeric_limits<double>::min();
        	  (*theta_sample)[i]=centre;		
    	  }
    	  else
    	  {	    
    	      Lrow (&L, i, &Lr);
    	      mean=(1.0/(1.0 + psi_draw[i]))*(dataWy[i]- VVprod(&Lr,&beta_draw));
    	      (*theta_sample)[i]=sample_truncated_normal_gen(mean, sqrt(Temp/(1.0+psi_draw[i])/lambda_draw), botlim,toplim, rng);
    	  }	
    	  prodbuf=B[i];
    	  Vconstprod(&prodbuf,(*theta_sample)[i], &prodbuf);
    	  VVdif(&diff, &prodbuf,&diff);				   
    } 
  }
  else
  {	
    vector<double> Lr(L.size(),0);
    for (int i = 0; i< pars.n; i++)
	{   
	    Lrow (&L, i, &Lr);
	    mean=(1.0/(1.0 + psi_draw[i]))*(datay[i] - VVprod(&Lr,&beta_draw));
	    (*theta_sample)[i]=sample_truncated_normal_gen(mean, sqrt(Temp/(1.0+psi_draw[i])/lambda_draw), tau_draw[i],Inf,rng);
	}
  }
}

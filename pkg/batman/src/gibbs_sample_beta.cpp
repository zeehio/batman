// written by Dr. Jie Hao, Dr William Astle
#include "chain_template.h"

void chain_template::gibbs_sample_beta(vector<double> *beta_sample)
{
     vector<double> v(spectra[0].L[0].size());
     vector<double> mvpd(v);
     vector<double> mvpd2(v);
     double prec = 0;
     double var = 0;
     double meanvec = 0;
     vector<int> temp_indices(pars.l-1,0);
     (*beta_sample)= beta_draw;
     for (int t = 0; t < pars.l; t++)
     {
        prec = pars.r;
        for (int sit = 0; sit < pars.s; sit++)
        {
            v = spectra[sit].L[t];
            prec=prec+ spectra[sit].lambda_draw*VVprod(&v, &v); 
        }
        var=1/prec;       
        assign_ind(t, &temp_indices);
        meanvec = 0;
        if(!temp_indices.empty())
        { 
           for (int sit = 0; sit < pars.s; sit++)
           {
               MVprodrange(&spectra[sit].L, beta_sample, &temp_indices ,&mvpd);
    
    	       for (unsigned int i = 0; i < mvpd.size(); i++)
                    mvpd2[i] = spectra[sit].dataWy[i]- spectra[sit].theta_draw[i]-mvpd[i];
               meanvec = meanvec + spectra[sit].lambda_draw*(VVprod(&spectra[sit].L[t], &mvpd2));
           }
        }
        else
        { 
            for (int sit = 0; sit < pars.s; sit++)
            {
                VVdif(&spectra[sit].dataWy, &spectra[sit].theta_draw, &mvpd2);     
                meanvec = meanvec + spectra[sit].lambda_draw*(VVprod(&spectra[sit].L[t], &mvpd2));
            }
        }
        
        meanvec=meanvec*var;
        (*beta_sample)[t]=sample_truncated_normal_below(meanvec, sqrt(Temp*var), 0.0, rng);
    }
  
}

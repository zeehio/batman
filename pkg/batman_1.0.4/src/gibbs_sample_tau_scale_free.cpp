#include "spec_template.h"

void spec_template::gibbs_sample_tau_scale_free(vector<double> *tau_sample)
{
    (*tau_sample).assign(pars.n, 0);
    double MInf=-numeric_limits<double>::infinity();
    vector<double> recon(pars.n);
    idwt(theta_draw, pars.levsize, pars.h_vec, recon); 
    
    for (int i = 0; i<pars.n; i++)
    {
     
    if(WAVELETS)
    	(*tau_sample)[i]=sample_truncated_normal_gen(pars.thresh, sqrt(Temp/lambda_draw/pars.steep),MInf, min(pars.thresh,recon[i]),rng);

    else
    	(*tau_sample)[i]=sample_truncated_normal_gen(pars.thresh, sqrt(Temp/lambda_draw/pars.steep),MInf, min(pars.thresh,theta_draw[i]),rng);
    }
}

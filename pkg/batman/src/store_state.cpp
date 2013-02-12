// written by Dr. Jie Hao, Dr William Astle
#include "spec_template.h"
void spec_template::store_state()
{
    if(largestore)
    {
        theta_sam.push_back(theta_draw); 
        lambda_sam.push_back(lambda_draw); 
        psi_sam.push_back(psi_draw);
        beta_sam.push_back(beta_draw);
    }
    if(!splitshape)
    { 
	    log_fwhh_sam.push_back(log_fwhh_draw);
    }
    for (int t = 0; t < pars.l; t++)
    {  
        FTems[t].store_state();
    }

    vector<double> V(datay.size());
    matrix reconL(L.size());

    for (unsigned int t2 = 0; t2<reconL.size(); t2++)
    {   
        idwt(L[t2], pars.levsize, pars.h_vec, reconL[t2]);
      	meta_sam.push_back(reconL[t2]);
    }
    MVprod( &reconL, &(beta_draw), &V);
		
	sfit.push_back(V);
}

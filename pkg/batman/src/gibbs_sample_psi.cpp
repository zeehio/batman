// written by Dr. Jie Hao, Dr William Astle
#include "spec_template.h"


void spec_template::gibbs_sample_psi(vector<double> *psi_sample)
{
    double scalepar;

    for (int i = 0; i <pars.p; i++)
	{
	  scalepar=2.0*Temp/(lambda_draw*pow(theta_draw[i],2.0)+(*(pars.pc))[(*(pars.pwlevels))[i]]);
	  (*psi_sample)[i]=my_gamrnd(((*(pars.pd))[(*(pars.pwlevels))[i]]-0.5)/Temp+1, scalepar, rng);
    }    
}



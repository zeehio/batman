#include "spec_template.h"

double spec_template::calculate_metropolis_ratio_eta_theta(matrix *Lprop, vector<double> *thetaprop)
{
    // written by Dr. Jie Hao, Dr William Astle
	//prior
    vector<double> pdd((*thetaprop).size(),0);
    vector<double> pdd2(theta_draw.size(),0);
    vector<double> mvpd(L[0].size(),0);
    vector<double> mvpd2(L[0].size(),0);
    
    vector<double> vec(pars.p,0);
    vector<double> vecprop(vec);
    
    VVproddot(thetaprop, &psi_draw, &pdd);
    VVproddot(&theta_draw, &psi_draw, &pdd2);
    
    double logratio =-0.5*lambda_draw*(VVprod(&pdd,thetaprop) -VVprod(&pdd2,&theta_draw));
    
    //likelihood
    MVprod(&L, &beta_draw, &mvpd);
    MVprod(Lprop, &beta_draw, &mvpd2);

    for (unsigned int i = 0; i < L[0].size(); i++)
    {
      vec[i] = dataWy[i]-theta_draw[i]-mvpd[i]; 
      vecprop[i] = dataWy[i]-(*thetaprop)[i]-mvpd2[i];
    }
    
    logratio = logratio-0.5*(VVprod(&vecprop, &vecprop)-VVprod(&vec, &vec))*lambda_draw;
    logratio = logratio/Temp;
    return logratio;
}

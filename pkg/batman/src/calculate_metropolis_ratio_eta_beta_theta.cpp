// written by Dr. Jie Hao, Dr William Astle
#include "chain_template.h"

double chain_template::calculate_metropolis_ratio_eta_beta_theta(int s, matrix *Lprop, vector<double> *betaprop, matrix * thetaprop)
{
	//prior
    double logratio = 0;
    vector<double> vvdf((*betaprop).size(),0);
    
    vector<double> vvdf2(spectra[0].pars.p,0);
    
    vector<double> vvdf3(vvdf2);
    vector<double> vsq(vvdf);
    vector<double> vsq2(vvdf);
    vector<double> mvpd(vvdf2);
    vector<double> mvpd2(vvdf2);
    
    vector<double> vec(vvdf2);
    vector<double> vecprop(vvdf2);
    vector<double> vvpdd(spectra[0].psi_draw.size(),0);
    vector<double> vvpdd2(spectra[0].psi_draw.size(),0);
    
    
    int sit = 0;
    
    Vpow(betaprop, 2.0, &vsq);
    Vpow(&beta_draw, 2.0, &vsq2);
    VVdif(&vsq, &vsq2, &vvdf);
    double sumv = accumulate( vvdf.begin(), vvdf.end(), 0.0);
    logratio = -0.5 * pars.r * sumv;
    
    for (sit = 0; sit < pars.s; sit++)
    {
        VVproddot (&((*thetaprop)[sit]), &spectra[sit].psi_draw, &vvpdd);
        VVproddot (&spectra[sit].theta_draw, &spectra[sit].psi_draw, &vvpdd2);
        logratio = logratio-0.5*spectra[sit].lambda_draw*(VVprod(&vvpdd, &((*thetaprop)[sit]))-VVprod(&vvpdd2, &(spectra[sit].theta_draw)));
    }
    
    vector<int> spec_indices(pars.s-1);
    assign_ind(s, &spec_indices);
    
    //likelihood
    VVdif(&(spectra[s].dataWy), &spectra[s].theta_draw, &vvdf2);
    MVprod(&spectra[s].L,&beta_draw,&mvpd);
    VVdif(&vvdf2, &mvpd, &vec);
    
    VVdif(&(spectra[s].dataWy), &((*thetaprop)[s]), &vvdf3);
    MVprod(Lprop,betaprop,&mvpd2);
    VVdif(&vvdf3, &mvpd2, &vecprop);
    
    logratio = logratio-0.5*(VVprod(&vecprop, &vecprop)-VVprod(&vec, &vec))* spectra[s].lambda_draw;
    for (unsigned int st = 0; st < spec_indices.size(); st++)
    {          
        sit = spec_indices[st]; 

        VVdif(&(spectra[sit].dataWy), &spectra[sit].theta_draw, &vvdf2);
        MVprod(&spectra[sit].L,&beta_draw,&mvpd);
        VVdif(&vvdf2, &mvpd, &vec);
        
        VVdif(&(spectra[sit].dataWy), &((*thetaprop)[sit]), &vvdf3);
        MVprod(&spectra[sit].L,betaprop,&mvpd2);
        VVdif(&vvdf3, &mvpd2, &vecprop);
    
        logratio = logratio-0.5*(VVprod(&vecprop, &vecprop)-VVprod(&vec, &vec))* spectra[sit].lambda_draw;
    }
    
    logratio=logratio / Temp;
    return logratio;
}

#include "chain_template.h"

double chain_template::calculate_metropolis_ratio_beta_theta(vector<double> *betaprop, matrix *thetaprop)
{
    //prior
    vector<double> vvdf((*betaprop).size(),0);
    //vector<double> vvdf2(spectra[0].datay.size(),0);
    vector<double> vvdf2(spectra[0].pars.p,0);
    vector<double> vvdf3(vvdf2);
    vector<double> vsq((*betaprop).size(),0);
    vector<double> vsq2((*betaprop).size(),0);
    
    vector<double> vvpdd(spectra[0].psi_draw.size(),0);
    vector<double> vvpdd2(spectra[0].psi_draw.size(),0);
    
    vector<double> mvpd(vvdf2);
    vector<double> mvpd2(vvdf2);
    
    vector<double> vec(vvdf2);
    vector<double> vecprop(vvdf2);
    
    
    Vpow(betaprop, 2.0, &vsq);
    Vpow(&beta_draw, 2.0, &vsq2);
    VVdif(&vsq, &vsq2, &vvdf);
    
    double logratio = -0.5 * pars.r * accumulate( vvdf.begin(), vvdf.end(), 0.0 );
    
    for (int sit = 0; sit < pars.s; sit++)
    {
        VVproddot(&((*thetaprop)[sit]), &spectra[sit].psi_draw, &vvpdd);
        VVproddot(&spectra[sit].theta_draw, &spectra[sit].psi_draw, &vvpdd2);
        logratio = logratio-0.5*spectra[sit].lambda_draw*(VVprod(&vvpdd, &((*thetaprop)[sit]))-VVprod(&vvpdd2, &spectra[sit].theta_draw));
    }
    
    
    for (int sit = 0; sit < pars.s; sit++)
    {   
        //VVdif(&(spectra[sit].datay), &spectra[sit].theta_draw, &vvdf2);
        VVdif(&(spectra[sit].dataWy), &spectra[sit].theta_draw, &vvdf2);
        MVprod(&spectra[sit].L,&beta_draw,&mvpd);
        VVdif(&vvdf2, &mvpd, &vec);
    
        //VVdif(&(spectra[sit].datay), &((*thetaprop)[sit]), &vvdf3);
        VVdif(&(spectra[sit].dataWy), &((*thetaprop)[sit]), &vvdf3);
        MVprod(&spectra[sit].L,betaprop,&mvpd2);
        VVdif(&vvdf3, &mvpd2, &vecprop);
        
        logratio = logratio-0.5*(VVprod(&vecprop, &vecprop)-VVprod(&vec, &vec))* spectra[sit].lambda_draw;
    }
    
    logratio = logratio/ Temp;
    return logratio;
}

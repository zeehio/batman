#include "chain_template.h"

void chain_template::sample_beta_noshift_conditional_prop(int t, double* beta_sample, 
double* logprob)
{
    vector<double> v(spectra[0].L[t].size());
    double meanvec=0;
    double prec = pars.r;
    for (int sit = 0; sit < pars.s; sit++)
    {   
        v = spectra[sit].L[t];
        prec = prec+spectra[sit].lambda_draw*VVprod(&v, &v);
    }

    double var = 1/prec;
    vector<int> temp_indices(pars.l-1,0);
    assign_ind(t, &temp_indices);
    
    meanvec=0;
    vector<double> b(pars.s*spectra[0].pars.n,0);
    vector<double> A(b.size(),0);

    vector<double> mvpd(spectra[0].pars.p,0);
    vector<double> vvdf(mvpd);
    vector<double> vvdf2(spectra[0].pars.n);

    if(!temp_indices.empty())
    {
        for (int sit = 0; sit < pars.s; sit++)
        {          
            MVprodrange(&(spectra[sit].L), &beta_draw, &temp_indices, &mvpd);
    
            VVdif(&(spectra[sit].dataWy), &mvpd, &vvdf);
            meanvec = meanvec + spectra[sit].lambda_draw*(VVprod(&spectra[sit].L[t], &vvdf));
    
    	    if(WAVELETS)
            { 
        		vector<double> recon(spectra[sit].pars.n);
        		idwt(vvdf, spectra[sit].pars.levsize, spectra[sit].pars.h_vec, recon); 
           		VVdif(&recon, &spectra[sit].tau_draw, &vvdf2);	
        		
        		idwt(spectra[sit].L[t], spectra[sit].pars.levsize, spectra[sit].pars.h_vec, recon); 
        		Vpartassign(&A, ((sit)*((spectra[0].pars).n)), &recon);
    	    }
    	    else
    	    {
        		VVdif(&vvdf, &spectra[sit].tau_draw, &vvdf2);
        		Vpartassign(&A, ((sit)*((spectra[0].pars).n)), &spectra[sit].L[t]);
    	    }
    	    Vpartassign(&b, ((sit)*((spectra[0].pars).n)), &vvdf2);                
        }
    }     
    else
    {
        for (int sit = 0; sit < pars.s; sit++)
        {    
    	    meanvec = meanvec+spectra[sit].lambda_draw*VVprod(&spectra[sit].L[t], &(spectra[sit].dataWy));
    	    vvdf=spectra[sit].dataWy;

    	    if(WAVELETS)
            { 
        		vector<double> recon(spectra[sit].pars.n);
        		idwt(vvdf, spectra[sit].pars.levsize, spectra[sit].pars.h_vec, recon); 
           		VVdif(&recon, &spectra[sit].tau_draw, &vvdf2);	
        		
        		idwt(spectra[sit].L[t], spectra[sit].pars.levsize, spectra[sit].pars.h_vec, recon); 
        		Vpartassign(&A, ((sit)*((spectra[0].pars).n)), &recon);
    	    }
    	    else
    	    {
        		VVdif(&vvdf, &spectra[sit].tau_draw, &vvdf2);
        		Vpartassign(&A, ((sit)*((spectra[0].pars).n)), &spectra[sit].L[t]);
    	    }
    	    Vpartassign(&b, ((sit)*((spectra[0].pars).n)), &vvdf2);     
        }
    }    
    meanvec = meanvec*var;
    vector<double> dd(A.size(),0);
    for(unsigned int i = 0; i < (A).size(); i++)
    (dd)[i] = (b)[i]/(A)[i];

    double mode = max(0.0,min((*std::min_element((dd).begin(),(dd).end())), meanvec));  
    double Inf = numeric_limits<double>::infinity();

    *beta_sample=sample_truncated_t_twosided(mode, sqrt(var), 1.0, 0.0, Inf, rng); 
    *logprob =log(truncated_t_twosided_pdf(beta_draw[t],mode, sqrt(Temp*var), 1.0, 0.0, Inf));
}

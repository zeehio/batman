#include "chain_template.h"

void chain_template::sample_beta_conditional_prop(int s, int t, matrix *Lprop, 
double * beta_sample, double *logprob)
{
    int sit = 0;
    vector<int> spec_indices(pars.s-1,0);
    
    assign_ind(s, &spec_indices);

    vector<double> v((*Lprop)[t]);
    v = (*Lprop)[t];
    double prec = 0;
    prec = pars.r+spectra[s].lambda_draw*(VVprod(&v,&v));

    for (unsigned int st = 0; st < spec_indices.size(); st++)
    {   
        sit = spec_indices[st];
        v = spectra[sit].L[t];
        prec=prec+spectra[sit].lambda_draw*(VVprod(&v, &v));
    }    
    
    double var = 1/prec;
    vector<int> temp_indices(pars.l-1,0);
    assign_ind(t, &temp_indices);
    
    vector<double> b((pars.s)*(spectra[0].pars.n),10.0);
    vector<double> A(b.size(),0.1);

    vector<double> mvpd(spectra[0].pars.p,0);
    vector<double> mvpd2(mvpd);
    vector<double> tem(spectra[0].pars.n,0);
    vector<double> tem2(mvpd);
    vector<double> vvdf(mvpd);
    vector<double> vvdf2(spectra[0].pars.n,0);

    double meanvec = 0;

    vector<double> recon(spectra[0].pars.n);
  
    if(!temp_indices.empty())
    {   
        MVprodrange(&spectra[s].L, &beta_draw, &temp_indices, &mvpd);
    
        VVdif(&(spectra[s].dataWy), &mvpd, &tem2);
        meanvec=spectra[s].lambda_draw*(VVprod(&((*Lprop)[t]), &tem2));
    
        if(WAVELETS)
    	{ 
    	  idwt(tem2, spectra[s].pars.levsize, spectra[s].pars.h_vec, recon); 
    	  VVdif(&recon, &spectra[s].tau_draw, &tem);
    	  
    	  idwt(spectra[s].L[t], spectra[s].pars.levsize, spectra[s].pars.h_vec, recon); 
    	  Vpartassign(&A, (s*(spectra[0].pars.n)), &recon);  
    	}
        else
    	{
    	  VVdif(&tem2, &spectra[s].tau_draw, &tem);
    	  Vpartassign(&A, (s*(spectra[0].pars.n)), &spectra[s].L[t]);
    	}    
        Vpartassign(&b, (s*(spectra[0].pars.n)), &tem);
    
        for (unsigned int st = 0; st <spec_indices.size(); st++)
        {   
            sit = spec_indices[st];
            MVprodrange(&spectra[sit].L, &beta_draw, &temp_indices, &mvpd2);
    
    	    VVdif(&(spectra[sit].dataWy), &mvpd2, &vvdf);
            meanvec=meanvec+spectra[sit].lambda_draw*(VVprod(&spectra[sit].L[t],&vvdf));
    
            if(WAVELETS)
    	    { 
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
        meanvec = spectra[s].lambda_draw*(VVprod(&((*Lprop)[t]), &(spectra[s].dataWy)));
        
        if(WAVELETS)
        { 
            idwt(spectra[s].dataWy, spectra[s].pars.levsize, spectra[s].pars.h_vec, recon); 
            VVdif(&recon, &spectra[s].tau_draw, &vvdf2);
            idwt(((*Lprop)[t]), spectra[s].pars.levsize, spectra[s].pars.h_vec, recon); 
            Vpartassign(&A, 0, &recon);	  
        }
        else
        { 
            VVdif(&(spectra[s].dataWy), &spectra[s].tau_draw, &vvdf2);
            Vpartassign(&A, 0, &((*Lprop)[t]));
        }
        Vpartassign(&b, 0, &vvdf2);
        
        for (unsigned int st = 0; st <spec_indices.size(); st++)
        {   
            sit = spec_indices[st];
            meanvec = meanvec + spectra[sit].lambda_draw*VVprod(&(spectra[sit].L[t]), &(spectra[sit].dataWy));
            
            if(WAVELETS)
            { 
              idwt(spectra[sit].dataWy, spectra[s].pars.levsize, spectra[sit].pars.h_vec, recon); 
              VVdif(&recon, &spectra[sit].tau_draw, &vvdf2);
              idwt(spectra[sit].L[t], spectra[sit].pars.levsize, spectra[sit].pars.h_vec, recon); 
              Vpartassign(&A, (sit*(spectra[0].pars.n)), &recon);
            }
            else
            {
              VVdif(&(spectra[sit].dataWy), &spectra[sit].tau_draw, &vvdf2);
              Vpartassign(&A, (sit*(spectra[0].pars.n)), &spectra[sit].L[t]);
            }
            Vpartassign(&b, (sit*(spectra[0].pars.n)), &vvdf2);
    	}  
     }

  meanvec=meanvec*var;
  vector<double> dd(A.size(),0);
  for(unsigned int i = 0; i < (A).size(); i++)
    (dd)[i] = (b)[i]/(A)[i];
  
  double mode = max(0.0,min((*std::min_element((dd).begin(),(dd).end())), meanvec));

  double Inf = numeric_limits<double>::infinity();
  *beta_sample=sample_truncated_t_twosided(mode, sqrt(Temp*var), 1.0, 0.0, Inf, rng);
  *logprob=log(truncated_t_twosided_pdf(beta_draw[t], mode, sqrt(Temp*var), 1.0, 0.0, Inf));
  
}

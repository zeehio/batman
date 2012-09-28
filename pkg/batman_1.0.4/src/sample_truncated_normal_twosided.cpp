#include "myheader.h"

double sample_truncated_normal_twosided(double mu, double sigma, double bot_lim,
double top_lim, rngType * rng)
// C. Roberts importance sampling for truncated normals (there is another slice-sampling paper not implemented here)
{
    double qb, qt, bot_lim_std, top_lim_std, z=0.0, sample=0.0, rho;
    bool qbc, qtc;
    int f;
    bot_lim_std=(bot_lim-mu)/sigma;
    top_lim_std=(top_lim-mu)/sigma;
    qb=sqrt(pow(bot_lim_std,2.0)+4.0);
    qt=sqrt(pow(top_lim_std,2.0)+4.0);
    qbc=(top_lim_std<bot_lim_std+2*exp(0.5)*exp((pow(bot_lim_std,2.0)-bot_lim_std*qb)/4.0)/(bot_lim_std+qb));
    qtc=(-bot_lim_std<-top_lim_std+2*exp(0.5)*exp((pow(top_lim_std,2.0)+top_lim_std*qt)/4.0)/(-top_lim_std+qt));
     
    if(((bot_lim_std*top_lim_std<=0)&&(top_lim_std-bot_lim_std<sqrt(2*M_PI)))||((bot_lim_std>0)&&qbc)||((top_lim_std<0)&&qtc))
    {
        rho=-1;
        f=0;
        // 'here'
        double mr = my_rand(rng);
        while(rho<mr)
        {
            f = f+1;
            mr = my_rand(rng);
            z=(top_lim_std-bot_lim_std)*mr+bot_lim_std;
            if (bot_lim_std>0)
            rho=exp((pow(bot_lim_std,2.0)-pow(z,2.0))/2.0);
            else
            {
                if(top_lim_std<0)
                    rho=exp((pow(top_lim_std,2.0)-pow(z,2.0))/2.0);
                else
                    rho=exp(-pow(z,2.0)/2.0);
            }
            mr = my_rand(rng);
        }
    sample=mu+sigma*z;
    return sample;
    }

    f=0;
    if((bot_lim_std*top_lim_std<=0.0))
    {  
        f=f+1;
        double t = my_rand(rng);
        sample=mu+sigma*my_norminv(my_normcdf(top_lim, mu, sigma)*t+(1-t)*my_normcdf(bot_lim, mu,sigma));
        if(f>1)
            cout<<"f "<< f <<endl;
        return sample;
    }
    f=0;
    if (bot_lim_std>0)
    {   f = f+1;
        z=top_lim_std+1;
        while(z>top_lim_std)
            z=sample_truncated_normal_below(0.0, 1.0, bot_lim_std,rng);  
        sample=mu+sigma*z;
        
        if(f>1)
            cout<<"f "<< f <<endl;
        return sample;
    }
    f=0;
    if (top_lim_std<0)
    {   f=f+1;
        z=-bot_lim_std+1;
        while(z>-bot_lim_std)
            z=sample_truncated_normal_below(0.0, 1.0, -top_lim_std,rng);
        
        sample=mu-sigma*z;
        return sample;
    }
    return sample;
}

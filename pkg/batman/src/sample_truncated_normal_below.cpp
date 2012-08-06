#include "myheader.h"
double sample_truncated_normal_below(double mu, double sigma, double bot_lim, 
rngType * rng)
// C. Roberts importance sampling for truncated normals (there is another slice-sampling paper not implemented here)
// check the paper to find out when it is optimal to use the two sided
// importance sampling alg
{
    double bot_lim_std, alpha_star, z, rho, sample;    
    bot_lim_std=(bot_lim-mu)/sigma;

    if (bot_lim_std > 0.0)
    {    
        alpha_star = 0.5 * (bot_lim_std + sqrt(pow(bot_lim_std,2.0) + 4.0));
        
        z = my_exprnd(1.0/alpha_star,rng);
        z =z + bot_lim_std;
        rho = exp(-0.5 * pow((z - alpha_star),2.0));
        
        double mr = my_rand(rng);
        while(rho < mr)
        {    
             //    'reject'
            z = my_exprnd(1.0/alpha_star, rng);
            z = z + bot_lim_std; 
            rho = exp(-0.5 * pow((z - alpha_star),2.0));
            mr = my_rand(rng);
        }
        sample = z * sigma + mu;
    }
    else
    {    
        double t = my_rand(rng);
        sample = mu + sigma*my_norminv( t + (1.0 - t) * my_normcdf(bot_lim, mu,sigma));
    }
    return sample;
}

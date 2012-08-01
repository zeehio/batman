#include "myheader.h"

double sample_truncated_normal_gen(double mu, double sigma, double bot_lim, 
double top_lim,rngType * rng)
/*C. Roberts importance sampling for truncated normals (there is another slice-sampling paper not implemented here)
check the paper to find out when it is optimal to use the two sided
importance sampling alg*/
{
    double sample = 0;
    if(bot_lim>top_lim)
    {
	  std::cerr<<"truncated norm limit error"<<std::endl;
	  exit(0);
    }
    if(bot_lim==top_lim)
    {
	  return bot_lim;
    }
    if((isinf(bot_lim))&&(isinf(top_lim)))
    {
      sample = my_normrnd(mu,sigma, rng);
      return sample;
    }
    if((!isinf(bot_lim))&&(isinf(top_lim)))
    {   
	  sample=sample_truncated_normal_below(mu, sigma, bot_lim,rng);
      return sample;
    }
    if((isinf(bot_lim))&&(!isinf(top_lim)))
    {   
      sample=sample_truncated_normal_above(mu, sigma, top_lim,rng);
      return sample;
    }
    sample=sample_truncated_normal_twosided(mu, sigma, bot_lim, top_lim, rng);
    return sample;
}

#include "myheader.h"
double sample_truncated_t_twosided(double mu, double scale, double df, 
double bot_lim, double top_lim, rngType *rng)
{
    double std_bot, std_top,z,sample;
    std_bot=(bot_lim-mu)/scale;
    std_top=(top_lim-mu)/scale;
    z = my_unifrnd(0,1, rng);
    sample = mu+scale*my_tinv(z*my_tcdf(std_top,df)+(1-z)*my_tcdf(std_bot,df),df);
    return sample;
}


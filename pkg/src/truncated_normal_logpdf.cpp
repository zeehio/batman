#include "myheader.h"
double truncated_normal_logpdf(double x, double mu, double sigma, double bot_lim)
{
    double log2=0.6931472;
    double sqrt2=1.414214;
    double logdensity = 0;
    double t = 0;
    double Inf = numeric_limits<double>::infinity();
    if(x<bot_lim)
        logdensity = -Inf;
    else
    {
        t = (bot_lim-mu)/sigma;
        if(mu<bot_lim)
            logdensity=normlogpdf(x, mu, sigma)-log(my_erfcx(t/sqrt2))+pow(t,2.0)/2.0+log2;
        else
            logdensity=normlogpdf(x,mu,sigma)-log(erfc(t/sqrt2))+log2;   
    }
    return logdensity;
}


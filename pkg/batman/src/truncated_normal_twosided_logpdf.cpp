// written by Dr. Jie Hao, Dr William Astle
#include "myheader.h"

double truncated_normal_twosided_logpdf(double x, double mu, double sigma, double bot_lim, double top_lim)
{
     double t, s;
     double logdensity =0;
     double log2=0.6931472;
     double sqrt2=1.414214;  
     double Inf = numeric_limits<double>::infinity();
     if((x<bot_lim)||(x>top_lim))
        logdensity = -Inf;
     else
     {
        t=(bot_lim-mu)/sigma;
        s=(top_lim-mu)/sigma;
     
        if(mu<bot_lim)
             logdensity=normlogpdf(x, mu, sigma)+log2-log(my_erfcx(t/sqrt2))+ pow(t,2.0)/2.0-log(1.0-(my_erfcx(s/sqrt2)/my_erfcx(t/sqrt2))*exp((pow(t,2.0)-pow(s,2.0))/2.0));
        else
        {   if(mu>top_lim)
                logdensity=truncated_normal_twosided_logpdf(-x,-mu, sigma, -top_lim, -bot_lim);
            else
                logdensity=normlogpdf(x, mu, sigma)-log(erfc(t/sqrt2)-erfc(s/sqrt2))+log2;
        }
     }
     return logdensity ;
}



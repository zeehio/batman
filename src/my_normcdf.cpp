#include "myheader.h"
double my_normcdf(double x, double mu, double sigma)
{
    double p;
    boost::math::normal myDist(mu, sigma);
    p = cdf(myDist, x);
    return p;
}

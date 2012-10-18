#include "myheader.h"

double my_unifrnd(double a, double b,rngType *rng) 
{
    double r = 0;
    boost::uniform_real<> unifr(a,b);      
    boost::variate_generator<boost::lagged_fibonacci607&, boost::uniform_real<> >
    var_u(*rng, unifr);             
    r = var_u();  
    return r;
}

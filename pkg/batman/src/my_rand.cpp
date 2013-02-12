// written by Dr. Jie Hao, Dr William Astle
#include "myheader.h"

double my_rand(rngType *rng)   
{   
    double r = 0;
    boost::uniform_real<> unifr(0.0,1.0);      
    boost::variate_generator<boost::lagged_fibonacci607&, boost::uniform_real<> >
             var_u(*rng, unifr);             
    r = var_u();  
    return r;

}

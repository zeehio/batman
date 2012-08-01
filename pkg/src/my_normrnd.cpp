#include "myheader.h"

double my_normrnd(double mu, double sigma, rngType *rng)   
{        
    //boost::mt19937 rng();    
    boost::normal_distribution<double> nd(mu, sigma);    
    boost::variate_generator<boost::lagged_fibonacci607&,
                             boost::normal_distribution<> > var_nor((*rng), nd);  
    return var_nor();  
}

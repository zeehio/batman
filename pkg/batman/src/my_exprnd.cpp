#include "myheader.h"

double my_exprnd( double mu,  rngType *rng) 
{   
    boost::exponential_distribution<> gd( 1/mu );  
    boost::variate_generator<boost::lagged_fibonacci607&,boost::exponential_distribution<> > var_exprnd( (*rng), gd );    
    return var_exprnd(); 
}

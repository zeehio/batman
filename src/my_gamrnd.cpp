#include "myheader.h"

double my_gamrnd( double shape, double scale,  rngType *rng) 
{   
    boost::gamma_distribution<> gd( shape );  
    boost::variate_generator<boost::lagged_fibonacci607&,
                             boost::gamma_distribution<> > var_gamma( (*rng), gd );    
    return scale*var_gamma(); 
}

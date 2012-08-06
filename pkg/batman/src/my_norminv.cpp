#include "myheader.h"

double my_norminv (double p)
{
    return boost::math::quantile(boost::math::normal(0.0, 1.0), p);  
}

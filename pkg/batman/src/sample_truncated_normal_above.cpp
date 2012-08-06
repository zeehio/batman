#include "myheader.h"

double sample_truncated_normal_above(double mu, double sigma, double top_lim, rngType * rng)
{
    return -sample_truncated_normal_below(-mu, sigma, -top_lim, rng);
}

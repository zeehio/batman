#include "myheader.h"
double MInf=-numeric_limits<double>::infinity();

double logmynchoosekfrac(vector<unsigned int>& n, vector<unsigned int>& k)
{
    double prodout=0.0;
    for (size_t it=0; it<n.size();it++)
      prodout+=lognchoosekfrac(n[it], k[it]);
    return prodout;
}

// maybe not the best way to calculate but ok for small integers
double lognchoosekfrac(unsigned int n, unsigned int k)
{
  unsigned int i;
  double ret=0.0;
  if (n<k )
    return MInf;
  else
    {
      for(i=0;i<k;i++)
        ret+=log((double)(n-k+i+1))-log((double)(i+1));
      ret-=(double)n*log(2.0);
    }
  return ret;
}

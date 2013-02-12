// written by Dr. Jie Hao, Dr William Astle
#include "myheader.h"

double my_betapdf(double x, double a, double b)
{
  double p;
  boost::math::beta_distribution<> myDist(a, b);
  p =pdf(myDist, x);
  return p;
}

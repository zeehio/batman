// written by Dr. Jie Hao, Dr William Astle
#include "myheader.h"

double my_tcdf(double t,double dof)
{
  double p;
  boost::math::students_t myDist(dof);
  p = cdf(myDist, t);	
  return p;
}



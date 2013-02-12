// written by Dr. Jie Hao, Dr William Astle
#include "myheader.h"

double my_tinv(double t, double dof)
{
  return boost::math::quantile(boost::math::students_t(dof), t);  
}


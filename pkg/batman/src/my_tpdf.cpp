#include <boost/math/distributions/students_t.hpp>

double my_tpdf(double t,double dof)
{
  double p;
  boost::math::students_t myDist(dof);
  p = pdf(myDist, t);	
  return p;
}



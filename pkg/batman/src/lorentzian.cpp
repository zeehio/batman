#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include "myheader.h"
using namespace std;

// Add a lorentzian to curve curve
void lorentzian(double mu, double fwhh, double height, vector<double> & x, vector<double> &y)
{
  vector<double> xtemp(x);
  vector<double>::iterator it;
  for (it=xtemp.begin();it<xtemp.end();it++)
      *it-=mu;
  
  for (it=y.begin();it<y.end();it++)
    *it+=height*(fwhh/M_PI)/(xtemp[it-y.begin()]*xtemp[it-y.begin()]+(fwhh*fwhh));   
 
}


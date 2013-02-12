// written by Dr. Jie Hao, Dr William Astle
#include <cmath>
double normlogpdf(double x, double mu, double sigma)
{
    double x_std,logpdf;
    x_std=(x-mu)/sigma;
    logpdf=-0.918938533-log(sigma)-0.5*pow(x_std,2.0); 
    return logpdf;
}

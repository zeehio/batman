#ifndef MYHEADER
#define MYHEADER
#define WAVELETS 1

#define BOOST_MATH_OVERFLOW_ERROR_POLICY errno_on_error //
#include <vector>
#define _USE_MATH_DEFINES
#include <cmath>
#include <ctime>
#include <iostream>
#include <string>
#include <numeric>
#include <iomanip>
#include <cstdlib>
#include <algorithm>
#include <cstdio>
#include <fstream>
#include <limits>
#include "opt.h"
#include "dwt.h"
//#include <omp.h>

#include "VVdif.h"

#include <boost/random.hpp> 
#include <boost/random/variate_generator.hpp>
#include <boost/random/exponential_distribution.hpp> 
#include <boost/random/normal_distribution.hpp>  
#include <boost/random/gamma_distribution.hpp>  
#include <boost/random/uniform_real.hpp>
 
#include <boost/math/special_functions/erf.hpp>
#include <boost/math/distributions/students_t.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/beta.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/type_traits/is_float.hpp>


//typedef boost::mt19937 rngType; 

typedef boost::lagged_fibonacci607 rngType; 
using namespace std;

typedef vector<double> row;
typedef vector<row> matrix;

typedef vector<int> rowI;
typedef vector<rowI> matrixI;

void read_txt_batmanoptions(opt *data, char filename[]);

void lorentzian(double mu, double fwhh, double height, vector<double> & x, vector<double> &y);
double lognchoosekfrac(unsigned int n, unsigned int k);
double logmynchoosekfrac(vector<unsigned int>& n, vector<unsigned int>& k);

void datamod(matrix * Data, double div, double lowbound, int step, int dppm, vector<double>* st, vector<double>* ed, matrix * data, int ps);
void ppm_ranges(int dppm, vector<double> * x, vector<double>* start, vector<double>* end, matrixI * range);


int read_datf(vector<string> * name, matrix *c1, matrix *c2, matrix *c3, matrix *c4, 
matrix *c5, matrix *c6, matrix *c7, matrixI *c8, char filename[]);
void read_txtf(matrix *data, char filename[]);
void read_txt_metalist(vector<string> *data, char filename[]);
 

void vecftoi(vector<double>& invec, vector<unsigned int>& outvec);
void Lrow (matrix * L, int ind, vector<double> * Lr);

void datamod(matrix * Data, double div, double lowbound, double step, double start, double end, matrix * data);

void assign_ind(int s, vector<int> *spec_indices);

void MVprod (matrix * M, vector<double> * V, vector<double> * P);
void MVprodrange (matrix * tem3, vector<double> * tem4, vector<int> *temp_indices, vector<double> * P);
void MMdif(matrix* M, matrix* V, matrix * P);
void Msq(matrix* M, matrix*);

void Vpow (vector<double> * V, double s, vector<double> * P);
void VVproddot (vector<double> * V, vector<double> * V2, vector<double> *P );
double VVprod (vector<double> * V, vector<double> * V2 );
void Vpartassign (vector<double> *b, int start, vector<double> *V);
double deltaadapt(double noadapt);
void VVsum (vector<double> * M, vector<double> * V, vector<double> * P);
void Vconstprod (vector<double> * V, double consta, vector<double> *P );
void VVdivdot (vector<double> * V, vector<double> * V2, vector<double> *P );
void VSdiv (vector<double> * V, double V2, vector<double> *P );

double sample_truncated_normal_below(double mu, double sigma, double bot_lim);
double truncated_normal_logpdf(double x, double mu, double sigma, double bot_lim);

double my_rand(rngType *rng);
double my_gamrnd( double shape, double scale,rngType *rng);
double my_betapdf(double x, double a, double b);

double my_erfcx(double x);
double my_tcdf(double t,double dof);
double my_tinv(double t, double dof);
double my_unifrnd(double a, double b,rngType *rng);
double my_tpdf(double t,double dof);
double my_erfcinv(double y);
double my_normcdf(double x, double mu, double sigma);
double my_norminv (double p);
double my_normrnd(double mu, double sigma,rngType *rng);
double my_exprnd(double mu, rngType * rng );   
bool my_isinf(double x);

double sample_truncated_normal_gen(double mu, double sigma, double bot_lim, double top_lim,
rngType * rng);
double sample_truncated_normal_above(double mu, double sigma, double top_lim, rngType * rng);
double sample_truncated_normal_below(double mu, double sigma, double bot_lim, rngType * rng);
double sample_truncated_normal_twosided(double mu, double sigma, double bot_lim, double top_lim, 
rngType * rng);

double normlogpdf(double x, double mu, double sigma);
double truncated_normal_twosided_logpdf(double x, double mu, double sigma, double bot_lim, double top_lim);
double truncated_normal_logpdf(double x, double mu, double sigma, double bot_lim);
double truncated_t_twosided_pdf(double x, double mu, double scale, double df, double bot_lim, double top_lim);
double sample_truncated_t_twosided(double mu, double scale, double df, double bot_lim,
double top_lim, rngType *rng);


double deltaadapt(double noadapt);

void write_txtf(vector<double> * V, char opfile[]);
void write_txtf2(vector<double> * V, vector<double> * V2, vector<double> *V3, char opfile[]);
void write_txtf_M(matrix *V, char opfile[]);


template < class T >
inline std::ostream& operator << (std::ostream& os, const std::vector<T>& v) 
{
  for (typename std::vector<T>::const_iterator ii = v.begin(); ii != v.end(); ++ii)
  {
      if((boost::is_float<T>::value)&&(ii!=v.begin()))
    	os<<",";
        os<< *ii;
  }
  os<<std::endl;
  return os;
}

//void printvector(vector<double> *ps);
#endif

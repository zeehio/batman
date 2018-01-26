// written by Dr. Jie Hao, Dr William Astle
#ifndef OPT_H
#define OPT_H

#include <vector>
using namespace std;

class opt
{
 public:
    string metalist;

    vector<double> st;
    vector<double> ed;
    vector<int> nospec;
    double lowerlimit;
    int downsample;
    int ito;
    int stop_burnin;
	int post_burnin;
  int paraProc;
  int tempF;
    int savefreq;
    double spec_freq;
    double div;
    double a;
    double b;
    //double c;
    //double d;
    double log_fwhh_mean;
    //double centre_shape;
    double rshape;
    double log_fwhh_prior_mean;
    double log_fwhh_prior_var;
    double log_fwhh_prop_var;
    double log_fwhh_re_prior_var;
    double log_fwhh_re_prop_var;
    double thresh;
    double steep;
    double rdelta;
    int fix;
    int rerun;
    int saveHR;
    int usechemshift;
    double tempe;
    long seed;

public:
    opt()
    {
        lowerlimit = -0.5;
        div = 200000;
        nospec.assign(1,1);
        downsample = 10;
        ito = 5000;
        stop_burnin = 200;
		post_burnin = 100;
    paraProc = 1;
    tempF = 2;
        savefreq = 50;
        spec_freq = 600;
        a=0.00001;
        b=0.000000001;
        //c=0.0001;
        //d=100;
        fix = 2;
        rerun = 5000;
        tempe = 1000;
        saveHR = 1;
	usechemshift = 0;

	    // these correspond to the fwhh
	    log_fwhh_prior_mean=0; 
	    log_fwhh_prior_var=0.1;
	    log_fwhh_prop_var=0.002;
	    log_fwhh_re_prior_var=0.0025;
	    log_fwhh_re_prop_var=0.0001;
	 
            
	    thresh=-0.03;
	    steep=300;
	    rdelta = 0.030;
	    seed = 123983400;
     }
} ;

#endif

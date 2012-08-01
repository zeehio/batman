#ifndef MULTIPLET_SITE_H
#define MULTIPLET_SITE_H

#include "myheader.h"
#define _USE_MATH_DEFINES
#include <cmath>

// multiplet site type
typedef int multi_type;
#define MULTI_PARAM 0
#define EXTRA_PARAM 1
#define MULTI_RASTER 2

class multiplet_site 
{ 
public:
	multi_type type;
	double no_protons; 
	double pos_ppm; 
	bool curve_built; 

	vector<unsigned int> couple_code;
	vector<double> Jconst;
	vector<double> weights;


	double delta_draw;

	double delta_draw_mean;
	int delta_mean_count;
	int dmean;

	vector <double> delta_sam;

	int dscount;
	double delta_uplogsd;
	int delta_accepcount;


	vector<double> curve;
	vector<double> curve_prop;

	double fpos;
	double fvar;
	double frdelta;
	int varcount;
	int varcount_tb;
	int varcount_etbb;

	//raster
	bool raster;
	double raster_width;
	vector<double> raster_curve;
	//vector<double> pos;


public:
	multiplet_site(vector<double> * pc_code,  double pos_ppm1,  vector<double> * pJconst, 
			double no_protons, vector<double>* x, double pfpos, double pfvar, double pfrdelta, opt* options) 
	: type(MULTI_PARAM), no_protons(no_protons),  curve_built(false), delta_draw(0.0)
	{
		delta_draw = 0.0;

		pos_ppm = pos_ppm1;

		fpos = pfpos;
		fvar = pfvar;
		frdelta = pfrdelta;
		if (!((fpos+50)<0.0000001))
		{
			pos_ppm = fpos;
			//cout<<"name "<< pos_ppm<< "-50 "<<(fpos!=50)<<endl;
		}
		varcount = 0;
		varcount_tb = 0;
		varcount_etbb = 0;

		delta_mean_count = 0;
		delta_draw_mean = 0;     
		dmean = 0;

		dscount = 0;
		delta_uplogsd=0.0;
		delta_accepcount=0;  
	}

	void setup_param(vector<unsigned int>& couple_code_local, vector<double>&  Jconst_local)
	{
		couple_code.assign(couple_code_local.begin(), couple_code_local.end());
		Jconst.assign(Jconst_local.begin(), Jconst_local.end()); 
	}
	void setup_param_extra(vector<double> Jconst_local, vector<double> weights_local)
	{
		type=EXTRA_PARAM;
		Jconst.assign(Jconst_local.begin(), Jconst_local.end()); 
		weights.assign(weights_local.begin(), weights_local.end());
	}
	~multiplet_site(){}

	void build_curve_param_nest(vector<double> &x, vector<double> &local_curve, double local_delta, double log_fwhh, double freq, size_t depth,  vector<unsigned int> it_vec)
	{
		if(type!=MULTI_PARAM)
		{
			std::cerr<<"I should be throwing an exception but instead you just got this message";
			exit(1);
		} 

		if(((depth==0)&&(std::accumulate(it_vec.begin(), it_vec.end(),0)))||(depth==Jconst.size()))
		{
			vector<double> diffcode(couple_code.begin(),couple_code.end());
			vector<double>::iterator it;
			for(it=diffcode.begin();it<diffcode.end();it++)
			{
				(*it)*=0.5;
				(*it)-=(double)it_vec[it-diffcode.begin()];
			}   
			double mu=pos_ppm+local_delta+VVprod(&Jconst, &diffcode)/freq;
			double height=(no_protons)*exp(logmynchoosekfrac(couple_code,it_vec));
			lorentzian(mu,exp(log_fwhh)/freq, height, x, local_curve);
		}   
		else
		{  
			for(unsigned int it=0;it<couple_code[depth]+1;it++)
			{
				vector<unsigned int> itcp_vec(it_vec);
				itcp_vec.push_back(it);
				build_curve_param_nest(x, local_curve, local_delta, log_fwhh, freq, depth+1, itcp_vec);
			}
		}							  
	}

	void build_curve(vector<double>& x, vector<double>& local_curve, double log_fwhh, double freq)
	{  
		build_curve_delta_prop(x, local_curve, delta_draw, log_fwhh, freq);
	}

	void build_curve_delta_prop(vector<double>& x, vector<double>& local_curve, double delta_prop, double log_fwhh, double freq)
	{ 
		if((!curve_built)||(x.size()!=curve_prop.size()))
		{  
			curve.assign(x.size(),0.0);
			curve_prop.assign(x.size(),0.0);
			curve_built=true;
		}	  
		if(type==MULTI_PARAM)
		{
			vector<unsigned int> empty;
			curve_prop.assign(x.size(),0.0);
			build_curve_param_nest(x, curve_prop, delta_prop, log_fwhh, freq, 0, empty);
		}
		if(type==EXTRA_PARAM)
		{
			double mu;
			double height;
			curve_prop.assign(x.size(),0.0);
			for(unsigned int it=0;it<Jconst.size();it++)
			{
				mu=pos_ppm+delta_prop+Jconst[it]/freq;
				height=no_protons*weights[it];
				lorentzian(mu,exp(log_fwhh)/freq, height, x, curve_prop);
			}
		}
		if (type == MULTI_RASTER)
		{
			build_curve_raster(&x, delta_prop, &curve_prop);
		}

		VVsum(&local_curve, &curve_prop,  &local_curve);
	}

	void prop_accepted_delta(double delta_prop, int inc)
	{            
		curve = curve_prop;
		delta_draw = delta_prop;
		delta_accepcount= delta_accepcount+inc;
	}

	void set_delta_uplogsd(double val)
	{          
		delta_uplogsd=val;
	}

	void zero_delta_accepcount()
	{       
		delta_accepcount=0;
	}

	void store_state()
	{        
		delta_sam.push_back(delta_draw);
		dscount++;
	}
	// for raster
	void raster_setup(double draster_width, vector<double>* praster_curve)
	{
		//raster = true;
		type = MULTI_RASTER;
		raster_width = draster_width;
		raster_curve.assign(praster_curve->size(),0.0);
		raster_curve = *praster_curve;
		Vconstprod(&raster_curve, ((double)raster_curve.size()/std::accumulate(raster_curve.begin(),raster_curve.end(),0.0)/raster_width), &raster_curve);
	}


	void build_curve_raster(vector<double>* x, double delta, vector<double>* thecurve)
	{
		thecurve->erase(thecurve->begin(),thecurve->end());
		thecurve->assign((*x).size(),0.0);
		double centre=pos_ppm+delta;
		double rast_len=(double) raster_curve.size();
		for(unsigned int it = 0; it < x->size();it++)
		{
			if(((*x)[it]>centre-raster_width/2.0)&&((*x)[it]<centre+raster_width/2.0))
			{
				//interpolation
				double rast_pos=rast_len*((*x)[it]-centre+raster_width/2.0)/raster_width;
				size_t pos_ceil=(size_t) ceil(rast_pos);
				size_t pos_floor=(size_t) floor(rast_pos);
				if((pos_floor>0)&&(pos_ceil<rast_len))
				{
					if(pos_ceil==pos_floor)
						(*thecurve)[it]=raster_curve[pos_floor];
					else
						(*thecurve)[it]=(pos_ceil-rast_pos)*raster_curve[pos_floor-1]+(rast_pos-pos_floor)*raster_curve[pos_ceil-1];
				}
			}
		}
		Vconstprod(thecurve, no_protons, thecurve);

		//end
	}

};

#endif

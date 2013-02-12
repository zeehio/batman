// written by Dr. Jie Hao, Dr William Astle
#ifndef METAB_TEMPLATE_H
#define METAB_TEMPLATE_H

#include "multiplet_site.h"
#include "myheader.h"

class metab_template
{
    public:
    string name;
    bool curve_built; 
    
    vector<multiplet_site> multiplet_sites;
    size_t no_mults;
    
    vector<double> curve;
    vector<double> curve_prop;
    
    double log_fwhh_re_draw;
    vector<double> log_fwhh_re_sam;
  
    public:
    metab_template(string name, int xsize)
       : name(name), curve_built(false), no_mults(0), log_fwhh_re_draw(0.0)
    {
        curve.assign(xsize,0);
        curve_prop.assign(xsize,0);
        multiplet_sites.erase(multiplet_sites.begin(),multiplet_sites.end());
    }
  
    void build_curve(vector<double>& x, double log_fwhh, double freq)
    { 	
        check_curve_build(x);
        vector<double> local_curve(x.size(), 0.0);
        for(unsigned int it=0;it<no_mults;it++)
            multiplet_sites[it].build_curve(x, local_curve, log_fwhh+log_fwhh_re_draw, freq);
        curve_prop=local_curve; 
    }
  
  
    void build_curve_prop_delta(vector<double> &x, unsigned int multi, double delta_prop,  double log_fwhh, double freq)
    {   
        check_curve_build(x);
        vector<double> local_curve(x.size(), 0.0);
        for(unsigned int it=0;it<no_mults;it++)
        {
    	if(it==multi)
    	  multiplet_sites[it].build_curve_delta_prop(x, local_curve, delta_prop, log_fwhh+log_fwhh_re_draw, freq);
    	else
    	  VVsum(&local_curve, &multiplet_sites[it].curve, &local_curve);
        }
        curve_prop=local_curve;
    }
  
    void build_curve_fwhh_re_prop(vector<double> &x, double log_fwhh_re_prop, double log_fwhh, double freq)
    { 
        check_curve_build(x);
        vector<double> local_curve(x.size(), 0.0);
        for(unsigned int it = 0; it < no_mults; it++)
          multiplet_sites[it].build_curve(x, local_curve, log_fwhh_re_prop+log_fwhh, freq);
        curve_prop=local_curve;  
    }
  
    void prop_accept_delta(int multi, double delta_prop, int inc)
    {
        multiplet_sites[multi].prop_accepted_delta(delta_prop, inc);
        curve = curve_prop;
    }
        
    void prop_accept_log_fwhh()
    {    
        for(unsigned int multi= 0; multi < no_mults; multi++)
          multiplet_sites[multi].curve=multiplet_sites[multi].curve_prop;
        curve = curve_prop;
    }

    void prop_accept_log_fwhh_re(double log_fwhh_re_prop)
    {    
        log_fwhh_re_draw=log_fwhh_re_prop;
        for(unsigned int multi= 0; multi < no_mults; multi++)
          multiplet_sites[multi].curve=multiplet_sites[multi].curve_prop;
        curve = curve_prop;
    }
    
    void add_multiplet(multiplet_site multiplet)
    {    
        multiplet_sites.push_back(multiplet);
        no_mults++;
    }
  
    void store_state()
    {    
        for(unsigned int it = 0; it < no_mults; it++)
            multiplet_sites[it].store_state();
        log_fwhh_re_sam.push_back(log_fwhh_re_draw);
    }
    private:
    void check_curve_build(vector<double> x)
    {
        if((!curve_built)||(x.size()!=curve_prop.size()))
        {
        	curve.assign(x.size(),0);
        	curve_prop.assign(x.size(),0);
        }
        curve_built=true;
    }
};
#endif

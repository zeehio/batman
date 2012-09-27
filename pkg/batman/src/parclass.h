#ifndef PARCLASS_H
#define PARCLASS_H

#include "metab_template.h"
#include "myheader.h"

class parclass
{
    public:
        int s;
        int l;
        double r;
        double rdelta;

    	vector<double> *pd, *pc;
    public:
        parclass& operator =( const parclass& A )
        {
             s = A.s;
             l = A.l;
             r = A.r;
             rdelta = A.rdelta;     
             return * this;
        }                 
};
#endif

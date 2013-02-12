// written by Dr. Jie Hao, Dr William Astle
#include "chain_template.h"

void chain_template::set_repen(double repen2)
{
    repen=repen2;
    for (int lit = 0; lit < pars.s; lit++)
        spectra[lit].repen=repen2;
}

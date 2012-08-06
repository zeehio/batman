#include "chain_template.h"

void chain_template::set_temp(double Temp1 )
{
    Temp=Temp1;
    for (int lit = 0; lit < pars.s; lit++)
      spectra[lit].Temp=Temp1;
}

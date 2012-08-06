#include "chain_template.h"

void chain_template::stop_burn()
{
    BurnIn=false;
    for (int lit = 0; lit < pars.s; lit++)
        spectra[lit].BurnIn=false;
}

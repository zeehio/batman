#include "chain_template.h"

void chain_template::set_adapt(double adapt)
{
    Adapt=adapt;
    for (int lit = 0; lit < pars.s; lit++)
        spectra[lit].Adapt=adapt;
}

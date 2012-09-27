#include "chain_template.h"

void chain_template::up_date_children(char par[])
{

    if(!strcmp(par, "beta"))
    {
        for( int s = 0; s<pars.s; s++)
            spectra[s].beta_draw=beta_draw;
    }
    if(!strcmp(par, "it"))
    {
        for (int s = 0; s<pars.s; s++)
            spectra[s].it = it;
    }

}

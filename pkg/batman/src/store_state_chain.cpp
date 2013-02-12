// written by Dr. Jie Hao, Dr William Astle
#include "chain_template.h"

void chain_template::store_state()
{ 
    beta_sam.push_back(beta_draw);
    vector<spec_template> locspectra(spectra);
    
    for (int lit = 0; lit < pars.s; lit++)
        locspectra[lit].store_state();
        
    spectra = locspectra;

}

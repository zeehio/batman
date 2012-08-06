#include <vector>
using namespace std;
#include "VVdif.h"

void MMdif (vector<vector<double> > * M, vector<vector<double> >* V, vector<vector<double> > * P)
{   
    for(unsigned int j = 0; j<(*M).size(); j++)
    {
      VVdif(&((*M)[j]),&((*V)[j]),&((*P)[j]));
    }
}

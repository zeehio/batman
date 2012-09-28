#include <vector>
using namespace std;

void VVproddot (vector<double> * V, vector<double> * V2, vector<double> *P )
{ 
    for (unsigned int i = 0; i<(*V).size(); i++)
        (*P)[i] =(*V)[i]*(*V2)[i];      
}


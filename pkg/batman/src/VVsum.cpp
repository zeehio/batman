#include <vector>
using namespace std;
void VVsum (vector<double> * M, vector<double> * V, vector<double> * P)
{   
    for(unsigned int j = 0; j<(*M).size(); j++)
    {
        (*P)[j] = (*M)[j]+(*V)[j];
    }
}

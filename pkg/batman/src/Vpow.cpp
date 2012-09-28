#include <vector>
#include <cmath>

using namespace std;

void Vpow (vector<double> * V, double s, vector<double> * P)
{   
    for (unsigned int i = 0; i<(*V).size(); i++)
        (*P)[i] = pow((*V)[i],s);
}


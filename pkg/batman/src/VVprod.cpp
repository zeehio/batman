// written by Dr. Jie Hao, Dr William Astle
#include <vector>
using namespace std;

double VVprod (vector<double> * V, vector<double> * V2 )
{   
    double P = 0;
  
    for (unsigned int i = 0; i<(*V).size(); i++)
    {
         P = P +  (*V)[i]*(*V2)[i];  
    }
    return P;
}


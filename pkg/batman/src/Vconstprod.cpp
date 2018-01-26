// written by Dr. Jie Hao, Dr William Astle
#include <vector>
using namespace std;

void Vconstprod (vector<double> * V, double consta, vector<double> *P )
{ 
    for (unsigned int i = 0; i<(*V).size(); i++)
      (*P)[i] =(*V)[i]*consta;      
}


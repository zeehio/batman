// written by Dr. Jie Hao, Dr William Astle
#include <vector>

using namespace std;

void Vpartassign (vector<double> *b, int start, vector<double> *V) 
{
    for (unsigned int i = 0; i < (*V).size(); i++)
    {
        (*b)[start+i] = (*V)[i];
    }    
}


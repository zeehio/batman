// written by Dr. Jie Hao, Dr William Astle
#include <vector>

using namespace std;

typedef vector<double> row;
typedef vector<row> matrix;

void MVprod (matrix * M, vector<double> * V, vector<double> * P)
{   
    double temp1 = 0;
    if((*M).size()!=(*V).size())
    {
        for(unsigned int j = 0; j<(*M).size(); j++)
        {
            temp1 = 0;
            for (unsigned int i = 0; i<(*V).size(); i++)
            {
                temp1 = temp1 +  (*M)[j][i]*(*V)[i];  
            }
            (*P)[j] = temp1;
        }
    }
    else
    {
        for(unsigned int j = 0; j<(*M)[0].size(); j++)
        {
            temp1 = 0;
            for (unsigned int i = 0; i<(*V).size(); i++)
            {
                temp1 = temp1 +  (*M)[i][j]*(*V)[i];  
            }
            (*P)[j] = temp1;
        }
    }
}

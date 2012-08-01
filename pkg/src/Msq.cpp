#include <vector>
#include <iostream>

using namespace std;

typedef vector<double> row;
typedef vector<row> matrix;

void Msq (matrix *M, matrix *P)
{   
    unsigned int it1, it2, it3;
    
    for(it1=0;it1<M->size();it1++)
    { 
        for(it2=0;it2<M->size();it2++)
        { 
            (*P)[it1][it2]=0.0;
            for(it3=0;it3<(*M)[0].size(); it3++)
            {
                (*P)[it1][it2]+=(*M)[it1][it3]*(*M)[it2][it3];
            }
        }
    }
}
	  

		
	  

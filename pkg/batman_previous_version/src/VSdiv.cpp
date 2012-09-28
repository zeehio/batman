#include <vector>
using namespace std;

void VSdiv (vector<double> * V, double V2, vector<double> *P )
{ 
    for (unsigned int i = 0; i<(*V).size(); i++)
       (*P)[i] =(*V)[i]/V2;      
}


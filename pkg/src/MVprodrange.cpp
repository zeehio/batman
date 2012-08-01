#include <vector>

using namespace std;
typedef vector<double> row;
typedef vector<row> matrix;

void MVprodrange (matrix * tem3, vector<double> * tem4, vector<int> *temp_indices, vector<double> * P)
{   
    double temp1 = 0;
    int ind = 0;
    for (unsigned int i = 0 ; i< (*tem3)[0].size(); i++)
    {
         temp1 = 0;
         for (unsigned int j = 0; j < (*temp_indices).size(); j++)  
         {
            ind = (*temp_indices)[j];
            temp1 = temp1 + ((*tem3)[ind][i])*((*tem4)[ind]);  
         }
         (*P)[i] = temp1;  
    }
}
  

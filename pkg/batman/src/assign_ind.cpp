#include <vector>
#include <iostream>

using namespace std;

void assign_ind(int s, vector<int> *spec_indices)
{
// index vector base on s
    int m = 0;
    for (unsigned int i = 0; i<(*spec_indices).size()+1; i++)
    {
        if ((int)i!=s)
        {
        (*spec_indices)[m]=i;
        ++m;
        }
    }
}

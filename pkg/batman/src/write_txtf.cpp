#include <vector>
#include <cstdio>
using namespace std;

typedef vector<double> row;
typedef vector<row> matrix;


void write_txtf(vector<double> * V, char opfile[])
{
    FILE *out;

    out = fopen(opfile,"w");
    for (unsigned int t = 0; t<(*V).size(); t++)
    {
        fprintf(out, "%.17g", (*V)[t]);
        if (t<=(*V).size()-1)
            fprintf(out, "\n");
    }
    fclose(out);
} 


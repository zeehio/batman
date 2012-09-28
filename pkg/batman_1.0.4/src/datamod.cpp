#include "myheader.h"

void datamod(matrix * Data, double div, double lowbound, int step, int dppm, 
vector<double>* st, vector<double>* ed, matrix * data, int s)
{
    double buf, buf2;
    matrixI range(2);// start[0] and end[1] 
    
    ppm_ranges(dppm, &(*Data)[0], st, ed, &range);
    
    for (unsigned int j = 0; j < range[0].size(); j++)
    {
        for (int i= range[0][j]; i<=range[1][j]; i= i+step)
        {
            buf = ((*Data)[0][i]);
            (*data)[0].push_back(buf);
            
            if ((*data).size()==2)
            {
                buf2 = ((*Data)[s+1][i])/div; 
                if (buf2<lowbound)
                buf2 = lowbound;     
                (*data)[1].push_back(buf2);
            }
            else 
            {
                for (unsigned int ii = 1; ii<(*data).size(); ii++)
                {
                    buf2 = ((*Data)[ii][i])/div; 
                    if (buf2<lowbound)
                    buf2 = lowbound;     
                    (*data)[ii].push_back(buf2);
                }
            }
        }
    }   
}


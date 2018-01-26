// written by Dr. Jie Hao, Dr William Astle
#include "myheader.h"

void read_txt_metalist(vector<string> *data, char filename[])
{
    string s;
    ifstream fin(filename);  
    while(getline(fin, s))
    {
        if (s[0]!='%') 
        {
            (*data).push_back(s);
        }
    }
}


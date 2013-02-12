// written by Dr. Jie Hao, Dr William Astle
#include "myheader.h"

void Lrow (matrix * L, int ind, vector<double> * Lr)
{
 for(unsigned int i = 0; i <(*L).size(); i++)
 {
    (*Lr)[i] = (*L)[i][ind]; 
 }   
}


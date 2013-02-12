#include "myheader.h"

double deltaadapt(double noadapt)
{
	// written by Dr. Jie Hao, Dr William Astle
    return min(0.5,3/sqrt(noadapt));
}


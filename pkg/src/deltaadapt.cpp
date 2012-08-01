#include "myheader.h"

double deltaadapt(double noadapt)
{
    return min(0.5,3/sqrt(noadapt));
}


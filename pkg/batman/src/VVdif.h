#ifndef VVDIFF_H 
#define VVDIFF_H

template <typename T>
void VVdif(T* M, T* V, T* P)
{   
    for(unsigned int j = 0; j<(*M).size(); j++)
    {
        (*P)[j] = (*M)[j]-(*V)[j];
    }
}

#endif

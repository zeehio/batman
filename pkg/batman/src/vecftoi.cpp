#include "myheader.h"

using namespace std;

void vecftoi(vector<double>& invec, vector<unsigned int>& outvec)
{
  vector<double>::iterator it=invec.begin();
  while(it!=invec.end())
  {
      outvec[it-invec.begin()]=floor(*it+0.5);
      it++;
  }
}

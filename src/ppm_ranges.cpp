#include "myheader.h"

void ppm_ranges(int dppm, vector<double> * x, vector<double>* start, vector<double>* end, matrixI * range)
{
     vector<double> st;
     vector<double> ed;

     double df = (*x)[0] - (*x)[1];
     if (dppm == -1)
     {
         for (unsigned int n = 0; n < (*start).size(); n++)  
         {
             if ((*start)[n]<max((*x)[0],(*x)[(*x).size()-1]) && (*start)[n]>min((*x)[0],(*x)[(*x).size()-1]))
             {
                  if (n+1<=((*end).size()-1))
                  {
                  st.push_back((*end)[n+1]+abs(df/2));
                  ed.push_back((*start)[n]-abs(df/2));
                  }
             }
         }
     }
     else
     {
         st = *start;
         ed = *end;
     }
     if (df<0.0)
     {
         matrixI rg((*range).size());
         for (unsigned int p = 0; p < st.size(); p++)
         {
             for (unsigned int i = 0; i<(*x).size(); i++) 
             {
                 if (i == (*x).size()-1)
                 { 
                     if ((*x)[i]<=ed[p] && (*x)[i]>=st[p])
                     rg[1].push_back(i);
                     else if ((*x)[i]<st[p])
                     continue;         
                 }
                 else if (i==0)
                 {
                     if ((*x)[i]>=st[p] && (*x)[i]<=ed[p])
                     rg[0].push_back(i);
                     else if ((*x)[i]>ed[p])
                     continue;
                 }
                 else
                 {
                     if ((*x)[i]<=ed[p] && (*x)[i+1]>ed[p])
                     rg[1].push_back(i);
                     if ((*x)[i]<st[p] && (*x)[i+1]>=st[p])
                     rg[0].push_back(i+1);
                 }  
             }                         
         }
         for (unsigned int j = 0; j < rg[0].size(); j++)
         {
             for(unsigned int jj = 0; jj < rg.size(); jj++)
                (*range)[jj].push_back(rg[jj][rg[0].size()-j-1]);
         }
     }
     else
     {
         for (unsigned int p = 0; p < st.size(); p++)
         {
             for (unsigned int i = 0; i<(*x).size(); i++)
             {
                 if (i == 0)
                 { 
                   if ((*x)[i]<=ed[p] && (*x)[i]>=st[p])
                   (*range)[0].push_back(i);
                   else if ((*x)[i]<st[p])
                   continue;
                 }
                 else if (i==(*x).size()-1)
                 {
                  if ((*x)[i]>=st[p] && (*x)[i]<=ed[p])
                  (*range)[1].push_back(i);
                  else if ((*x)[i]>ed[p])
                  continue;
                 }
                 else
                 {
                  if ((*x)[i]>ed[p] && (*x)[i+1]<=ed[p])
                  (*range)[0].push_back(i+1);
                  if ((*x)[i]>=st[p] && (*x)[i+1]<st[p])
                  (*range)[1].push_back(i);
                 }         
             }
         }
     }         
}

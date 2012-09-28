#include "myheader.h"
#include "metab_template.h"

void read_multiplet_data(int lineno, char filename[], opt* opts, 
vector<string> *listname, vector<metab_template> *Tems, vector<double> *x)
// read in metabolite multiplet data, note at present this assumes the file
// is ordered in groups corresponding to different metabolites
{
  vector<string> names(lineno);// name characters the first line
  matrix c1(lineno);
  matrix c2(lineno);
  matrix c3(lineno);
  matrix c4(lineno);
  matrix c5(lineno);
  matrix c6(lineno);
  matrix c7(lineno);
  matrixI c8(lineno);
  int count = 0;
  double pst, ped;
  
  int nl = read_datf(&names,&c1, &c2, &c3, &c4, &c5, &c6, &c7, &c8, filename); 


  char prevname[80]=" ";
  char name[80] = {'\0'};
  char lis[80] = {'\0'};
  metab_template templa("",(*x).size());

  for (unsigned int i  = 0; i < (*listname).size(); i++)
  {
      strcpy(lis,(*listname)[i].c_str());

      for (int it = 0; it < nl; it++)
      {
    	  strcpy(name,names[it].c_str());
    	  if (!strcmp(lis,name) && c8[it][0]==1) // find a match in the list
          {   
    	      if (!strcmp(prevname, " ")&& c8[it][0]==1)//first prevname do not match prevname
              {
        		  metab_template templa1(name,(*x).size());
        		  templa = templa1;
              }
    	      
    	      if (strcmp(prevname,name) && strcmp(prevname, " ") && c8[it][0]==1) // make a new template
              {  
        		  metab_template templa1(name,(*x).size());
        		  templa = templa1;
              }
    
              strcpy(prevname,name);
              
              for (unsigned int n2 = 0; n2<(*opts).st.size(); n2++)
                {
                    if ((*opts).st[n2]>(*opts).ed[n2])
                    {
                    pst = (*opts).ed[n2];
                    ped = (*opts).st[n2];
                    } else {
                    ped = (*opts).ed[n2];
                    pst = (*opts).st[n2];
                    }
                    if ((c5[it][0]+50)<0.0000001)
                    {
                        if (!((c7[it][0]+50)<0.0000001))
                        {
                            if (c1[it][0]<ped+(15.0*(*opts).log_fwhh_prop_var)+c7[it][0]
                            && c1[it][0]>pst-(15.0*(*opts).log_fwhh_prop_var)-c7[it][0])
                            {
                            //cout<<c1[it][0];
                            count = 1;
                            }  
                        } else {
                            if (c1[it][0]<ped+(15.0*(*opts).log_fwhh_prop_var)+(*opts).rdelta
                            && c1[it][0]>pst-(15.0*(*opts).log_fwhh_prop_var)-(*opts).rdelta)
                            {
                            //cout<<c1[it][0];
                            count = 1;
                            }
                        }
                    } else {
                        if (!((c7[it][0]+50)<0.0000001))
                        {
                            if (c5[it][0]<ped+(15.0*(*opts).log_fwhh_prop_var)+c7[it][0]
                            && c5[it][0]>pst-(15.0*(*opts).log_fwhh_prop_var)-c7[it][0])
                            {
                            //cout<<c1[it][0];
                            count = 1;
                            }  
                        } else {
                            if (c5[it][0]<ped+(15.0*(*opts).log_fwhh_prop_var)+(*opts).rdelta
                            && c5[it][0]>pst-(15.0*(*opts).log_fwhh_prop_var)-(*opts).rdelta)
                            {
                            count = 1;
                            }
                        }
                    }
                }
                if (count == 1)
                {
                   if(c2[it][0]+1 <0.0000001)          
			       {
				double prot=0;
			  
			  for(unsigned int locit=0;locit<c4[it].size();locit++)
			    prot+=c4[it][locit];
			  vector<double> weights(c4[it]);
			  for(unsigned int locit=0;locit<c4[it].size();locit++)
			    weights[locit]/=prot;
			  
			  //multiplet_site ms(prot, c1[it][0]);
			  multiplet_site ms(&c2[it], c1[it][0], &c3[it], prot, 
        		    x, c5[it][0], c6[it][0], c7[it][0], opts);
        		
			 // vector<unsigned int> c2int(c2[it].size(),0);
			  //vecftoi(c2[it], c2int);
			  ms.setup_param_extra(c3[it],weights);
			  templa.add_multiplet(ms); 
			} else {		
		            multiplet_site ms(&c2[it], c1[it][0], &c3[it], c4[it][0], 
        		    x, c5[it][0], c6[it][0], c7[it][0], opts);
        		    vector<unsigned int> c2int(c2[it].size(),0);
        	        vecftoi(c2[it], c2int);
                    ms.setup_param(c2int, c3[it]);
        	        templa.add_multiplet(ms);
			}
			
                    count = 0;
                }
    		  /*multiplet_site ms(&c2[it], c1[it][0], &c3[it], c4[it][0], 
    		  x, c5[it][0], c6[it][0], c7[it][0], opts);
    
    	      vector<unsigned int> c2int(c2[it].size(),0);
    	      vecftoi(c2[it], c2int);
    	      ms.setup_param(c2int, c3[it]);
    	      templa.add_multiplet(ms); */
          }
    	  if (!strcmp(lis,prevname) && c8[it][0]==1) 
          {
    	      (*Tems)[i] = templa; 
    	      count = 0;
          }
        }
    }
}

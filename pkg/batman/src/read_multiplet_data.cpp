#include "myheader.h"
#include "metab_template.h"

int read_multiplet_data(int lineno, char filename[], opt* opts,
						vector<string> *listname, vector<metab_template> *Tems, vector<double> *x,
						char chemshift[], int s, char inputdir[])
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
	// for ph
	vector<string> names2(lineno);
	matrix c11(lineno);
	matrix c12(lineno);
	
	
	int count = 0;
	double pst, ped;
	
	int nl = read_datf(&names,&c1, &c2, &c3, &c4, &c5, &c6, &c7, &c8, filename);
	
	if (nl < 0)
	{
		return nl;
	}
	
	
	if ((*opts).usechemshift == 1)
	{
		int nl2 = read_dat_chemshift(&names2, &c11,&c12, chemshift, s);
		if (nl2 < 0)
		{
			return nl2;
		}
		
		if (nl != nl2)
		{
			printf("Different number of multiplets, exiting ...\n");
			system("PAUSE");
			exit(1);
			//return -999;
		}
		for (unsigned int j = 0; j < nl2; j++)
		{
			if (!((c12[j][0]+50)<0.0000001))
			{
				//	if(names[j].compare(names2[j]) == 0 && c11[j][0] == c1[j][0])
				//{
				c1[j][0] = c12[j][0];
				//}else{
				//	printf("something wrong with multi chemshift, exiting ...\n");
				//	system("PAUSE");
				//	exit(1);
				//}
			}
		}
	}
	
	
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
					if((c2[it][0]+1 <0.0000001) && (c2[it][0]+1 > -0.0000001))
					{
						if (c3[it].size()!=c4[it].size())
						{
							cout<<"\nNo. of protons do not match no. of J constant for metabolite "<<names[it]<<", exiting ...\n";
							//system("PAUSE");
							exit(1);
						}
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
					} else if ((c2[it][0]+2 <0.0000001) && (c2[it][0]+2 > -0.0000001)) {
						
						//multiplet_site new_mult2(&pos_vec,&nprot_vec, &x);
						//cout << "c2 "<<c2[it][0] << endl;
						
						vector<double> raster(0.0);
						double vec_el;
						
						multiplet_site ms(&c2[it], c1[it][0], &c3[it], c4[it][0],
										  x, c5[it][0], c6[it][0], c7[it][0], opts);
						
						if (c3[it].size() == 2)
						{
							// raster
							//printf("find input raster\n");
							raster.clear();
							char fdirR[3000]={'\0'};
							strcpy(fdirR,inputdir);
							strcat(fdirR,name);
							strcat(fdirR,".txt");
							ifstream inA3_str(fdirR);
							//cout<<"route "<< fdirR<<endl;
							//cout<<"file is "<< inA3_str<< endl;
							
							//if (!inA3_str)
							//cout<<"empty, no file "<< inA3_str<< endl;
							
							
							while(inA3_str.good())
							{
								inA3_str>>vec_el;
								inA3_str.ignore(1);
								if(vec_el<max(c3[it][1],c3[it][0]) && vec_el>min(c3[it][1],c3[it][0]))
								{
									//cout<<"vec_el ppm"<<vec_el<<endl;
									inA3_str.peek();
									inA3_str>>vec_el;
									inA3_str.ignore(1);
									//cout<<"vec_el"<<vec_el<<endl;
									raster.push_back(vec_el);
								}
								//cout<<"vec_el"<<vec_el<<endl;
								// need to peek to modify state flag
								inA3_str.peek();
							}
							//cout << "c3 "<<c3[it].size()<<endl;
							//for (int ii = 0; ii <raster.size(); ii++)
							//cout<<"raster "<< raster[ii] << " ";
							
							ms.raster_setup(abs(c3[it][1]-c3[it][0]), &raster);
						} else {
							printf("Wrong ppm ranges for raster (-2) in J_constant, exiting ...\n");
							system("PAUSE");
							exit(1);
						}
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

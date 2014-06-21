// written by Dr. Jie Hao, Dr William Astle
#include "myheader.h"

    void read_txt_batmanoptions(opt *data, char filename[])
    {
        char buffer[1000] = {'\0'};
        char *p;
        int i = 0;
        int count = 0;
        FILE *fp=fopen(filename,"r");
        if( !fp)
        {
            printf("Can't open file %s, exiting ...\n", filename);
            system("PAUSE");
            exit(1);
        }
        while(1)
        {
            count++;
            char *tp0 = fgets(buffer, 1000, fp);
            if(feof(fp))
                break;
        }
        fclose(fp);

        vector<string> s(count-1);

        fp=fopen(filename,"r");

        while(1)
        {
            char buffer[1000] = {'\0'};
            char *tp0 = fgets(buffer, 1000, fp);
            p = buffer;
            if(feof(fp))
                break;
            if(*p =='%')
                continue;
            while(1)
            {
                if(*p == ':')
                {
                    p++;
                    if (*p == ' ')
                        p++;
                    while (*p != '\n')
                    {
                        s[i]+= *p;
                        p++;
                    }
                    if (*p == '\n')
                    {
                        s[i]+='\0';
                        i++;
                        break;
                    }
                }
                else if (*p == '\n')
                    break;
                p++;
            }
        }
        fclose(fp);
        int lind = 0;

        //read in ppm ranges
        int m = 0;
        char *pch1, *pch2, *pch3, *pch4;
	//acppm asigned according to line size.
	char ac[80], ac1[10], ac2[10], acppm[s[lind].size()];
        while (s[lind].size()>5)
        {
            //char ac[80] = {'\0'};
	    memset (acppm,'\0',s[lind].size());
            strcpy(acppm, s[lind].c_str());

            pch1 = strrchr(acppm,'(');
            pch2 = strrchr(acppm,',');
            pch3 = strrchr(acppm,')');

            int start2 = pch2-acppm+1;
            int lend2 = pch3-pch2-1;
            int start1 = pch1-acppm+1;
            int lend1 = pch2-pch1-1;
            if (start1<0 || lend1<0 || start2<0 || lend2<0)
            {
                printf("Problem with ppm ranges in %s, exiting ...\n", filename);
                system("PAUSE");
                exit(1);
            }

            //char ac1[10]={'\0'};
            //char ac2[10]={'\0'};

	    memset(ac1, '\0', 10);
	    memset(ac2, '\0', 10);
            strncpy(ac1, acppm+start1,lend1);
            strncpy(ac2, acppm+start2,lend2);

            double stt = atof(ac1);
            double endt = atof(ac2);

            if (stt < endt)
            {
                (*data).st.push_back(stt);
                (*data).ed.push_back(endt);
            }
            else
            {
                (*data).ed.push_back(stt);
                (*data).st.push_back(endt);
            }
            m++;
            int sst = start1-1;
            s[lind].erase(sst,pch3-pch1+1+1);
        }

        double temp = 0.0;
        int ct = 0;
        if ((*data).st.size()>1)
        {
            while (1)
            {
                for (unsigned int i = 0; i <(*data).st.size()-1; i++)
                {
                    if ((*data).st[i]< (*data).st[i+1])
                    {
                        temp = (*data).st[i];
                        (*data).st[i] = (*data).st[i+1];
                        (*data).st[i+1] = temp;
                        temp = (*data).ed[i];
                        (*data).ed[i] = (*data).ed[i+1];
                        (*data).ed[i+1] = temp;
                        ct--;
                    }
                }
                if (ct == 0 )
                    break;
                else
                    ct = 0;
            }
        }
        lind++;

        //char ac[80] = {'\0'};
	memset (ac,'\0',80);
        //read in spectra ranges
        char *psno;
        strcpy(ac, s[lind].c_str());
        psno = ac;

        int j1 = 0;
        int j2 = 0;
        while (s[lind].size() > j2)
        {
            if(*psno == ',')
            {
                j1 ++;
            }
            j2++;
            psno++;
        }
        vector< vector<int> > c(j1+1);
        j2 = 0;
        j1 = 0;
        int j3 = 0;
        char buf[100] = {'\0'};
        psno = ac;
        while (s[lind].size() > j2)
        {
            if(*psno == ',')
            {
                buf[j3] = '\0';
                c[j1].push_back(atoi(buf));
                buf[100] = '\0';
                j3 = 0;
                j1 ++;

            } else if (*psno == '-') {
                buf[j3] = '\0';
                c[j1].push_back(atoi(buf));
                buf[100] = '\0';
                j3 = 0;
            } else {
                buf[j3] = *psno;
                j3++;
            }
            j2++;
            //cout<<(*psno)<<endl;
            psno++;
        }
        buf[j3] = '\0';
        c[j1].push_back(atoi(buf));
        int tmpsn = 0;
        (*data).nospec.clear();
        for (int m = 0; m <c.size(); m ++)
        {
            if (c[m].size() == 2)
            {
                if (c[m][0] > c[m][1])
                {
                    tmpsn = c[m][0];
                    c[m][0] = c[m][1];
                    c[m][1] = tmpsn;
                }
                for (int n = c[m][0]; n <= c[m][1]; n++)
                    (*data).nospec.push_back(n);
            } else if (c[m].size() == 1)
            {
                (*data).nospec.push_back(c[m][0]);
            }
        }
        lind ++;

        memset (ac,'\0',80);
        strcpy(ac,s[lind].c_str());
        (*data).paraProc = atoi(ac);
        lind ++;
        // all the other values from options
        memset (ac,'\0',80);
        strcpy(ac,s[lind].c_str());
        (*data).lowerlimit = atof(ac);
        lind ++;

        memset (ac,'\0',80);
        strcpy(ac,s[lind].c_str());
        (*data).div = atof(ac);
        lind ++;

        memset (ac,'\0',80);
        strcpy(ac,s[lind].c_str());
        (*data).downsample = atoi(ac);
        lind ++;

        memset (ac,'\0',80);
        strcpy(ac,s[lind].c_str());
        (*data).saveHR = atoi(ac);
        lind ++;

        memset (ac,'\0',80);
        strcpy(ac,s[lind].c_str());
        (*data).seed = atoi(ac);
        lind ++;

        memset (ac,'\0',80);
        strcpy(ac,s[lind].c_str());
        (*data).stop_burnin = atoi(ac);
        lind ++;
		
		    memset (ac,'\0',80);
        strcpy(ac,s[lind].c_str());
        (*data).post_burnin = atoi(ac);
        lind ++;
        
        memset (ac,'\0',80);
        strcpy(ac,s[lind].c_str());
        (*data).tempF = atoi(ac);
        lind ++;

        memset (ac,'\0',80);
        strcpy(ac,s[lind].c_str());
        (*data).savefreq = atoi(ac);
        lind ++;

        memset (ac,'\0',80);
        strcpy(ac,s[lind].c_str());
        (*data).fix = atoi(ac);
        lind ++;

        memset (ac,'\0',80);
        strcpy(ac,s[lind].c_str());
        (*data).rerun = atoi(ac);
        lind ++;

        memset (ac,'\0',80);
        strcpy(ac,s[lind].c_str());
        (*data).tempe = atof(ac);
        lind ++;

        memset (ac,'\0',80);
        strcpy(ac,s[lind].c_str());
        (*data).spec_freq = atof(ac);
        lind ++;

        memset (ac,'\0',80);
        strcpy(ac,s[lind].c_str());
        (*data).a = atof(ac);
        lind ++;

        memset (ac,'\0',80);
        strcpy(ac,s[lind].c_str());
        (*data).b = atof(ac);
        lind ++;

        memset (ac,'\0',80);
        strcpy(ac,s[lind].c_str());
        (*data).log_fwhh_prior_mean = atof(ac);
        lind ++;

        memset (ac,'\0',80);
        strcpy(ac,s[lind].c_str());
        (*data).log_fwhh_prior_var = atof(ac);
        lind ++;

        memset (ac,'\0',80);
        strcpy(ac,s[lind].c_str());
        (*data).log_fwhh_prop_var = atof(ac);
        lind ++;

        memset (ac,'\0',80);
        strcpy(ac,s[lind].c_str());
        (*data).log_fwhh_re_prior_var = atof(ac);
        lind ++;

        memset (ac,'\0',80);
        strcpy(ac,s[lind].c_str());
        (*data).log_fwhh_re_prop_var = atof(ac);
        lind ++;

        memset (ac,'\0',80);
        strcpy(ac,s[lind].c_str());
        (*data).thresh = atof(ac);
        lind ++;

        memset (ac,'\0',80);
        strcpy(ac,s[lind].c_str());
        (*data).steep = atof(ac);
        lind ++;

        memset (ac,'\0',80);
        strcpy(ac,s[lind].c_str());
        (*data).rdelta = atof(ac);
        lind ++;

        memset (ac,'\0',80);
        strcpy(ac,s[lind].c_str());
        (*data).usechemshift = atoi(ac);
        lind ++;
    }



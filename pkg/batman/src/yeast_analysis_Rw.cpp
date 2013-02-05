    #include "myheader.h"
    #include "parclass.h"
    #include "multiplet_site.h"
    #include "metab_template.h"
    #include "chain_template.h"
    #include "spec_template.h"

    void read_multiplet_data(int lineno, char filename[], opt* opts,
            vector<string> *listname, vector<metab_template> *Tems, vector<double> *x,
            char chemshift[], int s, char inputdir[]);
    void write_results_all(int s, int sit, int rr, chain_template * chain, matrix * data,
            char opfile[], int ito, int burnin, vector<double> *x, int saveHR, matrix * dataH);
    /* This M-file provides an example to demonstrate MCMC-based Bayesian wavelet estimation.
    Load the data and set values for parameters related to the data.
    choose periodic boundary conditions - i.e. we do the wavelet
    decomposition on a circle
     */

    extern "C" {
    #include <R.h>
    #include <Rdefines.h>
    #include <Rinternals.h>
    #include <R_ext/RS.h>
    #include <R_ext/Print.h>
    SEXP batman(SEXP filenames, SEXP bn_int, SEXP ito_int, SEXP rr_int, SEXP sno_int, SEXP pBar);

    SEXP batman(SEXP filenames, SEXP bn_int, SEXP ito_int, SEXP rr_int, SEXP sno_int, SEXP pBar)
    {
        char *Pfilenames[7];
        int bn,itR,rr, s;
        bn = INTEGER_VALUE(bn_int);
        itR = INTEGER_VALUE(ito_int);
        rr = INTEGER_VALUE(rr_int);
        s = INTEGER_VALUE(sno_int);

        PROTECT(filenames = AS_CHARACTER(filenames));
        // allocate memory:
        Pfilenames[0] = R_alloc(strlen(CHAR(STRING_ELT(filenames, 0))), sizeof(char));
        Pfilenames[1] = R_alloc(strlen(CHAR(STRING_ELT(filenames, 1))),  sizeof(char));
        Pfilenames[2] = R_alloc(strlen(CHAR(STRING_ELT(filenames, 2))),  sizeof(char));
        Pfilenames[3] = R_alloc(strlen(CHAR(STRING_ELT(filenames, 3))),  sizeof(char));
        Pfilenames[4] = R_alloc(strlen(CHAR(STRING_ELT(filenames, 4))),  sizeof(char));
        Pfilenames[5] = R_alloc(strlen(CHAR(STRING_ELT(filenames, 5))),  sizeof(char));
        Pfilenames[6] = R_alloc(strlen(CHAR(STRING_ELT(filenames, 6))),  sizeof(char));



        // ... and copy mychar to Pmychar:
        strcpy(Pfilenames[0], CHAR(STRING_ELT(filenames, 0)));
        strcpy(Pfilenames[1], CHAR(STRING_ELT(filenames, 1)));
        strcpy(Pfilenames[2], CHAR(STRING_ELT(filenames, 2)));
        strcpy(Pfilenames[3], CHAR(STRING_ELT(filenames, 3)));
        strcpy(Pfilenames[4], CHAR(STRING_ELT(filenames, 4)));
        strcpy(Pfilenames[5], CHAR(STRING_ELT(filenames, 5)));
        strcpy(Pfilenames[6], CHAR(STRING_ELT(filenames, 6)));

        UNPROTECT(1);


        opt options;
        read_txt_batmanoptions(&options, Pfilenames[0]);

        rngType rng(options.seed);

        int iteration = options.ito;
        iteration = itR;

        matrix Data;
        read_txtf(&Data,Pfilenames[1]);

        int dsize = 2;
        if (options.fix == 1)
        {  dsize = options.nospec.size()+ 1;}
        else
        { dsize = 2;}

        matrix data(dsize);

        int dppm = 1; //select 1, remove -1
        datamod(&Data,&options, options.downsample, dppm, &data, s);

        vector<double> x(data[0]);

        vector<string> listn;
        read_txt_metalist(&listn,Pfilenames[2]);

        vector<int> zi(1, 0);
        vector<double> zd(1,0.0);
        metab_template ta("",x.size());
        multiplet_site multi_site(&zd, 0.0, &zd, zi[0], &x, -50, -50, -50, &options);
        ta.add_multiplet(multi_site);

        vector<metab_template> TemsU(listn.size(),ta);

        time_t seconds5,seconds6;
        seconds5 = time (NULL);

		int lineno = bn + 1;//input from R
        read_multiplet_data(lineno, Pfilenames[3], &options, &listn, &TemsU, &x, Pfilenames[5], s, Pfilenames[6]);
        vector<int> metaUsed;
        vector<metab_template> Tems;
        // check for meta in ppm range
        for (unsigned int cv = 0; cv <TemsU.size(); cv++)
        {
            if (TemsU[cv].multiplet_sites.size() != 0)
            {
                Tems.push_back(TemsU[cv]);
            }
        }
        if (Tems.size()==0)
        {
        Rprintf("\nNo metabolites with resonances in the region for analysis, exiting ...\n");
			/*SEXP myint;
			int *p_myint;
			int len = 2;
			int tp = 0;
			// Allocating storage space:
			PROTECT(myint = NEW_INTEGER(len));
			p_myint = INTEGER_POINTER(myint);
			p_myint[0] = tp;
			p_myint[1] = 1;
			UNPROTECT(1);
			return myint;*/

            system("PAUSE");
            exit(1);
        }
        char fdirML[3000]={'\0'};
        strcpy(fdirML,Pfilenames[4]);
        strcat(fdirML,"metabolitesListUsed.txt");

        FILE *outM;
        outM = fopen(fdirML,"w");
        // wirte to file the metabolites in range for anaylsis
        for (unsigned int cv = 0; cv <Tems.size(); cv++)
        {
            fprintf(outM, "%s",Tems[cv].name.c_str());
            if (cv <=(Tems.size()-1))
			{
                fprintf(outM, "\n");
			}
        }
        fclose(outM);

        chain_template chain(&Tems, &data,&rng, &options);
        seconds6 = time (NULL);

        matrix reconL0(chain.spectra[0].L.size());
        for (unsigned int t2 = 0; t2<reconL0.size(); t2++)
            idwt((chain).spectra[0].L[t2], (chain).spectra[0].pars.levsize, (chain).spectra[0].pars.h_vec, reconL0[t2]);

        char sbuf [33] = {'\0'};
        string str = boost::lexical_cast <string>(s+1);
        strcpy(sbuf, str.c_str());

        char fdirL[3000]={'\0'};
        strcpy(fdirL,Pfilenames[4]);
        strcat(fdirL,"L_");
        strcat(fdirL,sbuf);
        strcat(fdirL,".txt");
        write_txtf_M(&reconL0,fdirL);

        time_t seconds, seconds2;
        seconds = time (NULL);
        char fdirrr[3000]={'\0'};
        strcpy(fdirrr,Pfilenames[4]);
        strcat(fdirrr,"delta_draw_mean_");
        strcat(fdirrr,sbuf);
        strcat(fdirrr,".txt");
        FILE *out;

        Rprintf("\nSize of each spectrum is %i.\n",x.size());
        Rprintf("Size of metabolite list is %i.\n",TemsU.size());
        Rprintf("of which %i have resonances in/near the specified region and will be fit.\n",Tems.size());
        Rprintf("Constructing chain data structure...\n");
        Rprintf("time used is %i seconds.\n", seconds6-seconds5);
        Rprintf("Running MCMC...\n");
        // progress bar calling R from c++
        int *rPercentComplete;
        SEXP utilsPackage, percentComplete;
        PROTECT(utilsPackage = eval(lang2(install("getNamespace"), ScalarString(mkChar("utils"))), R_GlobalEnv));
        PROTECT(percentComplete = NEW_INTEGER(1));
        rPercentComplete = INTEGER_POINTER(percentComplete);

        *rPercentComplete = 0;
        eval(lang4(install("setTxtProgressBar"), pBar, percentComplete, R_NilValue), utilsPackage);

        if ( rr == 0)
        {
            seconds5 = time (NULL);
            for(int ito = 0; ito<iteration; ito++)
            {
                if(ito<options.stop_burnin-1)
                {
                    (chain).set_repen(1.0);
                    (chain).set_temp(1+options.tempe*(1-my_normcdf(ito-double(options.stop_burnin)
                            -1000.0+1.0, 0.0, double(pow(options.stop_burnin,2.0))*0.00005)));
                }

                *rPercentComplete = ito+1;
                eval(lang4(install("setTxtProgressBar"), pBar, percentComplete, R_NilValue), utilsPackage);

                chain.sample_metropolis((chain.it % 10 ==0 ),rr, options.log_fwhh_re_prior_var);

                if(ito==(options.stop_burnin/2)-1)
                    chain.set_adapt(true);

                if(ito==options.stop_burnin-1)
                {
                    chain.set_temp(1.0);
                    chain.stop_burn();
                    time_t seconds7;
                    seconds7 = time (NULL);
                    Rprintf("\ntime used for burnin is %i seconds.\n", seconds7-seconds);
                }
                if (options.savefreq!=0)
                {
                    if(chain.it%options.savefreq == 0 || chain.it== iteration)
					{
                        chain.store_state();
					}
                }
            }
            // save ppm shift mean
            out = fopen(fdirrr,"w");
            for (unsigned int d1 = 0; d1 < chain.spectra.size();d1++)
            {
                for (unsigned int d2 = 0; d2 < chain.spectra[d1].FTems.size();d2++)
                {
                    for (unsigned int d3 = 0; d3 < chain.spectra[d1].FTems[d2].multiplet_sites.size(); d3++)
                    {
                        if (chain.spectra[d1].FTems[d2].multiplet_sites[d3].delta_draw_mean != 0)
						{
                            chain.spectra[d1].FTems[d2].multiplet_sites[d3].delta_draw_mean =
                                    chain.spectra[d1].FTems[d2].multiplet_sites[d3].delta_draw_mean/chain.spectra[d1].smeancount;
						}
                        fprintf(out, "%f", chain.spectra[d1].FTems[d2].multiplet_sites[d3].delta_draw_mean);

                        if (d1 != chain.spectra.size()-1 ||d2 != chain.spectra[d1].FTems.size()-1||d3 != chain.spectra[d1].FTems[d2].multiplet_sites.size()-1 )
						{
							fprintf(out, "\n");
						}
                    }
                }
            }
            fclose(out);
        } else {
            vector<double> rrdata;
            string strs;
            char buf[100] = {'\0'};
            ifstream fin(fdirrr);
            while( getline(fin, strs) )
            {
                strcpy(buf, strs.c_str());
                rrdata.push_back(atof(buf));
            }
            int d4 = 0;
            // reading back ppm shift mean
            for (unsigned int d1 = 0; d1 < chain.spectra.size();d1++)
            {
                for (unsigned int d2 = 0; d2 < chain.spectra[d1].FTems.size();d2++)
                {
                    for (unsigned int d3 = 0; d3 < chain.spectra[d1].FTems[d2].multiplet_sites.size(); d3++)
                    {
                        chain.spectra[d1].FTems[d2].multiplet_sites[d3].delta_draw = rrdata[d4];
                        d4 = d4 +1;
                    }
                }
            }
            chain.set_repen(1.0);
            chain.set_temp(1.0);
            for (int ito = 0; ito <options.rerun; ito++)
            {
                chain.sample_metropolis((chain.it % 10 ==0),rr, options.log_fwhh_re_prior_var);
                *rPercentComplete = ito+1;
                eval(lang4(install("setTxtProgressBar"), pBar, percentComplete, R_NilValue), utilsPackage);
                if (options.savefreq!=0)
                {
                    if(chain.it%options.savefreq == 0)
					{
                        chain.store_state();
					}
                }
            }
        }
        UNPROTECT(2);

        seconds2 = time (NULL);
        Rprintf("\ntime used is %i seconds.\n", seconds2-seconds);
        Rprintf("saving posteriors...\n ");
        // save results
        matrix dataH(dsize);
        datamod(&Data,&options, 1, dppm, &dataH, s);
        if (options.saveHR == 1)
        {
            char fdirdH[3000]={'\0'};
            strcpy(fdirdH,Pfilenames[4]);
            strcat(fdirdH,"NMRdata_mod_");
            strcat(fdirdH,sbuf);
            strcat(fdirdH,".txt");
            write_txtf_M(&dataH, fdirdH);
        }
        int sit = 0;
        if (options.fix == 1)
        {
            for (sit =  0; sit <= s; sit ++)
            {
                write_results_all(sit, sit, rr, &chain, &data, Pfilenames[4],iteration,options.stop_burnin, &x,
                        options.saveHR, &dataH);
            }
        }
        else
        {
            sit = 0;
            write_results_all(s, sit, rr, &chain, &data, Pfilenames[4],iteration,options.stop_burnin, &x,
                    options.saveHR, &dataH);
        }
//for (int nii = 0; nii <chain.spectra[0].FTems.size(); nii++)
//cout <<"chain "<<chain.spectra[0].FTems[nii].name<< "_"<<chain.spectra[0].FTems[nii].multiplet_sites[0].pos_ppm<<endl;

        SEXP myint;
        int *p_myint;
        int len = 2;
        int tp = seconds2-seconds;
        // Allocating storage space:
        PROTECT(myint = NEW_INTEGER(len));
        p_myint = INTEGER_POINTER(myint);
        p_myint[0] = tp;
		    p_myint[1] = 0;
        UNPROTECT(1);
        return myint;
    }
    }

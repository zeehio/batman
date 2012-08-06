#include "myheader.h"
#include "chain_template.h"
void write_results_all(int s, int sit, int rr, chain_template * chain, matrix * data, 
char opfile[], int ito, int burnin, vector<double> *x, int saveHR, matrix * dataH)
{
    FILE *out;
    char rrbuf [33] = {'\0'};
    string str = boost::lexical_cast <string>(rr);
    strcpy(rrbuf, str.c_str());
    
    char sbuf [33] = {'\0'};
    str = boost::lexical_cast <string>(s+1);
    strcpy(sbuf, str.c_str());
      
    vector<double> xH((*dataH)[0]); 
    matrix LH((*chain).spectra[sit].FTems.size(),xH);
    
    char fdirsf[3000]={'\0'};
    strcpy(fdirsf, opfile);
    strcat(fdirsf,"specFit_");
    strcat(fdirsf,sbuf);
    strcat(fdirsf,"_rr_");
    strcat(fdirsf,rrbuf);
    strcat(fdirsf,".txt");

    out = fopen(fdirsf,"w");
    fprintf(out, "%s \t %s \t %s \t %s \t %s\n","ppm","Original spectrum", "Metabolites fit", "Wavelet fit", "Overall fit");

    vector<double> V((*data)[0].size(),0.0);
    matrix reconL((*chain).spectra[0].L.size());

    vector<double> yu((*data)[0]);
    
    VSdiv(&((*chain).spectra[sit].beta_mean),(*chain).spectra[sit].smeancount, &((*chain).spectra[sit].beta_mean));
    VSdiv(&((*chain).beta_mean),(*chain).spectra[sit].smeancount, &((*chain).beta_mean));
    VSdiv(&((*chain).spectra[sit].theta_mean), (*chain).spectra[sit].smeancount, &((*chain).spectra[sit].theta_mean));
    idwt((*chain).spectra[sit].theta_mean, (*chain).spectra[sit].pars.levsize, (*chain).spectra[sit].pars.h_vec, yu); 

    if (rr == 0)
    {
        for (unsigned int d2 = 0; d2 < (*chain).spectra[sit].FTems.size();d2++)
        {
            for (unsigned int d3 = 0; d3 < (*chain).spectra[sit].FTems[d2].multiplet_sites.size(); d3++)
                (*chain).spectra[sit].FTems[d2].multiplet_sites[d3].delta_draw =(*chain).spectra[sit].FTems[d2].multiplet_sites[d3].delta_draw_mean;   
        }
    }
    vector<metab_template> TemsH((*chain).spectra[sit].FTems);  
    for (unsigned int cv = 0; cv <TemsH.size(); cv++)
    {
        TemsH[cv].build_curve(*x,(*chain).spectra[sit].log_fwhh_draw, (*chain).spectra[sit].pars.freq);  
        reconL[cv]= TemsH[cv].curve;                             
    } 
    MVprod( &reconL, &((*chain).spectra[sit].beta_mean), &V);

    
    for (unsigned int t = 0; t<V.size(); t++)
    {	
    	fprintf(out, "%f \t %f \t %f \t %f \t %f",(*data)[0][t], (*data)[sit+1][t], V[t], yu[t], (V[t]+yu[t]));
        if (t<=V.size()-1)
        	fprintf(out, "\n");
    }
    fclose(out);
    
    char fdirbet[3000]={'\0'};
    strcpy(fdirbet, opfile);
    strcat(fdirbet,"beta_");
    strcat(fdirbet,sbuf);
    strcat(fdirbet,"_rr_");
    strcat(fdirbet,rrbuf);
    strcat(fdirbet,".txt");
    write_txtf(&((*chain).spectra[sit].beta_mean), fdirbet);
    
    char fdirbetsam[3000]={'\0'};
    strcpy(fdirbetsam, opfile);
    strcat(fdirbetsam,"beta_sam_");
    strcat(fdirbetsam,sbuf);
    strcat(fdirbetsam,"_rr_");
    strcat(fdirbetsam,rrbuf);
    strcat(fdirbetsam,".txt");
    write_txtf_M(&((*chain).spectra[sit].beta_sam), fdirbetsam);
    
    char savef[3000]={'\0'};
    strcpy(savef,opfile);
    strcat(savef,"metaTemp_");
    strcat(savef,sbuf);
    strcat(savef,"_rr_");
    strcat(savef,rrbuf);
    strcat(savef,".txt");
    write_txtf_M(&reconL, savef);
    
    char fdirthe[3000]={'\0'};
    strcpy(fdirthe,opfile);
    strcat(fdirthe,"theta_sam_");
    strcat(fdirthe,sbuf);
    strcat(fdirthe,"_rr_");
    strcat(fdirthe,rrbuf);
    strcat(fdirthe,".txt");
    matrix reconT((*chain).spectra[sit].theta_sam.size());
    for (unsigned int n = 0; n<(*chain).spectra[sit].theta_sam.size(); n++)
    idwt((*chain).spectra[sit].theta_draw, (*chain).spectra[sit].pars.levsize, (*chain).spectra[sit].pars.h_vec, reconT[n]); 
    write_txtf_M(&reconT, fdirthe);  
    
    char fdirsfit[3000]={'\0'};
    strcpy(fdirsfit,opfile);
    strcat(fdirsfit,"metaFit_sam_");
    strcat(fdirsfit,sbuf);
    strcat(fdirsfit,"_rr_");
    strcat(fdirsfit,rrbuf);
    strcat(fdirsfit,".txt");
    write_txtf_M(&((*chain).spectra[sit].sfit), fdirsfit);  
    
    char fdirmetaS[3000]={'\0'};
    strcpy(fdirmetaS,opfile);
    strcat(fdirmetaS,"metaIndFit_sam_");
    strcat(fdirmetaS,sbuf);
    strcat(fdirmetaS,"_rr_");
    strcat(fdirmetaS,rrbuf);
    strcat(fdirmetaS,".txt");
    write_txtf_M(&((*chain).spectra[sit].meta_sam), fdirmetaS);  
  
    if (rr == 0)
    { 
        char fdirlam[3000]={'\0'};
    
        strcpy(fdirlam,opfile);
        strcat(fdirlam,"lambda_sam_");
        strcat(fdirlam,sbuf);
        strcat(fdirlam,"_rr_");
        strcat(fdirlam,rrbuf);
        strcat(fdirlam,".txt");
        write_txtf(&((*chain).spectra[sit].lambda_sam), fdirlam);

        char fdirdel[3000]={'\0'};
        strcpy(fdirdel,opfile);
        strcat(fdirdel,"delta_sam_");
        strcat(fdirdel,sbuf);
        //strcat(fdirdel,"_rr_");
        //strcat(fdirdel,rrbuf);
        strcat(fdirdel,".txt");
            
        out = fopen(fdirdel,"w");
     
        for (unsigned int t = 0; t<(*chain).spectra[sit].FTems.size(); t++)
        {
            for (unsigned int i = 0; i<(*chain).spectra[sit].FTems[t].multiplet_sites.size(); i++)
            {
                char metaname[3000]={'\0'};
                strcpy(metaname, (*chain).spectra[sit].FTems[t].name.c_str());
                fprintf(out, "%s_",metaname);
                fprintf(out, "%.3f",(*chain).spectra[sit].FTems[t].multiplet_sites[i].pos_ppm);          
        
                if ((t == ((*chain).spectra[sit].FTems.size()-1)) && ((i == (*chain).spectra[sit].FTems[t].multiplet_sites.size()-1)))	 
                    fprintf(out, "\n");
        	    else
                    fprintf(out, "\t");
            }
        }
        
        for (unsigned int j = 0; j<(*chain).spectra[sit].FTems[0].multiplet_sites[0].delta_sam.size(); j++)
        {
            for (unsigned int t = 0; t<(*chain).spectra[sit].FTems.size(); t++)
            {
                for (unsigned int i = 0; i<(*chain).spectra[sit].FTems[t].multiplet_sites.size(); i++)
                {
                    fprintf(out, "%.17g", (*chain).spectra[sit].FTems[t].multiplet_sites[i].delta_sam[j]);
                    if ((t == ((*chain).spectra[sit].FTems.size()-1)) && (i == ((*chain).spectra[sit].FTems[t].multiplet_sites.size()-1)))
                	    fprintf(out, "\n");
                    else
                        fprintf(out, "\t");
                }
            }
        }
        fclose(out);
    }
    // Higher resolution spectrum and metabolites fit
    if (saveHR == 1)
    {
        vector<metab_template> TemsHR(TemsH);  
        for (unsigned int cv = 0; cv <TemsHR.size(); cv++)
        { 
        	TemsHR[cv].curve.resize(xH.size(),0.0);
            TemsHR[cv].build_curve(xH,(*chain).spectra[sit].log_fwhh_draw, (*chain).spectra[sit].pars.freq); 
            LH[cv]= TemsHR[cv].curve_prop;                      
        }  

        
        char savefHR[3000]={'\0'};
        strcpy(savefHR,opfile);
        strcat(savefHR,"metaTempHR_");
        strcat(savefHR,sbuf);
        strcat(savefHR,"_rr_");
        strcat(savefHR,rrbuf);
        strcat(savefHR,".txt");
        write_txtf_M(&LH, savefHR);
    
        char fdirsfHR[3000]={'\0'};
        strcpy(fdirsfHR, opfile);
        strcat(fdirsfHR,"specFitHR_");
        strcat(fdirsfHR,sbuf);
        strcat(fdirsfHR,"_rr_");
        strcat(fdirsfHR,rrbuf);
        strcat(fdirsfHR,".txt");
        vector<double> VHR((*dataH)[0].size());
        MVprod( &LH, &((*chain).spectra[sit].beta_draw), &VHR);
        write_txtf(&VHR, fdirsfHR);
    }         
}

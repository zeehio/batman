// written by Dr. Jie Hao, Dr William Astle
#include "myheader.h"
double truncated_t_twosided_pdf(double x, double mu, double scale, double df, double bot_lim, double top_lim)
{
    double std_bot, std_top, std_x, pdf;
    std_bot = (bot_lim-mu)/scale;
    std_top = (top_lim-mu)/scale;
    std_x = (x-mu)/scale;
    pdf=my_tpdf(std_x,df)/(my_tcdf(std_top,df)-my_tcdf(std_bot,df));
    
    return pdf;
}


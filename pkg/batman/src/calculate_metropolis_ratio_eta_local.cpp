#include "spec_template.h"

double spec_template::calculate_metropolis_ratio_eta_local(matrix * Lprop)
{
    // written by Dr. Jie Hao, Dr William Astle
	double logratio = 0;

    vector<double> vec(dataWy.size(),0);
    vector<double> vecprop(vec);
    vector<double> mvpd(vec);
    vector<double> mvpd2(vec);
    MVprod(&L, &beta_draw,&mvpd);

    MVprod(Lprop, &beta_draw,&mvpd2);

    for (unsigned int j = 0; j<dataWy.size(); j++)
    {
        vec[j] = dataWy[j] - theta_draw[j] - mvpd[j]; 
        vecprop[j] = dataWy[j] - theta_draw[j] -mvpd2[j];
    }
    
    logratio=-0.5*(VVprod(&vecprop, &vecprop)-VVprod(&vec, &vec))*lambda_draw;
    logratio=logratio/Temp;
    return logratio;
}

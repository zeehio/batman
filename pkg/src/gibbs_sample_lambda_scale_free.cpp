#include "spec_template.h"

double spec_template::gibbs_sample_lambda_scale_free()
{

  vector<double> vec(theta_draw.size());
  vector<double> vec2(vec);
  vector<double> mvpd(vec);
  vector<double> vec3(tau_draw.size());
  MVprod(&L, &beta_draw, &mvpd);

  for (int i = 0; i < pars.p; i++)
  {
      vec[i]=dataWy[i]-mvpd[i]-theta_draw[i];
      vec2[i]=sqrt(psi_draw[i])*theta_draw[i];
  }  
  for (int i = 0; i < pars.n; i++)
  {
      vec3[i]=sqrt(pars.steep)*(pars.thresh-tau_draw[i]);
  }
  
  double scalepar=(2.0*Temp/((VVprod(&vec,&vec)+VVprod(&vec2,&vec2)+VVprod(&vec3,&vec3))+ pars.b));
  
  //deliberately dont temper this
  double lambda_sample= my_gamrnd((pars.a+pars.p+pars.n/2-1)+1, scalepar, rng);
  return lambda_sample;
}

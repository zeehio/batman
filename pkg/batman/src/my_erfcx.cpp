#include "myheader.h"
using namespace boost;
using namespace math;

double my_erfcx(double z)
{   
     double result = 0.0;
     if(z < 1.5)
     {
         result = exp(pow(z, 2.0))*erfc(z);
     }
     else if(z < 2.5)
     {
         static const double Y = 0.50672817230224609375f;
         static const double P[] = {    
            -0.024350047620769840217L,
            0.0343522687935671451309L,
            0.0505420824305544949541L,
            0.0257479325917757388209L,
            0.00669349844190354356118L,
            0.00090807914416099524444L,
            0.515917266698050027934e-4L,
         };
         static const double Q[] = {    
            1L,
            1.71657861671930336344L,
            1.26409634824280366218L,
            0.512371437838969015941L,
            0.120902623051120950935L,
            0.0158027197831887485261L,
            0.000897871370778031611439L,
         };
         result = Y + tools::evaluate_polynomial(P, z - 1.5f) / tools::evaluate_polynomial(Q, z - 1.5f);
         result = result/z;
     }
     else if(z < 4.5)
     {
         static const double Y  = 0.5405750274658203125f;
         static const double P[] = {    
            0.0029527671653097284033L,
            0.0141853245895495604051L,
            0.0104959584626432293901L,
            0.00343963795976100077626L,
            0.00059065441194877637899L,
            0.523435380636174008685e-4L,
            0.189896043050331257262e-5L,
         };
         static const double Q[] = {    
            1L,
            1.19352160185285642574L,
            0.603256964363454392857L,
            0.165411142458540585835L,
            0.0259729870946203166468L,
            0.00221657568292893699158L,
            0.804149464190309799804e-4L,
         };
         result = Y + tools::evaluate_polynomial(P, z - 3.5f) / tools::evaluate_polynomial(Q, z - 3.5f);
         result = result/z;
      }
      else
      {
         static const double Y = 0.55825519561767578125f;
         static const double P[] = {    
            0.00593438793008050214106L,
            0.0280666231009089713937L,
            -0.141597835204583050043L,
            -0.978088201154300548842L,
            -5.47351527796012049443L,
            -13.8677304660245326627L,
            -27.1274948720539821722L,
            -29.2545152747009461519L,
            -16.8865774499799676937L,
         };
         static const double Q[] = {    
            1L,
            4.72948911186645394541L,
            23.6750543147695749212L,
            60.0021517335693186785L,
            131.766251645149522868L,
            178.167924971283482513L,
            182.499390505915222699L,
            104.365251479578577989L,
            30.8365511891224291717L,
         };
         result = Y + tools::evaluate_polynomial(P, 1 / z) / tools::evaluate_polynomial(Q, 1 / z);
         result = result/z;
      }
    return result;
} 

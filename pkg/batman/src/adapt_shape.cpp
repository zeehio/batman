#include "spec_template.h"

void spec_template::adapt_shape(double timepast)
{
    //what acceptance rate should we target ?
    if((double)shape_accepcount/timepast<0.4)
        shape_uplogsd=shape_uplogsd-deltaadapt(noshadaptations);
    else
        shape_uplogsd=shape_uplogsd+deltaadapt(noshadaptations);
    
    shape_accepcount=0;
    noshadaptations=noshadaptations+1;

}

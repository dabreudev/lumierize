#include "modulos.h"

float RNDGauss(float s)

{
  float r;

  r=sqrt(2.0)*s;
  r *= sqrt(-log(1-((float)rand()/(RAND_MAX+1.0))));
  r *= cos(2*M_PI*((float)rand()/(RAND_MAX+1.0)));
  return(r);
}

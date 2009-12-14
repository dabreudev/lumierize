#include "modulos.h"


int ResolvEcu_deg2(struct Pol_deg2 pol, float *x1, float *x2) {

  float rad;

  rad=pol.b*pol.b-4*pol.a*pol.c;
  if(rad==0) {
    *x1=-pol.b/2/pol.a;
    *x2=-pol.b/2/pol.a;
    return(1);
  }
  else if(rad<0) {
    *x1=0;
    *x2=0;
    return(0);
  }
  else {
    *x1=(-pol.b+sqrt(rad))/2/pol.a;
    *x2=(-pol.b-sqrt(rad))/2/pol.a;
    return(2);
  }

}


#include "modulos.h"


void pgLimits(int n, float *x, float *wmin, float *wmax)

{
  float min,max,r;

  MinMax(n,x,&min,&max);
/*   //  printf("AAAAAAAAAAA min %f max  %f \n",min,max); */
  r=max-min;
  if(r == 0) {
    *wmin=min*0.8;
    *wmax=max*1.1;
  }
  else {
    *wmin=min-r*0.1;
    *wmax=max+r*0.1;
  }
/*   //  printf("BBBBBBB min %f max  %f r %f\n",*wmin,*wmax,r); */
}




#include "cpgdoble.h"
#include "minmax.h"


void pgLimits_d(int n, double *x, float *wmin, float *wmax)

{
  double min,max;
  float r;

  MinMax_d(n,x,&min,&max);
/*   //  printf("AAAAAAAAAAA min %f max  %f \n",min,max); */
  r=max-min;
  if(r == 0) return;
  *wmin=min-r*0.1;
  *wmax=max+r*0.1;
/*   //  printf("BBBBBBB min %f max  %f r %f\n",*wmin,*wmax,r); */
}




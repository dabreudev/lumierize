#include "modulos.h"

/* //De Javier Gorgas */

float StWeightMedia(int n, float *x, float *w,float *sigma)
{
  int i;
  float m=0,s=0;
  float sumw=0;


  for(i=0;i<n;i++) {
    m+=x[i]*w[i];
    sumw+=w[i];
  }


  m=m/sumw;
  for(i=0;i<n;i++) {
    s+=(x[i]-m)*(x[i]-m)*w[i];
  }
  *sigma=sqrt(n*s/((n-1.)*sumw));
  return(m);
}

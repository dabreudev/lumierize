#include "modulos.h" 


float StErrWeightMedia(int n, float *x, float *err,float *sigma)
{
  int i;
  float m=0,s=0;
  float sumerr=0;
  
  
  for(i=0;i<n;i++) {
    m+=x[i]/err[i]/err[i];
    sumerr+=1./err[i]/err[i];
    /*     printf(" TW m %g x_i %g w_i %g \n",m,x[i],w[i]); */
  }
  
  
  m=m/sumerr;
  for(i=0;i<n;i++) {
    s+=(x[i]-m)*(x[i]-m)/err[i]/err[i];
  }
  *sigma=sqrt(n*s/((n-1.)*sumerr));
  return(m);
}

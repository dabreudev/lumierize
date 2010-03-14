#include "stmedia.h"

/* //De Javier Gorgas */

double StWeightMedia_d(int n, double *x, double *w,double *sigma)
{
  int i;
  double m=0,s=0;
  double sumw=0;


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

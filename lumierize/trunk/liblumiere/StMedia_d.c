#include "stmedia.h"

double StMedia_d(int n, double *x, double *sigma)

{
  int i;
  float m,s=0;

  m=StSuma1_d(n,x,1)/n;
  for(i=0; i<n; i++) s += (x[i]-m)*(x[i]-m);
  *sigma = sqrt(s/(n-1));
  return(m);
}

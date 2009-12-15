#include "modulos.h"


float StCovar(int n, float *a, float *b, float *sigmaa, float *sigmab)

{

  float mediaa,mediab;
  float covar;
  float sum2=0;
  int i;

  mediaa=StMedia(n,a,sigmaa);
  mediab=StMedia(n,b,sigmab);

  for(i=0;i<n;i++) {
    sum2+=(a[i]-mediaa)*(b[i]-mediab);
/*      sum+=a[i]*b[i]; */
  }
/*    sum=sum/(n-1.); */
/*    covar=sum-mediaa*mediab; */
 
  covar=sum2/(n-1.);
  return(covar);
}

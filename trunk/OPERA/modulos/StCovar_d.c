#include "modulos.h"


double StCovar_d(int n, double  *a, double *b, double *sigmaa, double *sigmab)

{

  double mediaa,mediab;
  double covar;
  double sum2=0;
  int i;

  mediaa=StMedia_d(n,a,sigmaa);
  mediab=StMedia_d(n,b,sigmab);

  for(i=0;i<n;i++) {
    sum2+=(a[i]-mediaa)*(b[i]-mediab);
/*      sum+=a[i]*b[i]; */
  }
/*    sum=sum/(n-1.); */
/*    covar=sum-mediaa*mediab; */
 
  covar=sum2/(n-1.);
  return(covar);
}

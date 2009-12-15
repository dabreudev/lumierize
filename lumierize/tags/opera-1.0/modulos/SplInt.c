#include "modulos.h"

float SplInt(int n, float *xa, float *ya, float *y2a, float x)
{
  int k,klo,khi;
  float h,a,b,y;

  klo=0;
  khi=n-1;
  while(khi-klo > 1) {
    k=(khi+klo)/2;
    if(xa[k] > x) {
      khi=k;
      }
    else {
      klo=k;
      }
    }
  h=xa[khi]-xa[klo];
  if(h == 0) {
    printf("SplInt: ERROR. Los valores de xa deben ser distintos\n");
    exit(1);
    }
  a=(xa[khi]-x)/h;
  b=(x-xa[klo])/h;
  y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*h*h/6.;
  return(y);
}
    

#include "modulos.h"

void MCLine(int n,float *x,float *y,float *a,float *erra) {

  int i;
  float sumxy=0,sumx2=0,sumyax=0;
  
  for(i=0;i<n;i++) {
    sumxy=sumxy+x[i]*y[i];
    sumx2=sumx2+x[i]*x[i];
  }

  *a=sumxy/sumx2;
  
  for(i=0;i<n;i++)     sumyax=sumyax+(y[i]-*a*x[i])*(y[i]-*a*x[i]);
  
  *erra=sqrt(sumyax/((n-1)*sumx2));
}

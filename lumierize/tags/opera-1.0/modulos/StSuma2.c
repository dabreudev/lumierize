#include "modulos.h"


float StSuma2(int n, float *a, int i, float *b, int j)

{

  int k;
  float suma2=0;

  if(i == 0) return(StSuma1(n,b,j));
  if(j == 0) return(StSuma1(n,a,i));

  for(k=0; k<n; k++) suma2 += pow(a[k],(float)i)*pow(b[k],(float)j);
  return(suma2);
}

#include "modulos.h"


double StSuma2_d(int n, double *a, int i, double *b, int j)

{

  int k;
  double suma2=0;

  if(i == 0) return(StSuma1_d(n,b,j));
  if(j == 0) return(StSuma1_d(n,a,i));

  for(k=0; k<n; k++) suma2 += pow(a[k],(float)i)*pow(b[k],(float)j);
  return(suma2);
}

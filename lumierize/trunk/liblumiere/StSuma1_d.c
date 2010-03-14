#include <math.h>
#include "stmedia.h"


double StSuma1_d(int n, double *a, int i)

{
  int k;
  double suma1=0;

  if(i == 0) return ((double)n);

  for(k=0; k<n; k++) suma1 += pow(a[k],(float)i);
  return(suma1);
}

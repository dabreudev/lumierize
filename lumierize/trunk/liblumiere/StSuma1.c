#include <math.h>
#include "stmedia.h"


float StSuma1(int n, float *a, int i)

{
  int k;
  float suma1=0;

  if(i == 0) return ((float)n);

  for(k=0; k<n; k++) suma1 += pow(a[k],(float)i);
  return(suma1);
}

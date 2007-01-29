#include "modulos.h"

void MinMax_d(int n, double *a, double *min, double *max)

{
  int i;

  *min=a[0];
  *max=a[0];
  for(i=1; i<n; i++) {
    if(a[i] > *max) *max=a[i];
    if(a[i] < *min) *min=a[i];
    }
}

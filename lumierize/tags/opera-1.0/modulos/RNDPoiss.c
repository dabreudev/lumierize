#include "modulos.h"

int RNDPoiss(int l)

{

  float p=0.0;
  int r=0;

  while (p<l) {
    p -= logf((float)rand()/RAND_MAX);
    r++;
    }

  return(--r);
}

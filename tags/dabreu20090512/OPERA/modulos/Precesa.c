#include "modulos.h"


void Precesa(double *p, float a1, float d1, float *a2, float *d2)

{
  float x[3],y[3];
  int n,j;

  /* Paso AR y DEC a vector */
  x[0] = cosf(a1) * cosf(d1);
  x[1] = sinf(a1) * cosf(d1);
  x[2] = sinf(d1);

  /* Multiplico por la matriz de precesion */
  for (n=0; n<3; n++) {
    y[n] = 0.0;
    for (j=0; j<3; j++) y[n] = y[n] + p[n*3+j] * x[j];
    }

  /* Y vuelvo a pasar a AR y DEC */
  *a2 = atanf(y[1]/y[0]);
  if(y[0] < 0) *a2 = *a2 + M_PI;
  if( *a2 < 0) *a2 = *a2 + 2*M_PI;
  *d2 = asinf(y[2]);

}

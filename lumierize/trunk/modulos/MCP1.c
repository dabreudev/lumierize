#include "modulos.h"

float MCP1(int n, float *x, float *y, float *a, float *b)

{
  int i;
  float sx,sy,sxx,sxy,d,s;

  sx= StSuma1(n,x,1);
  sy= StSuma1(n,y,1);
  sxx=StSuma1(n,x,2);
  sxy=StSuma2(n,x,1,y,1);
  d=n*sxx-sx*sx;

  *a=(sy*sxx-sxy*sx)/d;
  *b=(n*sxy-sx*sy)/d;

/*  Calculo de la desviacion estandar */
/*  --------------------------------- */
  d=0;
  for(i=0; i<n; i++) {
    s = y[i]-(*a + *b * x[i]);
    d += s*s;
    }
  d = sqrtf(d/(n-2.0));

  return(d);
}

#include "modulos.h"

double MCP1_d(int n, double *x, double *y, double *a, double *b)

{
  int i;
  double sx,sy,sxx,sxy,d,s;

  sx= StSuma1_d(n,x,1);
  sy= StSuma1_d(n,y,1);
  sxx=StSuma1_d(n,x,2);
  sxy=StSuma2_d(n,x,1,y,1);
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

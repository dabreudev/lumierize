#include "modulos.h"

void Ecu2Plac(double a, double d, double ac, double dc,
              double *x, double *y)

{
  double h;

  h = sin(d)*sin(dc)+cos(d)*cos(dc)*cos(a-ac);
  *x = cos(d)*sin(a-ac)/h;
  *y = (sin(d)*cos(dc)-cos(d)*sin(dc)*cos(a-ac))/h;
  
  return;

}

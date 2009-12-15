#include <math.h>

void Plac2Ecu(double x, double y, double ac, double dc,
              double *a, double *d)

{
  *a = ac+atan(x/(cos(dc)-y*sin(dc)));
  *d = atan((sin(dc)+y*cos(dc))/(cos(dc)-y*sin(dc))*cos(*a-ac));
  if(*a > 2*M_PI) *a -= 2*M_PI;
  return;
}


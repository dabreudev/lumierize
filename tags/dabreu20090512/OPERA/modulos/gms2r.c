#include <math.h>

float gms2r(char dsig, int dg, int dm, float ds)

{
  float t;

  t = (dg + dm/60.0 + ds/3600.0) * M_PI / 180.0;
  if(dsig == '-') t = -t;
  return(t);
}

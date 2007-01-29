#include <math.h>

void Ecu2Cart(float ar, float dec, float *v)

{
  v[0] = cos(ar)*cos(dec);
  v[1] = sin(ar)*cos(dec);
  v[2] = sin(dec);
}

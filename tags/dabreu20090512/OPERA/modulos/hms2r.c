#include <math.h>

float hms2r(int ah, int am, float as)

{
  float r; 
  r = ((float)ah + (float)am/60.0 + (float)as/3600.0) * M_PI / 12.0;
  return(r);
}

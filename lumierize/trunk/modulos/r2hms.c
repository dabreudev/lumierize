#include <math.h>

void r2hms(double r, int *ah, int *am, float *as)

{
  double rr;

  rr = r*12/M_PI;
  *ah = (int) rr;
  *am = (int) ((rr - *ah)*60);
  *as = (rr - *ah - *am / 60.0)*3600.0;
  if(*as>59.99) {
    *as = 0.00;
    (*am) ++;
    if(*am > 59.99) {
      *am = 0;
      (*ah) ++;
      if(*ah == 24) *ah = 0;
      }
    }

}

#include "modulos.h"


void r2gms(double r, char *dsig, int *dg, int *dm, float *ds)

{
  double rr;

  rr = r*180/M_PI;
  *dsig = '+';
  if(rr < 0) { *dsig = '-' ; rr = -rr; };

  *dg = (int) rr;
  *dm = (int) ((rr - *dg)*60);
  *ds = (rr - *dg - *dm / 60.0)*3600.0;
  if(*ds>59.995) {
    *ds = 0.00;
    (*dm)++;
    if(*dm > 59.995) {
      *dm = 0;
      (*dg)++;
      }
    }

}

#include "modulos.h"

double *MPrec(float datei, float datef)

{
  float ti,tf;
  double gio,zeta,teta;
  double *p;

  if((p=malloc(9*sizeof(double))) == NULL) {
    fprintf(stderr,"MPrec: ERROR. No puedo dimensionar la matriz de %d bytes\n",
	    9*sizeof(double));
    return(NULL);
    }

  ti = (datei - 2000.0) / 100.0;
  tf = (datef - 2000.0 - 100*ti) / 100.0;


/* Calculo de los elementos precesionales */
/* -------------------------------------- */

  gio = ((2306.2181 + 1.39656 * ti - 0.000139 * ti * ti) * tf +
         (0.30188 - 0.000344 * ti) * tf * tf + 0.017998 * tf * tf) / 3600.0;
  zeta = gio + ((0.7928 + 0.00041 * ti) * tf * tf +
         0.000205 * tf * tf * tf) / 3600.0;
  teta = ((2004.3109 - 0.8533 * ti - 0.000217 * ti * ti) * tf -
          (0.42665 + 0.000217 * ti) * tf * tf -
          0.041833 * tf * tf * tf ) / 3600.0;

  gio = gio * M_PI / 180.0;
  zeta = zeta * M_PI / 180.0;
  teta = teta * M_PI / 180.0;


/* Calculo de la matriz de precesion */
/* --------------------------------- */

  p[0] = -sin(gio) * sin(zeta) + cos(gio) * cos(teta) * cos(zeta);
  p[1] = -cos(gio) * sin(zeta) - sin(gio) * cos(teta) * cos(zeta);
  p[2] = -sin(teta) * cos(zeta);

  p[3] =  sin(gio) * cos(zeta) + cos(gio) * cos(teta) * sin(zeta);
  p[4] =  cos(gio) * cos(zeta) - sin(gio) * cos(teta) * sin(zeta);
  p[5] = -sin(teta) * sin(zeta);

  p[6] =  cos(gio) * sin(teta);
  p[7] = -sin(gio) * sin(teta);
  p[8] = cos(teta);

  return(p);

}

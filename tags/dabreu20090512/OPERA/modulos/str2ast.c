#include "modulos.h"


void strtoast(char *c, char *ident, float *ar, float *dec)

{
  char  i[3],a[7],dsig;
  int   i1,i2;
  float i3;

  strncpy(ident,c+1,7); ident[7] = '\0';  /* Leo la Indentificacion */

  strncpy(i,c+ 9,2); i[2] = '\0'; i1 = atoi(i);  /* Leo la AR y paso a rad. */
  strncpy(i,c+12,2); i[2] = '\0'; i2 = atoi(i);
  strncpy(a,c+14,6); a[6] = '\0'; i3 =  atof(a);
  *ar = hms2r(i1,i2,i3);

  dsig = *(c+22);                                /* Leo la DEC y paso a rad.*/
  strncpy(i,c+23,2); i[2] = '\0'; i1 = atoi(i);
  strncpy(i,c+26,2); i[2] = '\0'; i2 = atoi(i);
  strncpy(a,c+28,6); a[6] = '\0'; i3 = atof(a);
  *dec = gms2r(dsig,i1,i2,i3);

}

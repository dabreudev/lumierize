#include "modulos.h"


char *r2hmsch(float r)

{
  int h,m;
  float s;
  char *ch;

  if((ch=malloc(13))== NULL) {
    fprintf(stderr,"r2hmsch: ERROR. No puedo dimensionar ch de 13 bytes\n");
    return(NULL);
    }
  r2hms(r,&h,&m,&s);
  sprintf(ch," %2d %2d %5.2f",h,m,s);
  return(ch);
}

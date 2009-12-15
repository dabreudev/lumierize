#include "modulos.h"


char *r2gmsch(float r)

{
  int g,m;
  float s;
  char sig;
  char *ch;

  if((ch=malloc(13))== NULL) {
    fprintf(stderr,"r2hmsch: ERROR. No puedo dimensionar ch de 13 bytes\n");
    return(NULL);
    }
  r2gms(r,&sig,&g,&m,&s);
  sprintf(ch,"%c%2d %2d %5.1f ",sig,g,m,s);
  return(ch);
}

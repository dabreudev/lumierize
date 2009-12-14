#include "modulos.h"

void f_clrl(char label[][71])

{  
  char blanco[71];
  int i;

  for(i=0; i<70; i++) blanco[i]=' ';
  blanco[70]='\0';

  for(i=0; i<MAX_KEYLAB; i++) sprintf(label[i],"%s",blanco);
}

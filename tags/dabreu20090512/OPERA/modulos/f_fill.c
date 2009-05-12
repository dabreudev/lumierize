#include "modulos.h"

int f_fill(FILE *fp)

{

  int i,nbl,nbf;

  if(fp == NULL) {
    fprintf(stderr,"f_fill: ERROR. El fichero no esta abierto\n");
    return(0);
    }

  nbl=ftell(fp)/2880;
  nbf=2880-(ftell(fp)-nbl*2880);

  for(i=0;i<nbf;i++) putc('\0',fp);
  return(1);
}

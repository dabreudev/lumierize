#include "modulos.h"

int f_kff(FILE *fp,char *pclave,float *var)

{
  double a;
  int i;

  if((i=f_kfd(fp,pclave,&a))) *var=(float)a;
  return(i);
}

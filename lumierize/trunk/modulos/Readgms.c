#include "modulos.h"


float Readgms(void)

{
  int dg,dm;
  float ds;
  char dsig,ch1[4],ch2[4]="+";

  printf("Introduce declination in format gg mm ss.ss ");
  scanf("%s %d %f",ch1,&dm,&ds);
  if(ch1[0] != '-') {
    strcat(ch2,ch1);
    strcpy(ch1,ch2);
    }
  sscanf(ch1,"%c%d",&dsig,&dg);
  return(gms2r(dsig,dg,dm,ds));
}

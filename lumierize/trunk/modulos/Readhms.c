#include "modulos.h"


float Readhms(void)

{
  int ah,am;
  float as;

  printf("Introduce RA in format hh mm ss.ss ");
  scanf("%d %d %f",&ah,&am,&as);
  return(hms2r(ah,am,as));

}

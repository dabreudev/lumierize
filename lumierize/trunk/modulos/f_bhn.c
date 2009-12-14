#include "modulos.h"

int f_bhn(FILE *fp)

{
  int nh=0;
  char bloq[2880];


  rewind(fp);
  do {
    if(fread(bloq,2880,1,fp) == 0) {
      fprintf(stderr,"f_bhn: WARNING. He llegado al final \n");
      fprintf(stderr,"del fichero sin encontrar END\n");
      return(0);
      }
    nh++;
    } while(strstr(bloq,"END        ") == NULL);
  rewind(fp);
  return(nh);
}

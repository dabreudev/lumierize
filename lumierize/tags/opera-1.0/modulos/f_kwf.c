#include "modulos.h"

int f_kwf(FILE *fp, char *pclave, float var, char *comment)

{

  if(strlen(pclave) != 9) {
    fprintf(stderr,"f_kwf: ERROR. La Keyword debe tener 9 caracteres\n");
    fprintf(stderr,"               Keyword = >%s<\n",pclave);
    return(0);
    }

  if(strlen(comment) > 47) {
    fprintf(stderr,"f_kwf: ERROR. El comentario debe tener menos de ");
    fprintf(stderr,"47 caracteres\n");
    fprintf(stderr,"               Comment = >%s<\n",comment);
    return(0);
    }

  fprintf(fp,"%s %20.8g / %-47s",pclave,var,comment);
  return(1);
}

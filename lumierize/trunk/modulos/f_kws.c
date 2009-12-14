#include "modulos.h"

int f_kws(FILE *fp, char *pclave, char *var, char *comment)

{

  char ch21[21];

  if(strlen(pclave) != 9) {
    fprintf(stderr,"f_kws: ERROR. La Keyword debe tener 9 caracteres\n");
    fprintf(stderr,"               Keyword = >%s<\n",pclave);
    return(0);
    }

  if(strlen(comment) > 47) {
    fprintf(stderr,"f_kws: ERROR. El comentario debe tener menos de ");
    fprintf(stderr,"47 caracteres\n");
    fprintf(stderr,"               Comment = >%s<\n",comment);
    return(0);
    }

  if(strlen(var) > 18) {
    fprintf(stderr,"f_kws: ERROR. la variable debe tener menos de ");
    fprintf(stderr,"18 caracteres\n");
    fprintf(stderr,"              %s: Variable = >%s<\n",pclave,var);
    return(0);
    }

  sprintf(ch21,"'%s'",var);
  fprintf(fp,"%s %-20s / %-47s",pclave,ch21,comment);
  return(1);
}

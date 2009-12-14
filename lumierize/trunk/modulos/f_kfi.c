#include "modulos.h"

int f_kfi(FILE *fp,char *pclave,int *var)

{
  char *p,bloq[2880];
  int i;

  rewind(fp);
  do {
    i=fread(bloq,1,2880,fp);
    if((p=strstr(bloq,pclave)) != NULL) *var=atoi(p+9);
    } while ((p == NULL) && (strstr(bloq,"END      ") == NULL) && (i != 0));
  rewind(fp);

  if(p!=NULL) return(1);
  else {
   printf("f_kfi: WARNING. No he encontrado la Keyword %s\n",pclave);
   return(0);
   }
}
